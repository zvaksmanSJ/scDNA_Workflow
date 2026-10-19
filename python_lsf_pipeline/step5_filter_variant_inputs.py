#!/usr/bin/env python3
"""Filter small variants and Manta SVs for BASEJ-SOMATIC/lineage input."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


PASS_VALUES = {"PASS"}
BND_ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


@dataclass(frozen=True)
class SvEvent:
    sample: str
    svtype: str
    orientation: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    ref_support: int
    alt_support: int
    baf: float


def run(command: list[str], **kwargs) -> None:
    subprocess.run(command, check=True, **kwargs)


def parse_info(value: str) -> dict[str, str]:
    result = {}
    for item in value.split(";"):
        key, separator, item_value = item.partition("=")
        result[key] = item_value if separator else ""
    return result


def parse_sample(format_value: str, sample_value: str) -> dict[str, str]:
    return dict(zip(format_value.split(":"), sample_value.split(":")))


def parse_number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def vep_allele_candidates(ref: str, alt: str) -> set[str]:
    candidates = {alt}
    trimmed_ref = ref
    trimmed_alt = alt
    while trimmed_ref and trimmed_alt and trimmed_ref[-1] == trimmed_alt[-1]:
        trimmed_ref = trimmed_ref[:-1]
        trimmed_alt = trimmed_alt[:-1]
    while trimmed_ref and trimmed_alt and trimmed_ref[0] == trimmed_alt[0]:
        trimmed_ref = trimmed_ref[1:]
        trimmed_alt = trimmed_alt[1:]
    candidates.add(trimmed_alt or "-")
    return candidates


def is_coding_csq(
    csq_value: str,
    csq_fields: list[str],
    ref: str,
    alt: str,
) -> bool:
    field_index = {name: index for index, name in enumerate(csq_fields)}
    required = ("Allele", "BIOTYPE", "CDS_position")
    if any(name not in field_index for name in required):
        raise ValueError(f"CSQ header is missing required fields: {required}")
    allele_candidates = vep_allele_candidates(ref, alt)
    for annotation in csq_value.split(","):
        values = annotation.split("|")
        values.extend([""] * (len(csq_fields) - len(values)))
        if (
            values[field_index["Allele"]] in allele_candidates
            and values[field_index["BIOTYPE"]] == "protein_coding"
            and values[field_index["CDS_position"]] not in {"", ".", "-"}
        ):
            return True
    return False


def write_bgzip_vcf(
    output_path: Path,
    headers: list[str],
    records: list[str],
    bcftools: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
        temp_path = handle.name
        handle.writelines(headers)
        handle.writelines(records)
    try:
        run([bcftools, "view", "-Oz", "-o", str(output_path), temp_path])
        run([bcftools, "index", "-f", "-t", str(output_path)])
    finally:
        os.unlink(temp_path)


def filter_small_vcf(
    sample: str,
    input_path: Path,
    output_path: Path,
    chroms: set[str],
    min_depth: int,
    min_vaf: float,
    max_vaf: float,
    bcftools: str,
) -> Counter:
    counts = Counter()
    with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as handle:
        normalized_path = handle.name
    try:
        run(
            [
                bcftools,
                "norm",
                "-m",
                "-any",
                "-Oz",
                "-o",
                normalized_path,
                str(input_path),
            ]
        )
        headers: list[str] = []
        records: list[str] = []
        csq_fields: list[str] = []
        with gzip.open(normalized_path, "rt", encoding="utf-8") as source:
            for line in source:
                if line.startswith("##INFO=<ID=CSQ"):
                    match = re.search(r"Format: ([^\">]+)", line)
                    if not match:
                        raise ValueError(f"Cannot parse CSQ format header in {input_path}")
                    csq_fields = match.group(1).split("|")
                if line.startswith("#"):
                    headers.append(line)
                    continue
                counts["input"] += 1
                fields = line.rstrip("\n").split("\t")
                if chroms and fields[0] not in chroms:
                    counts["non_target_chromosome"] += 1
                    continue
                if fields[6] not in PASS_VALUES:
                    counts["non_pass"] += 1
                    continue
                sample_data = parse_sample(fields[8], fields[9])
                depth = parse_number(sample_data.get("DP", ""))
                if depth is None or depth < min_depth:
                    counts["depth_below_minimum"] += 1
                    continue
                vaf_values = [
                    value
                    for item in sample_data.get("VAF", "").split(",")
                    if (value := parse_number(item)) is not None
                ]
                if not vaf_values:
                    ad = [
                        int(value)
                        for value in sample_data.get("AD", "").split(",")
                        if value.isdigit()
                    ]
                    if len(ad) >= 2 and sum(ad) > 0:
                        vaf_values = [ad[1] / sum(ad)]
                if not any(min_vaf <= value <= max_vaf for value in vaf_values):
                    counts["vaf_outside_range"] += 1
                    continue
                info = parse_info(fields[7])
                if not csq_fields or "CSQ" not in info:
                    counts["missing_csq"] += 1
                    continue
                if not is_coding_csq(info["CSQ"], csq_fields, fields[3], fields[4]):
                    counts["not_protein_coding_cds"] += 1
                    continue
                records.append(line)
                counts["retained"] += 1
        write_bgzip_vcf(output_path, headers, records, bcftools)
    finally:
        os.unlink(normalized_path)
    print(f"[OK] {sample}: retained {counts['retained']}/{counts['input']} small variants")
    return counts


def contig_key(contig: str) -> tuple[int, int | str]:
    clean = contig.removeprefix("chr")
    if clean.isdigit():
        return (0, int(clean))
    return (1, clean)


def canonical_endpoints(
    chrom1: str,
    pos1: int,
    chrom2: str,
    pos2: int,
) -> tuple[str, int, str, int]:
    left = (contig_key(chrom1), pos1, chrom1)
    right = (contig_key(chrom2), pos2, chrom2)
    if right < left:
        return chrom2, pos2, chrom1, pos1
    return chrom1, pos1, chrom2, pos2


def parse_sv_event(
    sample: str,
    fields: list[str],
    min_vaf: float,
    max_vaf: float,
) -> tuple[SvEvent | None, str]:
    if fields[6] not in PASS_VALUES:
        return None, "non_pass"
    info = parse_info(fields[7])
    sample_data = parse_sample(fields[8], fields[9])
    if sample_data.get("FT") not in PASS_VALUES:
        return None, "sample_non_pass"
    try:
        pr = [int(value) for value in sample_data["PR"].split(",")]
        sr = [int(value) for value in sample_data["SR"].split(",")]
    except (KeyError, ValueError):
        return None, "missing_support"
    if len(pr) < 2 or len(sr) < 2:
        return None, "missing_support"
    ref_support = pr[0] + sr[0]
    alt_support = pr[1] + sr[1]
    total_support = ref_support + alt_support
    if total_support == 0:
        return None, "zero_support"
    baf = alt_support / total_support
    if not min_vaf <= baf <= max_vaf:
        return None, "vaf_outside_range"
    svtype = info.get("SVTYPE", "BND" if BND_ALT.search(fields[4]) else "UNK")
    chrom1 = fields[0]
    pos1 = int(fields[1])
    if svtype == "BND":
        match = BND_ALT.search(fields[4])
        if not match:
            return None, "missing_partner"
        chrom2, pos2_text = match.groups()
        pos2 = int(pos2_text)
        side = "L" if fields[4][0] in "[]" else "R"
        bracket = "OPEN" if "[" in fields[4] else "CLOSE"
        orientation = f"{side}-{bracket}"
    else:
        chrom2 = info.get("CHR2", chrom1)
        try:
            pos2 = int(info.get("END", pos1))
        except ValueError:
            return None, "missing_partner"
        orientation = "."
    chrom1, pos1, chrom2, pos2 = canonical_endpoints(chrom1, pos1, chrom2, pos2)
    return (
        SvEvent(
            sample=sample,
            svtype=svtype,
            orientation=orientation,
            chrom1=chrom1,
            pos1=pos1,
            chrom2=chrom2,
            pos2=pos2,
            ref_support=ref_support,
            alt_support=alt_support,
            baf=baf,
        ),
        "retained",
    )


def filter_sv_vcf(
    sample: str,
    input_path: Path,
    output_path: Path,
    min_vaf: float,
    max_vaf: float,
    bcftools: str,
) -> tuple[list[SvEvent], Counter]:
    counts = Counter()
    events: list[SvEvent] = []
    headers: list[str] = []
    records: list[str] = []
    with gzip.open(input_path, "rt", encoding="utf-8") as source:
        for line in source:
            if line.startswith("#"):
                headers.append(line)
                continue
            counts["input"] += 1
            fields = line.rstrip("\n").split("\t")
            event, verdict = parse_sv_event(sample, fields, min_vaf, max_vaf)
            counts[verdict] += 1
            if event is None:
                continue
            events.append(event)
            records.append(line)
    write_bgzip_vcf(output_path, headers, records, bcftools)
    print(f"[OK] {sample}: retained {counts['retained']}/{counts['input']} SV records")
    return events, counts


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def cluster_sv_events(
    events: list[SvEvent],
    tolerance: int,
) -> list[list[SvEvent]]:
    grouped: dict[tuple[str, str, str, str], list[SvEvent]] = defaultdict(list)
    for event in events:
        grouped[
            (event.svtype, event.orientation, event.chrom1, event.chrom2)
        ].append(event)
    clusters: list[list[SvEvent]] = []
    for group_events in grouped.values():
        ordered = sorted(group_events, key=lambda event: (event.pos1, event.pos2))
        union_find = UnionFind(len(ordered))
        first_active = 0
        for right, event in enumerate(ordered):
            while (
                first_active < right
                and event.pos1 - ordered[first_active].pos1 > tolerance
            ):
                first_active += 1
            for left in range(first_active, right):
                candidate = ordered[left]
                if abs(event.pos2 - candidate.pos2) <= tolerance:
                    union_find.union(left, right)
        members: dict[int, list[SvEvent]] = defaultdict(list)
        for index, event in enumerate(ordered):
            members[union_find.find(index)].append(event)
        clusters.extend(members.values())
    return clusters


def write_sv_matrices(
    clusters: list[list[SvEvent]],
    samples: list[str],
    output_dir: Path,
    tolerance: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    nr_rows = []
    nv_rows = []
    binary_rows = []
    feature_rows = []
    cell_feature_rows = []
    for cluster in clusters:
        representative = cluster[0]
        pos1 = int(statistics.median(event.pos1 for event in cluster))
        pos2 = int(statistics.median(event.pos2 for event in cluster))
        subtype = representative.svtype
        if representative.orientation != ".":
            subtype = f"{subtype}:{representative.orientation}"
        variant_id = (
            f"{representative.chrom1}_{pos1}_N_"
            f"<{subtype}:{representative.chrom2}:{pos2}>"
        )
        by_sample: dict[str, SvEvent] = {}
        for event in cluster:
            current = by_sample.get(event.sample)
            if current is None or event.alt_support > current.alt_support:
                by_sample[event.sample] = event
        nr_rows.append(
            [variant_id]
            + [
                by_sample[sample].ref_support + by_sample[sample].alt_support
                if sample in by_sample
                else 0
                for sample in samples
            ]
        )
        nv_rows.append(
            [variant_id]
            + [
                by_sample[sample].alt_support if sample in by_sample else 0
                for sample in samples
            ]
        )
        binary_rows.append(
            [variant_id] + [1 if sample in by_sample else 0 for sample in samples]
        )
        feature_rows.append(
            [
                variant_id,
                representative.svtype,
                representative.orientation,
                representative.chrom1,
                pos1,
                representative.chrom2,
                pos2,
                len(cluster),
                len(by_sample),
                tolerance,
            ]
        )
        for sample, event in sorted(by_sample.items()):
            cell_feature_rows.append(
                [
                    variant_id,
                    sample,
                    event.chrom1,
                    event.pos1,
                    event.chrom2,
                    event.pos2,
                    event.ref_support,
                    event.alt_support,
                    f"{event.baf:.6f}",
                ]
            )
    matrices = {
        "sv_NR.tsv": nr_rows,
        "sv_NV.tsv": nv_rows,
        "sv_binary.tsv": binary_rows,
    }
    for filename, rows in matrices.items():
        with open(output_dir / filename, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["VariantId", *samples])
            writer.writerows(rows)
    with open(
        output_dir / "sv_features.tsv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "VariantId",
                "SVTYPE",
                "orientation",
                "CHROM1",
                "POS1",
                "CHROM2",
                "POS2",
                "supporting_records",
                "supporting_cells",
                "breakpoint_tolerance",
            ]
        )
        writer.writerows(feature_rows)
    with open(
        output_dir / "sv_cell_features.tsv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "VariantId",
                "sample",
                "CHROM1",
                "POS1",
                "CHROM2",
                "POS2",
                "ref_support",
                "alt_support",
                "BAF",
            ]
        )
        writer.writerows(cell_feature_rows)


def resolve_manifest_path(path: str, base_dir: Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else base_dir / candidate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--manifest-base-dir", required=True)
    parser.add_argument("--variants-dir", required=True)
    parser.add_argument("--sv-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--group", default="somatic_cells")
    parser.add_argument("--chroms", default="21")
    parser.add_argument("--small-variant-method", choices=("legacy", "reference"), default="legacy",
                        help="reference delegates genotype, normalization and cohort PASS gates to BaseJ; no coding/DP/VAF prefilter")
    parser.add_argument("--min-depth", type=int, default=5)
    parser.add_argument("--small-vaf-min", type=float, default=0.25)
    parser.add_argument("--small-vaf-max", type=float, default=0.85)
    parser.add_argument("--sv-vaf-min", type=float, default=0.25)
    parser.add_argument("--sv-vaf-max", type=float, default=0.85)
    parser.add_argument("--sv-breakpoint-tolerance", type=int, default=500)
    parser.add_argument("--bcftools", default="bcftools")
    args = parser.parse_args()

    for lower, upper, label in (
        (args.small_vaf_min, args.small_vaf_max, "small-variant VAF"),
        (args.sv_vaf_min, args.sv_vaf_max, "SV VAF"),
    ):
        if not 0 <= lower <= upper <= 1:
            parser.error(f"{label} bounds must satisfy 0 <= minimum <= maximum <= 1")
    if args.min_depth < 0 or args.sv_breakpoint_tolerance < 0:
        parser.error("Depth and breakpoint tolerance cannot be negative")

    manifest_base = Path(args.manifest_base_dir).resolve()
    variants_dir = Path(args.variants_dir).resolve()
    sv_dir = Path(args.sv_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    small_output_dir = output_dir / "small_variants"
    sv_output_dir = output_dir / "structural_variants"
    chroms = (
        set() if args.chroms.lower() == "all"
        else {value.strip().removeprefix("chr") for value in args.chroms.split(",")}
    )

    manifest_rows = []
    with open(args.manifest, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("pass_qc", "true").lower() in {"true", "1", "yes"}:
                manifest_rows.append(row)
    if not manifest_rows:
        raise SystemExit("[ERROR] No passing samples in manifest")

    samples = [row["sample_id"] for row in manifest_rows]
    csv_rows = []
    small_summaries = []
    sv_summaries = []
    all_sv_events: list[SvEvent] = []
    for row in manifest_rows:
        sample = row["sample_id"]
        bam = resolve_manifest_path(row["bam"], manifest_base).resolve()
        if not bam.is_file():
            raise SystemExit(f"[ERROR] Missing BAM for {sample}: {bam}")
        small_input = (
            variants_dir
            / sample
            / "support"
            / "merged"
            / (f"{sample}.deepvariant.raw.vcf.gz" if args.small_variant_method == "reference"
               else f"{sample}.deepvariant.raw.vep.vcf.gz")
        )
        sv_input = sv_dir / sample / "manta" / "results" / "variants" / "diploidSV.vcf.gz"
        if not small_input.is_file():
            raise SystemExit(f"[ERROR] Missing annotated VCF for {sample}: {small_input}")
        if not sv_input.is_file():
            raise SystemExit(f"[ERROR] Missing Manta VCF for {sample}: {sv_input}")

        small_output = small_output_dir / (
            f"{sample}.reference_candidates.vcf.gz" if args.small_variant_method == "reference"
            else f"{sample}.coding.filtered.vcf.gz"
        )
        if args.small_variant_method == "reference":
            small_output.parent.mkdir(parents=True, exist_ok=True)
            command = [args.bcftools, "view", "-Oz", "-o", str(small_output)]
            if chroms:
                command.extend(["-r", ",".join(sorted(chroms, key=contig_key))])
            run([*command, str(small_input)])
            run([args.bcftools, "index", "-f", "-t", str(small_output)])
            with gzip.open(small_input, "rt") as source:
                n_input = sum(not line.startswith("#") for line in source)
            with gzip.open(small_output, "rt") as source:
                n_retained = sum(not line.startswith("#") for line in source)
            small_counts = Counter(input=n_input, retained=n_retained)
            print(f"[OK] {sample}: {n_retained} candidates; biological gates delegated to BaseJ")
        else:
            small_counts = filter_small_vcf(
                sample, small_input, small_output, chroms, args.min_depth,
                args.small_vaf_min, args.small_vaf_max, args.bcftools,
            )
        small_summaries.append((sample, small_counts))

        sv_output = sv_output_dir / f"{sample}.filtered.sv.vcf.gz"
        sv_events, sv_counts = filter_sv_vcf(
            sample,
            sv_input,
            sv_output,
            args.sv_vaf_min,
            args.sv_vaf_max,
            args.bcftools,
        )
        all_sv_events.extend(sv_events)
        sv_summaries.append((sample, sv_counts))
        csv_rows.append(
            {
                "biosampleName": sample,
                "bam": str(bam),
                "group": args.group,
                "vcf": str(small_output),
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("biosampleName", "bam", "group", "vcf")
        )
        writer.writeheader()
        writer.writerows(csv_rows)

    clusters = cluster_sv_events(all_sv_events, args.sv_breakpoint_tolerance)
    write_sv_matrices(clusters, samples, sv_output_dir, args.sv_breakpoint_tolerance)

    all_small_keys = sorted({key for _, counts in small_summaries for key in counts})
    with open(
        output_dir / "small_variant_filter_summary.tsv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *all_small_keys])
        for sample, counts in small_summaries:
            writer.writerow([sample, *[counts[key] for key in all_small_keys]])

    all_sv_keys = sorted({key for _, counts in sv_summaries for key in counts})
    with open(
        output_dir / "sv_filter_summary.tsv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *all_sv_keys])
        for sample, counts in sv_summaries:
            writer.writerow([sample, *[counts[key] for key in all_sv_keys]])

    with open(output_dir / "filter_config.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "small_variants": {
                    "method": args.small_variant_method,
                    "chromosomes": sorted(chroms, key=contig_key),
                    "required_filter": "cohort N_PASS(MLV==0)>=1" if args.small_variant_method == "reference" else "PASS",
                    "coding_definition": (
                        None if args.small_variant_method == "reference"
                        else "VEP BIOTYPE=protein_coding with non-empty CDS_position"
                    ),
                    "minimum_depth": None if args.small_variant_method == "reference" else args.min_depth,
                    "minimum_vaf": None if args.small_variant_method == "reference" else args.small_vaf_min,
                    "maximum_vaf": None if args.small_variant_method == "reference" else args.small_vaf_max,
                },
                "structural_variants": {
                    "region_restriction": None,
                    "required_record_filter": "PASS",
                    "required_sample_filter": "PASS",
                    "minimum_vaf": args.sv_vaf_min,
                    "maximum_vaf": args.sv_vaf_max,
                    "vaf_formula": "(PR_alt + SR_alt) / sum(PR_ref, PR_alt, SR_ref, SR_alt)",
                    "breakpoint_tolerance_bp": args.sv_breakpoint_tolerance,
                    "matching": "same SVTYPE, BND orientation, and chromosome pair",
                },
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    print(f"[OK] Wrote filtered BASEJ input CSV: {args.output_csv}")
    print(f"[OK] Wrote {len(clusters)} clustered SV features to {sv_output_dir}")


if __name__ == "__main__":
    main()
