#!/usr/bin/env python3
"""Cross-cell recurrent-artifact panel built from our own cohort.

The reference paper uses a cross-sample panel built from *other, unrelated*
individuals (52 single cells + 18 bulks from 17 additional people, via
dbGaP) to flag recurrent PTA/technical artifacts: a variant seen in several
unrelated people's independent single-cell amplifications is far more
likely a systematic artifact (e.g. a mapping/homopolymer trap) than a
coincidence, because unrelated individuals don't share rare variants.

We have no matched bulk or other individuals for this cohort (all 15 cells
are the same clonal cell line), so we instead build a same-cohort
recurrence panel using exact-site-and-allele recurrence across our own 15
cells. Because all cells share one genotype, this cannot perfectly
separate "systematic artifact" from "clonal germline" the way a
cross-individual panel can -- so we report three explicit buckets instead
of silently merging them:

  - UNIVERSAL   (variant PASS-called in essentially all callable cells):
                most parsimoniously explained as clonal germline, not FP.
  - SINGLETON   (PASS-called in exactly one cell): the expected signature
                of either a true private mutation or ordinary independent
                PTA noise; cannot be resolved further by recurrence alone.
  - RECURRENT_SUBSET (PASS-called in more than one but not (nearly) all
                cells, at the identical site+allele): the pattern an
                unrelated-individual panel is specifically designed to
                catch -- independent single-cell amplifications coincidentally
                producing an identical false call at the identical position
                is improbable, so recurrence in an intermediate subset is the
                most actionable signal that this site is a candidate
                systematic artifact (e.g. a recurrent alignment/homopolymer
                trap) rather than true subclonal biology. Flagged for
                exclusion from private-somatic-candidate lists pending
                further review; NOT silently dropped.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


def pass_sites(vcf_path: Path):
    result = subprocess.run(
        ["bcftools", "view", "-f", "PASS", "-H", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    for line in result.stdout.splitlines():
        fields = line.split("\t", 5)
        chrom, pos, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        for a in alt.split(","):
            if a in (".", "<NON_REF>"):
                continue
            yield (chrom, pos, ref, a)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variants-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--final-vcf-suffix", default=".deepvariant.raw.vep.snpeff.vcf.gz")
    p.add_argument("--universal-fraction", type=float, default=0.90,
                    help="Fraction of cells above which a site is called UNIVERSAL")
    p.add_argument("--output-sites", required=True, type=Path)
    p.add_argument("--output-summary", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    site_to_samples: dict[tuple, list[str]] = defaultdict(list)
    per_sample_total = {}

    for sample in samples:
        vcf_path = args.variants_dir / sample / f"{sample}{args.final_vcf_suffix}"
        if not vcf_path.is_file():
            print(f"[WARN] Missing VCF for {sample}, skipping")
            continue
        count = 0
        for site in pass_sites(vcf_path):
            site_to_samples[site].append(sample)
            count += 1
        per_sample_total[sample] = count
        print(f"[INFO] {sample}: {count} PASS site-alleles")

    n_cells = len(per_sample_total)
    universal_threshold = max(1, round(n_cells * args.universal_fraction))

    args.output_sites.parent.mkdir(parents=True, exist_ok=True)
    per_sample_bucket_counts = {s: {"UNIVERSAL": 0, "SINGLETON": 0, "RECURRENT_SUBSET": 0} for s in samples}

    with args.output_sites.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "n_cells", "classification", "cells"])
        for (chrom, pos, ref, alt), cells in site_to_samples.items():
            n = len(cells)
            if n >= universal_threshold:
                cls = "UNIVERSAL"
            elif n == 1:
                cls = "SINGLETON"
            else:
                cls = "RECURRENT_SUBSET"
            for sample in cells:
                per_sample_bucket_counts[sample][cls] += 1
            writer.writerow([chrom, pos, ref, alt, n, cls, ",".join(cells)])

    with args.output_summary.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample_id", "n_pass_total", "n_universal", "n_singleton",
                          "n_recurrent_subset_flagged_artifact"])
        for sample in samples:
            b = per_sample_bucket_counts.get(sample, {"UNIVERSAL": 0, "SINGLETON": 0, "RECURRENT_SUBSET": 0})
            writer.writerow([sample, per_sample_total.get(sample, ""),
                              b["UNIVERSAL"], b["SINGLETON"], b["RECURRENT_SUBSET"]])

    print(f"[OK] {len(site_to_samples)} distinct PASS sites across {n_cells} cells "
          f"-> {args.output_sites}")
    print(f"[OK] Per-sample summary -> {args.output_summary}")


if __name__ == "__main__":
    main()
