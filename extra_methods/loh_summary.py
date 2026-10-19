#!/usr/bin/env python3
"""Summarize Control-FREEC LOH calls into a per-cell cohort table.

Control-FREEC's "_CNVs" output reports, per segment: chromosome, start,
end, copy number, status (gain/loss/neutral), genotype, and uncertainty.
LOH is reported natively via the genotype column -- a single repeated
parental-allele letter (e.g. "A", "B", "AA", "BB") means only one parental
allele is present at that copy-number state (i.e. LOH), while a mixed
genotype (e.g. "AB", "AAB") means both parental alleles are retained, and
"-" means FREEC could not assign a genotype (typically too few informative
SNPs in the segment). No custom LOH classifier is implemented here -- this
script only aggregates FREEC's own native genotype field.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

LOH_GENOTYPE_RE = re.compile(r"^(A+|B+)$")


def is_loh_genotype(genotype: str) -> bool:
    return bool(LOH_GENOTYPE_RE.match(genotype))


def summarize_cnvs(cnvs_path: Path) -> dict:
    n_regions = 0
    n_loh = 0
    loh_bp = 0
    total_bp = 0
    with cnvs_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            start, end, genotype = int(parts[1]), int(parts[2]), parts[5]
            span = end - start
            n_regions += 1
            total_bp += span
            if is_loh_genotype(genotype):
                n_loh += 1
                loh_bp += span
    return {
        "n_regions": n_regions,
        "n_loh_regions": n_loh,
        "loh_region_fraction": (n_loh / n_regions) if n_regions else float("nan"),
        "loh_bp": loh_bp,
        "total_bp": total_bp,
        "loh_genome_fraction": (loh_bp / total_bp) if total_bp else float("nan"),
    }


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freec-dir", required=True, type=Path,
                    help="Output dir from loh_freec_lsf.py, containing <sample>/<sample>.dedup.bam_CNVs")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    fieldnames = ["sample_id", "n_regions", "n_loh_regions", "loh_region_fraction",
                  "loh_bp", "total_bp", "loh_genome_fraction", "status"]
    rows = []
    for sample in samples:
        cnvs_path = args.freec_dir / sample / f"{sample}.dedup.bam_CNVs"
        if not cnvs_path.is_file():
            rows.append({"sample_id": sample, "status": "MISSING_CNVS"})
            continue
        try:
            result = summarize_cnvs(cnvs_path)
            row = {"sample_id": sample, "status": "OK"}
            row.update({
                "n_regions": result["n_regions"],
                "n_loh_regions": result["n_loh_regions"],
                "loh_region_fraction": f"{result['loh_region_fraction']:.4f}",
                "loh_bp": result["loh_bp"],
                "total_bp": result["total_bp"],
                "loh_genome_fraction": f"{result['loh_genome_fraction']:.4f}",
            })
            rows.append(row)
        except Exception as exc:  # noqa: BLE001 - report per-sample, keep going
            rows.append({"sample_id": sample, "status": f"ERROR: {exc}"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    ok = [r for r in rows if r["status"] == "OK"]
    print(f"[OK] LOH summarized for {len(ok)}/{len(rows)} samples -> {args.output}")
    if ok:
        fracs = [float(r["loh_genome_fraction"]) for r in ok]
        print(f"[INFO] Cohort LOH genome fraction: "
              f"median={sorted(fracs)[len(fracs)//2]:.4f} "
              f"min={min(fracs):.4f} max={max(fracs):.4f}")


if __name__ == "__main__":
    main()
