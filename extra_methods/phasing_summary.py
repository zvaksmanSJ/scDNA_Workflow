#!/usr/bin/env python3
"""Summarize per-cell Beagle phasing results into a cohort table.

Reports, per cell: whether the phased VCF exists, total phased genotype
count, and the fraction of heterozygous genotypes that received a real
(non-"|.") phase call -- i.e. how much of the het-SNP backbone needed for
downstream LOH/haplotype-block work is actually phased.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def summarize_phased_vcf(vcf_path: Path) -> dict:
    n_records = 0
    n_het = 0
    n_het_phased = 0
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            gt = parts[9].split(":", 1)[0]
            n_records += 1
            sep = "|" if "|" in gt else ("/" if "/" in gt else None)
            if sep is None:
                continue
            alleles = gt.split(sep)
            if len(alleles) == 2 and alleles[0] != alleles[1] and "." not in alleles:
                n_het += 1
                if sep == "|":
                    n_het_phased += 1
    return {
        "n_records": n_records,
        "n_het": n_het,
        "n_het_phased": n_het_phased,
        "het_phased_fraction": (n_het_phased / n_het) if n_het else float("nan"),
    }


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phasing-dir", required=True, type=Path,
                    help="Output dir from phase_variants_lsf.py, containing <sample>/<sample>.phased.vcf.gz")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    fieldnames = ["sample_id", "n_records", "n_het", "n_het_phased",
                  "het_phased_fraction", "status"]
    rows = []
    for sample in samples:
        vcf_path = args.phasing_dir / sample / f"{sample}.phased.vcf.gz"
        if not vcf_path.is_file():
            rows.append({"sample_id": sample, "status": "MISSING"})
            continue
        try:
            result = summarize_phased_vcf(vcf_path)
            row = {"sample_id": sample, "status": "OK"}
            row.update({
                "n_records": result["n_records"],
                "n_het": result["n_het"],
                "n_het_phased": result["n_het_phased"],
                "het_phased_fraction": f"{result['het_phased_fraction']:.4f}",
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
    print(f"[OK] Phasing summarized for {len(ok)}/{len(rows)} samples -> {args.output}")


if __name__ == "__main__":
    main()
