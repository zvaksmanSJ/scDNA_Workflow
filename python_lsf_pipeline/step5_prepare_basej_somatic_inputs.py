#!/usr/bin/env python3
"""Create the VCF-start BASEJ-SOMATIC input CSV from local outputs."""

import argparse
import csv
import glob
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--variants-dir", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--group", default="somatic_cells")
    args = parser.parse_args()
    variants = {}
    for pattern in ("**/*.deepvariant.raw.vcf.gz", "**/*.raw.vcf.gz"):
        for path in glob.glob(os.path.join(args.variants_dir, pattern), recursive=True):
            name = Path(path).name
            sample = name.split(".deepvariant.raw.vcf.gz")[0].split(".raw.vcf.gz")[0]
            variants[sample] = os.path.abspath(path)
    if not variants:
        raise SystemExit("[ERROR] No merged DeepVariant VCFs found")
    rows = []
    with open(args.manifest, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("pass_qc", "true").lower() not in ("true", "1", "yes"):
                continue
            sample = row["sample_id"]
            bam = row["bam"]
            if not os.path.isabs(bam):
                bam = os.path.abspath(bam)
            vcf = variants.get(sample)
            if not vcf:
                raise SystemExit(f"[ERROR] Missing merged VCF for {sample}")
            rows.append({"biosampleName": sample, "bam": bam, "group": args.group, "vcf": vcf})
    if not rows:
        raise SystemExit("[ERROR] No passing cells found in manifest")
    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("biosampleName", "bam", "group", "vcf"))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] Wrote {len(rows)} BASEJ-SOMATIC inputs to {args.output_csv}")


if __name__ == "__main__":
    main()
