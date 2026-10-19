#!/usr/bin/env python3
"""Summarize per-cell Control-FREEC ploidy/purity estimates into a cohort table.

Reads each cell's own `_info.txt` output. When Control-FREEC is given a
single candidate ploidy (e.g. "ploidy=2") Output_Ploidy simply echoes that
fixed assumption back; when given multiple candidates (e.g.
"ploidy=2,3,4,5,6,7,8") Control-FREEC selects whichever value best explains
the observed copy-number-alteration profile (official feature -- see the
"ploidy" parameter in the Control-FREEC manual) and reports the selected
value here. No custom ploidy-estimation logic is implemented -- this script
only aggregates FREEC's own native output.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_info_txt(info_path: Path) -> dict:
    fields = {}
    with info_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                fields[parts[0]] = parts[1]
    return fields


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freec-dir", required=True, type=Path,
                    help="Output dir from loh_freec_lsf.py, containing <sample>/<sample>.dedup.bam_info.txt")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    fieldnames = ["sample_id", "output_ploidy", "sample_purity",
                  "good_polynomial_fit", "status"]
    rows = []
    for sample in samples:
        info_path = args.freec_dir / sample / f"{sample}.dedup.bam_info.txt"
        if not info_path.is_file():
            rows.append({"sample_id": sample, "status": "MISSING"})
            continue
        try:
            fields = parse_info_txt(info_path)
            rows.append({
                "sample_id": sample,
                "output_ploidy": fields.get("Output_Ploidy", ""),
                "sample_purity": fields.get("Sample_Purity", ""),
                "good_polynomial_fit": fields.get("Good_Polynomial_Fit", ""),
                "status": "OK",
            })
        except Exception as exc:  # noqa: BLE001 - report per-sample, keep going
            rows.append({"sample_id": sample, "status": f"ERROR: {exc}"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    ok = [r for r in rows if r["status"] == "OK"]
    print(f"[OK] Ploidy summarized for {len(ok)}/{len(rows)} samples -> {args.output}")
    if ok:
        ploidies = [r["output_ploidy"] for r in ok]
        from collections import Counter
        counts = Counter(ploidies)
        print(f"[INFO] Ploidy distribution across cohort: {dict(counts)}")


if __name__ == "__main__":
    main()
