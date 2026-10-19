#!/usr/bin/env python3
"""MAPD (Median Absolute Pairwise Difference) amplification-evenness QC.

MAPD is the single-cell WGA/PTA QC metric used by SCAN2 (and widely in
single-cell CNV tools such as 10x Genomics' CNV pipeline): for equal-sized
genomic bins, take the log2 read-depth ratio per bin, then compute the
median of the absolute difference between *consecutive* bins along each
chromosome. Because real copy-number signal is piecewise-constant across
many adjacent bins, this adjacent-bin-difference measure is dominated by
amplification noise (allelic dropout / uneven WGA) rather than by true
CNV structure -- lower MAPD means more even, higher-quality amplification.

This reads the per-cell CNVkit .cnr file (chromosome, start, end, log2,
weight) that our Step 3 variants/CNV stage already produces, so no new
upstream tool is required -- this is purely a QC calculation layered on
existing output.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd


def compute_mapd(cnr_path: Path, min_weight: float = 0.0) -> dict:
    df = pd.read_csv(cnr_path, sep="\t", dtype={"chromosome": str})
    required = {"chromosome", "start", "log2", "weight"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{cnr_path}: missing expected columns {missing}")

    df = df[np.isfinite(df["log2"])]
    if min_weight > 0:
        df = df[df["weight"] >= min_weight]

    diffs = []
    for _, chrom_df in df.groupby("chromosome", sort=False):
        chrom_df = chrom_df.sort_values("start")
        if len(chrom_df) < 2:
            continue
        diffs.append(np.abs(np.diff(chrom_df["log2"].to_numpy())))

    if not diffs:
        return {"mapd": float("nan"), "n_bins": len(df), "n_adjacent_pairs": 0}

    all_diffs = np.concatenate(diffs)
    return {
        "mapd": float(np.median(all_diffs)),
        "n_bins": int(len(df)),
        "n_adjacent_pairs": int(len(all_diffs)),
    }


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cnv-dir", required=True, type=Path,
                    help="Step 3 cnv/ directory containing <sample>/<sample>.dedup.cnr")
    p.add_argument("--manifest", required=True, type=Path,
                    help="passed_cells_only_manifest.tsv (sample_id column)")
    p.add_argument("--min-weight", type=float, default=0.0,
                    help="Exclude CNVkit bins below this weight (low-confidence bins)")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    rows = []
    for sample in samples:
        cnr_path = args.cnv_dir / sample / f"{sample}.dedup.cnr"
        if not cnr_path.is_file():
            rows.append({"sample_id": sample, "mapd": "", "n_bins": "",
                         "n_adjacent_pairs": "", "status": "MISSING_CNR"})
            continue
        try:
            result = compute_mapd(cnr_path, args.min_weight)
            rows.append({"sample_id": sample, "mapd": f"{result['mapd']:.6f}",
                         "n_bins": result["n_bins"],
                         "n_adjacent_pairs": result["n_adjacent_pairs"],
                         "status": "OK"})
        except Exception as exc:  # noqa: BLE001 - report per-sample, keep going
            rows.append({"sample_id": sample, "mapd": "", "n_bins": "",
                         "n_adjacent_pairs": "", "status": f"ERROR: {exc}"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "mapd", "n_bins",
                                                "n_adjacent_pairs", "status"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    ok = [r for r in rows if r["status"] == "OK"]
    print(f"[OK] MAPD computed for {len(ok)}/{len(rows)} samples -> {args.output}")
    if ok:
        values = [float(r["mapd"]) for r in ok]
        print(f"[INFO] Cohort MAPD: median={np.median(values):.4f} "
              f"min={min(values):.4f} max={max(values):.4f}")


if __name__ == "__main__":
    main()
