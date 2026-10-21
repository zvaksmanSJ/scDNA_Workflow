#!/usr/bin/env python3
"""Per-cell CNV burden summary, computed directly from CNVkit's own
`.dedup.call.cns` segment-level integer copy-number calls (no new upstream
tool required -- CNVkit already produces this file during the existing
variants_cnv_sv stage).

For each cell: counts CNVkit's called segments, and reports what fraction
of the analyzed genome is gained (cn > 2), lost (cn < 2), or altered overall
(cn != 2) -- a genome-wide CNV *burden* metric, complementary to (not a
replacement for) the per-cell CNVkit scatter/diagram figures themselves,
which the dashboard embeds separately alongside this table.

cn=0 segments are counted under "lost" (homozygous deletion is still a
loss), matching common convention; a cell with few/no segments departing
from cn=2 is reported as largely diploid.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def summarize_call_cns(path: Path) -> dict:
    total_bp = 0
    gain_bp = 0
    loss_bp = 0
    n_segments = 0
    max_cn = None
    min_cn = None
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            start, end = int(row["start"]), int(row["end"])
            length = end - start
            if length <= 0:
                continue
            cn = int(row["cn"])
            total_bp += length
            n_segments += 1
            if cn > 2:
                gain_bp += length
            elif cn < 2:
                loss_bp += length
            max_cn = cn if max_cn is None else max(max_cn, cn)
            min_cn = cn if min_cn is None else min(min_cn, cn)
    if total_bp == 0:
        return {"status": "EMPTY"}
    return {
        "n_segments": n_segments,
        "genome_bp_analyzed": total_bp,
        "pct_genome_gain": round(100 * gain_bp / total_bp, 3),
        "pct_genome_loss": round(100 * loss_bp / total_bp, 3),
        "pct_genome_altered": round(100 * (gain_bp + loss_bp) / total_bp, 3),
        "max_cn": max_cn,
        "min_cn": min_cn,
        "status": "OK",
    }


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cnv-dir", required=True, type=Path,
                    help="Step3's CNV output root, containing <sample>/<sample>.dedup.call.cns "
                         "and <sample>/<sample>.dedup-scatter.png")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    fieldnames = ["sample_id", "n_segments", "genome_bp_analyzed", "pct_genome_gain",
                  "pct_genome_loss", "pct_genome_altered", "max_cn", "min_cn",
                  "scatter_png", "status"]
    rows = []
    for sample in samples:
        call_cns = args.cnv_dir / sample / f"{sample}.dedup.call.cns"
        scatter_png = args.cnv_dir / sample / f"{sample}.dedup-scatter.png"
        if not call_cns.is_file():
            rows.append({"sample_id": sample, "status": "MISSING"})
            continue
        try:
            summary = summarize_call_cns(call_cns)
        except Exception as exc:  # noqa: BLE001 - report per-sample, keep going
            rows.append({"sample_id": sample, "status": f"ERROR: {exc}"})
            continue
        summary["sample_id"] = sample
        summary["scatter_png"] = str(scatter_png) if scatter_png.is_file() else ""
        rows.append(summary)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    ok = [r for r in rows if r.get("status") == "OK"]
    print(f"[OK] CNV burden summarized for {len(ok)}/{len(rows)} samples -> {args.output}")
    if ok:
        mean_altered = sum(r["pct_genome_altered"] for r in ok) / len(ok)
        print(f"[INFO] Cohort mean %% genome altered (non-diploid): {mean_altered:.2f}")


if __name__ == "__main__":
    main()
