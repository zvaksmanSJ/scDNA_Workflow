#!/usr/bin/env python3
"""Build a clean per-cell QC summary table + at-a-glance figures from the
existing MultiQC aggregate data (multiqc_general_stats.txt), which already
consolidates fastp, FastQC, Picard, and samtools metrics per sample. This
does not re-run or duplicate MultiQC -- it only re-presents its own output
table with readable column names and a few bar-chart figures (mirroring the
kind of per-sample bar plots MultiQC itself renders in its general stats
section) for a fast cohort overview; the full native MultiQC HTML report is
embedded separately in the dashboard for the complete/interactive view.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLUMN_MAP = {
    "Sample": "sample_id",
    "Picard: WgsMetrics_mqc-generalstats-picard_wgsmetrics-MEAN_COVERAGE": "mean_coverage",
    "Picard: WgsMetrics_mqc-generalstats-picard_wgsmetrics-MEDIAN_COVERAGE": "median_coverage",
    "Picard: WgsMetrics_mqc-generalstats-picard_wgsmetrics-PCT_30X": "pct_30x",
    "Samtools: flagstat_mqc-generalstats-samtools_flagstat-mapped_passed_pct": "pct_mapped",
    "fastp_mqc-generalstats-fastp-pct_duplication": "pct_duplication",
    "fastp_mqc-generalstats-fastp-after_filtering_q30_rate": "q30_rate",
    "fastp_mqc-generalstats-fastp-pct_surviving": "pct_surviving_trim",
    "FastQC_mqc-generalstats-fastqc-percent_gc": "pct_gc",
    "Picard: InsertSizeMetrics_mqc-generalstats-picard_insertsizemetrics-summed_mean": "mean_insert_size",
    "SnpEff_mqc-generalstats-snpeff-Ts_Tv_ratio": "snpeff_ts_tv_ratio",
    "SnpEff_mqc-generalstats-snpeff-Number_of_variants_before_filter": "snpeff_n_variants",
}

BAR_METRICS = [
    ("mean_coverage", "Mean genome coverage (Picard WgsMetrics)", "x"),
    ("pct_duplication", "Duplication rate (fastp)", "%"),
    ("q30_rate", "Q30 base rate after filtering (fastp)", "fraction"),
    ("pct_mapped", "Reads mapped (samtools flagstat)", "%"),
]


def short_name(sample_id: str) -> str:
    return sample_id.replace("SY5Y_", "").split("_DNA_")[0]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--multiqc-general-stats", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path,
                    help="Restrict to exactly these sample_id rows (MultiQC's general stats "
                         "table also has per-read-pair/pre-dedup/trim rows per cell that "
                         "would otherwise be double-counted).")
    p.add_argument("--output-tsv", required=True, type=Path)
    p.add_argument("--output-png", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        wanted_samples = {row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")}

    with args.multiqc_general_stats.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [r for r in reader if r.get("Sample") in wanted_samples]

    clean_rows = []
    for row in rows:
        clean = {}
        for src_col, dst_col in COLUMN_MAP.items():
            clean[dst_col] = row.get(src_col, "")
        clean_rows.append(clean)
    clean_rows.sort(key=lambda r: r["sample_id"])

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(COLUMN_MAP.values())
    with args.output_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(clean_rows)

    samples = [short_name(r["sample_id"]) for r in clean_rows]
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for ax, (metric, title, unit) in zip(axes.flat, BAR_METRICS):
        values = []
        for r in clean_rows:
            try:
                values.append(float(r[metric]))
            except (TypeError, ValueError):
                values.append(0.0)
        ax.bar(samples, values, color="#2c7fb8")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(unit, fontsize=9)
        ax.tick_params(axis="x", rotation=90, labelsize=7)
    fig.suptitle("Per-cell QC overview (from existing MultiQC aggregate data)", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=140)
    plt.close(fig)

    print(f"[OK] QC summary table -> {args.output_tsv} ({len(clean_rows)} samples)")
    print(f"[OK] QC overview figure -> {args.output_png}")


if __name__ == "__main__":
    main()
