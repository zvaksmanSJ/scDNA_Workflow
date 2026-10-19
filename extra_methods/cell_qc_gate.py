#!/usr/bin/env python3
"""Per-cell pass/Borderline/FAIL QC gate, ported from BaseJumper's own
composite-score formula (`basej-dnaqc`'s `function_plot_qc_dna.R` +
`cnvSummarizer.R`, both vendored for reference under
`../resolveome_align/basej-public-pipelines/basej-dnaqc/container/
custom_r_qcplots/scripts/`), so this pipeline gates cells the same way
BaseJumper does rather than the original pipeline's simple hard thresholds
on read count/mapping rate/proper-pair rate.

Five binary pass/fail indicators (BaseJumper's own default cutoffs, each
independently overridable below):
  1. preseq_count   > 3,500,000,000  (library-complexity extrapolation)
  2. PCT_CHIMERAS   <         0.20   (Picard CollectAlignmentSummaryMetrics)
  3. chrM fraction  <         0.20   (mitochondrial-read fraction, idxstats)
  4. MAPD_CNV_Log2  <         0.25   (CNV bin-to-bin noise)
  5. SKEW_CNV       <         0.25   (CNV directional-bias noise)

CompositeScore = sum of the 5 indicators (0-5). qc_status: score>=4 -> PASS,
score==3 -> Borderline, score<=2 -> FAIL; any metric that could not be
computed counts as NOT met (conservative: uncertain is treated as failing
that one indicator, not silently dropped from the denominator) and is
called out by name in the explanation text. A cell with zero total reads
is always FAIL regardless of score, matching BaseJumper's own override.

Metrics 1-3 reuse files this pipeline's AlignTrimQC task already produces
(Picard alignment_summary_metrics, samtools idxstats, preseq-bioskryb
extrapolation). Metrics 4-5 need a CNVkit `.cnr`/`.call.cns` pair -- the
WDL task wrapping this script runs a lightweight `cnvkit.py batch` purely
for this gate (the full-depth CNV/SV stage later in the pipeline still
runs its own, separate CNVkit call; this mirrors BaseJumper's own
architecture, where `basej-dnaqc` runs its own dedicated Ginkgo CNV call
purely for QC scoring, distinct from whatever deeper CNV/SV pipeline runs
afterwards).

SKEW_CNV adaptation note: BaseJumper computes it from Ginkgo's own
mixture-model copy-number assignment per bin. This pipeline uses CNVkit,
not Ginkgo, so SKEW_CNV here is the same formula (direction-of-deviation
balance around each called integer copy-number level, weighted by how
common that level is) applied to CNVkit's `.cnr` log2 ratios and
`.call.cns` integer segment calls instead of Ginkgo's bins -- not a
byte-for-byte reproduction, but the same metric definition and the same
cutoff.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_CUTOFF_PRESEQ = 3_500_000_000
DEFAULT_CUTOFF_CHIMERAS = 0.20
DEFAULT_CUTOFF_CHRM = 0.20
DEFAULT_CUTOFF_MAPD = 0.25
DEFAULT_CUTOFF_SKEW = 0.25

MITO_CONTIG_NAMES = {"chrM", "MT", "chrMT", "M"}


class MetricResult:
    """One of the 5 indicators: raw value (or None if not computable),
    the cutoff it was checked against, whether it passed, and why."""

    def __init__(self, name: str, value, cutoff: float, comparison: str, note: str = ""):
        self.name = name
        self.value = value
        self.cutoff = cutoff
        self.comparison = comparison  # "gt" (value must be > cutoff) or "lt" (value must be < cutoff)
        self.note = note
        if value is None:
            self.passed = False
            self.not_computable = True
        else:
            self.not_computable = False
            self.passed = (value > cutoff) if comparison == "gt" else (value < cutoff)

    def describe(self) -> str:
        op = ">" if self.comparison == "gt" else "<"
        if self.not_computable:
            return (f"{self.name}: NOT COMPUTABLE ({self.note}) -- "
                     f"treated as FAILED for this metric (conservative default)")
        verdict = "PASS" if self.passed else "FAIL"
        return f"{self.name}: {self.value:.6g} (needs {op} {self.cutoff:g}) -- {verdict}"


def parse_total_reads_from_idxstats(idxstats_path: Path) -> int:
    total = 0
    with idxstats_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                total += int(parts[2]) + int(parts[3])
            except ValueError:
                continue
    return total


def parse_chrm_fraction(idxstats_path: Path) -> float | None:
    mapped_total = 0
    mapped_chrm = 0
    with idxstats_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            contig = parts[0]
            try:
                mapped = int(parts[2])
            except ValueError:
                continue
            mapped_total += mapped
            if contig in MITO_CONTIG_NAMES:
                mapped_chrm += mapped
    if mapped_total == 0:
        return None
    return mapped_chrm / mapped_total


def parse_pct_chimeras(alignment_summary_metrics_path: Path) -> float | None:
    """Picard METRICS-CLASS text file: '#'-comment lines, then a header
    line, then one data row per CATEGORY (FIRST_OF_PAIR/SECOND_OF_PAIR/
    PAIR). We want the PAIR row's PCT_CHIMERAS column."""
    lines = alignment_summary_metrics_path.read_text().splitlines()
    header = None
    for i, line in enumerate(lines):
        if line.startswith("CATEGORY"):
            header = line.split("\t")
            data_lines = lines[i + 1:]
            break
    else:
        return None
    if "PCT_CHIMERAS" not in header:
        return None
    col = header.index("PCT_CHIMERAS")
    cat_col = header.index("CATEGORY")
    for data_line in data_lines:
        if not data_line.strip():
            continue
        fields = data_line.split("\t")
        if len(fields) <= max(col, cat_col):
            continue
        if fields[cat_col] == "PAIR":
            try:
                return float(fields[col])
            except ValueError:
                return None
    return None


def parse_preseq_count(preseq_path: Path) -> float | None:
    """Either a one-line extrapolated-count file (AlignTrimQC's normal
    output) or a 'status\\treason' NOT_COMPUTABLE file (AlignTrimQC's
    graceful-degradation output when the preseq model fit fails)."""
    text = preseq_path.read_text().strip()
    if not text:
        return None
    first_line = text.splitlines()[0]
    if first_line.startswith("status"):
        return None
    try:
        return float(first_line.split()[0])
    except (ValueError, IndexError):
        return None


def compute_mapd(cnr_path: Path, min_weight: float = 0.0) -> float | None:
    """Median Absolute Pairwise Difference between adjacent bins' log2
    ratios, within each chromosome -- identical logic to mapd_qc.py's
    compute_mapd (duplicated here, not imported, so this script has no
    cross-file dependency when WDL localizes it on its own)."""
    df = pd.read_csv(cnr_path, sep="\t", dtype={"chromosome": str})
    required = {"chromosome", "start", "log2", "weight"}
    if required - set(df.columns):
        return None
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
        return None
    return float(np.median(np.concatenate(diffs)))


def compute_skew(cnr_path: Path, call_cns_path: Path) -> float | None:
    """BaseJumper's SKEW_CNV, adapted from Ginkgo bins to CNVkit bins/
    segments (see module docstring). For each bin, assign the integer CN
    of the segment it falls in, compute whether the bin's own linear
    copy-number estimate sits above or below that integer, then measure
    how balanced (vs. one-sided) that above/below split is at each CN
    level, weighted by how many bins are at that level. 1 = perfectly
    balanced noise; 0 = i.e. SKEW_CNV near 0 would mean no skew, so
    SKEW_CNV = 1 - weighted_mean(balance_ratio)."""
    cnr = pd.read_csv(cnr_path, sep="\t", dtype={"chromosome": str})
    cns = pd.read_csv(call_cns_path, sep="\t", dtype={"chromosome": str})
    if "cn" not in cns.columns or not {"chromosome", "start", "end", "log2"}.issubset(cnr.columns):
        return None

    cnr = cnr[np.isfinite(cnr["log2"])].copy()
    if cnr.empty:
        return None

    # Assign each bin the CN of the segment it falls within (same chromosome,
    # bin start within [segment start, segment end)).
    cn_values = []
    for chrom, bins_in_chrom in cnr.groupby("chromosome", sort=False):
        segs = cns[cns["chromosome"] == chrom]
        if segs.empty:
            cn_values.extend([None] * len(bins_in_chrom))
            continue
        seg_starts = segs["start"].to_numpy()
        seg_ends = segs["end"].to_numpy()
        seg_cns = segs["cn"].to_numpy()
        for start in bins_in_chrom["start"]:
            idx = np.where((seg_starts <= start) & (start < seg_ends))[0]
            cn_values.append(int(seg_cns[idx[0]]) if len(idx) else None)
    cnr["cn"] = cn_values
    cnr = cnr.dropna(subset=["cn"])
    if cnr.empty:
        return None

    # y = linear copy-number estimate at ploidy 2 from the log2 ratio.
    cnr["y"] = 2.0 * np.power(2.0, cnr["log2"])
    cnr["diff_expected_cn"] = cnr["y"] - cnr["cn"]
    cnr["direction"] = np.where(cnr["diff_expected_cn"] < 0, "Negative", "Positive")

    counts = cnr.groupby(["cn", "direction"]).size().reset_index(name="n")
    total_bins = len(cnr)

    ratios, weights = [], []
    for cn_level, group in counts.groupby("cn"):
        if len(group) < 2:
            continue  # this CN level only has bins going one direction -- no balance to measure
        n_sorted = sorted(group["n"].tolist())
        ratios.append(n_sorted[0] / n_sorted[1])
        weights.append(sum(n_sorted) / total_bins)

    if not ratios:
        return None
    return float(1.0 - np.average(ratios, weights=weights))


def score_to_status(score: int | None) -> str:
    if score is None or score <= 2:
        return "FAIL"
    if score == 3:
        return "Borderline"
    return "PASS"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--idxstats", required=True, type=Path)
    p.add_argument("--alignment-summary-metrics", required=True, type=Path)
    p.add_argument("--preseq-file", required=True, type=Path)
    p.add_argument("--cnr", required=True, type=Path)
    p.add_argument("--call-cns", required=True, type=Path)
    p.add_argument("--cutoff-preseq", type=float, default=DEFAULT_CUTOFF_PRESEQ)
    p.add_argument("--cutoff-chimeras", type=float, default=DEFAULT_CUTOFF_CHIMERAS)
    p.add_argument("--cutoff-chrm", type=float, default=DEFAULT_CUTOFF_CHRM)
    p.add_argument("--cutoff-mapd", type=float, default=DEFAULT_CUTOFF_MAPD)
    p.add_argument("--cutoff-skew", type=float, default=DEFAULT_CUTOFF_SKEW)
    p.add_argument("--output-tsv", required=True, type=Path)
    p.add_argument("--output-explanation", required=True, type=Path,
                    help="Plain-text per-cell explanation (always written; most useful for FAIL/Borderline cells).")
    p.add_argument("--output-pass-flag", required=True, type=Path,
                    help="Single-word PASS/Borderline/FAIL file, read back by main.wdl's qc_pass boolean.")
    return p.parse_args()


def main():
    args = parse_args()

    total_reads = parse_total_reads_from_idxstats(args.idxstats)

    metrics = [
        MetricResult("preseq_count", parse_preseq_count(args.preseq_file),
                     args.cutoff_preseq, "gt",
                     note="preseq model fit failed or file missing"),
        MetricResult("PCT_CHIMERAS", parse_pct_chimeras(args.alignment_summary_metrics),
                     args.cutoff_chimeras, "lt",
                     note="no PAIR row in Picard alignment_summary_metrics"),
        MetricResult("chrM_fraction", parse_chrm_fraction(args.idxstats),
                     args.cutoff_chrm, "lt",
                     note="no mapped reads in idxstats"),
        MetricResult("MAPD_CNV_Log2", compute_mapd(args.cnr),
                     args.cutoff_mapd, "lt",
                     note="CNVkit .cnr missing required columns or too few adjacent-bin pairs"),
        MetricResult("SKEW_CNV", compute_skew(args.cnr, args.call_cns),
                     args.cutoff_skew, "lt",
                     note="CNVkit .cnr/.call.cns missing required columns or no balanced CN level found"),
    ]

    composite_score = sum(1 for m in metrics if m.passed)
    qc_status = score_to_status(composite_score)
    zero_reads_override = total_reads == 0
    if zero_reads_override:
        qc_status = "FAIL"
        composite_score = 0

    qc_pass = qc_status != "FAIL"

    # ---- machine-readable TSV row (one per cell; cohort-wide tasks concatenate these) ----
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "sample_id": args.sample_id,
        "total_reads": total_reads,
        "composite_score": composite_score,
        "qc_status": qc_status,
        "qc_pass": str(qc_pass).lower(),
    }
    for m in metrics:
        row[f"{m.name}_value"] = "" if m.value is None else f"{m.value:.6g}"
        row[f"{m.name}_cutoff"] = m.cutoff
        row[f"{m.name}_passed"] = str(m.passed).lower()
    with args.output_tsv.open("w") as fh:
        fh.write("\t".join(row.keys()) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

    # ---- human-readable explanation (always written; essential for FAIL/Borderline cells) ----
    args.output_explanation.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"QC gate result for {args.sample_id}: {qc_status} (composite score {composite_score}/5)",
        "",
    ]
    if zero_reads_override:
        lines.append("OVERRIDE: total_reads == 0 -- forced to FAIL regardless of composite score.")
        lines.append("")
    lines.append(f"total_reads = {total_reads}")
    lines.append("")
    lines.append("Per-metric breakdown (needs >=4/5 to PASS, 3/5 is Borderline, <=2/5 is FAIL):")
    for m in metrics:
        lines.append(f"  - {m.describe()}")
    lines.append("")
    if qc_status == "FAIL":
        lines.append("Downstream variant calling / CNV-SV / lineage stages for this cell "
                      "will be SKIPPED by main.wdl's per-cell qc_pass gate.")
    elif qc_status == "Borderline":
        lines.append("Downstream stages WILL still run for this cell (Borderline is not FAIL), "
                      "but review this cell's metrics before trusting its results.")
    else:
        lines.append("Downstream stages will run normally for this cell.")
    args.output_explanation.write_text("\n".join(lines) + "\n")

    # ---- single-word flag file main.wdl reads to build its Boolean gate ----
    args.output_pass_flag.write_text(qc_status + "\n")

    print(f"[{qc_status}] {args.sample_id}: composite score {composite_score}/5 "
          f"(total_reads={total_reads})")


if __name__ == "__main__":
    main()
