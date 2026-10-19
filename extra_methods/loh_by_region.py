#!/usr/bin/env python3
"""Genome-wide LOH recurrence across the cohort.

Bins the genome into fixed-size windows and, for each cell, marks a window
as "LOH" if any Control-FREEC segment overlapping that window was assigned
a single-repeated-parental-allele genotype (the same native FREEC genotype
field used by loh_summary.py -- no new/custom LOH classification logic).
Produces:
  - loh_by_region.tsv: one row per genomic bin with the count/fraction of
    cells showing LOH there, plus which cells.
  - loh_by_region_heatmap.png: cells x genome-position heatmap (rows=cells,
    columns=binned genomic position, colored by LOH presence), the
    established way to visualize cohort-wide CNA/LOH recurrence (GISTIC-
    style recurrence tracks).
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOH_GENOTYPE_RE = re.compile(r"^(A+|B+)$")
CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X"]


def is_loh_genotype(genotype: str) -> bool:
    return bool(LOH_GENOTYPE_RE.match(genotype))


def read_chrom_lengths(fai_path: Path) -> dict:
    lengths = {}
    with fai_path.open() as fh:
        for line in fh:
            parts = line.split("\t")
            if parts[0] in CHROM_ORDER:
                lengths[parts[0]] = int(parts[1])
    return lengths


def read_cnvs_segments(cnvs_path: Path) -> list[tuple]:
    segs = []
    with cnvs_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[0] not in CHROM_ORDER:
                continue
            chrom, start, end, genotype = parts[0], int(parts[1]), int(parts[2]), parts[5]
            segs.append((chrom, start, end, is_loh_genotype(genotype)))
    return segs


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freec-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--chrom-len-file", required=True, type=Path,
                    help="samtools .fai-format file restricted to chr1-22+X")
    p.add_argument("--bin-size", type=int, default=5_000_000)
    p.add_argument("--output-tsv", required=True, type=Path)
    p.add_argument("--output-png", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    chrom_lengths = read_chrom_lengths(args.chrom_len_file)

    # Build a flat, ordered list of (chrom, bin_start, bin_end) genomic bins.
    bins = []
    for chrom in CHROM_ORDER:
        length = chrom_lengths.get(chrom)
        if not length:
            continue
        start = 0
        while start < length:
            end = min(start + args.bin_size, length)
            bins.append((chrom, start, end))
            start = end
    n_bins = len(bins)
    bin_index = {(c, s, e): i for i, (c, s, e) in enumerate(bins)}

    loh_matrix = np.zeros((len(samples), n_bins), dtype=bool)
    missing_samples = []
    for row, sample in enumerate(samples):
        cnvs_path = args.freec_dir / sample / f"{sample}.dedup.bam_CNVs"
        if not cnvs_path.is_file():
            missing_samples.append(sample)
            continue
        segs = read_cnvs_segments(cnvs_path)
        for chrom, seg_start, seg_end, loh in segs:
            if not loh:
                continue
            bin_size = args.bin_size
            first_bin = seg_start // bin_size
            last_bin = (seg_end - 1) // bin_size
            for b in range(first_bin, last_bin + 1):
                bstart = b * bin_size
                bend = min(bstart + bin_size, chrom_lengths.get(chrom, bstart + bin_size))
                key = (chrom, bstart, bend)
                if key in bin_index:
                    loh_matrix[row, bin_index[key]] = True

    # --- table ---
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chrom", "start", "end", "n_cells_loh", "n_cells_total",
                          "fraction_cells_loh", "cells_with_loh"])
        for j, (chrom, start, end) in enumerate(bins):
            col = loh_matrix[:, j]
            cells_with = [samples[i] for i in range(len(samples)) if col[i]]
            n = int(col.sum())
            writer.writerow([chrom, start, end, n, len(samples),
                              f"{n/len(samples):.4f}" if samples else "0",
                              ";".join(cells_with)])

    # --- figure: cells x genome-position recurrence heatmap ---
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig_h = max(4, 0.35 * len(samples) + 1.5)
    fig, ax = plt.subplots(figsize=(16, fig_h))
    ax.imshow(loh_matrix, aspect="auto", cmap="Reds", interpolation="nearest", vmin=0, vmax=1)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels([s.replace("SY5Y_", "").split("_DNA_")[0] for s in samples], fontsize=8)

    chrom_boundaries = []
    pos = 0
    for chrom in CHROM_ORDER:
        length = chrom_lengths.get(chrom)
        if not length:
            continue
        n_chrom_bins = sum(1 for c, _, _ in bins if c == chrom)
        chrom_boundaries.append((pos, pos + n_chrom_bins, chrom))
        pos += n_chrom_bins
    tick_pos = [(a + b) / 2 for a, b, _ in chrom_boundaries]
    tick_lab = [c for _, _, c in chrom_boundaries]
    for a, _, _ in chrom_boundaries[1:]:
        ax.axvline(a - 0.5, color="black", linewidth=0.4)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, fontsize=8)
    ax.set_xlabel(f"Genomic position (chromosome; {args.bin_size // 1_000_000}Mb bins)")
    ax.set_title("Cohort-wide LOH recurrence (Control-FREEC native genotype field)")
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=140)
    plt.close(fig)

    print(f"[OK] LOH-by-region table -> {args.output_tsv} ({n_bins} bins x {len(samples)} cells)")
    print(f"[OK] LOH-by-region heatmap -> {args.output_png}")
    if missing_samples:
        print(f"[WARN] Missing FREEC output for: {', '.join(missing_samples)}")


if __name__ == "__main__":
    main()
