#!/usr/bin/env python3
"""Combine CNVkit segments across single cells into robust cohort CNV reports."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm


PRIMARY_CHROMS = [str(value) for value in range(1, 23)] + ["X", "Y"]
CELL_NUMBER = re.compile(r"_Cell(\d+)_")


@dataclass(frozen=True)
class Segment:
    start: int
    end: int
    log2: float


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.start + 1}-{self.end}"


def sample_sort_key(sample: str) -> tuple[int, str]:
    match = CELL_NUMBER.search(sample)
    return (int(match.group(1)) if match else 10**9, sample)


def contig_sort_key(chrom: str) -> tuple[int, int | str]:
    if chrom.isdigit():
        return (0, int(chrom))
    return (1, chrom)


def weighted_median(values: list[tuple[float, float]]) -> float:
    filtered = sorted(
        ((value, weight) for value, weight in values if weight > 0 and math.isfinite(value)),
        key=lambda item: item[0],
    )
    if not filtered:
        return math.nan
    midpoint = sum(weight for _, weight in filtered) / 2
    cumulative = 0.0
    for value, weight in filtered:
        cumulative += weight
        if cumulative >= midpoint:
            return value
    return filtered[-1][0]


def median_absolute_deviation(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return math.nan
    median = float(np.median(finite))
    return float(np.median(np.abs(finite - median)))


def finite_median(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else math.nan


def load_chrom_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0].removeprefix("chr")
            if chrom in PRIMARY_CHROMS:
                sizes[chrom] = int(fields[1])
    missing = [chrom for chrom in PRIMARY_CHROMS if chrom not in sizes]
    if missing:
        raise ValueError(f"Reference index is missing primary chromosomes: {missing}")
    return sizes


def load_segments(path: Path) -> dict[str, list[Segment]]:
    segments: dict[str, list[Segment]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "start", "end", "log2"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path} is missing columns: {sorted(required)}")
        for row in reader:
            chrom = row["chromosome"].removeprefix("chr")
            if chrom not in PRIMARY_CHROMS:
                continue
            start = int(row["start"])
            end = int(row["end"])
            value = float(row["log2"])
            if end > start and math.isfinite(value):
                segments[chrom].append(Segment(start, end, value))
    for chrom in segments:
        segments[chrom].sort(key=lambda segment: (segment.start, segment.end))
    return segments


def discover_samples(cnv_dir: Path) -> dict[str, Path]:
    samples = {}
    for path in cnv_dir.glob("*/*.dedup.cns"):
        sample = path.parent.name
        expected = f"{sample}.dedup.cns"
        if path.name == expected:
            samples[sample] = path
    if not samples:
        raise ValueError(f"No per-cell *.dedup.cns files found under {cnv_dir}")
    return dict(sorted(samples.items(), key=lambda item: sample_sort_key(item[0])))


def make_regions(chrom_sizes: dict[str, int], bin_size: int) -> list[Region]:
    regions = []
    for chrom in PRIMARY_CHROMS:
        length = chrom_sizes[chrom]
        for start in range(0, length, bin_size):
            regions.append(Region(chrom, start, min(start + bin_size, length)))
    return regions


def bin_segments(
    regions: list[Region],
    segments_by_chrom: dict[str, list[Segment]],
    min_covered_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.full(len(regions), np.nan, dtype=float)
    covered_fractions = np.zeros(len(regions), dtype=float)
    region_indices: dict[str, list[int]] = defaultdict(list)
    for index, region in enumerate(regions):
        region_indices[region.chrom].append(index)
    for chrom, indices in region_indices.items():
        segments = segments_by_chrom.get(chrom, [])
        segment_index = 0
        for region_index in indices:
            region = regions[region_index]
            while (
                segment_index < len(segments)
                and segments[segment_index].end <= region.start
            ):
                segment_index += 1
            cursor = segment_index
            weighted_sum = 0.0
            covered = 0
            while cursor < len(segments) and segments[cursor].start < region.end:
                segment = segments[cursor]
                overlap = max(
                    0,
                    min(region.end, segment.end) - max(region.start, segment.start),
                )
                if overlap:
                    weighted_sum += segment.log2 * overlap
                    covered += overlap
                cursor += 1
            region_length = region.end - region.start
            covered_fraction = covered / region_length
            covered_fractions[region_index] = covered_fraction
            if covered and covered_fraction >= min_covered_fraction:
                values[region_index] = weighted_sum / covered
    return values, covered_fractions


def write_matrix(
    path: Path,
    regions: list[Region],
    samples: list[str],
    matrix: np.ndarray,
    formatter,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chromosome", "start", "end", "region", *samples])
        for index, region in enumerate(regions):
            writer.writerow(
                [
                    region.chrom,
                    region.start,
                    region.end,
                    region.label,
                    *[formatter(value) for value in matrix[index]],
                ]
            )


def state_label(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return {-1: "loss", 0: "neutral", 1: "gain"}[int(value)]


def chromosome_boundaries(regions: list[Region]) -> tuple[list[int], list[float], list[str]]:
    boundaries = [0]
    centers = []
    labels = []
    start = 0
    for chrom in PRIMARY_CHROMS:
        count = sum(region.chrom == chrom for region in regions)
        if count:
            centers.append(start + count / 2)
            labels.append(chrom)
            start += count
            boundaries.append(start)
    return boundaries, centers, labels


def render_heatmap(
    matrix: np.ndarray,
    samples: list[str],
    regions: list[Region],
    title: str,
    bin_size: int,
    state: bool = False,
):
    width = min(40, max(18, len(regions) / 150))
    height = max(7, len(samples) * 0.38 + 2.5)
    figure, axis = plt.subplots(figsize=(width, height), constrained_layout=True)
    display = matrix.T
    if state:
        cmap = ListedColormap(["#2166ac", "#f7f7f7", "#b2182b"])
        norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
        image = axis.imshow(display, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        colorbar = figure.colorbar(image, ax=axis, ticks=[-1, 0, 1], fraction=0.015)
        colorbar.ax.set_yticklabels(["Loss", "Neutral", "Gain"])
    else:
        finite = np.abs(display[np.isfinite(display)])
        limit = max(0.6, float(np.percentile(finite, 98)) if finite.size else 1.0)
        limit = min(limit, 3.0)
        image = axis.imshow(
            display,
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        )
        colorbar = figure.colorbar(image, ax=axis, fraction=0.015)
        colorbar.set_label("Adjusted log2 ratio")
    boundaries, centers, labels = chromosome_boundaries(regions)
    axis.set_xticks(centers, labels, rotation=0)
    for boundary in boundaries:
        axis.axvline(boundary - 0.5, color="#555555", linewidth=0.35)
    axis.set_yticks(range(len(samples)), samples, fontsize=7)
    axis.set_xlabel(f"Chromosome ({bin_size:,} bp bins)")
    axis.set_title(title)
    return figure


def render_prevalence(
    summary_rows: list[dict[str, object]],
    chrom: str | None = None,
):
    selected = [
        row for row in summary_rows if chrom is None or row["chromosome"] == chrom
    ]
    if chrom is None:
        width = 20
        title = "Cohort CNV prevalence across all cells"
    else:
        width = 14
        title = f"Chromosome {chrom}: gain/loss prevalence across all cells"
    figure, axis = plt.subplots(figsize=(width, 4.8), constrained_layout=True)
    x = np.arange(len(selected))
    gains = np.array([float(row["gain_fraction"]) for row in selected])
    losses = -np.array([float(row["loss_fraction"]) for row in selected])
    axis.fill_between(x, gains, color="#b2182b", alpha=0.9, label="Gain")
    axis.fill_between(x, losses, color="#2166ac", alpha=0.9, label="Loss")
    axis.axhline(0, color="black", linewidth=0.6)
    axis.set_ylim(-1, 1)
    axis.set_ylabel("Fraction of cells\n(loss plotted negative)")
    axis.set_title(title)
    if chrom is None:
        regions = [
            Region(str(row["chromosome"]), int(row["start"]), int(row["end"]))
            for row in selected
        ]
        boundaries, centers, labels = chromosome_boundaries(regions)
        axis.set_xticks(centers, labels)
        for boundary in boundaries:
            axis.axvline(boundary - 0.5, color="#777777", linewidth=0.35)
    else:
        step = max(1, len(selected) // 10)
        ticks = np.arange(0, len(selected), step)
        axis.set_xticks(
            ticks,
            [f"{int(selected[index]['start']) / 1e6:.0f}" for index in ticks],
        )
        axis.set_xlabel("Position (Mb)")
    axis.legend(loc="upper right", ncols=2)
    return figure


def figure_data_uri(figure) -> str:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=130, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def write_html(
    path: Path,
    adjusted: np.ndarray,
    states: np.ndarray,
    samples: list[str],
    regions: list[Region],
    summary_rows: list[dict[str, object]],
    qc_rows: list[dict[str, object]],
    bin_size: int,
) -> None:
    adjusted_figure = render_heatmap(
        adjusted,
        samples,
        regions,
        "Single-cell-adjusted CNV log2 ratios",
        bin_size,
    )
    state_figure = render_heatmap(
        states,
        samples,
        regions,
        "Single-cell CNV gain/neutral/loss states",
        bin_size,
        state=True,
    )
    prevalence_figure = render_prevalence(summary_rows)
    images = [
        figure_data_uri(adjusted_figure),
        figure_data_uri(state_figure),
        figure_data_uri(prevalence_figure),
    ]
    plt.close(adjusted_figure)
    plt.close(state_figure)
    plt.close(prevalence_figure)
    recurrent = sorted(
        summary_rows,
        key=lambda row: max(float(row["gain_fraction"]), float(row["loss_fraction"])),
        reverse=True,
    )[:100]
    qc_table = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in (
                "sample",
                "autosomal_baseline_log2",
                "robust_sd",
                "gain_threshold",
                "loss_threshold",
                "callable_bins",
            )
        )
        + "</tr>"
        for row in qc_rows
    )
    recurrent_table = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in (
                "region",
                "gain_cells",
                "loss_cells",
                "gain_fraction",
                "loss_fraction",
                "combined_state",
            )
        )
        + "</tr>"
        for row in recurrent
    )
    path.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Combined single-cell CNV report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 24px; color: #222; }}
img {{ max-width: 100%; border: 1px solid #bbb; margin-bottom: 24px; }}
table {{ border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid #bbb; padding: 4px 7px; text-align: right; }}
th:first-child, td:first-child {{ text-align: left; }}
.note {{ background: #f2f5f8; padding: 12px; border-left: 5px solid #4c78a8; }}
</style></head><body>
<h1>Combined single-cell CNVkit report</h1>
<div class="note">
<p><strong>{len(samples)} cells; {bin_size:,} bp bins.</strong></p>
<p>Raw values are length-weighted CNVkit segmented log2 ratios. Adjusted values
subtract each cell's length-weighted autosomal median. Calls use cell-specific
thresholds of max(0.30, 2.5 × robust SD). This controls single-cell WGA noise
while preserving changes shared across cells.</p>
<p>Calls are relative to each cell's autosomal baseline; whole-genome ploidy
changes cannot be inferred without an external ploidy reference.</p>
</div>
<h2>Adjusted CNV heatmap</h2><img src="{images[0]}">
<h2>Discrete state heatmap</h2><img src="{images[1]}">
<h2>Gain/loss prevalence</h2><img src="{images[2]}">
<h2>Per-cell adjustment and noise thresholds</h2>
<table><thead><tr><th>Sample</th><th>Baseline log2</th><th>Robust SD</th>
<th>Gain threshold</th><th>Loss threshold</th><th>Callable bins</th></tr></thead>
<tbody>{qc_table}</tbody></table>
<h2>Top 100 recurrent regions</h2>
<table><thead><tr><th>Region</th><th>Gain cells</th><th>Loss cells</th>
<th>Gain fraction</th><th>Loss fraction</th><th>Combined state</th></tr></thead>
<tbody>{recurrent_table}</tbody></table>
</body></html>
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnv-dir", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", default="combined_single_cell_cnv")
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--min-covered-fraction", type=float, default=0.5)
    parser.add_argument("--minimum-log2-threshold", type=float, default=0.30)
    parser.add_argument("--noise-multiplier", type=float, default=2.5)
    parser.add_argument("--recurrent-min-fraction", type=float, default=0.20)
    args = parser.parse_args()
    if args.bin_size <= 0:
        parser.error("--bin-size must be positive")
    if not 0 <= args.min_covered_fraction <= 1:
        parser.error("--min-covered-fraction must be between 0 and 1")
    if not 0 <= args.recurrent_min_fraction <= 1:
        parser.error("--recurrent-min-fraction must be between 0 and 1")

    cnv_dir = Path(args.cnv_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / args.output_prefix
    chrom_sizes = load_chrom_sizes(Path(args.reference_fai))
    sample_paths = discover_samples(cnv_dir)
    samples = list(sample_paths)
    regions = make_regions(chrom_sizes, args.bin_size)
    raw = np.full((len(regions), len(samples)), np.nan, dtype=float)
    coverage = np.zeros_like(raw)
    segments_by_sample = {}
    for sample_index, (sample, path) in enumerate(sample_paths.items()):
        segments = load_segments(path)
        segments_by_sample[sample] = segments
        raw[:, sample_index], coverage[:, sample_index] = bin_segments(
            regions, segments, args.min_covered_fraction
        )

    autosomal_region_mask = np.array(
        [region.chrom.isdigit() for region in regions], dtype=bool
    )
    adjusted = np.full_like(raw, np.nan)
    states = np.full_like(raw, np.nan)
    qc_rows = []
    for sample_index, sample in enumerate(samples):
        autosomal_segments = [
            (segment.log2, segment.end - segment.start)
            for chrom, segments in segments_by_sample[sample].items()
            if chrom.isdigit()
            for segment in segments
        ]
        baseline = weighted_median(autosomal_segments)
        adjusted[:, sample_index] = raw[:, sample_index] - baseline
        autosomal_values = adjusted[autosomal_region_mask, sample_index]
        robust_sd = 1.4826 * median_absolute_deviation(autosomal_values)
        threshold = max(
            args.minimum_log2_threshold,
            args.noise_multiplier * robust_sd if math.isfinite(robust_sd) else 0,
        )
        finite = np.isfinite(adjusted[:, sample_index])
        states[finite, sample_index] = 0
        states[finite & (adjusted[:, sample_index] >= threshold), sample_index] = 1
        states[finite & (adjusted[:, sample_index] <= -threshold), sample_index] = -1
        qc_rows.append(
            {
                "sample": sample,
                "autosomal_baseline_log2": f"{baseline:.6f}",
                "robust_sd": f"{robust_sd:.6f}",
                "gain_threshold": f"{threshold:.6f}",
                "loss_threshold": f"{-threshold:.6f}",
                "callable_bins": int(finite.sum()),
                "median_covered_fraction": f"{np.median(coverage[finite, sample_index]):.6f}",
            }
        )

    summary_rows = []
    for region_index, region in enumerate(regions):
        row_states = states[region_index]
        callable_cells = int(np.isfinite(row_states).sum())
        gain_cells = int(np.sum(row_states == 1))
        loss_cells = int(np.sum(row_states == -1))
        neutral_cells = int(np.sum(row_states == 0))
        gain_fraction = gain_cells / callable_cells if callable_cells else 0
        loss_fraction = loss_cells / callable_cells if callable_cells else 0
        if (
            gain_fraction >= args.recurrent_min_fraction
            and loss_fraction >= args.recurrent_min_fraction
        ):
            combined_state = "mixed"
        elif gain_fraction >= args.recurrent_min_fraction:
            combined_state = "recurrent_gain"
        elif loss_fraction >= args.recurrent_min_fraction:
            combined_state = "recurrent_loss"
        else:
            combined_state = "no_recurrent_change"
        summary_rows.append(
            {
                "chromosome": region.chrom,
                "start": region.start,
                "end": region.end,
                "region": region.label,
                "callable_cells": callable_cells,
                "gain_cells": gain_cells,
                "neutral_cells": neutral_cells,
                "loss_cells": loss_cells,
                "gain_fraction": round(gain_fraction, 6),
                "loss_fraction": round(loss_fraction, 6),
                "combined_state": combined_state,
                "median_raw_log2": round(finite_median(raw[region_index]), 6),
                "median_adjusted_log2": round(
                    finite_median(adjusted[region_index]), 6
                ),
            }
        )

    write_matrix(
        Path(f"{prefix}.raw_log2_matrix.tsv"),
        regions,
        samples,
        raw,
        lambda value: "NA" if not math.isfinite(value) else f"{value:.6f}",
    )
    write_matrix(
        Path(f"{prefix}.adjusted_log2_matrix.tsv"),
        regions,
        samples,
        adjusted,
        lambda value: "NA" if not math.isfinite(value) else f"{value:.6f}",
    )
    write_matrix(
        Path(f"{prefix}.state_matrix.tsv"),
        regions,
        samples,
        states,
        state_label,
    )
    with Path(f"{prefix}.regional_summary.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    with Path(f"{prefix}.cell_qc.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(qc_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    pdf_path = Path(f"{prefix}.pdf")
    with PdfPages(pdf_path) as pdf:
        for figure in (
            render_heatmap(
                adjusted,
                samples,
                regions,
                "Single-cell-adjusted CNV log2 ratios",
                args.bin_size,
            ),
            render_heatmap(
                states,
                samples,
                regions,
                "Single-cell CNV gain/neutral/loss states",
                args.bin_size,
                state=True,
            ),
            render_prevalence(summary_rows),
        ):
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)
        for chrom in PRIMARY_CHROMS:
            figure = render_prevalence(summary_rows, chrom)
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)

    write_html(
        Path(f"{prefix}.html"),
        adjusted,
        states,
        samples,
        regions,
        summary_rows,
        qc_rows,
        args.bin_size,
    )
    print(f"[OK] Combined {len(samples)} cells across {len(regions)} regions")
    print(f"[OK] HTML: {prefix}.html")
    print(f"[OK] PDF: {prefix}.pdf")


if __name__ == "__main__":
    main()
