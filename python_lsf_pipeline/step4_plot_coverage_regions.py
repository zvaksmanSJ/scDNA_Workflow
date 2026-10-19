#!/usr/bin/env python3
"""Plot chromosome-wide >3X coverage from Step 4 annotated-region reports."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import io
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


def chromosome_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (int(value) if value.isdigit() else order.get(value, 100), value)


def read_reports(input_dir: Path):
    positions: dict[tuple[str, str], list[tuple[str, int, int]]] = defaultdict(list)
    with (input_dir / "annotation_features.tsv").open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 7 or row[1] != "gene_exon_union":
                continue
            positions[(row[5], row[6])].append((row[2], int(row[3]), int(row[4])))

    genes: dict[str, dict[tuple[str, str], dict[str, float | str]]] = {}
    for gene_path in sorted(input_dir.glob("*/*.gene_exon_coverage.tsv.gz")):
        gene_stats = {}
        with gzip.open(gene_path, "rt", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row["sample"]
                gene_stats[(row["gene_id"], row["gene_name"])] = {
                    "exon_bases": float(row["exon_bases"]),
                    "bases_gt3x": float(row["bases_gt3x"]),
                    "percent_gt3x": float(row["percent_exon_bases_gt3x"]),
                    "mean_gt3x": float(row["mean_exon_depth_gt3x"]),
                    "mean_all": float(row["mean_exon_depth_all"]),
                    "zero_bases": float(row["zero_exon_bases"]),
                }
        if not gene_stats:
            continue
        genes[sample] = gene_stats
    if not genes:
        raise SystemExit(f"No Step 4 gene coverage reports found under {input_dir}")
    return positions, genes


def build_chromosome_matrix(
    positions: dict[tuple[str, str], list[tuple[str, int, int]]],
    genes: dict[str, dict[tuple[str, str], dict[str, float | str]]],
    chrom: str,
    bin_size: int,
):
    selected_positions = {
        key: [(start, end) for row_chrom, start, end in intervals if row_chrom == chrom]
        for key, intervals in positions.items()
    }
    selected_positions = {
        key: intervals for key, intervals in selected_positions.items() if intervals
    }
    if not selected_positions:
        raise ValueError(f"No annotated gene-exon intervals found for chromosome {chrom}")
    chrom_end = max(end for intervals in selected_positions.values() for _, end in intervals)
    bin_count = (chrom_end + bin_size - 1) // bin_size
    matrix = np.full((len(genes), bin_count), np.nan)
    sample_names = sorted(genes)
    gene_positions: dict[tuple[str, str], tuple[int, int]] = {}

    for sample_index, sample in enumerate(sample_names):
        weighted = np.zeros(bin_count)
        observed = np.zeros(bin_count)
        for key, intervals in selected_positions.items():
            stats = genes[sample].get(key)
            if not stats:
                continue
            fraction = float(stats["percent_gt3x"]) / 100.0
            for start, end in intervals:
                if key not in gene_positions:
                    gene_positions[key] = (start, end)
                else:
                    old_start, old_end = gene_positions[key]
                    gene_positions[key] = (min(old_start, start), max(old_end, end))
                first, last = (start - 1) // bin_size, (end - 1) // bin_size
                for bin_index in range(first, last + 1):
                    bin_start = bin_index * bin_size + 1
                    bin_end = (bin_index + 1) * bin_size
                    overlap = max(0, min(end, bin_end) - max(start, bin_start) + 1)
                    weighted[bin_index] += fraction * overlap
                    observed[bin_index] += overlap
        valid = observed > 0
        matrix[sample_index, valid] = weighted[valid] / observed[valid]
    return sample_names, matrix, chrom_end, gene_positions


def low_coverage_genes(
    sample_names,
    genes,
    gene_positions,
    max_labels: int,
):
    rows = []
    for key, (start, end) in gene_positions.items():
        values = [genes[sample].get(key) for sample in sample_names]
        values = [value for value in values if value]
        if not values:
            continue
        mean_percent = sum(float(value["percent_gt3x"]) for value in values) / len(values)
        zero_samples = sum(
            float(value["zero_bases"]) >= float(value["exon_bases"])
            for value in values
        )
        no_gt3x_samples = sum(float(value["bases_gt3x"]) == 0 for value in values)
        rows.append({
            "gene_id": key[0],
            "gene_name": key[1],
            "start": start,
            "end": end,
            "mean_percent": mean_percent,
            "zero_samples": zero_samples,
            "no_gt3x_samples": no_gt3x_samples,
        })
    rows.sort(key=lambda row: (
        -row["zero_samples"],
        -row["no_gt3x_samples"],
        row["mean_percent"],
        row["start"],
    ))
    return rows[:max_labels], rows


def make_figure(
    chrom,
    sample_names,
    matrix,
    chrom_end,
    gene_rows,
    bin_size,
):
    width = max(18.0, chrom_end / 4_000_000)
    height = max(10.0, 0.48 * len(sample_names) + 5.0)
    fig, (ax, gene_ax) = plt.subplots(
        2, 1, figsize=(width, height),
        gridspec_kw={"height_ratios": [max(5, len(sample_names)), 2.4]},
        constrained_layout=True,
    )
    cmap = ListedColormap(["#b2182b", "#ef8a62", "#fee08b", "#1a9850"])
    cmap.set_bad("#d9d9d9")
    norm = BoundaryNorm([0, 0.001, 0.5, 0.9, 1.001], cmap.N)
    extent = [0, chrom_end / 1e6, len(sample_names), 0]
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm,
              extent=extent)
    ax.set_yticks(np.arange(len(sample_names)) + 0.5)
    ax.set_yticklabels(sample_names, fontsize=7)
    ax.set_xlabel(f"Chromosome {chrom} position (Mb); {bin_size // 1000:,} kb bins")
    ax.set_title(
        f"Chromosome {chrom}: fraction of annotated gene-exon bases with depth >3X\n"
        "Gray = no annotated exon bases in the bin"
    )
    ax.grid(axis="x", color="white", alpha=0.3, linewidth=0.5)
    ax.legend(
        handles=[
            Patch(color="#b2182b", label="0% of bases >3X"),
            Patch(color="#ef8a62", label="<50% of bases >3X"),
            Patch(color="#fee08b", label="50–<90% of bases >3X"),
            Patch(color="#1a9850", label="≥90% of bases >3X"),
            Patch(color="#d9d9d9", label="No annotated exon interval"),
        ],
        loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=5, fontsize=8,
    )

    gene_ax.set_xlim(0, chrom_end / 1e6)
    gene_ax.set_ylim(0, 1)
    gene_ax.set_yticks([])
    gene_ax.set_xlabel(f"Chromosome {chrom} position (Mb)")
    gene_ax.set_title("Genes prioritized in regions with poor or absent >3X coverage")
    for index, row in enumerate(sorted(gene_rows, key=lambda value: value["start"])):
        x = (row["start"] + row["end"]) / 2 / 1e6
        y = 0.18 + 0.2 * (index % 4)
        color = "#8b0000" if row["zero_samples"] else "#d95f02"
        gene_ax.vlines(x, 0.02, y - 0.03, color=color, linewidth=0.7)
        label = (
            f"{row['gene_name']} "
            f"({row['zero_samples']}/{len(sample_names)} at 0X; "
            f"{row['no_gt3x_samples']}/{len(sample_names)} no >3X)"
        )
        gene_ax.text(x, y, label, rotation=55, ha="left", va="bottom",
                     fontsize=7, color=color)
    gene_ax.spines[["left", "right", "top"]].set_visible(False)
    return fig


def write_html(path: Path, figures, low_gene_tables, threshold: int):
    sections = []
    for chrom, svg in figures:
        rows = "".join(
            "<tr>"
            f"<td>{html.escape(str(row['gene_name']))}</td>"
            f"<td>{html.escape(str(row['gene_id']))}</td>"
            f"<td>{row['start']:,}–{row['end']:,}</td>"
            f"<td>{row['zero_samples']}</td>"
            f"<td>{row['no_gt3x_samples']}</td>"
            f"<td>{row['mean_percent']:.2f}%</td>"
            "</tr>"
            for row in low_gene_tables[chrom]
        )
        sections.append(
            f"<section><h2>Chromosome {html.escape(chrom)}</h2>{svg}"
            "<h3>Lowest-coverage genes</h3>"
            "<table><thead><tr><th>Gene</th><th>Gene ID</th><th>Region</th>"
            "<th>Samples at exact 0X</th><th>Samples with no bases &gt;3X</th>"
            "<th>Cohort mean % &gt;3X</th>"
            f"</tr></thead><tbody>{rows}</tbody></table></section>"
        )
    path.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>ResolveOME coverage map</title><style>"
        "body{font-family:Arial,sans-serif;margin:24px;color:#222}"
        "section{margin-bottom:48px}svg{width:100%;height:auto;border:1px solid #ddd}"
        "table{border-collapse:collapse;width:100%;font-size:12px}"
        "th,td{border:1px solid #ccc;padding:5px;text-align:left}"
        "th{background:#eee}.note{background:#fff8dc;padding:12px;border-left:5px solid #d8a600}"
        "</style></head><body><h1>ResolveOME annotated-region coverage map</h1>"
        f"<p class='note'>Coverage is summarized over GENCODE gene-exon union intervals. "
        f'A base is considered adequately covered only when depth is greater than {threshold}X. '
        "Gray genomic bins contain no annotated exon interval and are not evidence of zero coverage."
        "</p>" + "".join(sections) + "</body></html>",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--chroms", default="")
    parser.add_argument("--bin-size", type=int, default=250_000)
    parser.add_argument("--max-gene-labels", type=int, default=28)
    parser.add_argument("--threshold", type=int, default=3)
    args = parser.parse_args()

    positions, genes = read_reports(args.input_dir)
    available = sorted(
        {chrom for intervals in positions.values() for chrom, _, _ in intervals},
        key=chromosome_key,
    )
    chroms = (
        [value.strip().removeprefix("chr") for value in args.chroms.split(",") if value.strip()]
        if args.chroms else available
    )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    html_figures = []
    tables = {}
    with PdfPages(args.output_prefix.with_suffix(".pdf")) as pdf:
        for chrom in chroms:
            sample_names, matrix, chrom_end, gene_positions = build_chromosome_matrix(
                positions, genes, chrom, args.bin_size
            )
            labels, all_low = low_coverage_genes(
                sample_names, genes, gene_positions, args.max_gene_labels
            )
            fig = make_figure(
                chrom, sample_names, matrix, chrom_end, labels, args.bin_size
            )
            pdf.savefig(fig, bbox_inches="tight")
            stream = io.StringIO()
            fig.savefig(stream, format="svg", bbox_inches="tight")
            html_figures.append((chrom, stream.getvalue()))
            tables[chrom] = all_low[:100]
            plt.close(fig)
    write_html(args.output_prefix.with_suffix(".html"), html_figures, tables, args.threshold)
    print(f"[OK] PDF: {args.output_prefix.with_suffix('.pdf')}")
    print(f"[OK] HTML: {args.output_prefix.with_suffix('.html')}")


if __name__ == "__main__":
    main()
