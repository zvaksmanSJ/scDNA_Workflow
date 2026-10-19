#!/usr/bin/env python3
"""Add GENCODE region and gene-level views to combined single-cell CNV matrices."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
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


PRIMARY_CHROMS = [str(value) for value in range(1, 23)] + ["X", "Y"]


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    label: str


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    gene_name: str
    gene_id: str
    gene_type: str

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start + 1}-{self.end}"


def contig_key(chrom: str) -> tuple[int, int | str]:
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)


def parse_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for item in value.rstrip(";").split(";"):
        key, separator, field_value = item.strip().partition(" ")
        if separator:
            attributes[key] = field_value.strip().strip('"')
    return attributes


def load_genes(path: Path) -> list[Gene]:
    opener = gzip.open if path.suffix == ".gz" else open
    genes = []
    seen = set()
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom = fields[0].removeprefix("chr")
            if chrom not in PRIMARY_CHROMS:
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id", "").split(".")[0]
            if not gene_id or gene_id in seen:
                continue
            seen.add(gene_id)
            genes.append(
                Gene(
                    chrom=chrom,
                    start=int(fields[3]) - 1,
                    end=int(fields[4]),
                    gene_name=attributes.get("gene_name", gene_id),
                    gene_id=gene_id,
                    gene_type=attributes.get("gene_type", "unknown"),
                )
            )
    return sorted(
        genes,
        key=lambda gene: (contig_key(gene.chrom), gene.start, gene.end, gene.gene_id),
    )


def read_matrix(path: Path) -> tuple[list[str], list[Region], np.ndarray]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:4] != ["chromosome", "start", "end", "region"]:
            raise ValueError(f"Unexpected matrix header in {path}: {header[:4]}")
        samples = header[4:]
        regions = []
        rows = []
        for row in reader:
            regions.append(Region(row[0], int(row[1]), int(row[2]), row[3]))
            values = []
            for value in row[4:]:
                if value == "NA":
                    values.append(math.nan)
                elif value in {"loss", "neutral", "gain"}:
                    values.append({"loss": -1.0, "neutral": 0.0, "gain": 1.0}[value])
                else:
                    values.append(float(value))
            rows.append(values)
    return samples, regions, np.asarray(rows, dtype=float)


def load_thresholds(path: Path, samples: list[str]) -> np.ndarray:
    with path.open(encoding="utf-8") as handle:
        rows = {
            row["sample"]: float(row["gain_threshold"])
            for row in csv.DictReader(handle, delimiter="\t")
        }
    missing = [sample for sample in samples if sample not in rows]
    if missing:
        raise ValueError(f"Cell QC is missing samples: {missing}")
    return np.asarray([rows[sample] for sample in samples], dtype=float)


def map_region_genes(
    regions: list[Region],
    genes: list[Gene],
) -> tuple[list[list[Gene]], dict[str, list[int]]]:
    genes_by_chrom: dict[str, list[Gene]] = defaultdict(list)
    for gene in genes:
        genes_by_chrom[gene.chrom].append(gene)
    region_genes = []
    gene_regions: dict[str, list[int]] = defaultdict(list)
    cursors: dict[str, int] = defaultdict(int)
    for region_index, region in enumerate(regions):
        chromosome_genes = genes_by_chrom.get(region.chrom, [])
        cursor = cursors[region.chrom]
        while cursor < len(chromosome_genes) and chromosome_genes[cursor].end <= region.start:
            cursor += 1
        cursors[region.chrom] = cursor
        overlaps = []
        scan = cursor
        while scan < len(chromosome_genes) and chromosome_genes[scan].start < region.end:
            gene = chromosome_genes[scan]
            if gene.end > region.start:
                overlaps.append(gene)
                gene_regions[gene.gene_id].append(region_index)
            scan += 1
        region_genes.append(overlaps)
    return region_genes, gene_regions


def write_annotated_region_matrix(
    path: Path,
    samples: list[str],
    regions: list[Region],
    region_genes: list[list[Gene]],
    matrix: np.ndarray,
    state: bool,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chromosome",
                "start",
                "end",
                "region",
                "gene_names",
                "gene_ids",
                "gene_types",
                *samples,
            ]
        )
        for index, region in enumerate(regions):
            genes = region_genes[index]
            values = []
            for value in matrix[index]:
                if not math.isfinite(value):
                    values.append("NA")
                elif state:
                    values.append({-1: "loss", 0: "neutral", 1: "gain"}[int(value)])
                else:
                    values.append(f"{value:.6f}")
            writer.writerow(
                [
                    region.chrom,
                    region.start,
                    region.end,
                    region.label,
                    ";".join(sorted({gene.gene_name for gene in genes})),
                    ";".join(sorted({gene.gene_id for gene in genes})),
                    ";".join(sorted({gene.gene_type for gene in genes})),
                    *values,
                ]
            )


def write_annotated_regional_summary(
    input_path: Path,
    output_path: Path,
    region_genes: list[list[Gene]],
) -> None:
    with input_path.open(encoding="utf-8") as source, output_path.open(
        "w", newline="", encoding="utf-8"
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        insert_at = fieldnames.index("region") + 1
        fieldnames[insert_at:insert_at] = ["gene_names", "gene_ids", "gene_types"]
        writer = csv.DictWriter(
            destination, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, row in enumerate(reader):
            genes = region_genes[index]
            row["gene_names"] = ";".join(
                sorted({gene.gene_name for gene in genes})
            )
            row["gene_ids"] = ";".join(sorted({gene.gene_id for gene in genes}))
            row["gene_types"] = ";".join(
                sorted({gene.gene_type for gene in genes})
            )
            writer.writerow(row)


def aggregate_genes(
    genes: list[Gene],
    regions: list[Region],
    gene_regions: dict[str, list[int]],
    raw: np.ndarray,
    adjusted: np.ndarray,
    thresholds: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gene_raw = np.full((len(genes), raw.shape[1]), np.nan)
    gene_adjusted = np.full_like(gene_raw, np.nan)
    gene_states = np.full_like(gene_raw, np.nan)
    for gene_index, gene in enumerate(genes):
        indices = gene_regions.get(gene.gene_id, [])
        if not indices:
            continue
        weights = np.asarray(
            [
                max(
                    0,
                    min(gene.end, regions[index].end)
                    - max(gene.start, regions[index].start),
                )
                for index in indices
            ],
            dtype=float,
        )
        for sample_index in range(raw.shape[1]):
            for source, target in (
                (raw, gene_raw),
                (adjusted, gene_adjusted),
            ):
                values = source[indices, sample_index]
                valid = np.isfinite(values) & (weights > 0)
                if valid.any():
                    target[gene_index, sample_index] = np.average(
                        values[valid], weights=weights[valid]
                    )
            value = gene_adjusted[gene_index, sample_index]
            if math.isfinite(value):
                threshold = thresholds[sample_index]
                gene_states[gene_index, sample_index] = (
                    1 if value >= threshold else -1 if value <= -threshold else 0
                )
    return gene_raw, gene_adjusted, gene_states


def write_gene_matrix(
    path: Path,
    genes: list[Gene],
    samples: list[str],
    matrix: np.ndarray,
    state: bool,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chromosome",
                "start",
                "end",
                "gene_name",
                "gene_id",
                "gene_type",
                *samples,
            ]
        )
        for index, gene in enumerate(genes):
            values = []
            for value in matrix[index]:
                if not math.isfinite(value):
                    values.append("NA")
                elif state:
                    values.append({-1: "loss", 0: "neutral", 1: "gain"}[int(value)])
                else:
                    values.append(f"{value:.6f}")
            writer.writerow(
                [
                    gene.chrom,
                    gene.start,
                    gene.end,
                    gene.gene_name,
                    gene.gene_id,
                    gene.gene_type,
                    *values,
                ]
            )


def finite_median(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else math.nan


def build_gene_summary(
    genes: list[Gene],
    raw: np.ndarray,
    adjusted: np.ndarray,
    states: np.ndarray,
    recurrent_fraction: float,
) -> list[dict[str, object]]:
    rows = []
    for index, gene in enumerate(genes):
        row_states = states[index]
        callable_cells = int(np.isfinite(row_states).sum())
        gains = int(np.sum(row_states == 1))
        losses = int(np.sum(row_states == -1))
        neutral = int(np.sum(row_states == 0))
        gain_fraction = gains / callable_cells if callable_cells else 0
        loss_fraction = losses / callable_cells if callable_cells else 0
        if gain_fraction >= recurrent_fraction and loss_fraction >= recurrent_fraction:
            combined = "mixed"
        elif gain_fraction >= recurrent_fraction:
            combined = "recurrent_gain"
        elif loss_fraction >= recurrent_fraction:
            combined = "recurrent_loss"
        else:
            combined = "no_recurrent_change"
        rows.append(
            {
                "chromosome": gene.chrom,
                "start": gene.start,
                "end": gene.end,
                "region": gene.region,
                "gene_name": gene.gene_name,
                "gene_id": gene.gene_id,
                "gene_type": gene.gene_type,
                "callable_cells": callable_cells,
                "gain_cells": gains,
                "neutral_cells": neutral,
                "loss_cells": losses,
                "gain_fraction": round(gain_fraction, 6),
                "loss_fraction": round(loss_fraction, 6),
                "combined_state": combined,
                "median_raw_log2": round(finite_median(raw[index]), 6),
                "median_adjusted_log2": round(finite_median(adjusted[index]), 6),
            }
        )
    return rows


def render_top_genes(rows: list[dict[str, object]]):
    ranked = sorted(
        rows,
        key=lambda row: max(float(row["gain_fraction"]), float(row["loss_fraction"])),
        reverse=True,
    )[:80]
    ranked.reverse()
    figure, axis = plt.subplots(figsize=(12, 18), constrained_layout=True)
    y = np.arange(len(ranked))
    axis.barh(y, [float(row["gain_fraction"]) for row in ranked], color="#b2182b", label="Gain")
    axis.barh(y, [-float(row["loss_fraction"]) for row in ranked], color="#2166ac", label="Loss")
    axis.set_yticks(
        y,
        [f"{row['gene_name']} ({row['chromosome']})" for row in ranked],
        fontsize=7,
    )
    axis.axvline(0, color="black", linewidth=0.6)
    axis.set_xlim(-1, 1)
    axis.set_xlabel("Fraction of cells (loss plotted negative)")
    axis.set_title("Top recurrent gene-associated CNV calls")
    axis.legend()
    return figure


def render_chromosome_genes(rows: list[dict[str, object]], chrom: str):
    selected = [row for row in rows if row["chromosome"] == chrom]
    figure, axis = plt.subplots(figsize=(15, 5), constrained_layout=True)
    x = np.asarray([(int(row["start"]) + int(row["end"])) / 2e6 for row in selected])
    gains = np.asarray([float(row["gain_fraction"]) for row in selected])
    losses = -np.asarray([float(row["loss_fraction"]) for row in selected])
    axis.scatter(x, gains, s=9, color="#b2182b", alpha=0.65, label="Gain")
    axis.scatter(x, losses, s=9, color="#2166ac", alpha=0.65, label="Loss")
    candidates = sorted(
        selected,
        key=lambda row: max(float(row["gain_fraction"]), float(row["loss_fraction"])),
        reverse=True,
    )[:20]
    for row in candidates:
        y = (
            float(row["gain_fraction"])
            if float(row["gain_fraction"]) >= float(row["loss_fraction"])
            else -float(row["loss_fraction"])
        )
        axis.annotate(
            str(row["gene_name"]),
            ((int(row["start"]) + int(row["end"])) / 2e6, y),
            fontsize=6,
            rotation=35,
        )
    axis.axhline(0, color="black", linewidth=0.6)
    axis.set_ylim(-1.05, 1.05)
    axis.set_xlabel("Position (Mb)")
    axis.set_ylabel("Fraction of cells\n(loss plotted negative)")
    axis.set_title(f"Chromosome {chrom}: gene-associated gain/loss prevalence")
    axis.legend(loc="upper right")
    return figure


def data_uri(figure) -> str:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=130, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def write_report(
    prefix: Path,
    rows: list[dict[str, object]],
    samples: list[str],
    bin_size: int,
) -> None:
    top_figure = render_top_genes(rows)
    image = data_uri(top_figure)
    plt.close(top_figure)
    ranked = sorted(
        rows,
        key=lambda row: max(float(row["gain_fraction"]), float(row["loss_fraction"])),
        reverse=True,
    )[:200]
    table_rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in (
                "gene_name",
                "gene_id",
                "gene_type",
                "region",
                "gain_cells",
                "loss_cells",
                "gain_fraction",
                "loss_fraction",
                "combined_state",
            )
        )
        + "</tr>"
        for row in ranked
    )
    Path(f"{prefix}.gene_report.html").write_text(
        f"""<!doctype html><html><head><meta charset="utf-8">
<title>Gene-annotated single-cell CNV report</title>
<style>body{{font-family:Arial,sans-serif;margin:24px}} img{{max-width:100%}}
table{{border-collapse:collapse;font-size:12px}} th,td{{border:1px solid #bbb;padding:4px 7px}}
.note{{background:#f2f5f8;padding:12px;border-left:5px solid #4c78a8}}</style>
</head><body><h1>Gene-annotated combined single-cell CNV report</h1>
<div class="note"><strong>{len(samples)} cells; {bin_size:,} bp bins; GENCODE genes.</strong>
Gene values are overlap-length-weighted from the underlying bins. Gain/loss calls
use the same cell-specific thresholds as the coordinate-only analysis. A gene may
overlap multiple bins.</div>
<h2>Top recurrent genes</h2><img src="{image}">
<h2>Top 200 gene-associated calls</h2><table><thead><tr>
<th>Gene</th><th>Gene ID</th><th>Type</th><th>Region</th><th>Gain cells</th>
<th>Loss cells</th><th>Gain fraction</th><th>Loss fraction</th><th>State</th>
</tr></thead><tbody>{table_rows}</tbody></table></body></html>
""",
        encoding="utf-8",
    )
    with PdfPages(f"{prefix}.gene_report.pdf") as pdf:
        figure = render_top_genes(rows)
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)
        for chrom in PRIMARY_CHROMS:
            figure = render_chromosome_genes(rows, chrom)
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--recurrent-min-fraction", type=float, default=0.20)
    args = parser.parse_args()
    prefix = Path(args.input_prefix).resolve()

    samples, regions, raw = read_matrix(Path(f"{prefix}.raw_log2_matrix.tsv"))
    adjusted_samples, adjusted_regions, adjusted = read_matrix(
        Path(f"{prefix}.adjusted_log2_matrix.tsv")
    )
    state_samples, state_regions, states = read_matrix(
        Path(f"{prefix}.state_matrix.tsv")
    )
    if samples != adjusted_samples or samples != state_samples:
        raise ValueError("Sample order differs across input matrices")
    if regions != adjusted_regions or regions != state_regions:
        raise ValueError("Region order differs across input matrices")
    thresholds = load_thresholds(Path(f"{prefix}.cell_qc.tsv"), samples)
    genes = load_genes(Path(args.gtf))
    region_genes, gene_regions = map_region_genes(regions, genes)
    write_annotated_regional_summary(
        Path(f"{prefix}.regional_summary.tsv"),
        Path(f"{prefix}.gene_annotated.regional_summary.tsv"),
        region_genes,
    )

    for suffix, matrix, state in (
        ("raw_log2_matrix", raw, False),
        ("adjusted_log2_matrix", adjusted, False),
        ("state_matrix", states, True),
    ):
        write_annotated_region_matrix(
            Path(f"{prefix}.gene_annotated.{suffix}.tsv"),
            samples,
            regions,
            region_genes,
            matrix,
            state,
        )

    gene_raw, gene_adjusted, gene_states = aggregate_genes(
        genes, regions, gene_regions, raw, adjusted, thresholds
    )
    for suffix, matrix, state in (
        ("raw_log2_matrix", gene_raw, False),
        ("adjusted_log2_matrix", gene_adjusted, False),
        ("state_matrix", gene_states, True),
    ):
        write_gene_matrix(
            Path(f"{prefix}.gene.{suffix}.tsv"),
            genes,
            samples,
            matrix,
            state,
        )
    summary_rows = build_gene_summary(
        genes,
        gene_raw,
        gene_adjusted,
        gene_states,
        args.recurrent_min_fraction,
    )
    with Path(f"{prefix}.gene_summary.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    write_report(prefix, summary_rows, samples, args.bin_size)
    print(f"[OK] Annotated {len(regions)} regions with {len(genes)} genes")
    print(f"[OK] Gene HTML: {prefix}.gene_report.html")
    print(f"[OK] Gene PDF: {prefix}.gene_report.pdf")


if __name__ == "__main__":
    main()
