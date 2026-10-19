#!/usr/bin/env python3
"""Exploratory clustering of detected SV features, not rooted phylogenetic inference.

Consume the filtered Step 5 SV matrices without changing their biological filters.
Use Jaccard distances and average linkage; a zero call means not detected, not
confirmed absence. Symbolic SV alleles are never treated as SNVs or small indels.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import importlib.metadata
import io
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from numpy.typing import NDArray
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform


LIMITATIONS = (
    "This is exploratory clustering of detected SV features, not a rooted "
    "phylogeny or evidence of ancestor-descendant direction. Zero means not "
    "detected, not confirmed biological absence. Detection sensitivity, "
    "coverage, SV calling artifacts, breakpoint clustering, and correlated "
    "features can affect similarity. Branch lengths are Jaccard merge "
    "distances, not evolutionary time or mutation counts. No bootstrap "
    "confidence or statistical clade support is inferred."
)


@dataclass(frozen=True)
class Matrix:
    samples: list[str]
    features: list[str]
    values: NDArray[np.int64]


def load_matrix(path: Path) -> Matrix:
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 3 or header[0] != "VariantId":
            raise ValueError(f"{path}: expected VariantId and at least two sample columns")
        samples = header[1:]
        if any(not sample for sample in samples) or len(set(samples)) != len(samples):
            raise ValueError(f"{path}: empty or duplicate sample names")
        features, values = [], []
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[0]:
                raise ValueError(f"{path}:{line_number}: invalid matrix row width or ID")
            features.append(row[0])
            try:
                values.append([int(value) for value in row[1:]])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: expected integer values") from error
    if not features or len(set(features)) != len(features):
        raise ValueError(f"{path}: no features or duplicate feature IDs")
    matrix = np.asarray(values, dtype=np.int64)
    if np.any(matrix < 0):
        raise ValueError(f"{path}: negative matrix value")
    return Matrix(samples, features, matrix)


def load_cohort(path: Path) -> tuple[list[str], list[str]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample_id", "cell_label", "pass_qc"}.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected sample_id, cell_label, pass_qc columns")
        rows = list(reader)
    if len(rows) < 2 or any(row["pass_qc"].lower() != "true" for row in rows):
        raise ValueError(f"{path}: expected at least two passed cells only")
    samples = [row["sample_id"] for row in rows]
    labels = [row["cell_label"] for row in rows]
    if any(not value for value in samples + labels) or len(set(samples)) != len(samples) or len(set(labels)) != len(labels):
        raise ValueError(f"{path}: empty or duplicate sample IDs/cell labels")
    return samples, labels


def align_matrix(matrix: Matrix, samples: list[str], features: list[str]) -> NDArray[np.int64]:
    if set(matrix.samples) != set(samples) or set(matrix.features) != set(features):
        raise ValueError("SV matrices and manifest must have identical sample/feature membership")
    sample_index = {sample: index for index, sample in enumerate(matrix.samples)}
    feature_index = {feature: index for index, feature in enumerate(matrix.features)}
    return matrix.values[np.ix_(
        [feature_index[feature] for feature in features],
        [sample_index[sample] for sample in samples],
    )]


def jaccard_counts(binary: NDArray[np.int64]) -> tuple[NDArray[np.float64], NDArray[np.int64], NDArray[np.int64]]:
    if binary.ndim != 2 or binary.shape[1] < 2 or not np.all((binary == 0) | (binary == 1)):
        raise ValueError("Expected a feature-by-sample 0/1 matrix with at least two samples")
    if np.any(binary.sum(axis=1) == 0):
        raise ValueError("SV feature has no detected cells")
    counts = binary.sum(axis=0)
    if np.any(counts == 0):
        raise ValueError("A cell has no detected SVs; Jaccard clustering would be uninformative")
    shared = binary.T @ binary
    union = counts[:, None] + counts[None, :] - shared
    distance = 1.0 - shared / union
    np.fill_diagonal(distance, 0.0)
    return distance, shared, union


def clustering_newick(tree: NDArray[np.float64], samples: list[str]) -> str:
    def render(node, parent_distance: float | None) -> str:
        if node.is_leaf():
            text = "'" + samples[node.id].replace("'", "''") + "'"
        else:
            text = f"({render(node.left, node.dist)},{render(node.right, node.dist)})"
        if parent_distance is not None:
            text += f":{max(0.0, parent_distance - node.dist):.10g}"
        return text

    return render(to_tree(tree), None) + ";\n"


def write_table(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def image_uri(figure) -> str:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=130, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-matrix-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-prevalence", type=int, default=1,
                        help="Minimum detected-cell count; use 2 for shared-feature sensitivity analysis")
    args = parser.parse_args()
    if args.min_prevalence < 1:
        parser.error("--min-prevalence must be positive")
    sources = {kind: args.sv_matrix_dir / f"sv_{kind}.tsv" for kind in ("binary", "NR", "NV")}
    sources["features"] = args.sv_matrix_dir / "sv_features.tsv"
    sources["manifest"] = args.manifest
    binary_matrix = load_matrix(sources["binary"])
    samples, labels = load_cohort(args.manifest)
    features = binary_matrix.features
    binary = align_matrix(binary_matrix, samples, features)
    nr = align_matrix(load_matrix(sources["NR"]), samples, features)
    nv = align_matrix(load_matrix(sources["NV"]), samples, features)
    if np.any(nv > nr) or not np.array_equal(binary, (nv > 0).astype(np.int64)):
        raise ValueError("SV read support and binary calls are inconsistent")
    with sources["features"].open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"VariantId", "SVTYPE", "supporting_cells", "breakpoint_tolerance"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("SV feature metadata is missing required columns")
        feature_rows = list(reader)
    metadata = {row["VariantId"]: row for row in feature_rows}
    if len(metadata) != len(feature_rows) or set(metadata) != set(features):
        raise ValueError("SV metadata and matrices have inconsistent or duplicate feature IDs")
    prevalence = binary.sum(axis=1)
    if any(int(metadata[feature]["supporting_cells"]) != count for feature, count in zip(features, prevalence)):
        raise ValueError("SV feature prevalence disagrees with metadata supporting_cells")
    source_feature_count = len(features)
    keep = prevalence >= args.min_prevalence
    features = [feature for feature, included in zip(features, keep) if included]
    binary, nr, nv, prevalence = binary[keep], nr[keep], nv[keep], prevalence[keep]
    distance, shared, union = jaccard_counts(binary)
    tree = linkage(squareform(distance, checks=True), method="average", optimal_ordering=True)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "jaccard_distance.tsv", ["sample_id", *samples],
                ([sample, *[f"{value:.10f}" for value in row]] for sample, row in zip(samples, distance)))
    write_table(out / "pairwise_sv_similarity.tsv",
                ["sample_a", "sample_b", "detected_a", "detected_b", "shared_features", "union_features", "jaccard_similarity", "jaccard_distance"],
                ([samples[i], samples[j], int(shared[i, i]), int(shared[j, j]), int(shared[i, j]), int(union[i, j]), f"{1-distance[i,j]:.10f}", f"{distance[i,j]:.10f}"]
                 for i in range(len(samples)) for j in range(i+1, len(samples))))
    sv_types = sorted({row["SVTYPE"] for row in metadata.values()})
    type_vector = np.asarray([metadata[feature]["SVTYPE"] for feature in features])
    counts = binary.sum(axis=0)
    sample_rows = []
    for index, (sample, label) in enumerate(zip(samples, labels)):
        detected = binary[:, index] == 1
        baf = nv[detected, index] / nr[detected, index]
        nearest = min((j for j in range(len(samples)) if j != index), key=lambda j: (distance[index, j], j))
        sample_rows.append([
            sample, label, int(counts[index]), int(np.sum(detected & (prevalence == 1))),
            int(np.sum(detected & (prevalence >= 2))), f"{np.median(baf):.6f}",
            samples[nearest], f"{distance[index,nearest]:.10f}",
            *[int(np.sum(detected & (type_vector == kind))) for kind in sv_types],
        ])
    write_table(out / "sample_sv_summary.tsv",
                ["sample_id", "cell_label", "detected_features", "private_features", "shared_features", "median_detected_baf", "nearest_cell", "nearest_jaccard_distance", *sv_types],
                sample_rows)
    write_table(out / "sv_feature_prevalence.tsv",
                ["VariantId", "SVTYPE", "detected_cells", "breakpoint_tolerance", "detected_samples"],
                ([feature, metadata[feature]["SVTYPE"], int(prevalence[i]), metadata[feature]["breakpoint_tolerance"],
                  ",".join(sample for sample, present in zip(samples, binary[i]) if present)]
                 for i, feature in enumerate(features)))
    write_table(out / "average_linkage.tsv", ["left_node", "right_node", "jaccard_merge_distance", "cell_count"],
                ([int(left), int(right), f"{height:.10f}", int(count)] for left, right, height, count in tree))
    (out / "sv_similarity.dendrogram.newick").write_text(clustering_newick(tree, samples), encoding="utf-8")

    figure, axis = plt.subplots(figsize=(12, 5))
    dendro = dendrogram(tree, labels=labels, ax=axis, leaf_rotation=45, color_threshold=0, above_threshold_color="#35618d")
    axis.set_title("Detected SV similarity: average-linkage clustering (not ancestry)")
    axis.set_ylabel("Jaccard merge distance")
    axis.set_ylim(0, 1.05)
    figure.tight_layout()
    order = dendro["leaves"]
    heatmap, heat_axis = plt.subplots(figsize=(10, 9))
    image = heat_axis.imshow(distance[np.ix_(order, order)], vmin=0, vmax=1, cmap="viridis")
    ordered_labels = [labels[index] for index in order]
    heat_axis.set_xticks(range(len(samples)), ordered_labels, rotation=90)
    heat_axis.set_yticks(range(len(samples)), ordered_labels)
    heat_axis.set_title("Detected SV Jaccard distances (0 = identical detected sets)")
    heatmap.colorbar(image, ax=heat_axis, label="Jaccard distance")
    heatmap.tight_layout()
    figure.savefig(out / "sv_similarity.dendrogram.png", dpi=160, bbox_inches="tight")
    heatmap.savefig(out / "sv_similarity.distance_heatmap.png", dpi=160, bbox_inches="tight")
    with PdfPages(out / "sv_similarity.report.pdf") as pdf:
        pdf.infodict()["Title"] = "Exploratory detected-SV similarity; not rooted ancestry"
        pdf.infodict()["Subject"] = LIMITATIONS
        for plot in (figure, heatmap):
            pdf.savefig(plot, bbox_inches="tight")
    rows_html = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row[:8]) + "</tr>" for row in sample_rows)
    headers = ["Sample", "Cell", "Detected", "Private", "Shared", "Median detected BAF", "Nearest cell", "Jaccard distance"]
    (out / "sv_similarity.report.html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>Exploratory SV similarity</title>"
        "<style>body{font-family:sans-serif;max-width:1400px;margin:24px auto;padding:16px}"
        "td,th{padding:6px;border:1px solid #ddd}table{border-collapse:collapse;font-size:12px}img{max-width:100%}</style>"
        "</head><body><h1>Exploratory detected-SV similarity</h1>"
        f"<p>{len(samples)} passed cells; {len(features):,} filtered SV features; "
        f"{int(np.sum(prevalence >= 2)):,} features detected in multiple cells.</p>"
        "<p>Jaccard distance = 1 - shared detected features / union of detected features. "
        f"Features detected in at least {args.min_prevalence} cell(s) are included with equal weight. "
        "Average linkage joins clusters by their mean pairwise Jaccard distance. "
        "Newick is a clustering serialization, not a biologically rooted tree.</p>"
        f"<p><strong>Limitations:</strong> {html.escape(LIMITATIONS)}</p>"
        f"<img alt='Exploratory SV dendrogram' src='{image_uri(figure)}'>"
        f"<img alt='Jaccard distance heatmap' src='{image_uri(heatmap)}'>"
        "<h2>Per-cell detected-feature summary</h2><table><thead><tr>"
        + "".join(f"<th>{header}</th>" for header in headers)
        + f"</tr></thead><tbody>{rows_html}</tbody></table></body></html>",
        encoding="utf-8",
    )
    plt.close(figure)
    plt.close(heatmap)
    provenance = {
        "analysis": "exploratory_detected_sv_similarity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "distance": "Jaccard", "linkage": "average", "rooted_phylogeny": False,
        "zero_call_semantics": "not_detected", "limitations": LIMITATIONS,
        "samples": samples, "feature_count": len(features),
        "source_feature_count": source_feature_count, "minimum_prevalence": args.min_prevalence,
        "multi_cell_feature_count": int(np.sum(prevalence >= 2)),
        "detected_cell_feature_pairs": int(binary.sum()),
        "breakpoint_tolerances": sorted({int(row["breakpoint_tolerance"]) for row in metadata.values()}),
        "inherited_filters": f"Step 5 filtered SV matrices; additional minimum detected-cell count={args.min_prevalence}",
        "inputs": {key: {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for key, path in sources.items()},
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")},
    }
    (out / "analysis_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Exploratory SV similarity: {len(samples)} cells, {len(features)} features")
    print(f"[OK] Reports and tables: {out.resolve()}")
    print("[NOTE] Dendrogram is similarity clustering, not rooted ancestry.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"[ERROR] {error}") from error
