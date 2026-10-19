#!/usr/bin/env python3
"""Append filtered SV features to BaseJ lineage NR/NV/genotype matrices."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_matrix(path: Path) -> tuple[list[str], list[list[str]]]:
    with open(path, encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty matrix: {path}") from error
        return header, list(reader)


def align_rows(
    source_header: list[str],
    source_rows: list[list[str]],
    target_header: list[str],
    *,
    fill_missing_samples: bool = False,
) -> list[list[str]]:
    if source_header[0] != "VariantId" or target_header[0] != "VariantId":
        raise ValueError("Matrices must use VariantId as their first column")
    source_index = {sample: index for index, sample in enumerate(source_header)}
    missing = [sample for sample in target_header[1:] if sample not in source_index]
    if missing and not fill_missing_samples:
        raise ValueError(f"SV matrix is missing BaseJ samples: {', '.join(missing)}")
    return [
        [
            row[0],
            *[
                row[source_index[sample]] if sample in source_index else "0"
                for sample in target_header[1:]
            ],
        ]
        for row in source_rows
    ]


def append_matrix(base_path: Path, sv_path: Path, output_path: Path) -> None:
    base_header, base_rows = read_matrix(base_path)
    sv_header, sv_rows = read_matrix(sv_path)
    base_samples_missing_from_sv = [
        sample for sample in base_header[1:] if sample not in sv_header[1:]
    ]
    if base_samples_missing_from_sv:
        raise ValueError(
            "SV matrix is missing BaseJ samples: "
            + ", ".join(base_samples_missing_from_sv)
        )

    # The PlusMandatoryNonEmpty BaseJ bundle can omit cells with no mandatory SNVs.
    # The SV matrices retain the full passed-cell cohort, so keep their sample header
    # and zero-fill those omitted BaseJ columns.
    target_header = sv_header
    aligned_base_rows = align_rows(
        base_header,
        base_rows,
        target_header,
        fill_missing_samples=True,
    )
    aligned_sv_rows = align_rows(sv_header, sv_rows, target_header)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(target_header)
        writer.writerows(aligned_base_rows)
        writer.writerows(aligned_sv_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage-inputs", required=True)
    parser.add_argument("--sv-matrix-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--group", default="somatic_cells")
    args = parser.parse_args()

    lineage_inputs = Path(args.lineage_inputs).resolve()
    sv_matrix_dir = Path(args.sv_matrix_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    manifest_rows = []
    with open(lineage_inputs, encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle))
    group_rows = {
        row["param"]: row
        for row in manifest_rows
        if row["group"] == args.group
    }
    required = ("nr_matrix", "nv_matrix", "binary_matrix")
    missing = [name for name in required if name not in group_rows]
    if missing:
        raise SystemExit(
            f"[ERROR] Lineage manifest is missing for {args.group}: {', '.join(missing)}"
        )

    base_nr = Path(group_rows["nr_matrix"]["path"])
    base_nv = Path(group_rows["nv_matrix"]["path"])
    base_binary = Path(group_rows["binary_matrix"]["path"])
    for path in (base_nr, base_nv, base_binary):
        if not path.is_file():
            raise SystemExit(f"[ERROR] BaseJ matrix does not exist: {path}")

    augmented_nr = output_dir / f"{args.group}_NR_with_SV.tsv"
    augmented_nv = output_dir / f"{args.group}_NV_with_SV.tsv"
    augmented_genotype = output_dir / f"{args.group}_genotype_bin_with_SV.tsv"
    append_matrix(base_nr, sv_matrix_dir / "sv_NR.tsv", augmented_nr)
    append_matrix(base_nv, sv_matrix_dir / "sv_NV.tsv", augmented_nv)
    append_matrix(
        base_binary,
        sv_matrix_dir / "sv_binary.tsv",
        augmented_genotype,
    )

    updated_rows = []
    genotype_added = False
    for row in manifest_rows:
        updated = dict(row)
        if row["group"] == args.group and row["param"] == "nr_matrix":
            updated["path"] = str(augmented_nr)
        elif row["group"] == args.group and row["param"] == "nv_matrix":
            updated["path"] = str(augmented_nv)
        elif row["group"] == args.group and row["param"] == "genotype_bin":
            updated["path"] = str(augmented_genotype)
            genotype_added = True
        updated_rows.append(updated)
    if not genotype_added:
        updated_rows.append(
            {
                "group": args.group,
                "param": "genotype_bin",
                "path": str(augmented_genotype),
            }
        )

    output_manifest = output_dir / "lineage_inputs_with_sv.csv"
    with open(output_manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("group", "param", "path"))
        writer.writeheader()
        writer.writerows(updated_rows)

    print(f"[OK] Wrote SV-integrated lineage manifest: {output_manifest}")
    print("[OK] Kept binary_matrix unchanged for SNP/indel mutational signatures")


if __name__ == "__main__":
    main()
