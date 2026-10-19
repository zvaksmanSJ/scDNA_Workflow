#!/usr/bin/env python3
"""Reorganize a WDL workflow run's outputs into a clean, browsable folder tree.

Consumes either of the two common engine output formats for ``wdl/main.wdl``:
  - miniwdl's ``miniwdl run -o outputs.json`` flat format:
    ``{"ResolveOmeDnaWorkflow.<output_name>": <value>, ...}``
  - Cromwell's ``cromwell run ... -m metadata.json`` format, which nests the
    same flat mapping one level deeper under an ``"outputs"`` key:
    ``{"outputs": {"ResolveOmeDnaWorkflow.<output_name>": <value>, ...}, ...}``
    (Cromwell's plain stdout result from `run`, `{"outputs": {...}, "id": ...}`,
    also works since it has the same shape.)

and lays every output
file out under a destination directory organized by pipeline component so the
final results for each stage are easy to find and "stitch together", and so
all the supporting run machinery (logs, scripts, JSON configs) is isolated
into its own folders that can be deleted independently of the actual results.

Destination layout:

    <dest>/
      per_cell/<sample_id>/{alignment,qc,variants,cnv,sv,phasing,loh_ploidy}/
      cohort/{qc,phasing,loh_ploidy,mutational_signatures,lineage}/
      dashboards/dashboard.html
      logs/<task_name>/...              (only if --run-dir is given)
      run_scripts/                       (inputs.json + any .sh/.json you pass in)

By default files are **symlinked** (not copied) to avoid duplicating
potentially large BAM/VCF data; pass --copy to make real copies instead.
This mirrors the same non-destructive, re-runnable design used by
extra_methods/organize_final_results.py for the real cohort run.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

# Workflow output name -> (per-cell dest subdir, is per-cell-indexed-array)
PER_CELL_OUTPUTS = {
    "dedup_bams": "alignment",
    "dedup_bam_bais": "alignment",
    "flagstats": "alignment",
    "idxstats": "alignment",
    "wgs_metrics": "qc",
    "deepvariant_vep_snpeff_vcfs": "variants",
    "deepvariant_vep_snpeff_vcf_csis": "variants",
    "snpeff_htmls": "variants",
    "snpeff_genes_csvs": "variants",
    "cnvkit_cns": "cnv",
    "cnvkit_cnr": "cnv",
    "cnvkit_call_cns": "cnv",
    "cnvkit_genemetrics": "cnv",
    "manta_vcfs": "sv",
    "phased_vcfs": "phasing",
    "phased_vcf_csis": "phasing",
    "freec_cnvs_fixed_ploidy": "loh_ploidy/fixed_ploidy",
    "freec_ratio_fixed_ploidy": "loh_ploidy/fixed_ploidy",
    "freec_info_fixed_ploidy": "loh_ploidy/fixed_ploidy",
    "freec_cnvs_candidate_ploidy": "loh_ploidy/candidate_ploidy",
    "freec_ratio_candidate_ploidy": "loh_ploidy/candidate_ploidy",
    "freec_info_candidate_ploidy": "loh_ploidy/candidate_ploidy",
}
# Per-cell outputs that are themselves nested arrays (Array[Array[File]]):
# one inner array of files per cell, flattened into that cell's QC folder.
PER_CELL_NESTED_OUTPUTS = {
    "fastp_out": "qc/fastp",
    "fastqc_raw_out": "qc/fastqc_raw",
    "fastqc_trimmed_out": "qc/fastqc_trimmed",
    "picard_out": "qc/picard",
}
# Cohort-level outputs (single file or flat array, not indexed by sample).
COHORT_OUTPUTS = {
    "multiqc_report": "qc/multiqc",
    "multiqc_data": "qc/multiqc",
    "qc_summary_tsv": "qc",
    "qc_overview_png": "qc",
    "mapd_qc_tsv": "qc",
    "fdr_estimate_tsv": "qc",
    "artifact_panel_summary_tsv": "qc",
    "artifact_panel_sites_tsv": "qc",
    "phasing_summary_tsv": "phasing",
    "loh_summary_tsv": "loh_ploidy",
    "ploidy_summary_tsv": "loh_ploidy",
    "loh_by_region_tsv": "loh_ploidy",
    "loh_by_region_png": "loh_ploidy",
    "signature_contributions_tsv": "mutational_signatures",
    "signature_cosine_similarity_tsv": "mutational_signatures",
    "mutation_spectrum_png": "mutational_signatures",
    "signature_barplot_png": "mutational_signatures",
    "lineage_outputs": "lineage",
    "lineage_nv_matrix": "lineage",
    "lineage_nr_matrix": "lineage",
    "bootstrap_support_report": "lineage",
    "verification_summary_json": "lineage",
}
DASHBOARD_OUTPUTS = {"dashboard_html": "dashboards"}

WORKFLOW_NAME = "ResolveOmeDnaWorkflow"


def strip_prefix(key: str) -> str:
    return key[len(WORKFLOW_NAME) + 1 :] if key.startswith(WORKFLOW_NAME + ".") else key


def place(src: str, dest_dir: Path, link: bool) -> None:
    src_path = Path(src)
    if not src_path.exists():
        print(f"[WARN] missing source file, skipping: {src}", file=sys.stderr)
        return
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src_path.name
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if link:
        dest.symlink_to(src_path.resolve())
    else:
        shutil.copy2(src_path, dest)


def organize_logs(run_dir: Path, dest: Path, link: bool) -> int:
    """Copy/symlink every *.log, stdout, stderr found under an engine run dir
    into logs/<task_name>/<call_dir_name>/, keeping logs completely separate
    from result files."""
    n = 0
    for log_file in run_dir.rglob("*"):
        if not log_file.is_file():
            continue
        if log_file.name not in ("stdout", "stderr") and not log_file.name.endswith(".log"):
            continue
        # miniwdl call dirs look like .../call-<task_name>/... or
        # .../call-<task_name>/shard-<i>/...; use the nearest call-* ancestor
        # (and its shard, if any) as the grouping key.
        parts = [p for p in log_file.relative_to(run_dir).parts if p.startswith("call-") or p.startswith("shard-")]
        group = "/".join(p.split("-", 1)[1] if "-" in p else p for p in parts) or "workflow"
        place(str(log_file), dest / "logs" / group, link)
        n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-json", required=True, type=Path,
                     help="outputs.json from `miniwdl run -o outputs.json main.wdl ...`")
    ap.add_argument("--sample-ids", required=True, nargs="+",
                     help="sample IDs in the SAME order used for sample_ids/r1_fastqs/r2_fastqs in the run")
    ap.add_argument("--dest", required=True, type=Path, help="destination directory for the organized results")
    ap.add_argument("--run-dir", type=Path, default=None,
                     help="optional: the engine's run directory, to also collect logs into dest/logs/")
    ap.add_argument("--inputs-json", type=Path, default=None,
                     help="optional: the inputs.json used for this run, copied into dest/run_scripts/")
    ap.add_argument("--extra-script", action="append", default=[],
                     help="optional: additional .sh/.json file(s) to copy into dest/run_scripts/ (repeatable)")
    ap.add_argument("--copy", action="store_true", help="copy files instead of symlinking (uses more disk)")
    args = ap.parse_args()

    link = not args.copy
    raw = json.loads(args.outputs_json.read_text())
    # Cromwell nests the flat output mapping one level deeper, under "outputs"
    # (both its -m metadata.json and its plain `run` stdout have this shape);
    # miniwdl's -o outputs.json is already the flat mapping itself.
    outputs = raw["outputs"] if isinstance(raw, dict) and "outputs" in raw else raw
    outputs = {strip_prefix(k): v for k, v in outputs.items()}
    sample_ids = args.sample_ids
    dest = args.dest

    placed, missing_keys = 0, []

    for key, subdir in PER_CELL_OUTPUTS.items():
        if key not in outputs:
            missing_keys.append(key)
            continue
        values = outputs[key]
        if len(values) != len(sample_ids):
            print(f"[WARN] {key}: {len(values)} values but {len(sample_ids)} sample_ids -- skipping", file=sys.stderr)
            continue
        for sid, v in zip(sample_ids, values):
            place(v, dest / "per_cell" / sid / subdir, link)
            placed += 1

    for key, subdir in PER_CELL_NESTED_OUTPUTS.items():
        if key not in outputs:
            missing_keys.append(key)
            continue
        values = outputs[key]
        if len(values) != len(sample_ids):
            print(f"[WARN] {key}: {len(values)} inner arrays but {len(sample_ids)} sample_ids -- skipping", file=sys.stderr)
            continue
        for sid, file_list in zip(sample_ids, values):
            for v in file_list:
                place(v, dest / "per_cell" / sid / subdir, link)
                placed += 1

    for key, subdir in {**COHORT_OUTPUTS, **DASHBOARD_OUTPUTS}.items():
        if key not in outputs:
            missing_keys.append(key)
            continue
        value = outputs[key]
        values = value if isinstance(value, list) else [value]
        for v in values:
            place(v, dest / subdir, link)
            placed += 1

    if missing_keys:
        print(f"[INFO] {len(missing_keys)} expected output key(s) not present in outputs.json "
              f"(fine if this was a partial/subset run): {', '.join(missing_keys)}", file=sys.stderr)

    n_logs = 0
    if args.run_dir:
        n_logs = organize_logs(args.run_dir, dest, link)

    run_scripts_dir = dest / "run_scripts"
    run_scripts_dir.mkdir(parents=True, exist_ok=True)
    if args.inputs_json:
        place(str(args.inputs_json), run_scripts_dir, link)
    for extra in args.extra_script:
        place(extra, run_scripts_dir, link)
    place(str(args.outputs_json), run_scripts_dir, link)

    readme = dest / "README.md"
    readme.write_text(
        "# ResolveOME DNA workflow -- organized results\n\n"
        "- `per_cell/<sample_id>/` -- every per-cell result, grouped by pipeline stage\n"
        "  (alignment, qc, variants, cnv, sv, phasing, loh_ploidy)\n"
        "- `cohort/` -- cohort-wide aggregates (qc, phasing, loh_ploidy, mutational_signatures, lineage)\n"
        "- `dashboards/dashboard.html` -- the final tabbed dashboard; open this first\n"
        "- `logs/` -- **logs only**; safe to delete entirely without losing any result\n"
        "- `run_scripts/` -- the inputs.json/outputs.json (and any extra .sh files) used for this\n"
        "  run; safe to delete once you no longer need to reproduce/audit the run\n\n"
        f"Generated by wdl/scripts/organize_outputs.py from {args.outputs_json}.\n"
    )

    print(f"[DONE] placed {placed} result file(s), {n_logs} log file(s) under {dest}")
    print(f"[DONE] wrote {readme}")


if __name__ == "__main__":
    main()
