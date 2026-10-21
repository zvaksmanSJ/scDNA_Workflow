#!/usr/bin/env python3
"""Build (or refresh) a single, browsable 'final results' tree for this DNA
cohort, without moving, copying, or duplicating any large existing output.

Every entry under final-root/components/<name>/ is a symlink to the real,
already-existing directory or file produced by the corresponding pipeline
stage (core BaseJumper-parity stage or an extra_methods addition). Only the
dashboard HTML itself and small aggregate summary tables/figures are
physically written under final-root/dashboards/, since those ARE the final
artifact (not a pointer to something bigger elsewhere).

Safe to re-run any time a stage produces new/updated output -- symlinks are
recreated each run, nothing is deleted outside the final-results tree.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def relink(link_path: Path, target: Path) -> None:
    """(Re)point a symlink at `target`, creating parent dirs as needed."""
    if not target.exists():
        print(f"[SKIP] target does not exist yet: {target}")
        return
    link_path.parent.mkdir(parents=True, exist_ok=True)
    if link_path.is_symlink() or link_path.exists():
        link_path.unlink()
    link_path.symlink_to(target)
    print(f"[LINK] {link_path} -> {target}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dna-results-dir", required=True, type=Path)
    p.add_argument("--final-dir", required=True, type=Path,
                    help="Root of the organized final-results tree, e.g. DNA_results/final")
    return p.parse_args()


def main():
    args = parse_args()
    D = args.dna_results_dir
    F = args.final_dir
    EM = D / "extra_methods"
    comp = F / "components"

    # --- core BaseJumper-parity pipeline stages ---
    relink(comp / "variants_per_cell", D / "variants")
    relink(comp / "cnv_cnvkit_per_cell", D / "cnv")
    relink(comp / "sv_manta_per_cell", D / "sv")
    relink(comp / "bam_per_cell", D / "bam")
    relink(comp / "qc_raw" / "fastp", D / "qc" / "fastp")
    relink(comp / "qc_raw" / "fastqc_raw", D / "qc" / "fastqc_raw")
    relink(comp / "qc_raw" / "fastqc_trimmed", D / "qc" / "fastqc_trimmed")
    relink(comp / "qc_raw" / "picard", D / "qc" / "picard")
    relink(comp / "qc_raw" / "preseq", D / "qc" / "preseq")
    relink(comp / "multiqc", D / "multiqc_step3")

    # --- extra_methods additions (this session) ---
    relink(comp / "mapd_qc.tsv", EM / "mapd_qc.tsv")
    relink(comp / "cnv_overview.tsv", EM / "cnv_overview.tsv")
    relink(comp / "fdr_estimate.tsv", EM / "fdr_estimate.tsv")
    relink(comp / "artifact_panel", EM)  # summary + sites tsv live at EM root; see README
    relink(comp / "phasing_per_cell", D / "extra_methods" / "phasing")
    relink(comp / "phasing_summary.tsv", EM / "phasing_summary.tsv")
    relink(comp / "loh_per_cell", D / "extra_methods" / "freec_loh")
    relink(comp / "loh_ploidy_per_cell", D / "extra_methods" / "freec_loh_ploidy")
    relink(comp / "loh_summary.tsv", EM / "loh_summary.tsv")
    relink(comp / "loh_by_region.tsv", EM / "loh_by_region.tsv")
    relink(comp / "ploidy_summary.tsv", EM / "ploidy_summary.tsv")
    relink(comp / "ploidy_candidates.tsv", EM / "ploidy_candidates.tsv")
    relink(comp / "ploidy_narrative.tsv", EM / "ploidy_narrative.tsv")
    relink(comp / "lineage", EM / "lineage")
    relink(comp / "lineage_verification", EM / "lineage_verification")
    relink(comp / "mutational_signatures", EM / "mutational_signatures")
    relink(comp / "qc_summary.tsv", EM / "qc_summary.tsv")

    print(f"\n[OK] Final-results tree refreshed under {F}")


if __name__ == "__main__":
    main()
