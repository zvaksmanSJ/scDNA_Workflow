#!/usr/bin/env python3
"""Orchestrate lineage verification:
  1. Re-run lineage_track.py with an independent random seed (different
     subsample of lineage-informative sites), producing a second,
     independently-built tree.
  2. Invoke verify_lineage.R (ape::read.tree + ape::dist.topo, an
     established phylogenetics package already used by the Sequoia step)
     to report bootstrap support per branch and the Robinson-Foulds
     topology distance between the primary and resampled trees.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def find_tree_file(directory: Path, donor_id: str) -> Path | None:
    candidate = directory / f"{donor_id}_snv_tree_with_branch_length_selectedscheme.tree"
    if candidate.is_file():
        return candidate
    matches = list(directory.glob(f"{donor_id}_*tree_with_branch_length_selectedscheme.tree"))
    return matches[0] if matches else None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lineage-track-script", required=True, type=Path)
    p.add_argument("--variants-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--final-vcf-suffix", default=".deepvariant.raw.vep.snpeff.vcf.gz",
                    help="Must match lineage_track.py's own default/the suffix used for the "
                         "primary tree build, since this VCF lives directly in each cell's "
                         "variants directory (unlike the un-annotated support/merged/ VCF).")
    p.add_argument("--primary-lineage-dir", required=True, type=Path)
    p.add_argument("--donor-id", default="SY5Y_cohort")
    p.add_argument("--gender", default="male")
    p.add_argument("--mpboot-bootstrap", type=int, default=1000)
    p.add_argument("--sequoia-sif", type=Path, default=None,
                    help="Apptainer/Singularity SIF bundling Sequoia + MPBoot + R/ape. Mutually "
                         "exclusive with --sequoia-script/--r-env.")
    p.add_argument("--sequoia-script", type=Path, default=None,
                    help="Direct path to Sequoia's build_phylogeny.R (standalone GPLv3, no "
                         "container). Mutually exclusive with --sequoia-sif.")
    p.add_argument("--mpboot-path", default="/usr/local/bin/scripts/")
    p.add_argument("--max-shared-sites", type=int, default=20000)
    p.add_argument("--rerun-seed", type=int, default=123)
    p.add_argument("--bootstrap-threshold", type=float, default=70.0)
    p.add_argument("--verify-r-script", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--skip-rerun", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.sequoia_sif is None and args.sequoia_script is None:
        raise SystemExit("Exactly one of --sequoia-sif or --sequoia-script is required")
    if args.sequoia_sif is not None and args.sequoia_script is not None:
        raise SystemExit("--sequoia-sif and --sequoia-script are mutually exclusive")

    primary_tree = find_tree_file(args.primary_lineage_dir, args.donor_id)
    if primary_tree is None:
        raise SystemExit(f"Primary tree not found in {args.primary_lineage_dir}")
    print(f"[INFO] Primary tree: {primary_tree}")

    resample_tree = None
    if not args.skip_rerun:
        resample_dir = args.output_dir / "resample_rerun"
        resample_donor = args.donor_id + "_resample"
        cmd = [sys.executable, str(args.lineage_track_script),
               "--variants-dir", str(args.variants_dir),
               "--manifest", str(args.manifest),
               "--final-vcf-suffix", args.final_vcf_suffix,
               "--output-dir", str(resample_dir),
               "--donor-id", resample_donor,
               "--gender", args.gender,
               "--mpboot-bootstrap", str(args.mpboot_bootstrap)]
        if args.sequoia_sif is not None:
            cmd += ["--sequoia-sif", str(args.sequoia_sif)]
        else:
            cmd += ["--sequoia-script", str(args.sequoia_script)]
        cmd += ["--mpboot-path", args.mpboot_path,
               "--max-shared-sites", str(args.max_shared_sites),
               "--random-seed", str(args.rerun_seed)]
        print("[CMD]", " ".join(cmd))
        subprocess.run(cmd, check=True)
        resample_tree = find_tree_file(resample_dir, resample_donor)
        if resample_tree is None:
            print("[WARN] Resample tree not produced; skipping topology comparison")

    bootstrap_out = args.output_dir / "bootstrap_support_report.tsv"
    summary_out = args.output_dir / "verification_summary.json"
    r_args = [str(primary_tree), str(args.bootstrap_threshold),
              str(bootstrap_out), str(summary_out),
              str(resample_tree) if resample_tree else ""]
    if args.sequoia_sif is not None:
        r_cmd = ["apptainer", "exec", "--cleanenv", "--no-mount", "/programs",
                  "-B", "/research_jude:/research_jude", str(args.sequoia_sif),
                  "Rscript", str(args.verify_r_script)] + r_args
    else:
        # verify_lineage.R only needs the `ape` package, already present in
        # any R environment that can run Sequoia itself (MPBoot's R wrapper
        # depends on ape too), so a plain Rscript call suffices here.
        r_cmd = ["Rscript", str(args.verify_r_script)] + r_args
    print("[CMD]", " ".join(r_cmd))
    subprocess.run(r_cmd, check=True)
    print("[OK] Lineage verification complete")


if __name__ == "__main__":
    main()
