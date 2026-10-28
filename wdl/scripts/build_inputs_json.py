#!/usr/bin/env python3
"""Build a complete main.wdl inputs.json from a simple FASTQ manifest.

Takes (1) a manifest TSV with one row per cell -- sample_id, r1_fastq,
r2_fastq -- and (2) a "base" JSON file already containing every fixed
infrastructure input (reference, containers, Beagle panel, extra_methods
script paths, etc. -- see wdl/example_inputs.json for a starting point you
customize ONCE per deployment), and merges the manifest's FASTQ list into
the base, producing a ready-to-run inputs.json. This avoids hand-editing
JSON arrays for every cell in a real cohort.

Manifest TSV format (tab-separated, no header required but recommended):
    sample_id	r1_fastq	r2_fastq
    cell1	/path/to/cell1_R1.fastq.gz	/path/to/cell1_R2.fastq.gz
    cell2	/path/to/cell2_R1.fastq.gz	/path/to/cell2_R2.fastq.gz
    ...

Usage:
    python3 build_inputs_json.py \\
      --manifest cells.tsv \\
      --base wdl/example_inputs.json \\
      --donor-id my_cohort \\
      --output inputs.json

Every non-FASTQ key already present in --base is copied through unchanged;
any "_comment*" keys (used for documentation in example_inputs.json) are
stripped. --donor-id overrides the base file's donor_id if given.

IMPORTANT -- autosome maps are REGENERATED, not merely copied from --base:
main.wdl's phasing stage hardcodes its own chromosome list (chr1-22; see
main.wdl's `Array[String] autosomes` declaration) and looks up each one in
the 3 Map[String, File] inputs below by that exact key. wdl/example_inputs.
json intentionally only ILLUSTRATES this with 2 example chromosomes (1, 2)
to keep the file readable -- copying that abbreviated example through
unchanged (which this script used to do, and a literal reading of this
project's own quick-start instructions would otherwise produce) silently
produces an inputs.json that runs cleanly through alignment/variant-calling
for every cell, then crashes with an undefined-Map-key error the first time
the phasing scatter reaches chromosome 3 -- after the expensive stages have
already run. To prevent that, this script ALWAYS regenerates the full
chr1-22 set for these 3 keys by globbing the real per-chromosome files
setup.sh actually downloads, rather than trusting whatever --base happens
to contain for them:
    phasing_ref_panel_by_chrom       <- refs/phasing_panel/CCDG_14151_B01_GRM_WGS_2020-08-05_chr<N>.filtered.shapeit2-duohmm-phased.vcf.gz
    phasing_ref_panel_tbi_by_chrom   <- same file + ".tbi"
    genetic_map_by_chrom             <- refs/genetic_maps/chr_in_chrom_field/plink.chrchr<N>.GRCh38.map
--wdl-runtime-dir controls where these are searched for (default:
"wdl_runtime", matching setup.sh's own default install location); use
--skip-autosome-regen to restore the old copy-through-unchanged behavior
(e.g. if you've installed these files somewhere setup.sh doesn't expect and
already hand-built a complete --base with every chromosome present).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

WORKFLOW_NAME = "ResolveOmeDnaWorkflow"
AUTOSOMES = [str(i) for i in range(1, 23)]  # matches main.wdl's own hardcoded `Array[String] autosomes`


def regenerate_autosome_maps(wdl_runtime_dir: Path) -> dict[str, dict[str, str]]:
    """Glob the real per-chromosome files setup.sh downloads (NOT whatever
    --base happens to contain) and build complete chr1-22 Map[String, File]
    values for phasing_ref_panel_by_chrom / phasing_ref_panel_tbi_by_chrom /
    genetic_map_by_chrom. Raises with a clear, actionable message listing
    every missing chromosome/file if setup.sh hasn't been run (or was run
    with SKIP_PHASING=1) rather than silently producing a partial map."""
    panel_dir = wdl_runtime_dir / "refs/phasing_panel"
    map_dir = wdl_runtime_dir / "refs/genetic_maps/chr_in_chrom_field"
    panel_by_chrom, panel_tbi_by_chrom, map_by_chrom = {}, {}, {}
    missing: list[str] = []
    for chrom in AUTOSOMES:
        panel = panel_dir / f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
        panel_tbi = panel.with_name(panel.name + ".tbi")
        genetic_map = map_dir / f"plink.chrchr{chrom}.GRCh38.map"
        for path, label in ((panel, f"chr{chrom} phasing panel VCF"), (panel_tbi, f"chr{chrom} phasing panel .tbi"),
                            (genetic_map, f"chr{chrom} genetic map")):
            if not path.is_file():
                missing.append(f"{label}: {path}")
        panel_by_chrom[chrom] = str(panel)
        panel_tbi_by_chrom[chrom] = str(panel_tbi)
        map_by_chrom[chrom] = str(genetic_map)
    if missing:
        joined = "\n  - ".join(missing)
        raise SystemExit(
            f"[ERROR] {len(missing)} file(s) needed for the full chr1-22 phasing/genetic-map maps are "
            f"missing under {wdl_runtime_dir}:\n  - {joined}\n"
            "Run setup/setup.sh (without SKIP_PHASING=1) first, or pass --skip-autosome-regen to fall "
            "back to copying --base's own (likely incomplete) maps through unchanged, or --wdl-runtime-dir "
            "to point this at wherever setup.sh actually installed to."
        )
    return {
        f"{WORKFLOW_NAME}.phasing_ref_panel_by_chrom": panel_by_chrom,
        f"{WORKFLOW_NAME}.phasing_ref_panel_tbi_by_chrom": panel_tbi_by_chrom,
        f"{WORKFLOW_NAME}.genetic_map_by_chrom": map_by_chrom,
    }


def read_manifest(path: Path) -> tuple[list[str], list[str], list[str]]:
    sample_ids, r1s, r2s = [], [], []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 3:
                raise SystemExit(f"[ERROR] manifest row has fewer than 3 columns: {row!r}")
            sample_id, r1, r2 = row[0].strip(), row[1].strip(), row[2].strip()
            if sample_id.lower() in ("sample_id", "sample", "id") and r1.lower() in ("r1_fastq", "r1", "fastq_r1"):
                continue  # skip an optional header row
            for p, label in ((r1, "R1"), (r2, "R2")):
                if not Path(p).is_file():
                    print(f"[WARN] {label} for {sample_id} does not exist (yet?): {p}", file=sys.stderr)
            sample_ids.append(sample_id)
            r1s.append(r1)
            r2s.append(r2)
    if not sample_ids:
        raise SystemExit(f"[ERROR] no cells found in manifest: {path}")
    dupes = {s for s in sample_ids if sample_ids.count(s) > 1}
    if dupes:
        raise SystemExit(f"[ERROR] duplicate sample_id(s) in manifest: {sorted(dupes)}")
    return sample_ids, r1s, r2s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, type=Path, help="TSV: sample_id, r1_fastq, r2_fastq")
    ap.add_argument("--base", required=True, type=Path, help="base inputs.json with every fixed/infra input filled in")
    ap.add_argument("--donor-id", default=None, help="override the base file's donor_id")
    ap.add_argument("--wdl-runtime-dir", type=Path, default=Path("wdl_runtime"),
                     help="Where setup.sh installed phasing resources (default: wdl_runtime, matching "
                          "setup.sh's own default BASE_DIR) -- used to regenerate the full chr1-22 "
                          "phasing/genetic-map maps; see --skip-autosome-regen to disable.")
    ap.add_argument("--skip-autosome-regen", action="store_true",
                     help="Copy --base's own phasing_ref_panel_by_chrom/phasing_ref_panel_tbi_by_chrom/"
                          "genetic_map_by_chrom through unchanged instead of regenerating them from "
                          "--wdl-runtime-dir. Only safe if --base already has a complete, correct chr1-22 "
                          "set for all three -- the bundled wdl/example_inputs.json deliberately does NOT "
                          "(see this script's own module docstring for why).")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    base = json.loads(args.base.read_text())
    base = {k: v for k, v in base.items() if not k.startswith("_comment")}

    sample_ids, r1s, r2s = read_manifest(args.manifest)
    base[f"{WORKFLOW_NAME}.sample_ids"] = sample_ids
    base[f"{WORKFLOW_NAME}.r1_fastqs"] = r1s
    base[f"{WORKFLOW_NAME}.r2_fastqs"] = r2s
    if args.donor_id:
        base[f"{WORKFLOW_NAME}.donor_id"] = args.donor_id

    if args.skip_autosome_regen:
        print("[WARN] --skip-autosome-regen set: phasing_ref_panel_by_chrom/"
              "phasing_ref_panel_tbi_by_chrom/genetic_map_by_chrom are being copied through from --base "
              "UNCHANGED. If --base is wdl/example_inputs.json (or derived from it without further "
              "editing), these only cover chr1-2 and the resulting inputs.json WILL fail partway through "
              "a real run, at the phasing stage's chr3 scatter iteration.", file=sys.stderr)
    else:
        autosome_maps = regenerate_autosome_maps(args.wdl_runtime_dir)
        base.update(autosome_maps)
        print(f"[OK] Regenerated complete chr1-22 phasing/genetic-map maps from {args.wdl_runtime_dir}")

    args.output.write_text(json.dumps(base, indent=2) + "\n")
    print(f"[DONE] wrote {args.output} with {len(sample_ids)} cell(s): {', '.join(sample_ids)}")


if __name__ == "__main__":
    main()
