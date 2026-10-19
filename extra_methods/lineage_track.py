#!/usr/bin/env python3
"""Lightweight, same-cohort cell-lineage reconstruction.

This reconstructs a phylogenetic tree relating our 15 single cells directly
from their own per-cell merged DeepVariant calls, using the same Sequoia
tool (via the existing basej-somatic Sequoia container) that the full
official BaseJumper somatic+lineage Nextflow pipeline will eventually use
several stages downstream. It does NOT wait for that full pipeline (variant
cohort gating -> BaseJ somatic Nextflow -> integrated matrices -> basej-
lineage Nextflow), which is still many stages and likely days away on this
shared cluster; instead it builds NV (variant-supporting reads) and NR
(total reads) matrices directly from each cell's own merged annotated VCF,
matching the (variant x cell) matrix approach described in the reference
paper's Methods, with one documented simplification:

  The reference paper retrieves NV/NR for every cell at every site via a
  full pileup across all BAMs (so even a cell with zero DeepVariant
  evidence at a site still gets an accurate read-depth count). We instead
  use each cell's own DeepVariant FORMAT/AD field directly: a site absent
  from a given cell's own VCF (i.e. DeepVariant proposed no candidate
  allele there in that cell) is recorded as NR=0, NV=0 for that cell. This
  is a reasonable approximation for SNV sites (DeepVariant's non-gVCF
  "calling" mode already emits RefCall records at most positions where
  any cell showed even weak candidate evidence), but is a lighter-weight
  method than a true uniform pileup and is documented here explicitly.

Only SNVs are used for the tree topology (matching Sequoia's own
`only_snvs` default), restricted to PASS-called sites in at least one
cell, unioned across the cohort.
"""
from __future__ import annotations

import argparse
import csv
import random
import subprocess
from pathlib import Path


def biallelic_pass_snv_sites(vcf_path: Path) -> set[tuple]:
    result = subprocess.run(
        ["bcftools", "view", "-f", "PASS", "-m2", "-M2", "-v", "snps", "-H", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    sites = set()
    for line in result.stdout.splitlines():
        fields = line.split("\t", 5)
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        sites.add((chrom, pos, ref, alt))
    return sites


def read_sample_ad(vcf_path: Path, sites: set[tuple]) -> dict[tuple, tuple[int, int]]:
    """Return {site: (NR, NV)} for every site present in this sample's own
    VCF (any FILTER), limited to biallelic SNV sites matching `sites`."""
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    out = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 5:
            continue
        chrom, pos, ref, alt, ad = parts
        if "," in alt:
            continue  # multiallelic in this cell's own call; skip for matrix simplicity
        key = (chrom, pos, ref, alt)
        if key not in sites:
            continue
        ad_parts = ad.split(",")
        if len(ad_parts) != 2 or ad_parts[0] in (".", "") or ad_parts[1] in (".", ""):
            continue
        try:
            ref_depth, alt_depth = int(ad_parts[0]), int(ad_parts[1])
        except ValueError:
            continue
        out[key] = (ref_depth + alt_depth, alt_depth)
    return out


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variants-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--final-vcf-suffix", default=".deepvariant.raw.vep.snpeff.vcf.gz")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--donor-id", default="SY5Y_cohort")
    p.add_argument("--gender", default="male", choices=["male", "female"])
    p.add_argument("--mpboot-bootstrap", type=int, default=1000)
    p.add_argument("--sequoia-sif", type=Path, default=None,
                    help="Apptainer/Singularity SIF bundling Sequoia + MPBoot (used for this "
                         "session's real cohort run). Mutually exclusive with --sequoia-script.")
    p.add_argument("--sequoia-script", type=Path, default=None,
                    help="Direct path to Sequoia's build_phylogeny.R (standalone GPLv3 script, "
                         "github.com/TimCoorens/Sequoia), for a plain conda-env Rscript "
                         "invocation with no container -- used by the general-purpose WDL "
                         "workflow. Mutually exclusive with --sequoia-sif.")
    p.add_argument("--mpboot-path", default="/usr/local/bin/scripts/")
    p.add_argument("--max-shared-sites", type=int, default=20000,
                    help="Cap on non-singleton informative sites (all singletons always kept)")
    p.add_argument("--random-seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    vcf_paths = {}
    union_sites: set[tuple] = set()
    for sample in samples:
        vcf_path = args.variants_dir / sample / f"{sample}{args.final_vcf_suffix}"
        if not vcf_path.is_file():
            print(f"[WARN] Missing VCF for {sample}, excluding from lineage")
            continue
        vcf_paths[sample] = vcf_path
        sample_sites = biallelic_pass_snv_sites(vcf_path)
        union_sites |= sample_sites
        print(f"[INFO] {sample}: {len(sample_sites)} PASS biallelic SNV sites")

    print(f"[INFO] Union of PASS SNV sites across cohort: {len(union_sites)}")
    if not union_sites:
        raise SystemExit("No PASS SNV sites found across cohort; cannot build lineage")

    per_sample_ad = {}
    for sample, vcf_path in vcf_paths.items():
        per_sample_ad[sample] = read_sample_ad(vcf_path, union_sites)
        print(f"[INFO] {sample}: AD retrieved for "
              f"{len(per_sample_ad[sample])}/{len(union_sites)} union sites")

    used_samples = sorted(per_sample_ad)

    # Pre-filter to lineage-informative (non-invariant, non-near-universal)
    # sites. A site that discretizes to "present" (VAF >= present threshold)
    # in literally every cell -- or "absent" in every cell -- carries zero
    # information for parsimony tree topology (standard phylogenetic
    # practice: constant/invariant characters don't inform tree search).
    # Because this is a clonal cell line, true germline variants are shared
    # by essentially all 15 cells but, empirically, PTA allele-dropout noise
    # produces a smooth recurrence gradient rather than a clean bimodal
    # split (observed directly on this cohort: site counts increase
    # monotonically from singletons up through near-universal rather than
    # separating cleanly) -- so a near-universal presence threshold (90%,
    # matching the artifact-panel script) is used to exclude the
    # germline+dropout-noise bulk, and the remaining informative set is
    # capped via reproducible random subsampling (all true singletons are
    # always kept, since those directly define terminal/private branches)
    # to keep the alignment length tractable for MPBoot's 1000-replicate
    # bootstrap search, which scales with alignment length as well as taxa
    # count. This subsampling is an explicit, documented simplification
    # for this lightweight track, not a claim of exhaustive mutation
    # recovery.
    vaf_absent, vaf_present = 0.10, 0.30
    universal_threshold = max(1, round(len(used_samples) * 0.90))
    informative_sites = []
    for chrom, pos, ref, alt in union_sites:
        xy = chrom.upper() in ("X", "Y", "CHRX", "CHRY")
        present_t = vaf_present * (2 if (xy and args.gender == "male") else 1)
        present_count = 0
        for sample in used_samples:
            nr, nv = per_sample_ad[sample].get((chrom, pos, ref, alt), (0, 0))
            vaf = (nv / nr) if nr > 0 else 0.0
            if vaf >= present_t:
                present_count += 1
        if 0 < present_count < universal_threshold:
            informative_sites.append((chrom, pos, ref, alt, present_count))

    print(f"[INFO] Lineage-informative (non-invariant, <{universal_threshold}/"
          f"{len(used_samples)} cells) sites: {len(informative_sites)} "
          f"of {len(union_sites)} union sites")
    if not informative_sites:
        raise SystemExit("No lineage-informative (variable) sites found across cohort")

    singletons = [s for s in informative_sites if s[4] == 1]
    shared = [s for s in informative_sites if s[4] > 1]
    kept = list(singletons)
    budget = max(0, args.max_shared_sites - len(singletons))
    if len(shared) > budget:
        rng = random.Random(args.random_seed)
        kept += rng.sample(shared, budget)
        print(f"[INFO] Subsampled shared (non-singleton) informative sites: "
              f"{budget} of {len(shared)} kept (seed={args.random_seed})")
    else:
        kept += shared
    print(f"[INFO] Final lineage matrix sites: {len(kept)} "
          f"({len(singletons)} singleton + {len(kept) - len(singletons)} shared)")
    sorted_sites = sorted((c, p, r, a) for c, p, r, a, _ in kept)

    nr_path = args.output_dir / f"{args.donor_id}_NR.tsv"
    nv_path = args.output_dir / f"{args.donor_id}_NV.tsv"
    with nr_path.open("w", newline="") as nr_fh, nv_path.open("w", newline="") as nv_fh:
        nr_writer = csv.writer(nr_fh, delimiter="\t")
        nv_writer = csv.writer(nv_fh, delimiter="\t")
        nr_writer.writerow(["Mut"] + used_samples)
        nv_writer.writerow(["Mut"] + used_samples)
        for chrom, pos, ref, alt in sorted_sites:
            mut_id = f"{chrom}_{pos}_{ref}_{alt}"
            nr_row, nv_row = [mut_id], [mut_id]
            for sample in used_samples:
                nr, nv = per_sample_ad[sample].get((chrom, pos, ref, alt), (0, 0))
                nr_row.append(nr)
                nv_row.append(nv)
            nr_writer.writerow(nr_row)
            nv_writer.writerow(nv_row)

    print(f"[OK] NR matrix -> {nr_path} ({len(sorted_sites)} variants x {len(used_samples)} cells)")
    print(f"[OK] NV matrix -> {nv_path}")

    if args.sequoia_sif is None and args.sequoia_script is None:
        raise SystemExit("Exactly one of --sequoia-sif or --sequoia-script is required")
    if args.sequoia_sif is not None and args.sequoia_script is not None:
        raise SystemExit("--sequoia-sif and --sequoia-script are mutually exclusive")

    if args.sequoia_sif is not None:
        cmd = [
            "apptainer", "exec", "--cleanenv", "--no-mount", "/programs",
            "-B", "/research_jude:/research_jude",
            str(args.sequoia_sif),
            "Rscript", "/usr/local/bin/scripts/rscript_sequoia_build_phylogeny_only.R",
            "--donor_id", args.donor_id,
            "--input_nv", str(nv_path),
            "--input_nr", str(nr_path),
            "--output_dir", str(args.output_dir) + "/",
            "--gender", args.gender,
            "--mpboot_path", args.mpboot_path,
            "--mpboot_bootstrap", str(args.mpboot_bootstrap),
        ]
    else:
        cmd = [
            "Rscript", str(args.sequoia_script),
            "--donor_id", args.donor_id,
            "--input_nv", str(nv_path),
            "--input_nr", str(nr_path),
            "--output_dir", str(args.output_dir) + "/",
            "--gender", args.gender,
            "--mpboot_path", args.mpboot_path,
            "--mpboot_bootstrap", str(args.mpboot_bootstrap),
        ]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("[OK] Sequoia lineage reconstruction complete")


if __name__ == "__main__":
    main()
