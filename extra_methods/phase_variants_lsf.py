#!/usr/bin/env python3
"""Submit one per-cell LSF phasing job using Beagle 5.5 (GPLv3, fully free)
and the official 1000 Genomes Project GRCh38 high-coverage phased callset
(public/open IGSR resource, no license restriction) as the reference panel.

This is a new addition relative to both BaseJumper and our prior workflow:
neither currently performs statistical phasing. Phasing is run per
chromosome (autosomes 1-22 only -- see CHROMS note below) against each
chromosome's own reference panel VCF and genetic map, then concatenated
back into one phased VCF per cell.

The 1000 Genomes reference panel and Beagle's genetic maps use "chr"-
prefixed contig names (chr1, chr2, ...); our own pipeline's reference and
VCFs use bare contig names (1, 2, ...). Contigs are renamed to "chr"-prefixed
before phasing and renamed back afterwards so the phased output matches our
existing contig-naming convention.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import textwrap
from pathlib import Path

# Autosomes only. chrX is intentionally excluded: this cohort is male, so
# chrX is hemizygous outside the small pseudoautosomal regions. The NYGC
# Eagle2-phased chrX reference panel VCF stores some reference samples as
# diploid (PAR) and others as haploid (non-PAR) within the same file, which
# Beagle 5.5's parser rejects outright (IllegalArgumentException: sample has
# an inconsistent number of alleles). This was confirmed to fail identically
# for all 15 cells. Phased het-SNP + CNV-segment LOH (our actual downstream
# use of phasing) is an autosomal-diploid concept that doesn't apply the same
# way to a hemizygous chromosome, so excluding chrX here is a correct
# scientific choice, not just a technical workaround.
CHROMS = [str(i) for i in range(1, 23)]


def panel_filename(c: str) -> str:
    # chrX in this NYGC release is phased with Eagle2 (not SHAPEIT2-duohmm,
    # which needs pedigree/duo-hmm info unavailable for the non-autosomal
    # chromosome) and is published under a different filename pattern.
    if c == "X":
        return "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
    return f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{c}.filtered.shapeit2-duohmm-phased.vcf.gz"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variants-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--raw-vcf-subdir", default="support/merged",
                    help="Subdirectory under each cell's variants dir holding the merged pre-annotation VCF")
    p.add_argument("--raw-vcf-suffix", default=".deepvariant.raw.vcf.gz",
                    help="Un-annotated merged DeepVariant VCF (GT-only input to Beagle)")
    p.add_argument("--beagle-jar", required=True, type=Path)
    p.add_argument("--ref-panel-dir", required=True, type=Path)
    p.add_argument("--genetic-map-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--logs-dir", required=True, type=Path)
    p.add_argument("--queue", default="standard")
    p.add_argument("--ncpus", type=int, default=4)
    p.add_argument("--mem-mb", type=int, default=49152)
    p.add_argument("--wall", default="145:00")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def chrom_rename_file(output_dir: Path, to_chr: bool) -> Path:
    name = "rename_to_chr.tsv" if to_chr else "rename_from_chr.tsv"
    path = output_dir / name
    with path.open("w") as fh:
        for c in CHROMS:
            if to_chr:
                fh.write(f"{c}\tchr{c}\n")
            else:
                fh.write(f"chr{c}\t{c}\n")
    return path


def build_job_script(sample: str, args) -> str:
    sample_out = args.output_dir / sample
    chrom_lines = []
    for c in CHROMS:
        ref_vcf = args.ref_panel_dir / panel_filename(c)
        gmap = args.genetic_map_dir / "chr_in_chrom_field" / f"plink.chrchr{c}.GRCh38.map"
        chrom_lines.append(textwrap.dedent(f"""
            echo "[INFO] Phasing chr{c}"
            bcftools view -r {c} "$RAW_VCF" -Oz -o "$WORK/input_chr{c}.vcf.gz"
            bcftools annotate --rename-chrs "$RENAME_TO_CHR" -Oz -o "$WORK/input_renamed_chr{c}.vcf.gz" "$WORK/input_chr{c}.vcf.gz"
            bcftools index -f "$WORK/input_renamed_chr{c}.vcf.gz"
            java -Xmx{max(8, args.mem_mb // 1024 - 4)}g -jar {args.beagle_jar} \\
              gt="$WORK/input_renamed_chr{c}.vcf.gz" \\
              ref={ref_vcf} \\
              map={gmap} \\
              chrom=chr{c} \\
              out="$WORK/phased_chr{c}" \\
              nthreads={args.ncpus} \\
              impute=false
            test -s "$WORK/phased_chr{c}.vcf.gz"
            bcftools index -f "$WORK/phased_chr{c}.vcf.gz"
            """).strip())
    chrom_block = "\n\n".join(chrom_lines)

    return textwrap.dedent(f"""
        #!/bin/bash
        #BSUB -J phase_{sample}
        #BSUB -q {args.queue}
        #BSUB -n {args.ncpus}
        #BSUB -M {args.mem_mb}
        #BSUB -R "span[hosts=1] rusage[mem={args.mem_mb}]"
        #BSUB -W {args.wall}
        #BSUB -oo {args.logs_dir}/phase_{sample}.%J.out
        #BSUB -eo {args.logs_dir}/phase_{sample}.%J.err

        set -euo pipefail
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate resolveome-align

        RAW_VCF={args.variants_dir}/{sample}/{args.raw_vcf_subdir}/{sample}{args.raw_vcf_suffix}
        test -s "$RAW_VCF"
        RENAME_TO_CHR={chrom_rename_file(args.output_dir, True)}
        RENAME_FROM_CHR={chrom_rename_file(args.output_dir, False)}
        WORK=$(mktemp -d {args.output_dir}/{sample}.phasing_work.XXXXXX)
        mkdir -p {sample_out}
        trap 'rm -rf "$WORK"' EXIT

        {chrom_block}

        bcftools concat -Oz -o "$WORK/phased_allchrom.vcf.gz" \\
          {" ".join(f'"$WORK/phased_chr{c}.vcf.gz"' for c in CHROMS)}
        bcftools annotate --rename-chrs "$RENAME_FROM_CHR" \\
          -Oz -o {sample_out}/{sample}.phased.vcf.gz "$WORK/phased_allchrom.vcf.gz"
        bcftools index -f {sample_out}/{sample}.phased.vcf.gz

        echo "[OK] Phasing completed for {sample}"
        echo "[OK] Phased VCF: {sample_out}/{sample}.phased.vcf.gz"
        """).strip() + "\n"


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir = args.output_dir / "scripts"
    scripts_dir.mkdir(exist_ok=True)

    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    for c in CHROMS:
        ref_vcf = args.ref_panel_dir / panel_filename(c)
        if not ref_vcf.is_file():
            raise SystemExit(f"Missing reference panel file for chr{c}: {ref_vcf}")

    job_ids = []
    for sample in samples:
        script_path = scripts_dir / f"phase_{sample}.sh"
        script_path.write_text(build_job_script(sample, args))
        script_path.chmod(0o755)
        cmd = ["bsub", "-J", f"phase_{sample}", "-q", args.queue, "-n", str(args.ncpus),
               "-M", str(args.mem_mb), "-R", f"span[hosts=1] rusage[mem={args.mem_mb}]",
               "-W", args.wall, "-oo", str(args.logs_dir / f"phase_{sample}.%J.out"),
               "-eo", str(args.logs_dir / f"phase_{sample}.%J.err"), str(script_path)]
        if args.dry_run:
            print(" ".join(cmd))
            continue
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        jid = result.stdout.strip().split("<")[1].split(">")[0]
        job_ids.append(jid)
        print(f"[SUBMITTED] phase_{sample} -> {jid}")

    if job_ids:
        (args.output_dir / "phase_job_ids.txt").write_text("\n".join(job_ids) + "\n")
        print(f"[OK] Submitted {len(job_ids)} phasing jobs")


if __name__ == "__main__":
    main()
