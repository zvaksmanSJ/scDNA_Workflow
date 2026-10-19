#!/usr/bin/env python3
"""Submit one per-cell LSF Control-FREEC job for CNA + LOH detection.

Control-FREEC (Boeva et al., Bioinformatics 2011/2012; GPL, bioconda) is a
widely used, peer-reviewed tool that natively computes B-allele-frequency
(BAF) profiles from a BAM plus known-SNP positions and calls LOH directly
(reported in its own _CNVs output's genotype/status columns) -- this is
the established method used instead of inventing a custom BAF/LOH
classifier. No matched-normal/bulk sample is required for WGS data.

Each cell is run tumor-only (no [control] section): ploidy=2, sex=XY
(cohort is male), window=50000 (standard WGS default). The known-SNP
positions for BAF calculation reuse the same 1000 Genomes GRCh38 panel
already downloaded for phasing (stripped to site-only VCFs by
build_snp_positions.sh), so no separate SNP database download is needed.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import textwrap
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam-dir", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--ref-fasta", required=True, type=Path)
    p.add_argument("--chrom-fasta-dir", required=True, type=Path)
    p.add_argument("--chr-len-file", type=Path, default=None,
                    help="Chromosome-length file (samtools .fai format) restricted to the "
                         "chromosomes actually present in --chrom-fasta-dir. FREEC requires a "
                         "FASTA entry for every contig listed here, so this must NOT be the full "
                         "genome .fai (which also lists ~170 decoy/alt/unplaced scaffolds plus "
                         "Y/MT that are not split out). Defaults to {ref_fasta}.fai if omitted.")
    p.add_argument("--snp-positions-vcf", required=True, type=Path)
    p.add_argument("--freec-bin", required=True, type=Path)
    p.add_argument("--sex", default="XY")
    p.add_argument("--ploidy", default="2",
                    help="Comma-separated candidate ploidies, e.g. '2,3,4'. Control-FREEC will "
                         "select whichever value best explains the observed CNA profile (official "
                         "feature, see 'ploidy' parameter in the FREEC manual) rather than forcing "
                         "a fixed assumption.")
    p.add_argument("--window", type=int, default=50000)
    p.add_argument("--min-coverage-per-position", type=int, default=5)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--logs-dir", required=True, type=Path)
    p.add_argument("--queue", default="standard")
    p.add_argument("--ncpus", type=int, default=4)
    p.add_argument("--mem-mb", type=int, default=32768)
    p.add_argument("--wall", default="48:00")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def write_config(sample: str, bam: Path, sample_out: Path, args) -> Path:
    config_path = sample_out / f"{sample}.freec_config.txt"
    config = textwrap.dedent(f"""
        [general]
        chrLenFile = {args.chr_len_file or (str(args.ref_fasta) + '.fai')}
        chrFiles = {args.chrom_fasta_dir}
        ploidy = {args.ploidy}
        sex = {args.sex}
        window = {args.window}
        breakPointType = 2
        outputDir = {sample_out}
        maxThreads = {args.ncpus}

        [sample]
        mateFile = {bam}
        inputFormat = BAM
        mateOrientation = FR

        [BAF]
        makePileup = {args.snp_positions_vcf}
        fastaFile = {args.ref_fasta}
        SNPfile = {args.snp_positions_vcf}
        minimalCoveragePerPosition = {args.min_coverage_per_position}
        shiftInQuality = 33
        """).strip() + "\n"
    config_path.write_text(config)
    return config_path


def build_job_script(sample: str, args) -> str:
    bam = args.bam_dir / sample / f"{sample}.dedup.bam"
    sample_out = args.output_dir / sample
    sample_out.mkdir(parents=True, exist_ok=True)
    config_path = write_config(sample, bam, sample_out, args)

    return textwrap.dedent(f"""
        #!/bin/bash
        #BSUB -J freec_{sample}
        #BSUB -q {args.queue}
        #BSUB -n {args.ncpus}
        #BSUB -M {args.mem_mb}
        #BSUB -R "span[hosts=1] rusage[mem={args.mem_mb}]"
        #BSUB -W {args.wall}
        #BSUB -oo {args.logs_dir}/freec_{sample}.%J.out
        #BSUB -eo {args.logs_dir}/freec_{sample}.%J.err

        set -euo pipefail
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate freec

        test -s {bam}
        {args.freec_bin} -conf {config_path}

        CNVS={sample_out}/{sample}.dedup.bam_CNVs
        test -s "$CNVS"
        echo "[OK] Control-FREEC completed for {sample}"
        echo "[OK] CNVs/LOH output: $CNVS"
        """).strip() + "\n"


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir = args.output_dir / "scripts"
    scripts_dir.mkdir(exist_ok=True)

    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    job_ids = []
    for sample in samples:
        script_path = scripts_dir / f"freec_{sample}.sh"
        script_path.write_text(build_job_script(sample, args))
        script_path.chmod(0o755)
        cmd = ["bsub", "-J", f"freec_{sample}", "-q", args.queue, "-n", str(args.ncpus),
               "-M", str(args.mem_mb), "-R", f"span[hosts=1] rusage[mem={args.mem_mb}]",
               "-W", args.wall, "-oo", str(args.logs_dir / f"freec_{sample}.%J.out"),
               "-eo", str(args.logs_dir / f"freec_{sample}.%J.err"), str(script_path)]
        if args.dry_run:
            print(" ".join(cmd))
            continue
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        jid = result.stdout.strip().split("<")[1].split(">")[0]
        job_ids.append(jid)
        print(f"[SUBMITTED] freec_{sample} -> {jid}")

    if job_ids:
        (args.output_dir / "freec_job_ids.txt").write_text("\n".join(job_ids) + "\n")
        print(f"[OK] Submitted {len(job_ids)} Control-FREEC jobs")


if __name__ == "__main__":
    main()
