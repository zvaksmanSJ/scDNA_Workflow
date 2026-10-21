#!/usr/bin/env python3
"""
Submit Step3 jobs to LSF for:

  1. DeepVariant SNV calling by genomic chunk
  2. VEP annotation of each DeepVariant chunk
  3. Merge VEP chunks per sample/cell
  4. SnpEff annotation once on the merged VEP VCF
  5. CNVkit CNV calling
  6. Manta SV calling
  7. Optional MultiQC

Annotation order:

    DeepVariant chunk
        -> VEP chunk
            -> VEP merge per sample
                -> SnpEff on merged VEP VCF

Intermediate files are placed under:

    <outdir>/variants/<sample>/support/

Only final SnpEff outputs are placed directly under:

    <outdir>/variants/<sample>/

The VCF FILTER field from DeepVariant is preserved. If
--snv-filter-pass-only is supplied, each DeepVariant chunk is filtered to
PASS and "." records before being sent to VEP.

Manifest:

    <outdir>/logs_step3/<--job-manifest-name, default step3_job_manifest.tsv>
"""

import argparse
import csv
import glob
import os
import re
import shlex
import subprocess
import textwrap
from datetime import datetime
from typing import Dict, List, Optional, Tuple


# ============================================================================
# Defaults
# ============================================================================

DEFAULT_OUTDIR = (
    "/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/"
    "Zalman/Ultima_test/ResolveOME_illumina/results/DNA_Ilumina_output_1"
)

DEFAULT_REF = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/refs/genome/"
    "GRCh38.primary.fa"
)

DEFAULT_ENV = "resolveome-align"
DEFAULT_QUEUE = "standard"

DEFAULT_DV_SIF = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/containers/"
    "deepvariant.sif"
)

DEFAULT_DV_MODEL_DIR = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/pipeline_resources/genomes/Homo_sapiens/NCBI/"
    "GRCh38/Annotation/DeepVariant/Illumina/"
    "bioskryb-af-20241102/model"
)

DEFAULT_DV_CKPT = "sc_af.model.ckpt"

DEFAULT_DV_POP_VCFS = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/pipeline_resources/genomes/Homo_sapiens/NCBI/"
    "GRCh38/Annotation/DeepVariant/Illumina/"
    "bioskryb-af-20241102/population_vcfs_nochr"
)

DEFAULT_MANTA_SIF = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/containers/manta.sif"
)

DEFAULT_CNV_GENE_BED = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/refs/annotation/"
    "genes.bed"
)

DEFAULT_CNV_REFERENCE_CNN = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/step3_bundle/"
    "resources/cnvkit/reference.cnn"
)

DEFAULT_SNPEFF_DB = "GRCh38.105"

DEFAULT_VEP_SIF = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/workflow_DNA/resolveome_align/containers/"
    "vep.r112.sif"
)

DEFAULT_VEP_CACHE = (
    "/research_jude/rgs01_jude/groups/SJCloud/common/Zalman/Ultima_test/"
    "ResolveOME_illumina/pipeline_resources/genomes/Homo_sapiens/NCBI/"
    "GRCh38/Annotation/VEP"
)

DEFAULT_VEP_FASTA = DEFAULT_REF

CANONICAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


# ============================================================================
# Shell and validation helpers
# ============================================================================

def run(cmd: str) -> str:
    """Run a command and raise with stdout/stderr if it fails."""
    proc = subprocess.run(
        cmd,
        shell=True,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"\nCOMMAND FAILED:\n{cmd}\n"
            f"\nSTDOUT:\n{proc.stdout}\n"
            f"\nSTDERR:\n{proc.stderr}\n"
        )
    return proc.stdout.strip()


def run_no_fail(cmd: str) -> Tuple[int, str, str]:
    proc = subprocess.run(
        cmd,
        shell=True,
        text=True,
        capture_output=True,
    )
    return proc.returncode, proc.stdout.strip(), proc.stderr.strip()


def parse_jobid(text: str) -> str:
    match = re.search(r"Job\s+<(\d+)>", text)
    if not match:
        raise RuntimeError(f"Could not parse LSF job ID from: {text}")
    return match.group(1)


def csvset(value: str) -> set:
    return {x.strip().lower() for x in value.split(",") if x.strip()}


def ensure_file(path: str, message: str):
    if not path or not os.path.isfile(path):
        raise SystemExit(message)


def ensure_dir(path: str, message: str):
    if not path or not os.path.isdir(path):
        raise SystemExit(message)


def file_nonempty(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def q(path: str) -> str:
    return shlex.quote(path)


def dep_and(job_ids: List[str]) -> Optional[str]:
    if not job_ids:
        return None
    return " && ".join(f"done({job_id})" for job_id in job_ids)


# ============================================================================
# Chromosome and chunk helpers
# ============================================================================

def normalize_chrom_token(token: str) -> str:
    value = token.strip().upper()
    if value.startswith("CHR"):
        value = value[3:]
    if value == "M":
        value = "MT"
    return value


def chrom_matches_fai(contig_name: str, wanted: str) -> bool:
    contig = contig_name.upper()
    target = wanted.upper()

    if contig == target:
        return True

    if contig.startswith("CHR") and contig[3:] == target:
        return True

    if target == "MT" and contig in {"M", "CHRM", "CHRMT"}:
        return True

    return False


def parse_chroms(value: str) -> List[str]:
    if value.strip().lower() == "all":
        # MT is intentionally excluded from the automatic genome-wide sweep:
        # the DeepVariant allele-frequency population VCF bundle has no MT/chrM
        # entry (only chr1-22,X,Y), and requesting the AF channel for a contig
        # with no population reference produced non-deterministic channel
        # counts in make_examples/call_variants (observed failures with varying,
        # inconsistent expected-vs-actual tensor shapes). MT remains a valid
        # explicit --chroms value for anyone who later adds a proper MT
        # population reference and wants to call it deliberately.
        return [c for c in CANONICAL_CHROMS if c != "MT"]

    chromosomes = [
        normalize_chrom_token(item)
        for item in value.split(",")
        if item.strip()
    ]

    invalid = [
        chromosome
        for chromosome in chromosomes
        if chromosome not in CANONICAL_CHROMS
    ]

    if invalid:
        raise SystemExit(
            f"[ERROR] Invalid chromosomes: {invalid}. "
            "Allowed values are 1-22,X,Y,MT or all."
        )

    return chromosomes


def read_fai_contigs(
    fai_path: str,
    chroms_arg: str,
) -> List[Tuple[str, int, str]]:
    """Return (FAI contig name, length, canonical chromosome)."""
    all_contigs = []

    with open(fai_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                all_contigs.append((fields[0], int(fields[1])))

    requested = parse_chroms(chroms_arg)
    selected = []

    for chromosome in requested:
        match = None

        for contig_name, length in all_contigs:
            if chrom_matches_fai(contig_name, chromosome):
                match = (contig_name, length, chromosome)
                break

        if match is None:
            raise SystemExit(
                f"[ERROR] Chromosome {chromosome} was not found in {fai_path}"
            )

        selected.append(match)

    return selected


def make_chunks(
    contigs: List[Tuple[str, int, str]],
    chunk_size: int,
):
    """
    Return:
      fai_name, canonical_chromosome, start, end, region, chunk_id
    """
    chunks = []

    for fai_name, length, canonical in contigs:
        start = 1

        while start <= length:
            end = min(length, start + chunk_size - 1)
            region = f"{fai_name}:{start}-{end}"
            chunk_id = f"{fai_name}_{start}_{end}"

            chunks.append(
                (
                    fai_name,
                    canonical,
                    start,
                    end,
                    region,
                    chunk_id,
                )
            )

            start = end + 1

    return chunks


def population_vcf_for_chrom(
    population_path: str,
    canonical_chrom: str,
) -> str:
    if not population_path:
        return ""

    if os.path.isfile(population_path):
        return population_path

    if not os.path.isdir(population_path):
        return ""

    file_chrom = "M" if canonical_chrom == "MT" else canonical_chrom

    return os.path.join(
        population_path,
        f"cohort-chr{file_chrom}.release_missing2ref.no_calls.vcf.gz",
    )


# ============================================================================
# Manifest and LSF helpers
# ============================================================================

MANIFEST_COLUMNS = [
    "submission_ts",
    "sample",
    "stage",
    "job_name",
    "job_id",
    "queue",
    "requested_ncpus",
    "requested_mem_mb",
    "requested_wall",
    "out_log",
    "err_log",
    "script_path",
    "depends_on",
]


def init_manifest(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()


def append_manifest(path: str, row: Dict[str, str]):
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_COLUMNS,
            delimiter="\t",
        )
        writer.writerow(row)


def write_job_script(
    script_dir: str,
    job_name: str,
    body: str,
) -> str:
    os.makedirs(script_dir, exist_ok=True)

    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", job_name)
    script_path = os.path.join(script_dir, f"{safe_name}.sh")

    with open(script_path, "w") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n")
        handle.write(body.strip())
        handle.write("\n")

    os.chmod(script_path, 0o755)
    return script_path


def submit_job(
    *,
    manifest_path: str,
    sample: str,
    stage: str,
    job_name: str,
    queue: str,
    ncpus: int,
    memory_mb: int,
    wall: str,
    out_log: str,
    err_log: str,
    body: str,
    script_dir: str,
    dependency: Optional[str] = None,
) -> str:
    script_path = write_job_script(script_dir, job_name, body)

    dependency_arg = ""
    if dependency:
        # Linux execve() enforces a per-argument length limit (MAX_ARG_STRLEN,
        # 131072 bytes) independent of the much larger overall ARG_MAX; a
        # dependency expression built from many thousands of job IDs (e.g.
        # MultiQC depending on every per-chunk DeepVariant/VEP/CNVkit/Manta
        # job across a 31-cell cohort, ~10,447 jobs) can exceed it, which
        # crashes the submitting bsub call itself with `OSError: [Errno 7]
        # Argument list too long` -- a failure to even submit the job, not a
        # job failure. Confirmed: this exact 31-cell cohort's dependency
        # string was ~194 KB. Skip the inline dependency above a safe margin
        # rather than crash; the caller is then responsible for only
        # invoking this job once its real prerequisites have completed
        # (see run_all_resolveome.py's own "multiqc" stage, which relies on
        # its prior "variants" stage's wait_jobs() having already done so).
        MAX_SAFE_DEPENDENCY_LEN = 100_000  # comfortable margin under 131072
        if len(dependency) > MAX_SAFE_DEPENDENCY_LEN:
            print(
                f"[WARN] Dependency expression for {job_name} is "
                f"{len(dependency):,} characters (limit ~131,072); skipping "
                f"the inline LSF dependency to avoid a submission crash. "
                f"This job is submitted WITHOUT a -w condition."
            )
        else:
            dependency_arg = f'-w "{dependency}" '

    command = (
        f'bsub -J "{job_name}" '
        f"-q {q(queue)} "
        f"-n {ncpus} "
        f"-M {memory_mb} "
        f'-R "span[hosts=1] rusage[mem={memory_mb}]" '
        f"-W {wall} "
        f"{dependency_arg}"
        f"-oo {q(out_log)} "
        f"-eo {q(err_log)} "
        f"{q(script_path)}"
    )

    output = run(command)
    job_id = parse_jobid(output)

    print(f"[SUBMITTED] {job_name} -> {job_id}")

    append_manifest(
        manifest_path,
        {
            "submission_ts": datetime.now().isoformat(timespec="seconds"),
            "sample": sample,
            "stage": stage,
            "job_name": job_name,
            "job_id": job_id,
            "queue": queue,
            "requested_ncpus": str(ncpus),
            "requested_mem_mb": str(memory_mb),
            "requested_wall": wall,
            "out_log": out_log,
            "err_log": err_log,
            "script_path": script_path,
            "depends_on": dependency or "",
        },
    )

    return job_id


# ============================================================================
# Preflight
# ============================================================================

def preflight_checklist(
    args,
    snv_tools: set,
    cnv_tools: set,
    sv_tools: set,
    annotators: set,
):
    print("[CHECK] Running preflight checklist ...")

    ensure_dir(args.outdir, f"[ERROR] Output directory not found: {args.outdir}")
    ensure_file(args.ref, f"[ERROR] Reference FASTA not found: {args.ref}")
    ensure_file(
        args.ref + ".fai",
        f"[ERROR] Reference index not found: {args.ref}.fai",
    )

    for executable in ["samtools", "bcftools", "bsub"]:
        rc, _, _ = run_no_fail(
            f"bash -lc 'command -v {q(executable)} >/dev/null 2>&1'"
        )
        if rc != 0:
            raise SystemExit(
                f"[ERROR] Required executable was not found: {executable}"
            )

    if "deepvariant" in snv_tools:
        ensure_file(
            args.deepvariant_sif,
            f"[ERROR] DeepVariant SIF not found: {args.deepvariant_sif}",
        )
        ensure_dir(
            args.deepvariant_model_dir,
            "[ERROR] DeepVariant model directory not found: "
            f"{args.deepvariant_model_dir}",
        )

        checkpoint = os.path.join(
            args.deepvariant_model_dir,
            args.deepvariant_checkpoint_filename,
        )

        ensure_file(
            checkpoint + ".index",
            f"[ERROR] Missing checkpoint index: {checkpoint}.index",
        )
        ensure_file(
            checkpoint + ".data-00000-of-00001",
            "[ERROR] Missing checkpoint data file: "
            f"{checkpoint}.data-00000-of-00001",
        )
        ensure_file(
            os.path.join(args.deepvariant_model_dir, "example_info.json"),
            "[ERROR] Missing DeepVariant example_info.json in "
            f"{args.deepvariant_model_dir}",
        )

        if args.dv_use_allele_frequency:
            if not args.deepvariant_population_vcfs:
                raise SystemExit(
                    "[ERROR] --dv-use-allele-frequency requires "
                    "--deepvariant-population-vcfs"
                )

            if not (
                os.path.isdir(args.deepvariant_population_vcfs)
                or os.path.isfile(args.deepvariant_population_vcfs)
            ):
                raise SystemExit(
                    "[ERROR] Population VCF path not found: "
                    f"{args.deepvariant_population_vcfs}"
                )

    if "cnvkit" in cnv_tools:
        check = (
            "bash -lc "
            "'source \"$(conda info --base)/etc/profile.d/conda.sh\" && "
            f"conda activate {q(args.env)} && "
            "command -v cnvkit.py >/dev/null 2>&1'"
        )

        rc, _, _ = run_no_fail(check)

        if rc != 0:
            raise SystemExit(
                f"[ERROR] cnvkit.py was not found in conda environment {args.env}"
            )

        if args.cnv_gene_bed:
            ensure_file(
                args.cnv_gene_bed,
                f"[ERROR] CNV gene BED not found: {args.cnv_gene_bed}",
            )

        if args.cnv_reference_cnn:
            ensure_file(
                args.cnv_reference_cnn,
                "[ERROR] CNV reference CNN not found: "
                f"{args.cnv_reference_cnn}",
            )

    if "manta" in sv_tools:
        ensure_file(
            args.manta_sif,
            f"[ERROR] Manta SIF not found: {args.manta_sif}",
        )

    if "snpeff" in annotators:
        check = (
            "bash -lc "
            "'source \"$(conda info --base)/etc/profile.d/conda.sh\" && "
            f"conda activate {q(args.env)} && "
            "command -v snpEff >/dev/null 2>&1'"
        )

        rc, _, _ = run_no_fail(check)

        if rc != 0:
            raise SystemExit(
                f"[ERROR] snpEff was not found in conda environment {args.env}"
            )

    if "vep" in annotators:
        ensure_dir(
            args.vep_cache_dir,
            f"[ERROR] VEP cache directory not found: {args.vep_cache_dir}",
        )

        vep_fasta = args.vep_fasta or args.ref
        ensure_file(
            vep_fasta,
            f"[ERROR] VEP FASTA not found: {vep_fasta}",
        )

        if args.vep_sif:
            ensure_file(
                args.vep_sif,
                f"[ERROR] VEP SIF not found: {args.vep_sif}",
            )
        else:
            check = (
                "bash -lc "
                "'source \"$(conda info --base)/etc/profile.d/conda.sh\" && "
                f"conda activate {q(args.env)} && "
                "command -v vep >/dev/null 2>&1'"
            )

            rc, _, _ = run_no_fail(check)

            if rc != 0:
                raise SystemExit(
                    "[ERROR] VEP executable was not found and --vep-sif was "
                    "not provided"
                )

    print("[CHECK] Preflight checklist PASSED.")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Submit Step3 SNV/CNV/SV/annotation jobs to LSF"
    )

    # Core options
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--ref", default=DEFAULT_REF)
    parser.add_argument("--env", default=DEFAULT_ENV)
    parser.add_argument("--queue", default=DEFAULT_QUEUE)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--mem-mb", type=int, default=48000)
    parser.add_argument("--wall", default="145:00")
    parser.add_argument("--chroms", default="all")
    parser.add_argument("--chunk-size-bp", type=int, default=20_000_000)
    parser.add_argument(
        "--job-manifest-name",
        default="step3_job_manifest.tsv",
        help=(
            "Filename (not a path) for this invocation's job manifest, written "
            "under <outdir>/logs_step3/. Defaults to the original hardcoded name "
            "for backward compatibility. Give concurrently-running invocations "
            "against the same --outdir distinct names (e.g. when splitting the "
            "small-variant track from the CNV/SV track) -- the manifest file is "
            "truncated at init and appended to without locking, so two "
            "invocations sharing one name would corrupt each other's rows."
        ),
    )

    # Tool switches
    parser.add_argument("--snv-tools", default="deepvariant")
    parser.add_argument("--cnv-tools", default="cnvkit")
    parser.add_argument("--sv-tools", default="manta")
    parser.add_argument("--annotators", default="snpeff,vep")

    # DeepVariant
    parser.add_argument("--deepvariant-sif", default=DEFAULT_DV_SIF)
    parser.add_argument("--deepvariant-model-dir", default=DEFAULT_DV_MODEL_DIR)
    parser.add_argument(
        "--deepvariant-checkpoint-filename",
        default=DEFAULT_DV_CKPT,
    )
    parser.add_argument(
        "--deepvariant-population-vcfs",
        default=DEFAULT_DV_POP_VCFS,
    )
    parser.add_argument(
        "--dv-use-allele-frequency",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--dv-add-hp-channel",
        action="store_true",
        # Default restored to True (2026-10-20): the 2026-10-16 change that
        # flipped this to False was based on a misdiagnosis. The real cause
        # of the "confirmed via the MT chunk" failure referenced in that
        # change was MT's missing population VCF (the bioskryb-af-20241102
        # bundle has no MT/chrM entry -- see parse_chroms()'s own comment
        # above, which already excludes MT from "--chroms all" for exactly
        # this reason): with --use_allele_frequency on, an MT chunk silently
        # produces NO allele-frequency channel, so disabling the HP channel
        # "by coincidence" restored MT's channel count to match a 7-channel
        # expectation -- while breaking every real (non-MT) chromosome,
        # since this checkpoint (sc_af.model.ckpt) was trained on 8 channels
        # (6 base + allele-frequency + HP) and raises `ValueError: Received
        # incompatible tensor with shape (3, 3, 8, 32) when attempting to
        # restore variable with shape (3, 3, 7, 32)` for every chunk once
        # the HP channel is missing. Confirmed by diffing this project's own
        # real, successful Batch 1 run (15 cells, 9,993 real VCFs produced),
        # whose own job scripts have --add_hp_channel present, against a
        # 31-cell run with this flag off, which failed 100% of DeepVariant
        # chunks across every one of its 31 cells with exactly that error.
        default=True,
    )
    parser.add_argument("--dv-channels", default="")
    parser.add_argument(
        "--snv-filter-pass-only",
        action="store_true",
        help=(
            "Filter each DeepVariant chunk to FILTER=PASS or FILTER=. "
            "before VEP annotation."
        ),
    )
    parser.add_argument(
        "--merge-existing-chunks",
        action="store_true",
        help=(
            "Reuse existing DeepVariant and VEP chunk VCFs, then submit only "
            "the per-sample merges and downstream SnpEff jobs."
        ),
    )

    # Manta
    parser.add_argument("--manta-sif", default=DEFAULT_MANTA_SIF)

    # CNVkit
    parser.add_argument("--cnv-gene-bed", default=DEFAULT_CNV_GENE_BED)
    parser.add_argument(
        "--cnv-reference-cnn",
        default=DEFAULT_CNV_REFERENCE_CNN,
    )
    parser.add_argument("--cnv-metrics-threshold", default="0.2")

    # SnpEff
    parser.add_argument("--snpeff-db", default=DEFAULT_SNPEFF_DB)
    parser.add_argument(
        "--snpeff-java-heap-mb",
        type=int,
        default=40000,
        help="Java heap for SnpEff in MB.",
    )

    # VEP
    parser.add_argument("--vep-sif", default=DEFAULT_VEP_SIF)
    parser.add_argument("--vep-cache-dir", default=DEFAULT_VEP_CACHE)
    parser.add_argument("--vep-fasta", default=DEFAULT_VEP_FASTA)
    parser.add_argument("--vep-species", default="homo_sapiens")
    parser.add_argument("--vep-assembly", default="GRCh38")
    parser.add_argument("--vep-cache-version", default="111")
    parser.add_argument("--vep-extra-args", default="")

    # Optional
    parser.add_argument("--run-multiqc-if-missing", action="store_true")

    args = parser.parse_args()

    if args.threads < 1:
        raise SystemExit("[ERROR] --threads must be at least 1")

    if args.chunk_size_bp < 1:
        raise SystemExit("[ERROR] --chunk-size-bp must be positive")

    if args.mem_mb < 1024:
        raise SystemExit("[ERROR] --mem-mb must be at least 1024 MB")

    if args.snpeff_java_heap_mb < 1024:
        raise SystemExit(
            "[ERROR] --snpeff-java-heap-mb must be at least 1024 MB"
        )

    snv_tools = csvset(args.snv_tools)
    cnv_tools = csvset(args.cnv_tools)
    sv_tools = csvset(args.sv_tools)
    annotators = csvset(args.annotators)

    snv_tools.discard("none")
    cnv_tools.discard("none")
    sv_tools.discard("none")
    annotators.discard("none")

    preflight_checklist(
        args,
        snv_tools,
        cnv_tools,
        sv_tools,
        annotators,
    )

    bam_paths = sorted(
        glob.glob(
            os.path.join(
                args.outdir,
                "bam",
                "*",
                "*.dedup.bam",
            )
        )
    )

    if not bam_paths:
        raise SystemExit(
            "[ERROR] No BAM files found under "
            f"{args.outdir}/bam/*/*.dedup.bam"
        )

    logs_dir = os.path.join(args.outdir, "logs_step3")
    scripts_dir = os.path.join(logs_dir, "job_scripts")

    variants_root = os.path.join(args.outdir, "variants")
    cnv_root = os.path.join(args.outdir, "cnv")
    sv_root = os.path.join(args.outdir, "sv")

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(scripts_dir, exist_ok=True)
    os.makedirs(variants_root, exist_ok=True)
    os.makedirs(cnv_root, exist_ok=True)
    os.makedirs(sv_root, exist_ok=True)

    manifest_path = os.path.join(
        logs_dir,
        args.job_manifest_name,
    )
    init_manifest(manifest_path)

    contigs = read_fai_contigs(
        args.ref + ".fai",
        args.chroms,
    )
    chunks = make_chunks(
        contigs,
        args.chunk_size_bp,
    )

    all_job_ids: List[str] = []

    for bam_path in bam_paths:
        sample = os.path.basename(bam_path).replace(
            ".dedup.bam",
            "",
        )

        sample_variants_dir = os.path.join(
            variants_root,
            sample,
        )
        sample_cnv_dir = os.path.join(
            cnv_root,
            sample,
        )
        sample_sv_dir = os.path.join(
            sv_root,
            sample,
        )

        # Intermediate files are kept here.
        support_dir = os.path.join(
            sample_variants_dir,
            "support",
        )
        deepvariant_chunks_dir = os.path.join(
            support_dir,
            "deepvariant_chunks",
        )
        vep_chunks_dir = os.path.join(
            support_dir,
            "vep_chunks",
        )
        merged_dir = os.path.join(
            support_dir,
            "merged",
        )
        temporary_dir = os.path.join(
            support_dir,
            "tmp",
        )

        for directory in [
            sample_variants_dir,
            sample_cnv_dir,
            sample_sv_dir,
            support_dir,
            deepvariant_chunks_dir,
            vep_chunks_dir,
            merged_dir,
            temporary_dir,
        ]:
            os.makedirs(directory, exist_ok=True)

        # Make sure the BAM has an index before submitting jobs.
        if not file_nonempty(bam_path + ".bai"):
            run(
                "bash -lc "
                + q(
                    f'source "$(conda info --base)/etc/profile.d/conda.sh" && '
                    f"conda activate {q(args.env)} && "
                    f"samtools index -@ 4 {q(bam_path)}"
                )
            )

        sample_jobs: Dict[str, List[str]] = {
            "deepvariant_chunks": [],
            "deepvariant_merge": [],
            "pass_filter_chunks": [],
            "vep_chunks": [],
            "vep_merge": [],
            "snpeff": [],
            "cnv": [],
            "sv": [],
        }

        deepvariant_chunk_paths: List[str] = []
        vep_chunk_paths: List[str] = []

        # ------------------------------------------------------------------
        # DeepVariant chunks
        # ------------------------------------------------------------------
        if "deepvariant" in snv_tools:
            checkpoint = os.path.join(
                args.deepvariant_model_dir,
                args.deepvariant_checkpoint_filename,
            )

            deepvariant_tmp_root = os.path.join(
                temporary_dir,
                "deepvariant",
            )
            os.makedirs(deepvariant_tmp_root, exist_ok=True)

            hp_arg = (
                "--add_hp_channel"
                if args.dv_add_hp_channel
                else ""
            )

            channel_arg = (
                f"--channels {q(args.dv_channels)}"
                if args.dv_channels
                else ""
            )

            for (
                fai_name,
                canonical_chrom,
                start,
                end,
                region,
                chunk_id,
            ) in chunks:
                deepvariant_chunk_vcf = os.path.join(
                    deepvariant_chunks_dir,
                    f"{sample}.deepvariant.{chunk_id}.vcf.gz",
                )
                deepvariant_chunk_paths.append(
                    deepvariant_chunk_vcf
                )

                if args.merge_existing_chunks:
                    if not file_nonempty(deepvariant_chunk_vcf):
                        raise SystemExit(
                            "[ERROR] Existing DeepVariant chunk not found or "
                            f"empty: {deepvariant_chunk_vcf}"
                        )

                    if "vep" in annotators:
                        vep_chunk_vcf = os.path.join(
                            vep_chunks_dir,
                            f"{sample}.deepvariant.{chunk_id}.vep.vcf.gz",
                        )
                        if not file_nonempty(vep_chunk_vcf):
                            raise SystemExit(
                                "[ERROR] Existing VEP chunk not found or "
                                f"empty: {vep_chunk_vcf}"
                            )
                        vep_chunk_paths.append(vep_chunk_vcf)
                    continue

                population_arg = ""
                allele_frequency_arg = ""

                if args.dv_use_allele_frequency:
                    population_vcf = population_vcf_for_chrom(
                        args.deepvariant_population_vcfs,
                        canonical_chrom,
                    )

                    if not population_vcf:
                        raise SystemExit(
                            "[ERROR] Could not resolve population VCF for "
                            f"{canonical_chrom}"
                        )

                    if not os.path.isfile(population_vcf):
                        if canonical_chrom == "MT":
                            print(
                                "[WARN] No mitochondrial population VCF found; "
                                "running MT without allele-frequency input."
                            )
                        else:
                            raise SystemExit(
                                "[ERROR] Population VCF not found: "
                                f"{population_vcf}"
                            )
                    else:
                        population_arg = (
                            f"--population_vcfs {q(population_vcf)}"
                        )
                        allele_frequency_arg = (
                            "--use_allele_frequency"
                        )

                deepvariant_body = textwrap.dedent(
                    f"""
                    set -euo pipefail

                    APPTAINER_BIN=$(command -v apptainer || command -v singularity || true)
                    [ -n "$APPTAINER_BIN" ] || {{
                        echo "[ERROR] apptainer/singularity not found"
                        exit 2
                    }}

                    test -s {q(args.deepvariant_sif)} || {{
                        echo "[ERROR] DeepVariant SIF not found"
                        exit 3
                    }}

                    mkdir -p {q(deepvariant_tmp_root)}

                    TMPBASE=$(mktemp -d {q(os.path.join(deepvariant_tmp_root, chunk_id + "_XXXXXX"))})
                    MODELBASE=$(mktemp -d {q(os.path.join(deepvariant_tmp_root, chunk_id + "_model_XXXXXX"))})
                    trap 'rm -rf "$TMPBASE" "$MODELBASE"' EXIT

                    echo "[INFO] Staging private DeepVariant model: $MODELBASE"
                    cp -a {q(args.deepvariant_model_dir)}/. "$MODELBASE/"
                    test -s "$MODELBASE/{args.deepvariant_checkpoint_filename}.index"
                    test -s "$MODELBASE/{args.deepvariant_checkpoint_filename}.data-00000-of-00001"

                    echo "[INFO] DeepVariant temporary directory: $TMPBASE"

                    "$APPTAINER_BIN" exec \\
                      --cleanenv \\
                      --no-home \\
                      --no-mount /programs \\
                      -B /research_jude:/research_jude \\
                      --env TMPBASE="$TMPBASE" \\
                      --env MODELBASE="$MODELBASE" \\
                      {q(args.deepvariant_sif)} \\
                      bash --noprofile --norc -c '
                        set -euo pipefail

                        test -d "$TMPBASE"

                        /opt/deepvariant/bin/make_examples \\
                          --mode calling \\
                          --ref {q(args.ref)} \\
                          --reads {q(bam_path)} \\
                          --regions {q(region)} \\
                          --examples "$TMPBASE/make_examples.tfrecord@{args.threads}.gz" \\
                          {population_arg} \\
                          {allele_frequency_arg} \\
                          {hp_arg} \\
                          {channel_arg}

                        /opt/deepvariant/bin/call_variants \\
                          --examples "$TMPBASE/make_examples.tfrecord@{args.threads}.gz" \\
                          --outfile "$TMPBASE/call_variants_output.tfrecord.gz" \\
                          --checkpoint "$MODELBASE/{args.deepvariant_checkpoint_filename}"

                        /opt/deepvariant/bin/postprocess_variants \\
                          --ref {q(args.ref)} \\
                          --infile "$TMPBASE/call_variants_output.tfrecord.gz" \\
                          --outfile {q(deepvariant_chunk_vcf)}
                      '

                    source "$(conda info --base)/etc/profile.d/conda.sh"
                    conda activate {q(args.env)}

                    test -s {q(deepvariant_chunk_vcf)}
                    bcftools index -f {q(deepvariant_chunk_vcf)}

                    echo "[OK] DeepVariant chunk completed"
                    echo "[OK] Chunk: {chunk_id}"
                    echo "[OK] DeepVariant VCF: {deepvariant_chunk_vcf}"
                    """
                ).strip()

                job_name = f"snv_dv_{sample}_{chunk_id}"
                out_log = os.path.join(
                    logs_dir,
                    f"snv_dv_{sample}_{chunk_id}.%J.out",
                )
                err_log = os.path.join(
                    logs_dir,
                    f"snv_dv_{sample}_{chunk_id}.%J.err",
                )

                deepvariant_job_id = submit_job(
                    manifest_path=manifest_path,
                    sample=sample,
                    stage="SNV_DEEPVARIANT_CHUNK",
                    job_name=job_name,
                    queue=args.queue,
                    ncpus=max(2, args.threads // 2),
                    # 32 GB floor (up from 24 GB, 2026-10-21): no chunk has
                    # actually hit TERM_MEMLIMIT across either production
                    # cohort (2,962 + 10,381 chunks checked), but this bumps
                    # headroom proactively rather than reactively, since
                    # --dv-add-hp-channel is now always on (the HP/haplotype
                    # channel was restored to its required default=True on
                    # 2026-10-20 -- see that flag's own comment). An extra
                    # input channel means a larger per-example tensor for
                    # call_variants' CNN to hold in memory than the earlier,
                    # briefly-used 7-channel (no-HP) configuration needed.
                    memory_mb=max(32000, args.mem_mb // 2),
                    wall=args.wall,
                    out_log=out_log,
                    err_log=err_log,
                    body=deepvariant_body,
                    script_dir=scripts_dir,
                )

                sample_jobs["deepvariant_chunks"].append(
                    deepvariant_job_id
                )
                all_job_ids.append(deepvariant_job_id)

                # ----------------------------------------------------------
                # Optional PASS-only filtering for this chunk
                # ----------------------------------------------------------
                vep_input_vcf = deepvariant_chunk_vcf

                if args.snv_filter_pass_only:
                    pass_chunk_vcf = os.path.join(
                        deepvariant_chunks_dir,
                        f"{sample}.deepvariant.{chunk_id}.pass.vcf.gz",
                    )

                    filter_body = textwrap.dedent(
                        f"""
                        source "$(conda info --base)/etc/profile.d/conda.sh"
                        conda activate {q(args.env)}

                        FINAL={q(pass_chunk_vcf)}
                        TMP="${{FINAL}}.tmp.${{LSB_JOBID:-$$}}"

                        trap 'rm -f "$TMP" "$TMP.tbi" "$TMP.csi"' EXIT
                        rm -f "$TMP" "$TMP.tbi" "$TMP.csi"

                        bcftools view \\
                          -f PASS,. \\
                          -Oz \\
                          -o "$TMP" \\
                          {q(deepvariant_chunk_vcf)}

                        test -s "$TMP"
                        bcftools view -h "$TMP" >/dev/null

                        mv -f "$TMP" "$FINAL"
                        bcftools index -f "$FINAL"

                        trap - EXIT
                        echo "[OK] PASS-only DeepVariant chunk created"
                        echo "[OK] PASS-only VCF: $FINAL"
                        """
                    ).strip()

                    filter_job_name = (
                        f"snv_filterpass_{sample}_{chunk_id}"
                    )
                    filter_out_log = os.path.join(
                        logs_dir,
                        f"snv_filterpass_{sample}_{chunk_id}.%J.out",
                    )
                    filter_err_log = os.path.join(
                        logs_dir,
                        f"snv_filterpass_{sample}_{chunk_id}.%J.err",
                    )

                    filter_job_id = submit_job(
                        manifest_path=manifest_path,
                        sample=sample,
                        stage="SNV_FILTER_PASS_CHUNK",
                        job_name=filter_job_name,
                        queue=args.queue,
                        ncpus=1,
                        memory_mb=max(4000, args.mem_mb // 8),
                        wall=args.wall,
                        out_log=filter_out_log,
                        err_log=filter_err_log,
                        body=filter_body,
                        script_dir=scripts_dir,
                        dependency=f"done({deepvariant_job_id})",
                    )

                    sample_jobs["pass_filter_chunks"].append(
                        filter_job_id
                    )
                    all_job_ids.append(filter_job_id)
                    vep_input_vcf = pass_chunk_vcf

                # ----------------------------------------------------------
                # VEP for this chunk
                # ----------------------------------------------------------
                if "vep" in annotators:
                    vep_chunk_vcf = os.path.join(
                        vep_chunks_dir,
                        f"{sample}.deepvariant.{chunk_id}.vep.vcf.gz",
                    )
                    vep_chunk_paths.append(vep_chunk_vcf)

                    vep_fasta = args.vep_fasta or args.ref

                    cache_version_arg = (
                        f"--cache_version {q(args.vep_cache_version)}"
                        if args.vep_cache_version
                        else ""
                    )

                    extra_args = args.vep_extra_args.strip()

                    if args.vep_sif:
                        vep_command = textwrap.dedent(
                            f"""
                            "$APPTAINER_BIN" exec \\
                              --cleanenv \\
                              --no-home \\
                              --no-mount /programs \\
                              -B /research_jude:/research_jude \\
                              {q(args.vep_sif)} \\
                              bash --noprofile --norc -c '
                                set -euo pipefail

                                vep \\
                                  --input_file {q(vep_input_vcf)} \\
                                  --output_file STDOUT \\
                                  --no_stats \\
                                  --format vcf \\
                                  --vcf \\
                                  --compress_output bgzip \\
                                  --offline \\
                                  --cache \\
                                  --dir_cache {q(args.vep_cache_dir)} \\
                                  --species {q(args.vep_species)} \\
                                  --assembly {q(args.vep_assembly)} \\
                                  --fasta {q(vep_fasta)} \\
                                  --everything \\
                                  --fork {max(1, args.threads // 2)} \\
                                  {cache_version_arg} \\
                                  {extra_args}
                              ' > "$TMP"
                            """
                        ).strip()
                    else:
                        vep_command = textwrap.dedent(
                            f"""
                            vep \\
                              --input_file {q(vep_input_vcf)} \\
                              --output_file STDOUT \\
                              --no_stats \\
                              --format vcf \\
                              --vcf \\
                              --compress_output bgzip \\
                              --offline \\
                              --cache \\
                              --dir_cache {q(args.vep_cache_dir)} \\
                              --species {q(args.vep_species)} \\
                              --assembly {q(args.vep_assembly)} \\
                              --fasta {q(vep_fasta)} \\
                              --everything \\
                              --fork {max(1, args.threads // 2)} \\
                              {cache_version_arg} \\
                              {extra_args} \\
                              > "$TMP"
                            """
                        ).strip()

                    vep_dependency_job = (
                        sample_jobs["pass_filter_chunks"][-1]
                        if args.snv_filter_pass_only
                        else deepvariant_job_id
                    )

                    vep_body = textwrap.dedent(
                        f"""
                        set -euo pipefail

                        APPTAINER_BIN=$(command -v apptainer || command -v singularity || true)

                        FINAL={q(vep_chunk_vcf)}
                        TMP="${{FINAL}}.tmp.${{LSB_JOBID:-$$}}"

                        trap 'rm -f "$TMP" "$TMP.tbi" "$TMP.csi"' EXIT
                        rm -f "$TMP" "$TMP.tbi" "$TMP.csi"

                        {vep_command}

                        test -s "$TMP"
                        source "$(conda info --base)/etc/profile.d/conda.sh"
                        conda activate {q(args.env)}

                        bcftools view -h "$TMP" >/dev/null
                        bcftools view -h "$TMP" | grep -Eqi '##INFO=<ID=CSQ|VEP'

                        mv -f "$TMP" "$FINAL"
                        bcftools index -f "$FINAL"

                        trap - EXIT

                        echo "[OK] VEP chunk annotation completed"
                        echo "[OK] Chunk: {chunk_id}"
                        echo "[OK] VEP chunk VCF: $FINAL"
                        """
                    ).strip()

                    vep_job_name = f"ann_vep_{sample}_{chunk_id}"
                    vep_out_log = os.path.join(
                        logs_dir,
                        f"ann_vep_{sample}_{chunk_id}.%J.out",
                    )
                    vep_err_log = os.path.join(
                        logs_dir,
                        f"ann_vep_{sample}_{chunk_id}.%J.err",
                    )

                    vep_job_id = submit_job(
                        manifest_path=manifest_path,
                        sample=sample,
                        stage="ANN_VEP_CHUNK",
                        job_name=vep_job_name,
                        queue=args.queue,
                        ncpus=max(2, args.threads // 2),
                        memory_mb=max(16000, args.mem_mb // 3),
                        wall=args.wall,
                        out_log=vep_out_log,
                        err_log=vep_err_log,
                        body=vep_body,
                        script_dir=scripts_dir,
                        dependency=f"done({vep_dependency_job})",
                    )

                    sample_jobs["vep_chunks"].append(vep_job_id)
                    all_job_ids.append(vep_job_id)

        # ------------------------------------------------------------------
        # Merge raw DeepVariant chunks into support/merged
        # ------------------------------------------------------------------
        raw_deepvariant_merged_vcf = os.path.join(
            merged_dir,
            f"{sample}.deepvariant.raw.vcf.gz",
        )

        if deepvariant_chunk_paths:
            merge_parallel_jobs = max(4, args.threads)
            raw_merge_body = textwrap.dedent(
                f"""
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate {q(args.env)}

                FINAL={q(raw_deepvariant_merged_vcf)}
                TMP="${{FINAL}}.tmp.${{LSB_JOBID:-$$}}"
                MERGE_TMP=$(mktemp -d {q(os.path.join(temporary_dir, "merge_raw_XXXXXX"))})
                SAMPLE_NAME={q(sample)}
                INPUTS=(
                  {" ".join(q(path) for path in deepvariant_chunk_paths)}
                )
                NORMALIZED_INPUTS=()

                trap 'rm -rf "$MERGE_TMP"; rm -f "$TMP" "$TMP.tbi" "$TMP.csi"' EXIT
                rm -f "$TMP" "$TMP.tbi" "$TMP.csi"

                printf '%s\n' "$SAMPLE_NAME" > "$MERGE_TMP/sample_name.txt"

                # Reheader+index each chunk in parallel (I/O-bound NFS round
                # trips dominate here, not CPU, so this is a significant
                # speedup over a sequential loop across ~166 chunks/cell).
                # Order is reconstructed afterwards from INPUTS so bcftools
                # concat -a still receives files in genomic order regardless
                # of parallel completion order.
                normalize_one() {{
                  local INPUT="$1" MERGE_TMP="$2"
                  test -s "$INPUT"
                  SAMPLE_COUNT=$(bcftools query -l "$INPUT" | awk 'NF {{count++}} END {{print count + 0}}')
                  if [ "$SAMPLE_COUNT" -ne 1 ]; then
                    echo "[ERROR] Expected one sample in chunk, found $SAMPLE_COUNT: $INPUT" >&2
                    return 4
                  fi
                  local NORMALIZED_INPUT="$MERGE_TMP/$(basename "$INPUT")"
                  bcftools reheader \\
                    -s "$MERGE_TMP/sample_name.txt" \\
                    -o "$NORMALIZED_INPUT" \\
                    "$INPUT"
                  bcftools index -f "$NORMALIZED_INPUT"
                }}
                export -f normalize_one

                printf '%s\n' "${{INPUTS[@]}}" \\
                  | xargs -P {merge_parallel_jobs} -I{{}} bash -c 'normalize_one "$1" "$2"' _ {{}} "$MERGE_TMP"

                for INPUT in "${{INPUTS[@]}}"; do
                  NORMALIZED_INPUTS+=("$MERGE_TMP/$(basename "$INPUT")")
                done
                for NORMALIZED_INPUT in "${{NORMALIZED_INPUTS[@]}}"; do
                  test -s "$NORMALIZED_INPUT"
                done

                bcftools concat \\
                  -a \\
                  -Oz \\
                  -o "$TMP" \\
                  "${{NORMALIZED_INPUTS[@]}}"

                test -s "$TMP"
                bcftools view -h "$TMP" >/dev/null

                mv -f "$TMP" "$FINAL"
                bcftools index -f "$FINAL"

                trap - EXIT
                echo "[OK] Raw DeepVariant chunks merged"
                echo "[OK] Raw merged VCF: $FINAL"
                """
            ).strip()

            raw_merge_job_id = submit_job(
                manifest_path=manifest_path,
                sample=sample,
                stage="SNV_DEEPVARIANT_MERGE",
                job_name=f"snv_dv_merge_{sample}",
                queue=args.queue,
                ncpus=merge_parallel_jobs,
                memory_mb=max(12000, args.mem_mb // 4),
                wall=args.wall,
                out_log=os.path.join(
                    logs_dir,
                    f"snv_dv_merge_{sample}.%J.out",
                ),
                err_log=os.path.join(
                    logs_dir,
                    f"snv_dv_merge_{sample}.%J.err",
                ),
                body=raw_merge_body,
                script_dir=scripts_dir,
                dependency=dep_and(sample_jobs["deepvariant_chunks"]),
            )

            sample_jobs["deepvariant_merge"].append(raw_merge_job_id)
            all_job_ids.append(raw_merge_job_id)

        # ------------------------------------------------------------------
        # Merge VEP chunks per sample
        # ------------------------------------------------------------------
        merged_vep_vcf = os.path.join(
            merged_dir,
            f"{sample}.deepvariant.raw.vep.vcf.gz",
        )

        if vep_chunk_paths:
            merge_parallel_jobs = max(4, args.threads)
            vep_merge_body = textwrap.dedent(
                f"""
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate {q(args.env)}

                FINAL={q(merged_vep_vcf)}
                TMP="${{FINAL}}.tmp.${{LSB_JOBID:-$$}}"
                MERGE_TMP=$(mktemp -d {q(os.path.join(temporary_dir, "merge_vep_XXXXXX"))})
                SAMPLE_NAME={q(sample)}
                INPUTS=(
                  {" ".join(q(path) for path in vep_chunk_paths)}
                )
                NORMALIZED_INPUTS=()

                trap 'rm -rf "$MERGE_TMP"; rm -f "$TMP" "$TMP.tbi" "$TMP.csi"' EXIT
                rm -f "$TMP" "$TMP.tbi" "$TMP.csi"

                printf '%s\n' "$SAMPLE_NAME" > "$MERGE_TMP/sample_name.txt"

                # Reheader+index each chunk in parallel (I/O-bound NFS round
                # trips dominate here, not CPU, so this is a significant
                # speedup over a sequential loop across ~166 chunks/cell).
                # Order is reconstructed afterwards from INPUTS so bcftools
                # concat -a still receives files in genomic order regardless
                # of parallel completion order.
                normalize_one() {{
                  local INPUT="$1" MERGE_TMP="$2"
                  test -s "$INPUT"
                  SAMPLE_COUNT=$(bcftools query -l "$INPUT" | awk 'NF {{count++}} END {{print count + 0}}')
                  if [ "$SAMPLE_COUNT" -ne 1 ]; then
                    echo "[ERROR] Expected one sample in chunk, found $SAMPLE_COUNT: $INPUT" >&2
                    return 4
                  fi
                  local NORMALIZED_INPUT="$MERGE_TMP/$(basename "$INPUT")"
                  bcftools reheader \\
                    -s "$MERGE_TMP/sample_name.txt" \\
                    -o "$NORMALIZED_INPUT" \\
                    "$INPUT"
                  bcftools index -f "$NORMALIZED_INPUT"
                }}
                export -f normalize_one

                printf '%s\n' "${{INPUTS[@]}}" \\
                  | xargs -P {merge_parallel_jobs} -I{{}} bash -c 'normalize_one "$1" "$2"' _ {{}} "$MERGE_TMP"

                for INPUT in "${{INPUTS[@]}}"; do
                  NORMALIZED_INPUTS+=("$MERGE_TMP/$(basename "$INPUT")")
                done
                for NORMALIZED_INPUT in "${{NORMALIZED_INPUTS[@]}}"; do
                  test -s "$NORMALIZED_INPUT"
                done

                bcftools concat \\
                  -a \\
                  -Oz \\
                  -o "$TMP" \\
                  "${{NORMALIZED_INPUTS[@]}}"

                test -s "$TMP"
                bcftools view -h "$TMP" >/dev/null
                bcftools view -h "$TMP" | grep -Eqi '##INFO=<ID=CSQ|VEP'

                mv -f "$TMP" "$FINAL"
                bcftools index -f "$FINAL"

                trap - EXIT
                echo "[OK] VEP chunks merged for sample"
                echo "[OK] Merged VEP VCF: $FINAL"
                """
            ).strip()

            vep_merge_job_id = submit_job(
                manifest_path=manifest_path,
                sample=sample,
                stage="ANN_VEP_MERGE",
                job_name=f"ann_vep_merge_{sample}",
                queue=args.queue,
                ncpus=merge_parallel_jobs,
                memory_mb=max(12000, args.mem_mb // 4),
                wall=args.wall,
                out_log=os.path.join(
                    logs_dir,
                    f"ann_vep_merge_{sample}.%J.out",
                ),
                err_log=os.path.join(
                    logs_dir,
                    f"ann_vep_merge_{sample}.%J.err",
                ),
                body=vep_merge_body,
                script_dir=scripts_dir,
                dependency=dep_and(sample_jobs["vep_chunks"]),
            )

            sample_jobs["vep_merge"].append(vep_merge_job_id)
            all_job_ids.append(vep_merge_job_id)

            # --------------------------------------------------------------
            # SnpEff once on merged VEP VCF
            # --------------------------------------------------------------
            if "snpeff" in annotators:
                final_snpeff_vcf = os.path.join(
                    sample_variants_dir,
                    f"{sample}.deepvariant.raw.vep.snpeff.vcf.gz",
                )
                snpeff_html = os.path.join(
                    sample_variants_dir,
                    f"{sample}.snpeff.html",
                )
                snpeff_csv = os.path.join(
                    sample_variants_dir,
                    f"{sample}.snpeff.genes.csv",
                )

                snpeff_heap_mb = args.snpeff_java_heap_mb
                snpeff_memory_mb = max(
                    args.mem_mb,
                    snpeff_heap_mb + 8000,
                )

                snpeff_body = textwrap.dedent(
                    f"""
                    source "$(conda info --base)/etc/profile.d/conda.sh"
                    conda activate {q(args.env)}

                    FINAL={q(final_snpeff_vcf)}
                    TMP="${{FINAL}}.tmp.${{LSB_JOBID:-$$}}"

                    trap 'rm -f "$TMP" "$TMP.tbi" "$TMP.csi"' EXIT
                    rm -f "$TMP" "$TMP.tbi" "$TMP.csi"

                    snpEff \\
                      -Xmx{snpeff_heap_mb}m \\
                      -XX:+ExitOnOutOfMemoryError \\
                      -v \\
                      -stats {q(snpeff_html)} \\
                      -csvStats {q(snpeff_csv)} \\
                      {q(args.snpeff_db)} \\
                      {q(merged_vep_vcf)} \\
                      | bgzip -c > "$TMP"

                    test -s "$TMP"
                    bcftools view -h "$TMP" >/dev/null
                    bcftools view -h "$TMP" | grep -q '##INFO=<ID=ANN'

                    mv -f "$TMP" "$FINAL"
                    bcftools index -f "$FINAL"

                    trap - EXIT

                    echo "[OK] SnpEff annotation completed"
                    echo "[OK] Final annotated VCF: $FINAL"
                    echo "[OK] SnpEff report: {snpeff_html}"
                    echo "[OK] SnpEff gene statistics: {snpeff_csv}"
                    """
                ).strip()

                snpeff_job_id = submit_job(
                    manifest_path=manifest_path,
                    sample=sample,
                    stage="ANN_SNPEFF",
                    job_name=f"ann_snpeff_{sample}",
                    queue=args.queue,
                    ncpus=1,
                    memory_mb=snpeff_memory_mb,
                    wall=args.wall,
                    out_log=os.path.join(
                        logs_dir,
                        f"ann_snpeff_{sample}.%J.out",
                    ),
                    err_log=os.path.join(
                        logs_dir,
                        f"ann_snpeff_{sample}.%J.err",
                    ),
                    body=snpeff_body,
                    script_dir=scripts_dir,
                    dependency=f"done({vep_merge_job_id})",
                )

                sample_jobs["snpeff"].append(snpeff_job_id)
                all_job_ids.append(snpeff_job_id)

        # ------------------------------------------------------------------
        # CNVkit
        # ------------------------------------------------------------------
        if "cnvkit" in cnv_tools:
            cns = os.path.join(
                sample_cnv_dir,
                f"{sample}.dedup.cns",
            )
            cnr = os.path.join(
                sample_cnv_dir,
                f"{sample}.dedup.cnr",
            )
            called_cns = os.path.join(
                sample_cnv_dir,
                f"{sample}.dedup.call.cns",
            )
            diagram_pdf = os.path.join(
                sample_cnv_dir,
                f"{sample}.cnvkit.diagram.pdf",
            )
            scatter_pdf = os.path.join(
                sample_cnv_dir,
                f"{sample}.cnvkit.scatter.pdf",
            )
            gene_metrics = os.path.join(
                sample_cnv_dir,
                f"{sample}.cnvkit.genemetrics.tsv",
            )
            notes = os.path.join(
                sample_cnv_dir,
                f"{sample}.cnvkit.notes.txt",
            )

            # cnvkit.py batch rejects --fasta together with -r/--reference: -r
            # implies reusing an existing flat reference, so --fasta (used only
            # to build a new reference on the fly) must be omitted in that case.
            #
            # fasta_arg/reference_arg are combined onto ONE template line below
            # (not two, each with a conditional trailing "\"), and neither one
            # carries its own embedded backslash: whichever of the two is set,
            # the template supplies exactly one trailing "\" for that combined
            # line. A prior version had each arg on its own line with its own
            # conditional backslash; when cnv_reference_cnn was set (this
            # project's normal/required path), fasta_arg was "" with no
            # backslash, leaving a bare blank line between "--method wgs \"
            # and "-r ... \" -- bash does not continue a "\"-terminated line
            # onto a blank unescaped line, so the command silently truncated
            # to "cnvkit.py batch BAM --method wgs" (dropping -r and every
            # argument after it), which 100% of this cohort's cells hit:
            # "Option -n/--normal must be given to build a new reference if
            # -r/--reference is not used." Confirmed by diffing this
            # project's own real, successful Batch 1 run (15 cells, 15 real
            # .cnr files produced; its job scripts have no blank line) against
            # this run's job scripts (blank line present, 0 real .cnr files
            # produced for any of 31 cells).
            if args.cnv_reference_cnn:
                reference_arg = f"-r {q(args.cnv_reference_cnn)}"
                fasta_arg = ""
            else:
                reference_arg = ""
                fasta_arg = f"--fasta {q(args.ref)}"

            cnv_body = textwrap.dedent(
                f"""
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate {q(args.env)}

                mkdir -p {q(sample_cnv_dir)}

                cnvkit.py batch {q(bam_path)} \\
                  --method wgs \\
                  {fasta_arg}{reference_arg} \\
                  --output-dir {q(sample_cnv_dir)} \\
                  --diagram \\
                  --scatter \\
                  --processes {max(1, args.threads // 2)}

                generated_cns={q(os.path.join(sample_cnv_dir, sample + ".cns"))}
                generated_cnr={q(os.path.join(sample_cnv_dir, sample + ".cnr"))}

                if [ -f "$generated_cns" ]; then
                    mv -f "$generated_cns" {q(cns)}
                fi

                if [ -f "$generated_cnr" ]; then
                    mv -f "$generated_cnr" {q(cnr)}
                fi

                test -s {q(cns)}
                test -s {q(cnr)}

                cnvkit.py call {q(cns)} -o {q(called_cns)}

                cnvkit.py genemetrics \\
                  -s {q(called_cns)} \\
                  {q(cnr)} \\
                  -o {q(gene_metrics)} || true

                cnvkit.py diagram \\
                  {q(cnr)} \\
                  -s {q(called_cns)} \\
                  -o {q(diagram_pdf)} || true

                cnvkit.py scatter \\
                  {q(cnr)} \\
                  -s {q(called_cns)} \\
                  -o {q(scatter_pdf)} || true

                cat > {q(notes)} <<'NOTE'
CNVkit genemetrics was run without a gene BED argument.
Use downstream overlap of the called CNS file with the configured genes.bed
file for gene-level interpretation.
NOTE

                echo "[OK] CNVkit completed"
                echo "[OK] Called segments: {called_cns}"
                """
            ).strip()

            cnv_job_id = submit_job(
                manifest_path=manifest_path,
                sample=sample,
                stage="CNV_CNVKIT",
                job_name=f"cnv_cnvkit_{sample}",
                queue=args.queue,
                ncpus=max(2, args.threads // 2),
                memory_mb=max(12000, args.mem_mb // 4),
                wall=args.wall,
                out_log=os.path.join(
                    logs_dir,
                    f"cnv_cnvkit_{sample}.%J.out",
                ),
                err_log=os.path.join(
                    logs_dir,
                    f"cnv_cnvkit_{sample}.%J.err",
                ),
                body=cnv_body,
                script_dir=scripts_dir,
            )

            sample_jobs["cnv"].append(cnv_job_id)
            all_job_ids.append(cnv_job_id)

        # ------------------------------------------------------------------
        # Manta
        # ------------------------------------------------------------------
        if "manta" in sv_tools:
            manta_dir = os.path.join(
                sample_sv_dir,
                "manta",
            )
            manta_vcf = os.path.join(
                sample_sv_dir,
                f"{sample}.manta.sv.vcf.gz",
            )

            manta_body = textwrap.dedent(
                f"""
                APPTAINER_BIN=$(command -v apptainer || command -v singularity || true)
                [ -n "$APPTAINER_BIN" ] || {{
                    echo "[ERROR] apptainer/singularity not found"
                    exit 2
                }}

                test -s {q(args.manta_sif)} || {{
                    echo "[ERROR] Manta SIF not found"
                    exit 3
                }}

                rm -rf {q(manta_dir)}
                mkdir -p {q(manta_dir)}

                "$APPTAINER_BIN" exec \\
                  --cleanenv \\
                  --no-home \\
                  --no-mount /programs \\
                  -B /research_jude:/research_jude \\
                  {q(args.manta_sif)} \\
                  bash --noprofile --norc -c "
                    set -euo pipefail

                    command -v configManta.py >/dev/null 2>&1
                    command -v python2 >/dev/null 2>&1
                    command -v tclsh >/dev/null 2>&1

                    configManta.py \\
                      --bam {q(bam_path)} \\
                      --referenceFasta {q(args.ref)} \\
                      --runDir {q(manta_dir)}

                    python2 {q(os.path.join(manta_dir, 'runWorkflow.py'))} \\
                      -m local \\
                      -j {max(1, args.threads // 2)}
                  "

                manta_source={q(os.path.join(
                    manta_dir,
                    "results",
                    "variants",
                    "diploidSV.vcf.gz",
                ))}

                test -s "$manta_source"
                cp -f "$manta_source" {q(manta_vcf)}

                if [ -f "$manta_source.tbi" ]; then
                    cp -f "$manta_source.tbi" {q(manta_vcf + ".tbi")}
                else
                    source "$(conda info --base)/etc/profile.d/conda.sh"
                    conda activate {q(args.env)}
                    bcftools index -f {q(manta_vcf)}
                fi

                echo "[OK] Manta completed"
                echo "[OK] Manta VCF: {manta_vcf}"
                """
            ).strip()

            manta_job_id = submit_job(
                manifest_path=manifest_path,
                sample=sample,
                stage="SV_MANTA",
                job_name=f"sv_manta_{sample}",
                queue=args.queue,
                ncpus=max(2, args.threads // 2),
                memory_mb=max(16000, args.mem_mb // 3),
                wall=args.wall,
                out_log=os.path.join(
                    logs_dir,
                    f"sv_manta_{sample}.%J.out",
                ),
                err_log=os.path.join(
                    logs_dir,
                    f"sv_manta_{sample}.%J.err",
                ),
                body=manta_body,
                script_dir=scripts_dir,
            )

            sample_jobs["sv"].append(manta_job_id)
            all_job_ids.append(manta_job_id)

    # =========================================================================
    # Optional MultiQC
    # =========================================================================
    if args.run_multiqc_if_missing:
        multiqc_dir = os.path.join(
            args.outdir,
            "multiqc_step3",
        )

        if not os.path.isdir(multiqc_dir):
            multiqc_body = textwrap.dedent(
                f"""
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate {q(args.env)}

                mkdir -p {q(multiqc_dir)}
                multiqc -f \\
                  -o {q(multiqc_dir)} \\
                  {q(args.outdir)}

                echo "[OK] MultiQC completed"
                echo "[OK] MultiQC output: {multiqc_dir}"
                """
            ).strip()

            multiqc_job_id = submit_job(
                manifest_path=manifest_path,
                sample="ALL",
                stage="MULTIQC",
                job_name="multiqc_step3",
                queue=args.queue,
                ncpus=1,
                memory_mb=4000,
                wall="04:00",
                out_log=os.path.join(
                    logs_dir,
                    "multiqc_step3.%J.out",
                ),
                err_log=os.path.join(
                    logs_dir,
                    "multiqc_step3.%J.err",
                ),
                body=multiqc_body,
                script_dir=scripts_dir,
                dependency=dep_and(all_job_ids),
            )

            all_job_ids.append(multiqc_job_id)

    print(f"[DONE] Submitted {len(all_job_ids)} LSF jobs.")
    print(f"[LOGS] {logs_dir}")
    print(f"[SCRIPTS] {scripts_dir}")
    print(f"[MANIFEST] {manifest_path}")


if __name__ == "__main__":
    main()
