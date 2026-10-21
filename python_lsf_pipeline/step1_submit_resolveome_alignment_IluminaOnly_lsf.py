#!/usr/bin/env python3
import argparse
import csv
import json
import re
import shlex
import subprocess
from datetime import datetime
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Step1: Submit Illumina FASTQ->BAM alignment jobs to LSF (one job per cell)."
    )
    p.add_argument("--cells", required=True, help="TSV header: sample_id\\tr1\\tr2")
    p.add_argument("--outdir", required=True, help="Output root directory")
    p.add_argument("--ref", required=True, help="Reference FASTA")
    p.add_argument("--queue", default="standard")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--mem-mb", type=int, default=32000)
    p.add_argument("--wall", default="145:00")
    p.add_argument("--env", default="resolveome-align")
    p.add_argument("--assay", default="DNA")
    p.add_argument("--platform", default="ILLUMINA", choices=["ILLUMINA"])
    p.add_argument("--sort-mem", default="3G", help="samtools sort memory per thread (e.g. 2G, 3G)")
    p.add_argument("--skip-fastqc", action="store_true", help="Skip FastQC step entirely")
    p.add_argument("--metrics-intervals", help="BED used for source-style chromosome depth sums; defaults to full canonical chromosomes")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def ensure_dirs(outdir: Path):
    for d in ["logs", "bam", "qc", "manifests", "scripts", "tmp", "trimmed_fastq"]:
        (outdir / d).mkdir(parents=True, exist_ok=True)


def validate_inputs(args):
    if args.threads < 2:
        raise ValueError("--threads must be at least 2")
    if not Path(args.cells).is_file():
        raise FileNotFoundError(f"--cells not found: {args.cells}")
    if not Path(args.ref).is_file():
        raise FileNotFoundError(f"--ref not found: {args.ref}")
    if not Path(args.ref + ".fai").is_file():
        raise FileNotFoundError(f"Missing FASTA index: {args.ref}.fai")
    # bwa-mem2 index presence (either style)
    if not (Path(args.ref + ".bwt.2bit.64").exists() or Path(args.ref + ".0123").exists()):
        raise FileNotFoundError(
            f"Missing bwa-mem2 index for {args.ref}. Run: bwa-mem2 index {args.ref}"
        )
    if not Path(args.ref).with_suffix(".dict").is_file():
        raise FileNotFoundError(f"Missing Picard sequence dictionary for {args.ref}")


def write_worker_script(worker_path: Path):
    worker = r"""#!/usr/bin/env bash
set -euo pipefail

sample_id="$1"; r1="$2"; r2="$3"; outdir="$4"; ref="$5"; threads="$6"; env_name="$7"; platform="$8"; assay="$9"; manifest="${10}"; sort_mem="${11}"; skip_fastqc="${12}"; metrics_intervals="${13}"

mkdir -p "${outdir}"/{bam,qc,logs,tmp,trimmed_fastq}
sample_dir="${outdir}/bam/${sample_id}"
fastp_dir="${outdir}/qc/fastp/${sample_id}"
picard_dir="${outdir}/qc/picard/${sample_id}"
raw_qc="${outdir}/qc/fastqc_raw/${sample_id}"
trim_qc="${outdir}/qc/fastqc_trimmed/${sample_id}"
trim_dir="${outdir}/trimmed_fastq/${sample_id}"
mkdir -p "${sample_dir}" "${fastp_dir}" "${picard_dir}" "${raw_qc}" "${trim_qc}" "${trim_dir}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${env_name}"

for tool in fastp fastqc bwa-mem2 samtools picard bam2mr preseq-bioskryb; do
  command -v "${tool}" >/dev/null || { echo "[ERROR] Missing tool: ${tool}" >&2; exit 2; }
done
[[ -s "${r1}" ]] || { echo "Missing R1: ${r1}" >&2; exit 2; }
[[ -s "${r2}" ]] || { echo "Missing R2: ${r2}" >&2; exit 2; }
[[ -s "${ref}" ]] || { echo "Missing REF: ${ref}" >&2; exit 2; }
[[ -s "${ref}.fai" ]] || { echo "Missing REF fai: ${ref}.fai" >&2; exit 2; }

# FASTQ gzip integrity checks (fast fail)
gzip -t "${r1}" || { echo "Corrupt gzip R1: ${r1}" >&2; exit 2; }
gzip -t "${r2}" || { echo "Corrupt gzip R2: ${r2}" >&2; exit 2; }

# Thread split: leave some for samtools
if (( threads >= 8 )); then
  bwa_t=$((threads - 2))
  sam_t=2
elif (( threads >= 4 )); then
  bwa_t=$((threads - 1))
  sam_t=1
else
  bwa_t=1
  sam_t=1
fi

rg="@RG\tID:${sample_id}\tSM:${sample_id}\tPL:${platform}\tLB:${assay}_lib\tPU:${sample_id}.unit1"

name_bam="${sample_dir}/${sample_id}.name.bam"
fixmate_bam="${sample_dir}/${sample_id}.fixmate.bam"
sorted_bam="${sample_dir}/${sample_id}.sorted.bam"
dedup_bam="${sample_dir}/${sample_id}.dedup.bam"
flagstat="${sample_dir}/${sample_id}.dedup.flagstat.txt"
idxstats="${sample_dir}/${sample_id}.dedup.idxstats.txt"
pre_flagstat="${picard_dir}/${sample_id}.pre_dedup.flagstat.txt"
metrics_prefix="${picard_dir}/${sample_id}"
trim_r1="${trim_dir}/${sample_id}_R1_trim.fastq.gz"
trim_r2="${trim_dir}/${sample_id}_R2_trim.fastq.gz"
completion="${sample_dir}/alignment.complete"
rm -f "${completion}"

# Faster sort temp on local scratch if available
scratch_root="${LS_TMPDIR:-${TMPDIR:-/tmp}}"
sort_tmp=$(mktemp -d "${scratch_root}/${sample_id}.sort.XXXXXX")

cleanup_on_error() {
  echo "[ERROR] ${sample_id}: cleaning partial outputs" >&2
  rm -f "${name_bam}" "${fixmate_bam}" "${sorted_bam}" "${dedup_bam}" "${dedup_bam}.bai" "${flagstat}" "${idxstats}"
}
trap cleanup_on_error ERR

if [[ "${skip_fastqc}" != "1" ]]; then
  fastqc -t 2 -o "${raw_qc}" "${r1}" "${r2}"
fi

# Force poly-G trimming for all instrument headers; keep the reference fastp defaults.
fastp --thread "${threads}" \
  --in1 "${r1}" --in2 "${r2}" --out1 "${trim_r1}" --out2 "${trim_r2}" \
  --json "${fastp_dir}/${sample_id}.fastp.json" \
  --html "${fastp_dir}/${sample_id}.fastp.html" \
  --detect_adapter_for_pe --trim_poly_g
gzip -t "${trim_r1}" "${trim_r2}"
if [[ "${skip_fastqc}" != "1" ]]; then
  fastqc -t 2 -o "${trim_qc}" "${trim_r1}" "${trim_r2}"
fi

bwa-mem2 mem -M -Y -K 2500000000 -t "${bwa_t}" -R "${rg}" "${ref}" "${trim_r1}" "${trim_r2}" \
  | samtools sort -n -@ "${sam_t}" -m "${sort_mem}" -T "${sort_tmp}/${sample_id}.name" -o "${name_bam}" -

samtools fixmate -m -@ "${sam_t}" "${name_bam}" "${fixmate_bam}"

samtools sort -@ "${sam_t}" -m "${sort_mem}" -T "${sort_tmp}/${sample_id}.coord" -o "${sorted_bam}" "${fixmate_bam}"

samtools flagstat "${sorted_bam}" > "${pre_flagstat}"
samtools markdup -r -s -f "${metrics_prefix}.markdup_stats.txt" -@ "${sam_t}" "${sorted_bam}" "${dedup_bam}"
samtools index -@ "${sam_t}" "${dedup_bam}"
samtools flagstat "${dedup_bam}" > "${flagstat}"
samtools idxstats "${dedup_bam}" > "${idxstats}"
samtools quickcheck -v "${dedup_bam}"

picard -Xmx8g CollectMultipleMetrics \
  R="${ref}" I="${dedup_bam}" O="${metrics_prefix}" \
  PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics \
  PROGRAM=CollectGcBiasMetrics PROGRAM=MeanQualityByCycle PROGRAM=CollectQualityYieldMetrics

# The WGS collector complements the DNAQC library-metrics suite.
picard -Xmx8g CollectWgsMetrics R="${ref}" I="${dedup_bam}" O="${metrics_prefix}.wgs_metrics"

for suffix in alignment_summary_metrics gc_bias.summary_metrics quality_yield_metrics wgs_metrics; do
  test -s "${metrics_prefix}.${suffix}" || { echo "[ERROR] Missing Picard ${suffix}" >&2; exit 3; }
done
if [[ ! -s "${metrics_prefix}.insert_size_metrics" ]]; then
  echo "[WARN] No eligible pairs for insert-size metrics: ${sample_id}" >&2
  printf 'not_computable\tno_eligible_pairs\n' > "${metrics_prefix}.insert_size.status.tsv"
fi
samtools depth -b "${metrics_intervals}" "${dedup_bam}" \
  | awk 'BEGIN {OFS="\t"; print "Target","total_coverage"} {total[$1]+=$3} END {for (c in total) print c,total[c]}' \
  > "${metrics_prefix}.cov_sentieonmetrics.sample_interval_summary"

preseq_dir="${outdir}/qc/preseq/${sample_id}"
mkdir -p "${preseq_dir}"
bam2mr -seg_len 100000 -o "${preseq_dir}/${sample_id}.mr" "${dedup_bam}"
preseq_estimate="${preseq_dir}/${sample_id}_preseq.txt"
preseq_status="${preseq_dir}/${sample_id}.preseq.status.tsv"
rm -f "${preseq_estimate}" "${preseq_status}"
if preseq-bioskryb gc_extrap -o "${preseq_dir}/${sample_id}.curve" "${preseq_dir}/${sample_id}.mr" \
    > "${preseq_dir}/${sample_id}.preseq.log" 2>&1; then
  tail -n 1 "${preseq_dir}/${sample_id}.curve" | cut -f 2 \
    | awk '{if ($1 !~ /^[0-9.eE+-]+$/) exit 2; printf "%.0f\n", $1}' > "${preseq_estimate}"
  test -s "${preseq_estimate}"
else
  echo "[WARN] ${sample_id}: preseq model fit failed; source complexity QC is NOT_COMPUTABLE. See ${preseq_dir}/${sample_id}.preseq.log" >&2
  printf 'status\treason\nNOT_COMPUTABLE\tpreseq model fit failed\n' > "${preseq_status}"
fi

rm -f "${name_bam}" "${fixmate_bam}" "${sorted_bam}"
rm -rf "${sort_tmp}"

(
  flock 200
  printf '%s\t%s\t%s\t%s\t%s\n' "${sample_id}" "${dedup_bam}" "${dedup_bam}.bai" "${pre_flagstat}" "${idxstats}" >> "${manifest}"
) 200>"${manifest}.lock"

printf 'fastp_bwamem2_picard_v2\n' > "${completion}"
trap - ERR
"""
    worker_path.write_text(worker)
    worker_path.chmod(0o755)


def read_cells_tsv(cells_path: Path):
    rows = []
    with cells_path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "r1", "r2"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError("cells TSV must contain: sample_id, r1, r2")
        for row in r:
            sid = row["sample_id"].strip()
            if sid:
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid):
                    raise ValueError(f"Unsafe sample_id: {sid!r}")
                if any(item["sample_id"] == sid for item in rows):
                    raise ValueError(f"Duplicate sample_id: {sid}")
                for key in ("r1", "r2"):
                    if not Path(row[key]).is_file():
                        raise FileNotFoundError(f"Missing {key} for {sid}: {row[key]}")
                rows.append({"sample_id": sid, "r1": row["r1"].strip(), "r2": row["r2"].strip()})
    return rows


def extract_job_id(submit_output: str):
    m = re.search(r"Job <(\d+)>", submit_output)
    return m.group(1) if m else None


def write_job_script(path: Path, sample_id: str, r1: str, r2: str, args, outdir: Path, worker: Path, manifest: Path):
    skip_fastqc_flag = "1" if args.skip_fastqc else "0"
    command = shlex.join([str(worker), sample_id, r1, r2, str(outdir), args.ref,
                          str(args.threads), args.env, args.platform, args.assay,
                          str(manifest), args.sort_mem, skip_fastqc_flag, args.metrics_intervals])
    content = f"""#!/bin/bash
#BSUB -J align_{sample_id}
#BSUB -q {args.queue}
#BSUB -n {args.threads}
#BSUB -M {args.mem_mb}
#BSUB -R "span[hosts=1] rusage[mem={args.mem_mb}]"
#BSUB -W {args.wall}
#BSUB -oo {outdir}/logs/{sample_id}.%J.out
#BSUB -eo {outdir}/logs/{sample_id}.%J.err

set -euo pipefail
{command}
"""
    path.write_text(content)
    path.chmod(0o755)


def main():
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    ensure_dirs(outdir)
    validate_inputs(args)
    if any((outdir / "bam").glob("*/*.bam")):
        raise RuntimeError("BAMs already exist in --outdir; use a new output root to preserve prior results")
    if args.metrics_intervals:
        args.metrics_intervals = str(Path(args.metrics_intervals).resolve())
        if not Path(args.metrics_intervals).is_file() or Path(args.metrics_intervals).stat().st_size == 0:
            raise FileNotFoundError(f"Missing/empty metrics BED: {args.metrics_intervals}")
    else:
        bed = outdir / "manifests" / "library_qc_intervals.bed"
        with open(args.ref + ".fai") as source, bed.open("w") as target:
            for line in source:
                chrom, length, *_ = line.split("\t")
                if chrom.removeprefix("chr") in {*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"}:
                    target.write(f"{chrom}\t0\t{length}\n")
        args.metrics_intervals = str(bed)
        print("[INFO] Library-QC depth uses explicit full canonical-chromosome intervals generated from the current FASTA")
    (outdir / "manifests" / "alignment_protocol.json").write_text(json.dumps({
        "protocol": "fastp_bwamem2_picard_v2", "fastp_version": "0.20.1",
        "adapter_detection": True, "poly_g": "forced", "duplicate_removal": True,
        "bwa_flags": ["-M", "-Y", "-K", "2500000000"],
        "qc_flagstat_source": "trimmed_alignment_before_duplicate_removal",
        "reference": str(Path(args.ref).resolve()), "fastqc": not args.skip_fastqc,
        "metrics_intervals": args.metrics_intervals, "preseq_version": "2.0.3",
    }, indent=2) + "\n")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest = outdir / "manifests" / f"aligned_bams_{run_id}.tsv"
    manifest.write_text("sample_id\tbam\tbai\tflagstat\tidxstats\n")

    worker = outdir / "scripts" / "run_one_cell_alignment.sh"
    write_worker_script(worker)

    rows = read_cells_tsv(Path(args.cells))
    if not rows:
        raise RuntimeError("No valid rows in cells TSV")

    job_ids = []
    for row in rows:
        sample_id, r1, r2 = row["sample_id"], row["r1"], row["r2"]
        job_script = outdir / "scripts" / f"{sample_id}.lsf.sh"
        write_job_script(job_script, sample_id, r1, r2, args, outdir, worker, manifest)

        cmd = f"bsub < {shlex.quote(str(job_script))}"
        if args.dry_run:
            print(cmd)
            continue

        proc = subprocess.run(["bash", "-lc", cmd], check=True, capture_output=True, text=True)
        print(proc.stdout.strip())
        jid = extract_job_id(proc.stdout)
        if not jid:
            raise RuntimeError(f"Cannot parse LSF job ID: {proc.stdout}")
        job_ids.append(jid)

    latest_manifest = outdir / "manifests" / "aligned_bams_latest.tsv"
    latest_manifest.write_text(str(manifest) + "\n")

    jobids_file = outdir / "manifests" / f"alignment_job_ids_{run_id}.txt"
    jobids_file.write_text("\n".join(job_ids) + ("\n" if job_ids else ""))

    print(f"[INFO] Submitted {len(rows)} alignment jobs")
    print(f"[INFO] Manifest (fills as jobs finish): {manifest}")
    print(f"[INFO] Latest pointer: {latest_manifest}")
    print(f"[INFO] Job IDs file: {jobids_file}")
    print(run_id)


if __name__ == "__main__":
    main()
