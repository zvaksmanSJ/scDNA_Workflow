## Per-cell FASTQ -> deduplicated BAM alignment + QC.
## Sprocket port of ../../wdl/tasks/align_qc.wdl -- same pipeline logic,
## adapted for Sprocket's per-task container model (see ../README.md for
## the full rationale). Starts from already-prepared, already-merged
## per-cell FASTQs (no multi-lane FASTQ-combining step here by design).
##
## THIS FILE IS IMPORTED BY ../main.wdl, NOT RUN STANDALONE: its `input {}`
## fields are filled in by main.wdl's `call align_qc_t.AlignTrimQC as align
## { ... }` block (inside a `scatter` over every cell).
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: no `conda activate` -- every
## tool this task needs (fastp, fastqc, bwa-mem2, samtools, picard, bam2mr,
## preseq-bioskryb) must already be on PATH inside the container named by
## `resolveome_align_container`, because Sprocket's LSF+Apptainer backend
## runs the ENTIRE command block inside that one container automatically
## (see ../containers/resolveome-align.Dockerfile for exactly what that
## image bundles, and its own header comment for an important caveat about
## bam2mr/preseq-bioskryb specifically).
version 1.2

task AlignTrimQC {
  input {
    String sample_id
    File r1
    File r2
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    # The next 5 inputs are bwa-mem2's index files. They are intentionally
    # not referenced by name in the command block below -- bwa-mem2 finds
    # them next to ref_fasta by its own file-extension convention -- but
    # they must still be declared as task inputs so Sprocket actually
    # stages them alongside ref_fasta in the task's working directory.
    File ref_bwamem2_0123
    File ref_bwamem2_bwt2bit64
    File ref_bwamem2_ann
    File ref_bwamem2_amb
    File ref_bwamem2_pac
    File metrics_intervals_bed
    Int threads = 8
    String sort_mem = "3G"
    Boolean skip_fastqc = false
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
    Int mem_gb = 32
  }

  command <<<
    set -euo pipefail

    # ---- setup: output dirs + fail fast if any required tool is missing ----
    mkdir -p bam qc/fastp qc/fastqc_raw qc/fastqc_trimmed qc/picard qc/preseq trimmed_fastq

    for tool in fastp fastqc bwa-mem2 samtools picard bam2mr preseq-bioskryb; do
      command -v "$tool" >/dev/null || { echo "[ERROR] Missing tool: $tool" >&2; exit 2; }
    done

    # ---- sanity-check the input FASTQs are not truncated/corrupt before doing any real work ----
    gzip -t ~{r1}
    gzip -t ~{r2}

    # ---- thread-count split: reserve 1-2 threads for samtools sort/fixmate, rest for bwa-mem2/fastp ----
    threads=~{threads}
    if (( threads >= 8 )); then bwa_t=$((threads - 2)); sam_t=2
    elif (( threads >= 4 )); then bwa_t=$((threads - 1)); sam_t=1
    else bwa_t=1; sam_t=1
    fi

    rg="@RG\tID:~{sample_id}\tSM:~{sample_id}\tPL:ILLUMINA\tLB:DNA_lib\tPU:~{sample_id}.unit1"

    name_bam="bam/~{sample_id}.name.bam"
    fixmate_bam="bam/~{sample_id}.fixmate.bam"
    sorted_bam="bam/~{sample_id}.sorted.bam"
    dedup_bam="bam/~{sample_id}.dedup.bam"
    trim_r1="trimmed_fastq/~{sample_id}_R1_trim.fastq.gz"
    trim_r2="trimmed_fastq/~{sample_id}_R2_trim.fastq.gz"
    metrics_prefix="qc/picard/~{sample_id}"

    sort_tmp=$(mktemp -d "${TMPDIR:-/tmp}/~{sample_id}.sort.XXXXXX")

    # ---- step 1: FastQC on the raw (untrimmed) reads ----
    if [ "~{skip_fastqc}" != "true" ]; then
      fastqc -t 2 -o qc/fastqc_raw ~{r1} ~{r2}
    fi

    # ---- step 2: fastp adapter/quality/poly-G trimming ----
    fastp --thread "$threads" \
      --in1 ~{r1} --in2 ~{r2} --out1 "$trim_r1" --out2 "$trim_r2" \
      --json "qc/fastp/~{sample_id}.fastp.json" \
      --html "qc/fastp/~{sample_id}.fastp.html" \
      --detect_adapter_for_pe --trim_poly_g
    gzip -t "$trim_r1" "$trim_r2"
    # ---- step 3: FastQC again, now on the trimmed reads ----
    if [ "~{skip_fastqc}" != "true" ]; then
      fastqc -t 2 -o qc/fastqc_trimmed "$trim_r1" "$trim_r2"
    fi

    # ---- step 4: bwa-mem2 alignment, piped straight into a name-sort (needed by fixmate next) ----
    bwa-mem2 mem -M -Y -K 2500000000 -t "$bwa_t" -R "$rg" ~{ref_fasta} "$trim_r1" "$trim_r2" \
      | samtools sort -n -@ "$sam_t" -m ~{sort_mem} -T "$sort_tmp/~{sample_id}.name" -o "$name_bam" -

    # fixmate (mate-info, needed for markdup) -> coordinate-sort -> markdup (the actual dedup step)
    samtools fixmate -m -@ "$sam_t" "$name_bam" "$fixmate_bam"
    samtools sort -@ "$sam_t" -m ~{sort_mem} -T "$sort_tmp/~{sample_id}.coord" -o "$sorted_bam" "$fixmate_bam"
    samtools flagstat "$sorted_bam" > "qc/picard/~{sample_id}.pre_dedup.flagstat.txt"
    samtools markdup -r -s -f "$metrics_prefix.markdup_stats.txt" -@ "$sam_t" "$sorted_bam" "$dedup_bam"
    samtools index -@ "$sam_t" "$dedup_bam"
    samtools flagstat "$dedup_bam" > "bam/~{sample_id}.dedup.flagstat.txt"
    samtools idxstats "$dedup_bam" > "bam/~{sample_id}.dedup.idxstats.txt"
    samtools quickcheck -v "$dedup_bam"

    # ---- step 5: Picard alignment/insert-size/GC-bias/quality-yield + WGS coverage metrics ----
    picard -Xmx8g CollectMultipleMetrics \
      R=~{ref_fasta} I="$dedup_bam" O="$metrics_prefix" \
      PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics PROGRAM=MeanQualityByCycle PROGRAM=CollectQualityYieldMetrics
    picard -Xmx8g CollectWgsMetrics R=~{ref_fasta} I="$dedup_bam" O="$metrics_prefix.wgs_metrics"

    for suffix in alignment_summary_metrics gc_bias.summary_metrics quality_yield_metrics wgs_metrics; do
      test -s "$metrics_prefix.$suffix" || { echo "[ERROR] Missing Picard $suffix" >&2; exit 3; }
    done
    # insert_size_metrics can legitimately be empty for libraries with no properly-paired reads
    if [ ! -s "$metrics_prefix.insert_size_metrics" ]; then
      printf 'not_computable\tno_eligible_pairs\n' > "$metrics_prefix.insert_size.status.tsv"
    fi
    # ---- step 6: per-contig coverage summary over metrics_intervals_bed ----
    samtools depth -b ~{metrics_intervals_bed} "$dedup_bam" \
      | awk 'BEGIN {OFS="\t"; print "Target","total_coverage"} {total[$1]+=$3} END {for (c in total) print c,total[c]}' \
      > "$metrics_prefix.cov_sentieonmetrics.sample_interval_summary"

    # ---- step 7: preseq library-complexity extrapolation (gracefully degrades if the model can't fit) ----
    mkdir -p "qc/preseq/~{sample_id}"
    bam2mr -seg_len 100000 -o "qc/preseq/~{sample_id}/~{sample_id}.mr" "$dedup_bam"
    if preseq-bioskryb gc_extrap -o "qc/preseq/~{sample_id}/~{sample_id}.curve" "qc/preseq/~{sample_id}/~{sample_id}.mr" \
        > "qc/preseq/~{sample_id}/~{sample_id}.preseq.log" 2>&1; then
      tail -n 1 "qc/preseq/~{sample_id}/~{sample_id}.curve" | cut -f 2 \
        | awk '{if ($1 !~ /^[0-9.eE+-]+$/) exit 2; printf "%.0f\n", $1}' > "qc/preseq/~{sample_id}/~{sample_id}_preseq.txt"
    else
      echo "[WARN] ~{sample_id}: preseq model fit failed" >&2
      printf 'status\treason\nNOT_COMPUTABLE\tpreseq model fit failed\n' > "qc/preseq/~{sample_id}/~{sample_id}.preseq.status.tsv"
    fi
    # Normalize to one deterministic filename regardless of which branch above
    # ran, so CellQcGate (tasks/cell_qc_gate.wdl) has a single path to read --
    # it already handles both the plain-count and the NOT_COMPUTABLE formats.
    cp "qc/preseq/~{sample_id}/~{sample_id}_preseq.txt" "qc/preseq/~{sample_id}/~{sample_id}.preseq_result.txt" 2>/dev/null \
      || cp "qc/preseq/~{sample_id}/~{sample_id}.preseq.status.tsv" "qc/preseq/~{sample_id}/~{sample_id}.preseq_result.txt"

    rm -f "$name_bam" "$fixmate_bam" "$sorted_bam"
    rm -rf "$sort_tmp"
  >>>

  output {
    File dedup_bam = "bam/~{sample_id}.dedup.bam"
    File dedup_bam_bai = "bam/~{sample_id}.dedup.bam.bai"
    File flagstat = "bam/~{sample_id}.dedup.flagstat.txt"
    File idxstats = "bam/~{sample_id}.dedup.idxstats.txt"
    File wgs_metrics = "qc/picard/~{sample_id}.wgs_metrics"
    # Named (not just globbed) outputs the per-cell QC gate (tasks/cell_qc_gate.wdl,
    # called right after this task in main.wdl) reads directly.
    File alignment_summary_metrics = "qc/picard/~{sample_id}.alignment_summary_metrics"
    File preseq_result = "qc/preseq/~{sample_id}/~{sample_id}.preseq_result.txt"
    Array[File] fastp_out = glob("qc/fastp/*")
    Array[File] fastqc_raw_out = glob("qc/fastqc_raw/*")
    Array[File] fastqc_trimmed_out = glob("qc/fastqc_trimmed/*")
    Array[File] picard_out = glob("qc/picard/*")
    Array[File] preseq_out = glob("qc/preseq/~{sample_id}/*")
  }

  requirements {
    container: resolveome_align_container
    cpu: threads
    memory: "~{mem_gb} GB"
  }
}
