## Per-cell CNV (CNVkit) and SV (Manta) calling. Sprocket port of
## ../../wdl/tasks/cnv_sv.wdl.
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: the Cromwell version's
## MantaCall task manually `apptainer exec`s the Manta container AND THEN
## `conda activate`s a second environment (just for the bcftools-index
## fallback). Sprocket runs the whole command block inside exactly ONE
## container per task, so that pattern is split into two tasks: `MantaCall`
## (runs entirely inside the Manta container, producing a raw VCF -- and
## its own `.tbi` index, if Manta happened to make one) and the shared
## `IndexVcf` helper task in common.wdl (resolveome-align container,
## bcftools only), called from ../main.wdl ONLY when Manta's own index is
## missing -- see ../main.wdl's `if (!defined(manta.manta_vcf_tbi))` block.
## This preserves the exact original behavior (prefer Manta's own `.tbi`;
## fall back to a `.csi` from `bcftools index` only if needed) rather than
## silently always re-indexing.
version 1.2

# CNVkit (free, Apache-2.0, pypi/bioconda) -- copy-number calling straight
# from the dedup BAM (`--method wgs`, no target/antitarget BED needed). Two
# modes depending on whether a pooled reference .cnn was supplied: reuse it
# (fast, consistent across cells) or build a flat reference on the fly from
# just this cell's own BAM + the FASTA (slower, less ideal for a cohort, but
# works standalone). Followed by `call` (CN-state calling) and `genemetrics`
# (per-gene CN summary) on the resulting .cns/.cnr. Single container
# throughout -- no split needed.
task CnvkitCall {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    File? cnv_reference_cnn
    Int threads = 4
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p cnv_out

    if [ -n "~{cnv_reference_cnn}" ]; then
      cnvkit.py batch ~{bam} \
        --method wgs \
        -r ~{cnv_reference_cnn} \
        --output-dir cnv_out \
        --diagram --scatter \
        --processes ~{threads}
    else
      cnvkit.py batch ~{bam} \
        --method wgs \
        --fasta ~{ref_fasta} \
        --output-dir cnv_out \
        --diagram --scatter \
        --processes ~{threads}
    fi

    sample_base=$(basename ~{bam} .bam)
    [ -f "cnv_out/${sample_base}.cns" ] && mv -f "cnv_out/${sample_base}.cns" "cnv_out/~{sample_id}.dedup.cns" || true
    [ -f "cnv_out/${sample_base}.cnr" ] && mv -f "cnv_out/${sample_base}.cnr" "cnv_out/~{sample_id}.dedup.cnr" || true

    cnvkit.py call "cnv_out/~{sample_id}.dedup.cns" -o "cnv_out/~{sample_id}.dedup.call.cns"
    cnvkit.py genemetrics "cnv_out/~{sample_id}.dedup.cnr" -s "cnv_out/~{sample_id}.dedup.cns" \
      > "cnv_out/~{sample_id}.cnvkit.genemetrics.tsv" || true
  >>>
  output {
    File cns = "cnv_out/~{sample_id}.dedup.cns"
    File cnr = "cnv_out/~{sample_id}.dedup.cnr"
    File call_cns = "cnv_out/~{sample_id}.dedup.call.cns"
    Array[File] plots = glob("cnv_out/*.pdf")
    File genemetrics = "cnv_out/~{sample_id}.cnvkit.genemetrics.tsv"
  }
  requirements {
    container: resolveome_align_container
    cpu: threads
    memory: "16 GB"
  }
}

# Manta (free, GPLv3, Illumina/bioconda) -- germline/diploid structural-
# variant calling. Runs entirely inside the Manta container (no more manual
# apptainer-exec wrapping, since Sprocket's LSF+Apptainer backend already
# does that for the whole task based on `requirements.container` below).
# Outputs the raw VCF plus its own `.tbi`, IF Manta produced one -- the
# `manta_vcf_tbi` output is `File?` precisely so ../main.wdl can detect
# whether a fallback `IndexVcf` call is needed.
task MantaCall {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    Int threads = 4
    String manta_container = "docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6"
  }
  command <<<
    set -euo pipefail
    RUNDIR="manta_run"
    rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"

    configManta.py \
      --bam ~{bam} \
      --referenceFasta ~{ref_fasta} \
      --runDir "$RUNDIR"
    python2 "$RUNDIR/runWorkflow.py" -m local -j ~{threads}

    test -s "$RUNDIR/results/variants/diploidSV.vcf.gz"
    cp -f "$RUNDIR/results/variants/diploidSV.vcf.gz" "~{sample_id}.manta.vcf.gz"
    if [ -f "$RUNDIR/results/variants/diploidSV.vcf.gz.tbi" ]; then
      cp -f "$RUNDIR/results/variants/diploidSV.vcf.gz.tbi" "~{sample_id}.manta.vcf.gz.tbi"
    fi
  >>>
  output {
    File manta_vcf = "~{sample_id}.manta.vcf.gz"
    File? manta_vcf_tbi = "~{sample_id}.manta.vcf.gz.tbi"
  }
  requirements {
    container: manta_container
    cpu: threads
    memory: "16 GB"
  }
}

