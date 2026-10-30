## Per-cell CNV (CNVkit) and SV (Manta) calling. Faithful port of the
## corresponding section of step3_submit_variants_cnv_sv_annot_lsf.py, with
## the `--fasta`/`-r` mutual-exclusion bug fix baked in (CNVkit's `batch`
## subcommand rejects `--fasta` together with `-r/--reference`: `-r` implies
## reusing an existing flat reference, so `--fasta` -- only needed to build a
## new reference on the fly -- must be omitted whenever a reference .cnn is
## supplied) and the apptainer `--no-mount /programs` fix for Manta.
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: both tasks' inputs
## are filled in by main.wdl's `call cnvsv_t.CnvkitCall as cnvkit { ... }`
## and `call cnvsv_t.MantaCall as manta { ... }` blocks (each inside the
## per-cell scatter), using values from main.wdl's single inputs.json --
## there is no separate JSON file per cell/task in normal cohort-wide use.
version 1.0

# CNVkit (free, Apache-2.0, pypi/bioconda) -- copy-number calling straight
# from the dedup BAM (`--method wgs`, no target/antitarget BED needed). Two
# modes depending on whether a pooled reference .cnn was supplied: reuse it
# (fast, consistent across cells) or build a flat reference on the fly from
# just this cell's own BAM + the FASTA (slower, less ideal for a cohort, but
# works standalone). Followed by `call` (CN-state calling) and `genemetrics`
# (per-gene CN summary) on the resulting .cns/.cnr.
task CnvkitCall {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    File? cnv_reference_cnn
    Int threads = 4
    String conda_env = "wdl-resolveome-align"
  }
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}
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
  runtime {
    cpu: threads
    memory: "16 GB"
  }
}

# Manta (free, GPLv3, Illumina/bioconda) -- germline/diploid structural-
# variant calling, run inside its own apptainer container (configManta.py
# writes a run directory + Python2 workflow script, then runWorkflow.py
# actually executes the SV-calling graph with -j <threads> parallel jobs).
task MantaCall {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    File manta_sif
    Int threads = 4
    String conda_env = "wdl-resolveome-align"
  }
  command <<<
    set -euo pipefail
    RUNDIR="manta_run"
    rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"

    APPTAINER_BIN=$(command -v apptainer || command -v singularity)
    "$APPTAINER_BIN" exec \
      --cleanenv \
      --no-home \
      --no-mount /programs \
      -B "$(pwd):$(pwd)" \
      ~{manta_sif} \
      bash --noprofile --norc -c '
        set -euo pipefail
        configManta.py \
          --bam '~{bam}' \
          --referenceFasta '~{ref_fasta}' \
          --runDir '"$RUNDIR"'
        python2 '"$RUNDIR"'/runWorkflow.py -m local -j '~{threads}'
      '

    test -s "$RUNDIR/results/variants/diploidSV.vcf.gz"
    cp -f "$RUNDIR/results/variants/diploidSV.vcf.gz" "~{sample_id}.manta.vcf.gz"
    if [ -f "$RUNDIR/results/variants/diploidSV.vcf.gz.tbi" ]; then
      cp -f "$RUNDIR/results/variants/diploidSV.vcf.gz.tbi" "~{sample_id}.manta.vcf.gz.tbi"
    else
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate ~{conda_env}
      bcftools index -f "~{sample_id}.manta.vcf.gz"
    fi
  >>>
  output {
    File manta_vcf = "~{sample_id}.manta.vcf.gz"
  }
  runtime {
    cpu: threads
    memory: "16 GB"
  }
}
