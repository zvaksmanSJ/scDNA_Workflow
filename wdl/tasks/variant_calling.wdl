## Per-cell small-variant calling + annotation. Faithful port of the
## DeepVariant/VEP/SnpEff portion of step3_submit_variants_cnv_sv_annot_lsf.py,
## with every bug found and fixed this session baked in:
##   - apptainer `--no-mount /programs` (prevents a hang when the host has a
##     /programs mount that collides with the container's own filesystem)
##   - no `--add_hp_channel` (incompatible with the private sc_af checkpoint)
##   - MT excluded from the "all chromosomes" sweep (no MT entry in the
##     allele-frequency population VCF bundle -> non-deterministic tensor
##     shapes in make_examples/call_variants)
##   - merge step normalizes + indexes each chunk in parallel before
##     `bcftools concat` (was previously a serial bottleneck)
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: every task's inputs
## are filled in by main.wdl's `call variants_t.<Task> as <alias> { ... }`
## blocks, using values from main.wdl's single inputs.json; there is no
## separate JSON per cell/chunk/task in normal cohort-wide use. Five tasks
## in this file, in pipeline order: PlanChunks (once per cohort) ->
## DeepVariantChunk + VepAnnotateChunk (once per cell x chunk) ->
## MergeVcfChunks (twice per cell: raw + VEP-annotated) -> SnpEffAnnotate
## (once per cell).
version 1.0

# Runs ONCE for the whole cohort (not per cell). Reads the reference's .fai
# and splits the requested chromosome(s) into ~chunk_size_bp regions, so
# DeepVariant can be parallelized across (cell x chunk) pairs downstream --
# see ../scripts/make_region_chunks.py for the actual chunking logic.
task PlanChunks {
  input {
    File ref_fasta_fai
    String chroms = "all"
    Int chunk_size_bp = 20000000
    File make_region_chunks_script
    String conda_env = "resolveome-align"
  }
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}
    python3 ~{make_region_chunks_script} \
      --fai ~{ref_fasta_fai} \
      --chroms ~{chroms} \
      --chunk-size-bp ~{chunk_size_bp} \
      --output-tsv chunks.tsv
  >>>
  output {
    Array[Array[String]] chunks = read_tsv("chunks.tsv")
  }
  runtime { cpu: 1, memory: "2 GB" }
}

# Runs once per (cell, genomic chunk) pair -- small-variant calling via
# DeepVariant (free, BSD-3, Google). Two modes, chosen by whether
# deepvariant_model_dir_tar is set (see each branch's own comment below):
#   - unset (default, general-purpose): a single `run_deepvariant
#     --model_type=WGS` call using DeepVariant's own bundled model
#   - set (this project's own real-cohort run): the lower-level 3-step
#     make_examples/call_variants/postprocess_variants pipeline against a
#     custom checkpoint (e.g. BioSkryb's private single-cell sc_af model),
#     optionally with an allele-frequency population-VCF channel
task DeepVariantChunk {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    String region        # e.g. "1:1-20000000"
    String chunk_id       # e.g. "1_1_20000000"
    String canonical_chrom  # e.g. "1" -- used to look up the matching population VCF, if any
    File deepvariant_sif
    File? deepvariant_model_dir_tar   # optional tarball of a CUSTOM model dir (ckpt.index + ckpt.data-*).
                                       # Leave unset (default) for the general-purpose free path, which
                                       # uses DeepVariant's own bundled --model_type=WGS model and needs
                                       # no extra download at all. Only set this if you have your own
                                       # (e.g. BioSkryb's private single-cell sc_af) checkpoint -- see
                                       # setup.sh for how the real cohort run in this project used one.
    String deepvariant_checkpoint_filename = "sc_af.model.ckpt"
    File? population_vcfs_tar          # optional tar of per-chromosome AF population VCFs,
                                        # named cohort-chr<N>.release_missing2ref.no_calls.vcf.gz
                                        # (only meaningful together with deepvariant_model_dir_tar --
                                        # BioSkryb-private resource for our own real cohort -- see
                                        # setup.sh; leave unset for a general deployment)
    Boolean use_allele_frequency = false
    String conda_env = "resolveome-align"
    Int threads = 8
    Int mem_gb = 24
  }

  command <<<
    set -euo pipefail
    TMPBASE=$(mktemp -d)
    APPTAINER_BIN=$(command -v apptainer || command -v singularity)

    if [ -n "~{default='' deepvariant_model_dir_tar}" ]; then
      # ---- Custom-checkpoint path (e.g. BioSkryb's private single-cell sc_af model) ----
      MODELBASE=$(mktemp -d)
      tar -xf ~{deepvariant_model_dir_tar} -C "$MODELBASE"
      test -s "$MODELBASE/~{deepvariant_checkpoint_filename}.index"
      test -s "$MODELBASE/~{deepvariant_checkpoint_filename}.data-00000-of-00001"

      population_arg=""
      af_arg=""
      if ~{use_allele_frequency} && [ -n "~{default='' population_vcfs_tar}" ]; then
        POPDIR=$(mktemp -d)
        tar -xf ~{population_vcfs_tar} -C "$POPDIR"
        POP_VCF="$POPDIR/cohort-chr~{canonical_chrom}.release_missing2ref.no_calls.vcf.gz"
        if [ -s "$POP_VCF" ]; then
          population_arg="--population_vcfs $POP_VCF"
          af_arg="--use_allele_frequency"
        else
          echo "[WARN] No population VCF for chr~{canonical_chrom}; continuing without AF channel" >&2
        fi
      fi

      "$APPTAINER_BIN" exec \
        --cleanenv \
        --no-home \
        --no-mount /programs \
        -B "$(pwd):$(pwd)" \
        --env TMPBASE="$TMPBASE" \
        --env MODELBASE="$MODELBASE" \
        ~{deepvariant_sif} \
        bash --noprofile --norc -c '
          set -euo pipefail
          /opt/deepvariant/bin/make_examples \
            --mode calling \
            --ref '~{ref_fasta}' \
            --reads '~{bam}' \
            --regions '~{region}' \
            --examples "$TMPBASE/make_examples.tfrecord@'~{threads}'.gz" \
            '"$population_arg"' \
            '"$af_arg"'

          /opt/deepvariant/bin/call_variants \
            --examples "$TMPBASE/make_examples.tfrecord@'~{threads}'.gz" \
            --outfile "$TMPBASE/call_variants_output.tfrecord.gz" \
            --checkpoint "$MODELBASE/'~{deepvariant_checkpoint_filename}'"

          /opt/deepvariant/bin/postprocess_variants \
            --ref '~{ref_fasta}' \
            --infile "$TMPBASE/call_variants_output.tfrecord.gz" \
            --outfile '~{sample_id}'.deepvariant.'~{chunk_id}'.vcf.gz
        '
      rm -rf "$MODELBASE"
    else
      # ---- General-purpose free path: DeepVariant's own bundled WGS model, no extra download ----
      "$APPTAINER_BIN" exec \
        --cleanenv \
        --no-home \
        --no-mount /programs \
        -B "$(pwd):$(pwd)" \
        --env TMPBASE="$TMPBASE" \
        ~{deepvariant_sif} \
        /opt/deepvariant/bin/run_deepvariant \
          --model_type WGS \
          --ref ~{ref_fasta} \
          --reads ~{bam} \
          --regions ~{region} \
          --output_vcf ~{sample_id}.deepvariant.~{chunk_id}.vcf.gz \
          --intermediate_results_dir "$TMPBASE" \
          --num_shards ~{threads}
    fi

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}
    test -s ~{sample_id}.deepvariant.~{chunk_id}.vcf.gz
    bcftools index -f ~{sample_id}.deepvariant.~{chunk_id}.vcf.gz
    rm -rf "$TMPBASE"
  >>>

  output {
    File chunk_vcf = "~{sample_id}.deepvariant.~{chunk_id}.vcf.gz"
    File chunk_vcf_csi = "~{sample_id}.deepvariant.~{chunk_id}.vcf.gz.csi"
  }
  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
  }
}

task MergeVcfChunks {
  input {
    String sample_id
    Array[File] chunk_vcfs
    String out_suffix        # ".deepvariant.raw.vcf.gz" or ".deepvariant.raw.vep.vcf.gz"
    Int merge_parallel_jobs = 6
    String conda_env = "resolveome-align"
    Boolean require_vep_ann = false
  }
  # Runs once per cell, twice total (once for the raw DeepVariant chunk
  # VCFs via `merge_raw` in main.wdl, once for the VEP-annotated chunk VCFs
  # via `merge_vep`). Every chunk VCF is sorted+indexed IN PARALLEL first
  # (xargs -P) -- this was a serial bottleneck before this session's fix --
  # then `bcftools concat -a` stitches them back into one whole-genome VCF.
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}

    MERGE_TMP=$(mktemp -d)
    FINAL="~{sample_id}~{out_suffix}"

    normalize_one() {
      local input="$1" outdir="$2"
      bcftools sort -Oz -o "$outdir/$(basename "$input")" "$input"
      bcftools index -f "$outdir/$(basename "$input")"
    }
    export -f normalize_one

    printf '%s\n' ~{sep=" " chunk_vcfs} \
      | xargs -P ~{merge_parallel_jobs} -I{} bash -c 'normalize_one "$1" "$2"' _ {} "$MERGE_TMP"

    NORMALIZED=()
    for f in ~{sep=" " chunk_vcfs}; do
      NORMALIZED+=("$MERGE_TMP/$(basename "$f")")
    done

    bcftools concat -a -Oz -o "$FINAL" "${NORMALIZED[@]}"
    test -s "$FINAL"
    bcftools view -h "$FINAL" >/dev/null
    # when merging the VEP-annotated chunks, double-check the VEP CSQ
    # annotation header actually made it through before declaring success
    if ~{require_vep_ann}; then
      bcftools view -h "$FINAL" | grep -Eqi '##INFO=<ID=CSQ|VEP'
    fi
    bcftools index -f "$FINAL"
  >>>
  output {
    File merged_vcf = "~{sample_id}~{out_suffix}"
    File merged_vcf_csi = "~{sample_id}~{out_suffix}.csi"
  }
  runtime { cpu: merge_parallel_jobs, memory: "16 GB" }
}

# Runs once per (cell, genomic chunk) pair -- VEP annotation (free,
# Apache-2.0, Ensembl) of that chunk's DeepVariant VCF, via its own
# apptainer container + offline cache (untarred fresh per call from
# vep_cache_dir_tar, since VEP needs the cache files writable/local, not
# shared read-only across parallel tasks).
task VepAnnotateChunk {
  input {
    String sample_id
    File chunk_vcf
    String chunk_id
    File vep_sif
    File vep_cache_dir_tar
    File ref_fasta
    String vep_species = "homo_sapiens"
    String vep_assembly = "GRCh38"
    Int vep_cache_version = 111
    Int threads = 4
    String conda_env = "resolveome-align"
  }
  command <<<
    set -euo pipefail
    CACHE_DIR=$(mktemp -d)
    tar -xf ~{vep_cache_dir_tar} -C "$CACHE_DIR"

    FINAL="~{sample_id}.deepvariant.~{chunk_id}.vep.vcf.gz"
    TMP="$FINAL.tmp"

    APPTAINER_BIN=$(command -v apptainer || command -v singularity)
    "$APPTAINER_BIN" exec \
      --cleanenv \
      --no-home \
      --no-mount /programs \
      -B "$(pwd):$(pwd)" \
      ~{vep_sif} \
      bash --noprofile --norc -c '
        set -euo pipefail
        vep \
          --input_file '~{chunk_vcf}' \
          --output_file STDOUT \
          --no_stats \
          --format vcf \
          --vcf \
          --compress_output bgzip \
          --offline \
          --cache \
          --dir_cache "'"$CACHE_DIR"'" \
          --species '~{vep_species}' \
          --assembly '~{vep_assembly}' \
          --fasta '~{ref_fasta}' \
          --everything \
          --fork '~{threads}' \
          --cache_version '~{vep_cache_version}'
      ' > "$TMP"

    test -s "$TMP"
    mv -f "$TMP" "$FINAL"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}
    bcftools index -f "$FINAL"
    rm -rf "$CACHE_DIR"
  >>>
  output {
    File vep_chunk_vcf = "~{sample_id}.deepvariant.~{chunk_id}.vep.vcf.gz"
  }
  runtime { cpu: threads, memory: "16 GB" }
}

task SnpEffAnnotate {
  input {
    String sample_id
    File merged_vep_vcf
    String snpeff_db = "GRCh38.111"
    Int snpeff_heap_mb = 8000
    String conda_env = "resolveome-align"
  }
  # Runs once per cell -- a second, complementary annotation pass (snpEff,
  # free, LGPLv3, bioconda) layered on top of the already VEP-annotated,
  # whole-genome merged VCF from MergeVcfChunks(require_vep_ann=true); this
  # is the final per-cell VCF every downstream cohort-wide stage consumes.
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}

    FINAL="~{sample_id}.deepvariant.raw.vep.snpeff.vcf.gz"
    TMP="$FINAL.tmp"

    snpEff \
      -Xmx~{snpeff_heap_mb}m \
      -XX:+ExitOnOutOfMemoryError \
      -v \
      -stats ~{sample_id}.snpeff.html \
      -csvStats ~{sample_id}.snpeff.genes.csv \
      ~{snpeff_db} \
      ~{merged_vep_vcf} \
      | bgzip -c > "$TMP"

    test -s "$TMP"
    bcftools view -h "$TMP" >/dev/null
    bcftools view -h "$TMP" | grep -q '##INFO=<ID=ANN'
    mv -f "$TMP" "$FINAL"
    bcftools index -f "$FINAL"
  >>>
  output {
    File final_annotated_vcf = "~{sample_id}.deepvariant.raw.vep.snpeff.vcf.gz"
    File final_annotated_vcf_csi = "~{sample_id}.deepvariant.raw.vep.snpeff.vcf.gz.csi"
    File snpeff_html = "~{sample_id}.snpeff.html"
    File snpeff_genes_csv = "~{sample_id}.snpeff.genes.csv"
  }
  runtime { cpu: 1, memory: "~{snpeff_heap_mb / 1000 + 2} GB" }
}
