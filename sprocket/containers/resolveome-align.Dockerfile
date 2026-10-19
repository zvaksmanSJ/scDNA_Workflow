# ResolveOME DNA workflow -- "resolveome-align" container for Sprocket.
#
# Bundles the same tool stack as ../../setup/setup.sh's `resolveome-align`
# conda environment, as a container image, because Sprocket's execution
# model requires every task to run inside a container -- there is no
# "run bare on the host with a conda env activated" option (see
# ../README.md for the full explanation).
#
# Used by: align_qc.wdl, cnv_sv.wdl (CnvkitCall), phasing.wdl,
# qc_aggregate.wdl, variant_calling.wdl (PlanChunks, MergeVcfChunks,
# SnpEffAnnotate), extra_qc_dashboard.wdl (everything except
# MutationalSignatures), and tasks/common.wdl's IndexVcf helper.
#
# -----------------------------------------------------------------------
# !!! KNOWN GAP, NOT SOLVED BY THIS DOCKERFILE: bam2mr / preseq-bioskryb
# -----------------------------------------------------------------------
# align_qc.wdl's AlignTrimQC task calls `bam2mr` and `preseq-bioskryb`
# (BioSkryb's own reference-compatible preseq 2.0.3 build) for its library-
# complexity step. Neither tool has a public conda/pip/apt package -- they
# are historically sourced from a project-internal bundle (see
# Backup_code/setup_resolveome_alignment.sh), not from any public registry
# this Dockerfile could pull from. This gap ALREADY EXISTS in the Cromwell
# version of this workflow (../../setup/setup.sh never installs these two
# binaries either) -- it is not something this Sprocket port introduces.
#
# Before building this image for real use, obtain `bam2mr` and
# `preseq-bioskryb` yourself (e.g. copy them from this project's existing
# resolveome_align/tools/ installation) and place them in a local
# `./bin/` directory next to this Dockerfile, then uncomment the COPY line
# below. Without them, AlignTrimQC's preseq step will fail outright (not
# gracefully degrade -- `bam2mr` is called unconditionally, before the
# existing graceful NOT_COMPUTABLE fallback, which only covers
# preseq-bioskryb's own model-fit failures, not a missing binary).
# -----------------------------------------------------------------------
FROM condaforge/miniforge3:24.9.2-0

RUN mamba create -y -n resolveome-align -c conda-forge -c bioconda \
      python=3.11 \
      bwa=0.7.18 \
      bwa-mem2=2.2.1 \
      samtools=1.20 \
      bcftools=1.20 \
      htslib=1.20 \
      tabix=1.11 \
      picard=3.2.0 \
      fastp=0.20.1 \
      fastqc=0.12.1 \
      multiqc=1.21 \
      cnvkit=0.9.12 \
      snpeff=5.1 \
      gatk4=4.6.2.0 \
      bedtools \
      openjdk=17 \
      pandas \
      matplotlib \
      numpy \
    && mamba clean -afy

# Put the env's bin/ directly on PATH so tasks don't need `conda activate`
# (Sprocket tasks run as plain `bash` scripts, not login shells).
ENV PATH="/opt/conda/envs/resolveome-align/bin:${PATH}"

# ---- bam2mr / preseq-bioskryb: see the gap note above ----
# COPY bin/bam2mr bin/preseq-bioskryb /opt/conda/envs/resolveome-align/bin/
# RUN chmod +x /opt/conda/envs/resolveome-align/bin/bam2mr \
#               /opt/conda/envs/resolveome-align/bin/preseq-bioskryb

# Sanity check at build time that the environment's tools are reachable.
RUN bwa-mem2 version && samtools --version | head -1 && bcftools --version | head -1 \
    && picard CollectWgsMetrics --version 2>&1 | head -1 && fastp --version 2>&1 \
    && fastqc --version && multiqc --version && cnvkit.py version && snpEff -version \
    && gatk --version
