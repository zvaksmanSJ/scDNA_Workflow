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
# bam2mr / preseq-bioskryb: built from source here, not a private binary
# -----------------------------------------------------------------------
# align_qc.wdl's AlignTrimQC task calls `bam2mr` and `preseq-bioskryb` for
# its library-complexity step. Despite the pipeline-internal "-bioskryb"
# naming (used only to distinguish this exact build from the newer,
# CLI-incompatible bioconda `preseq` 3.x, which renamed `bam2mr` to
# `to-mr`), these are preseq v2.0.3's own `bam2mr`/`gc_extrap` commands
# (GPLv3, fully public -- github.com/smithlabcode/preseq) -- not a
# private/proprietary BioSkryb build, and not something you need to
# supply yourself. This mirrors exactly what ../../setup/setup_parity_
# tools.sh already does for the Cromwell path's `resolveome-align` conda
# environment.
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
      make \
      gxx \
      curl \
    && mamba clean -afy

# Put the env's bin/ directly on PATH so tasks don't need `conda activate`
# (Sprocket tasks run as plain `bash` scripts, not login shells).
ENV PATH="/opt/conda/envs/resolveome-align/bin:${PATH}"

# ---- bam2mr / preseq-bioskryb: built from the official public release ----
RUN curl --fail --location --retry 3 \
      https://github.com/smithlabcode/preseq/releases/download/v2.0.3/preseq_v2.0.3.tar.bz2 \
      --output /tmp/preseq_v2.0.3.tar.bz2 \
    && tar -xjf /tmp/preseq_v2.0.3.tar.bz2 -C /tmp \
    && cd /tmp/preseq \
    && make all CXXFLAGS="-O2 -Wall -fno-common" CFLAGS="-O2 -fcommon" \
    && install -m 755 preseq /opt/conda/envs/resolveome-align/bin/preseq-bioskryb \
    && install -m 755 bam2mr /opt/conda/envs/resolveome-align/bin/bam2mr \
    && cd / && rm -rf /tmp/preseq /tmp/preseq_v2.0.3.tar.bz2

# Sanity check at build time that the environment's tools are reachable.
RUN bwa-mem2 version && samtools --version | head -1 && bcftools --version | head -1 \
    && picard CollectWgsMetrics --version 2>&1 | head -1 && fastp --version 2>&1 \
    && fastqc --version && multiqc --version && cnvkit.py version && snpEff -version \
    && gatk --version && bam2mr 2>&1 | head -1 && preseq-bioskryb 2>&1 | head -1
