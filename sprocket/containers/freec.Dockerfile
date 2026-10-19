# ResolveOME DNA workflow -- "freec" container for Sprocket.
#
# Bundles the same tool stack as ../../setup/setup.sh's `freec` conda
# environment, as a container image (see resolveome-align.Dockerfile's
# header comment for why this is needed under Sprocket's execution model).
#
# Used by: loh_ploidy.wdl (ControlFreec).
FROM condaforge/miniforge3:24.9.2-0

RUN mamba create -y -n freec -c conda-forge -c bioconda \
      control-freec=11.6 \
      samtools=1.24 \
      bedtools=2.31 \
      bioconductor-rsamtools \
    && mamba clean -afy

ENV PATH="/opt/conda/envs/freec/bin:${PATH}"

RUN freec -version 2>&1 | head -1 && samtools --version | head -1 && bedtools --version
