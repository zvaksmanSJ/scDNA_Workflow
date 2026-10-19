# ResolveOME DNA workflow -- "mutsig" container for Sprocket.
#
# Bundles the same tool stack as ../../setup/setup.sh's `mutsig` conda
# environment, as a container image (see resolveome-align.Dockerfile's
# header comment for why this is needed under Sprocket's execution model).
#
# Used by: extra_qc_dashboard.wdl (MutationalSignatures) ONLY. For the
# combined python3+bcftools+R+Sequoia image used by lineage.wdl, see
# lineage.Dockerfile instead -- that task needs more than this image alone
# provides.
FROM condaforge/miniforge3:24.9.2-0

RUN mamba create -y -n mutsig -c conda-forge -c bioconda \
      r-base=4.5 \
      bioconductor-mutationalpatterns \
      bioconductor-bsgenome.hsapiens.ucsc.hg38 \
      bioconductor-variantannotation \
    && mamba clean -afy

ENV PATH="/opt/conda/envs/mutsig/bin:${PATH}"

RUN Rscript -e 'library(MutationalPatterns); library(BSgenome.Hsapiens.UCSC.hg38); library(VariantAnnotation)'
