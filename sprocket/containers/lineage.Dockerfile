# ResolveOME DNA workflow -- "lineage" container for Sprocket.
#
# Combines what the Cromwell version of this workflow runs as TWO separate
# conda environments in one task (`resolveome-align` for python3/bcftools,
# then `mutsig`'s bin/ prepended to PATH for R/Sequoia) into a single image,
# because Sprocket's execution model requires every task to run inside
# exactly one container (see ../README.md).
#
# Used by: lineage.wdl (BuildLineageTree, VerifyLineage).
#
# Pre-installs Sequoia's (github.com/TimCoorens/Sequoia, GPLv3)
# build_phylogeny.R dependencies at BUILD time, rather than relying on that
# script's own "install if missing" fallback at run time (which would need
# outbound internet access from the LSF compute node on every single task
# invocation, and would slow down -- and risk flaking -- every run).
# Dependency list taken directly from build_phylogeny.R's own
# `cran_packages <- c(...)` declaration plus its `require("treemut")` /
# `install_git(...)` call for the one GitHub-only package.
FROM condaforge/miniforge3:24.9.2-0

RUN mamba create -y -n lineage -c conda-forge -c bioconda \
      python=3.11 \
      bcftools=1.20 \
      htslib=1.20 \
      r-base=4.5 \
      r-ggplot2 \
      r-ape \
      r-seqinr \
      r-stringr \
      r-data.table \
      r-tidyr \
      r-dplyr \
      r-vgam \
      r-mass \
      r-devtools \
      r-optparse \
      bioconductor-biocmanager \
    && mamba clean -afy

ENV PATH="/opt/conda/envs/lineage/bin:${PATH}"

# Sequoia's one GitHub-only dependency (not on CRAN/Bioconductor): treemut.
# Installed here, at build time, so build_phylogeny.R's own
# `if (!require("treemut")) install_git(...)` check is already satisfied
# and never needs to reach out to github.com at task-run time.
RUN Rscript -e 'devtools::install_git("https://github.com/NickWilliamsSanger/treemut")'

RUN Rscript -e 'library(ggplot2); library(ape); library(seqinr); library(stringr); \
                library(data.table); library(tidyr); library(dplyr); library(VGAM); \
                library(MASS); library(devtools); library(optparse); library(treemut)' \
    && bcftools --version | head -1
