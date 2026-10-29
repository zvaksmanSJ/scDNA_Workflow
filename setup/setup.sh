#!/usr/bin/env bash
# ResolveOME single-cell DNA workflow -- general-purpose setup script.
#
# Installs every conda environment, reference file, and container this WDL
# workflow (../main.wdl) needs, using ONLY free/open tools and resources:
#   - stock Docker Hub containers (google/deepvariant, ensemblorg/ensembl-vep)
#     pulled with apptainer/singularity, NOT any proprietary/licensed image
#   - Ensembl's public GRCh38 FASTA + GATK's public legacy known-sites bundle
#   - Beagle 5.5 (Browning lab, GPL) + its public PLINK genetic maps
#   - the 1000 Genomes GRCh38 high-coverage phased panel (public, EBI FTP)
#   - Sequoia (TimCoorens/Sequoia, GPLv3, standalone GitHub script) + MPBoot
#     (GPLv2, standalone GitHub binary release) -- NOT the licensed BioSkryb
#     "basej-somatic" apptainer container used for this project's own run
#
# This script is idempotent (every step is skip-if-already-present) and safe
# to re-run. It does NOT touch, move, or re-derive anything under
# resolveome_align/refs or resolveome_align/env (those are left exactly as
# they are for this project's own real cohort run) -- it installs into
# ../wdl_runtime/ by default (a separate, WDL-only location at workflow_DNA
# root), or any other BASE_DIR of your choosing.
#
# Conda env names are prefixed "wdl-" (e.g. "wdl-resolveome-align") for the
# same reason: this project's OWN production pipeline separately creates an
# env literally named "resolveome-align" (via setup_resolveome_alignment.sh,
# not this script) with a different package set (e.g. it includes `delly`,
# which this script's env does not). Conda envs are a per-user/account
# resource, not per-project -- on any account where that production env
# already exists, this script's own update path
# (`conda env update -n <name> --prune`) would otherwise silently prune
# packages out of that OTHER, unrelated, already-in-use environment the
# first time this script ran. The "wdl-" prefix guarantees no collision
# regardless of what else is already installed on the account running this.
#
# Usage (from anywhere -- the default BASE_DIR is resolved relative to this
# script's own location, not your current directory):
#   ./setup.sh [BASE_DIR] [THREADS] [ENSEMBL_RELEASE] [VEP_CACHE_VERSION]
#
# Defaults:
#   BASE_DIR           = ../wdl_runtime   (i.e. workflow_DNA/wdl_runtime/)
#   THREADS            = 16
#   ENSEMBL_RELEASE     = 111
#   VEP_CACHE_VERSION   = 111
#
# Individual sections can be skipped with SKIP_<SECTION>=1 env vars, e.g.:
#   SKIP_CONDA=1 SKIP_CONTAINERS=1 ./setup.sh
# Recognized: SKIP_CONDA, SKIP_REFERENCE, SKIP_KNOWN_SITES, SKIP_CONTAINERS,
#             SKIP_PHASING, SKIP_FREEC_REFS, SKIP_PARITY_TOOLS, SKIP_LINEAGE
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${1:-${SCRIPT_DIR}/../wdl_runtime}"
THREADS="${2:-16}"
ENSEMBL_RELEASE="${3:-111}"
VEP_CACHE_VERSION="${4:-111}"

log(){ echo "[INFO] $*"; }
warn(){ echo "[WARN] $*" >&2; }
err(){ echo "[ERROR] $*" >&2; }
have(){ command -v "$1" >/dev/null 2>&1; }
download_if_missing() {
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then echo "[OK] exists: $out"; return 0; fi
  log "downloading: $url"
  mkdir -p "$(dirname "$out")"
  # Large downloads (notably the ~25GB VEP cache tarball) have been observed
  # to reliably drop mid-transfer with "curl: (56) OpenSSL SSL_read:
  # unexpected eof while reading" a few GB in, on this network path --
  # reproduced twice in a row at a similar point, so this is a real,
  # repeatable connection behavior here, not a one-off fluke. curl's own
  # --retry doesn't reliably cover a mid-transfer drop like this, and
  # without -C - (resume), a plain retry restarts the whole download from
  # byte 0 -- so a connection that reliably drops at the same multi-GB
  # point would never successfully finish. -C - resumes "${out}.part" from
  # where it left off on every attempt (curl creates it fresh if absent),
  # and the outer loop gives it repeated fresh attempts beyond curl's own
  # internal --retry, each one only needing to cover the remaining bytes.
  #
  # Separately: a .gz/.vcf.gz output that curl reports as fully downloaded
  # can still be silently corrupt (observed in practice -- a background
  # phasing-panel download was left truncated-but-non-empty by an LSF
  # walltime kill of the whole script, yet still ended up at the final
  # "$out" name and passed this function's own `[[ -s "$out" ]]` skip-check
  # on the next run; it only surfaced later as a tabix indexing crash).
  # For .gz outputs, verify with `gzip -t` before trusting curl's own
  # success and finalizing via mv -- on failure, discard "${out}.part"
  # entirely (not just resume) since corruption may be mid-stream, not
  # only at the tail, so resuming from the current byte offset would not
  # necessarily fix it.
  local attempt
  for attempt in 1 2 3 4 5 6 7 8; do
    if curl -fL --retry 3 --retry-all-errors -C - -o "${out}.part" "$url"; then
      if [[ "$out" == *.gz ]] && ! gzip -t "${out}.part" 2>/dev/null; then
        warn "downloaded file failed gzip integrity check (possible silent corruption) -- discarding and retrying from scratch: $url"
        rm -f "${out}.part"
        sleep 15
        continue
      fi
      mv "${out}.part" "$out"
      return 0
    fi
    warn "download attempt ${attempt}/8 failed (will resume from $(du -h "${out}.part" 2>/dev/null | cut -f1 || echo 0) if retried): $url"
    sleep 15
  done
  err "download failed after 8 attempts (even with resume): $url"
  return 1
}

log "BASE_DIR=${BASE_DIR}  THREADS=${THREADS}  ENSEMBL_RELEASE=${ENSEMBL_RELEASE}  VEP_CACHE_VERSION=${VEP_CACHE_VERSION}"

have conda || { err "conda not found on PATH"; exit 1; }
APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
[[ -n "$APPTAINER_BIN" ]] || { err "neither apptainer nor singularity found on PATH (required for DeepVariant/VEP/Manta containers)"; exit 1; }
log "Using container runtime: ${APPTAINER_BIN}"

mkdir -p "${BASE_DIR}"/{env,refs/{genome,annotation/vep_cache,gatk_known_sites,freec_chrom_fasta,phasing_panel,genetic_maps,freec_snp_positions},containers,tools/{beagle,sequoia,mpboot},logs}

source "$(conda info --base)/etc/profile.d/conda.sh"

# ===========================================================================
# 1. Conda environments
# ===========================================================================
if [[ "${SKIP_CONDA:-0}" != "1" ]]; then
  log "=== Conda environments ==="

  cat > "${BASE_DIR}/env/wdl-resolveome-align.yml" <<YAML
name: wdl-resolveome-align
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.11
  - bwa=0.7.18
  - bwa-mem2=2.2.1
  - samtools=1.20
  - bcftools=1.20
  - htslib=1.20
  - tabix=1.11
  - picard=3.2.0
  - fastqc=0.12.1
  - multiqc=1.21
  # cnvkit 0.9.12 requires pomegranate<=0.14.9, which cannot be installed
  # alongside python=3.11 (conda solver fails with
  # LibMambaUnsatisfiableError) -- 0.9.14 drops the pomegranate dependency
  # entirely, so it's pinned here instead, not just "latest".
  - cnvkit=0.9.14
  - snpeff=5.1
  - gatk4=4.6.2.0
  - bedtools
  - pandas
  - matplotlib
  - numpy
YAML

  cat > "${BASE_DIR}/env/freec.yml" <<YAML
name: freec
channels: [conda-forge, bioconda, defaults]
dependencies:
  - control-freec=11.6
  - samtools=1.24
  - bedtools=2.31
  - bioconductor-rsamtools
YAML

  cat > "${BASE_DIR}/env/manta_env.yml" <<YAML
name: manta_env
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=2.7
  - manta=1.6.0
YAML

  cat > "${BASE_DIR}/env/mutsig.yml" <<YAML
name: mutsig
channels: [conda-forge, bioconda, defaults]
dependencies:
  - r-base=4.5
  - bioconductor-mutationalpatterns
  - bioconductor-bsgenome.hsapiens.ucsc.hg38
  - bioconductor-variantannotation
YAML

  for envname in wdl-resolveome-align freec manta_env mutsig; do
    if conda env list | awk '{print $1}' | grep -qx "${envname}"; then
      log "conda env '${envname}' already exists, updating (--prune)"
      conda env update -n "${envname}" -f "${BASE_DIR}/env/${envname}.yml" --prune
    else
      log "creating conda env '${envname}'"
      conda env create -f "${BASE_DIR}/env/${envname}.yml"
    fi
  done
else
  log "SKIP_CONDA=1 -- skipping conda environment creation"
fi

# ===========================================================================
# 2. Reference genome (Ensembl GRCh38 primary assembly) + indices
# ===========================================================================
if [[ "${SKIP_REFERENCE:-0}" != "1" ]]; then
  log "=== Reference genome ==="
  conda activate wdl-resolveome-align

  GENOME_DIR="${BASE_DIR}/refs/genome"
  FA_GZ="${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
  FA="${GENOME_DIR}/GRCh38.primary.fa"

  download_if_missing \
    "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    "${FA_GZ}"

  if [[ ! -s "${FA}" ]]; then
    log "Decompressing reference FASTA..."
    gunzip -c "${FA_GZ}" > "${FA}"
  else
    echo "[OK] exists: ${FA}"
  fi

  [[ -s "${FA}.fai" ]] || samtools faidx "${FA}"
  [[ -s "${GENOME_DIR}/GRCh38.primary.dict" ]] || \
    picard CreateSequenceDictionary R="${FA}" O="${GENOME_DIR}/GRCh38.primary.dict"
  [[ -s "${FA}.0123" && -s "${FA}.bwt.2bit.64" ]] || bwa-mem2 index "${FA}"

  log "Reference genome ready: ${FA}"

  # ---- per-chromosome FASTA split (chr1-22+X ONLY) + filtered chrLenFile ----
  # Control-FREEC requires one FASTA file per contig listed in its chrLenFile,
  # and refuses decoy/alt/unplaced scaffolds -- so we split only the 23
  # canonical autosomes+X contigs, matching their names in the .fai exactly.
  CHROM_FASTA_DIR="${BASE_DIR}/refs/freec_chrom_fasta"
  CHR_LEN_FILE="${BASE_DIR}/refs/freec_chrLen.fai"
  if [[ ! -s "${CHR_LEN_FILE}" ]]; then
    log "Splitting per-chromosome FASTAs (chr1-22+X) for Control-FREEC..."
    : > "${CHR_LEN_FILE}.tmp"
    for c in $(seq 1 22) X; do
      out="${CHROM_FASTA_DIR}/${c}.fa"
      if [[ ! -s "${out}" ]]; then
        samtools faidx "${FA}" "${c}" > "${out}"
        samtools faidx "${out}"
      fi
      grep -P "^${c}\t" "${FA}.fai" >> "${CHR_LEN_FILE}.tmp"
    done
    mv "${CHR_LEN_FILE}.tmp" "${CHR_LEN_FILE}"
  else
    echo "[OK] exists: ${CHR_LEN_FILE}"
  fi

  # main.wdl takes chrom_fasta_dir_tar as a single tarball input -- repackage
  # (loh_ploidy.wdl untars this fresh per cell, so the tar must contain the
  # per-chromosome .fa + .fa.fai files directly, not the freec_chrom_fasta/
  # directory itself as a path component)
  CHROM_FASTA_TAR="${BASE_DIR}/refs/freec_chrom_fasta.tar"
  [[ -s "${CHROM_FASTA_TAR}" ]] || tar -cf "${CHROM_FASTA_TAR}" -C "${CHROM_FASTA_DIR}" .
else
  log "SKIP_REFERENCE=1 -- skipping reference genome download/indexing"
fi

# ===========================================================================
# 3. GATK known sites (dbSNP, Mills indels) -- used by the alignment/BQSR step
# ===========================================================================
if [[ "${SKIP_KNOWN_SITES:-0}" != "1" ]]; then
  log "=== GATK known sites ==="
  conda activate wdl-resolveome-align
  KNOWN="${BASE_DIR}/refs/gatk_known_sites"
  set +e
  download_if_missing "https://storage.googleapis.com/gatk-legacy-bundles/b38/dbsnp_146.hg38.vcf.gz" "${KNOWN}/dbsnp_146.hg38.vcf.gz"
  download_if_missing "https://storage.googleapis.com/gatk-legacy-bundles/b38/dbsnp_146.hg38.vcf.gz.tbi" "${KNOWN}/dbsnp_146.hg38.vcf.gz.tbi"
  download_if_missing "https://storage.googleapis.com/gatk-legacy-bundles/b38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "${KNOWN}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  download_if_missing "https://storage.googleapis.com/gatk-legacy-bundles/b38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" "${KNOWN}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  set -e
  for v in "${KNOWN}"/*.vcf.gz; do
    [[ -e "$v" ]] || continue
    [[ -s "${v}.tbi" ]] || tabix -f -p vcf "$v" || true
  done
else
  log "SKIP_KNOWN_SITES=1 -- skipping GATK known-sites download"
fi

# ===========================================================================
# 4. Containers (free, stock images only -- NOT BioSkryb's licensed images)
# ===========================================================================
if [[ "${SKIP_CONTAINERS:-0}" != "1" ]]; then
  log "=== Containers ==="
  DV_SIF="${BASE_DIR}/containers/deepvariant.sif"
  VEP_SIF="${BASE_DIR}/containers/vep.r${ENSEMBL_RELEASE}.sif"
  MANTA_SIF="${BASE_DIR}/containers/manta.sif"

  if [[ -s "${DV_SIF}" ]]; then
    echo "[OK] exists: ${DV_SIF}"
  else
    log "Pulling stock DeepVariant container -> ${DV_SIF}"
    "${APPTAINER_BIN}" pull "${DV_SIF}" docker://google/deepvariant:1.6.1
  fi

  if [[ -s "${VEP_SIF}" ]]; then
    echo "[OK] exists: ${VEP_SIF}"
  else
    log "Pulling stock VEP container -> ${VEP_SIF}"
    "${APPTAINER_BIN}" pull "${VEP_SIF}" "docker://ensemblorg/ensembl-vep:release_${ENSEMBL_RELEASE}.0"
  fi

  if [[ -s "${MANTA_SIF}" ]]; then
    echo "[OK] exists: ${MANTA_SIF}"
  else
    log "Pulling stock Manta container -> ${MANTA_SIF}"
    "${APPTAINER_BIN}" pull "${MANTA_SIF}" docker://halllab/manta:latest || \
      warn "Manta container pull failed -- you can instead use the manta_env conda env directly (see cnv_sv.wdl)"
  fi

  # VEP cache (needed by the VEP container at runtime; downloaded separately
  # since the container image itself does not bundle species caches)
  VEP_CACHE_DIR="${BASE_DIR}/refs/annotation/vep_cache"
  VEP_CACHE_TAR="${VEP_CACHE_DIR}/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz"
  if [[ -d "${VEP_CACHE_DIR}/homo_sapiens" ]]; then
    echo "[OK] exists: ${VEP_CACHE_DIR}/homo_sapiens"
  else
    download_if_missing \
      "https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz" \
      "${VEP_CACHE_TAR}"
    tar -xzf "${VEP_CACHE_TAR}" -C "${VEP_CACHE_DIR}"
  fi
  # main.wdl takes vep_cache_dir_tar as a single tarball input -- repackage
  VEP_CACHE_REPACK="${VEP_CACHE_DIR}/vep_cache.tar"
  [[ -s "${VEP_CACHE_REPACK}" ]] || tar -cf "${VEP_CACHE_REPACK}" -C "${VEP_CACHE_DIR}" homo_sapiens

  # SnpEff database (downloaded through the conda env's own snpEff binary)
  conda activate wdl-resolveome-align
  SNPEFF_DB_NAME="GRCh38.${ENSEMBL_RELEASE}"
  if ! snpEff databases 2>/dev/null | grep -q "${SNPEFF_DB_NAME}"; then
    log "Downloading snpEff DB: ${SNPEFF_DB_NAME} (best effort)"
    snpEff download -v "${SNPEFF_DB_NAME}" || warn "snpEff DB download failed; you may need a different snpeff_db value"
  fi
else
  log "SKIP_CONTAINERS=1 -- skipping container pulls / VEP cache / snpEff DB"
fi

# ===========================================================================
# 5. Phasing resources (Beagle 5.5 + PLINK genetic maps + 1000G GRCh38 panel)
# ===========================================================================
if [[ "${SKIP_PHASING:-0}" != "1" ]]; then
  log "=== Phasing resources (Beagle) ==="
  BEAGLE_JAR="${BASE_DIR}/tools/beagle/beagle.27Feb25.75f.jar"
  download_if_missing "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar" "${BEAGLE_JAR}"

  MAP_ZIP="${BASE_DIR}/refs/genetic_maps/plink.GRCh38.map.zip"
  download_if_missing "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" "${MAP_ZIP}"
  if [[ ! -d "${BASE_DIR}/refs/genetic_maps/chr_in_chrom_field" ]]; then
    (cd "${BASE_DIR}/refs/genetic_maps" && unzip -o -q "$(basename "${MAP_ZIP}")")
  fi

  PANEL_DIR="${BASE_DIR}/refs/phasing_panel"
  PANEL_BASE="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
  for c in $(seq 1 22) X; do
    if [[ "$c" == "X" ]]; then
      fname="CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
    else
      fname="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${c}.filtered.shapeit2-duohmm-phased.vcf.gz"
    fi
    download_if_missing "${PANEL_BASE}/${fname}" "${PANEL_DIR}/${fname}" &
    # throttle to 4 concurrent downloads
    if (( $(jobs -r | wc -l) >= 4 )); then wait -n; fi
  done
  wait
  conda activate wdl-resolveome-align
  for f in "${PANEL_DIR}"/*.vcf.gz; do
    [[ -s "${f}.tbi" ]] || tabix -f -p vcf "$f"
  done
  log "Phasing resources ready: ${BEAGLE_JAR}, genetic maps, ${PANEL_DIR}"
else
  log "SKIP_PHASING=1 -- skipping Beagle/genetic-map/1000G-panel download"
fi

# ===========================================================================
# 6. Control-FREEC SNP positions (site-only VCF derived from the same 1000G
#    panel downloaded above -- reused, not re-downloaded)
# ===========================================================================
if [[ "${SKIP_FREEC_REFS:-0}" != "1" ]]; then
  log "=== Control-FREEC SNP positions ==="
  conda activate wdl-resolveome-align
  PANEL_DIR="${BASE_DIR}/refs/phasing_panel"
  OUT="${BASE_DIR}/refs/freec_snp_positions"
  for c in $(seq 1 22) X; do
    out_vcf="${OUT}/chr${c}.sites.vcf.gz"
    [[ -s "${out_vcf}" && -s "${out_vcf}.tbi" ]] && { echo "[OK] exists: ${out_vcf}"; continue; }
    if [[ "$c" == "X" ]]; then
      fname="CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
    else
      fname="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${c}.filtered.shapeit2-duohmm-phased.vcf.gz"
    fi
    bcftools view -G -m2 -M2 -v snps "${PANEL_DIR}/${fname}" \
      | bcftools annotate --rename-chrs <(echo -e "chr${c}\t${c}") -Oz -o "${out_vcf}"
    bcftools index -t -f "${out_vcf}"
  done
  ALL_SITES="${OUT}/all_chrom.sites.vcf.gz"
  if [[ ! -s "${ALL_SITES}" ]]; then
    bcftools concat -Oz -o "${ALL_SITES}" $(for c in $(seq 1 22) X; do echo "${OUT}/chr${c}.sites.vcf.gz"; done)
    bcftools index -t -f "${ALL_SITES}"
  fi
  log "Control-FREEC SNP positions ready: ${ALL_SITES}"
else
  log "SKIP_FREEC_REFS=1 -- skipping Control-FREEC SNP-positions derivation"
fi

# ===========================================================================
# 7. Library-complexity QC: preseq v2.0.3's bam2mr/gc_extrap pair (GPLv3,
#    fully public -- https://github.com/smithlabcode/preseq -- NOT a
#    private/proprietary binary, despite the pipeline-internal "-bioskryb"
#    naming convention used to distinguish it from the newer, CLI-
#    incompatible bioconda preseq 3.x also installed in this same env).
#    align_qc.wdl's AlignTrimQC task calls these two commands directly.
# ===========================================================================
if [[ "${SKIP_PARITY_TOOLS:-0}" != "1" ]]; then
  log "=== Library-complexity QC tools (preseq 2.0.3 bam2mr/gc_extrap) ==="
  conda env update -n wdl-resolveome-align -f "${SCRIPT_DIR}/parity_tools.environment.yml"
  TARGET_DIR="${BASE_DIR}/tools/preseq-2.0.3" "${SCRIPT_DIR}/setup_parity_tools.sh"
else
  log "SKIP_PARITY_TOOLS=1 -- skipping preseq 2.0.3 bam2mr/gc_extrap install"
fi

# ===========================================================================
# 8. Lineage tools: Sequoia (standalone, GPLv3) + MPBoot (standalone, GPLv2)
#    -- NOT the licensed BioSkryb apptainer container used for the real run
# ===========================================================================
if [[ "${SKIP_LINEAGE:-0}" != "1" ]]; then
  log "=== Lineage tools (Sequoia + MPBoot) ==="
  SEQUOIA_DIR="${BASE_DIR}/tools/sequoia"
  if [[ ! -d "${SEQUOIA_DIR}/Sequoia" ]]; then
    git clone --depth 1 https://github.com/TimCoorens/Sequoia.git "${SEQUOIA_DIR}/Sequoia"
  else
    echo "[OK] exists: ${SEQUOIA_DIR}/Sequoia"
  fi

  MPBOOT_DIR="${BASE_DIR}/tools/mpboot"
  MPBOOT_TGZ="${MPBOOT_DIR}/mpboot-1.1.0-linux.tar.gz"
  if [[ ! -s "${MPBOOT_DIR}/mpboot-1.1.0-Linux/bin/mpboot" ]]; then
    download_if_missing \
      "https://github.com/diepthihoang/mpboot/releases/download/v1.1.0/mpboot-1.1.0-linux.tar.gz" \
      "${MPBOOT_TGZ}"
    tar -xzf "${MPBOOT_TGZ}" -C "${MPBOOT_DIR}"
  fi
  log "Lineage tools ready: ${SEQUOIA_DIR}/Sequoia, ${MPBOOT_DIR}"
  log "Use --sequoia-script \${SEQUOIA_DIR}/Sequoia/build_phylogeny.R --mpboot-path <mpboot bin dir> with lineage_track.py/verify_lineage.py"
else
  log "SKIP_LINEAGE=1 -- skipping Sequoia/MPBoot install"
fi

# ===========================================================================
# Final summary
# ===========================================================================
log "=== Setup summary ==="
conda env list | grep -E "wdl-resolveome-align|freec|manta_env|mutsig" || true
for f in "${BASE_DIR}/containers/deepvariant.sif" "${BASE_DIR}/containers/vep.r${ENSEMBL_RELEASE}.sif" \
         "${BASE_DIR}/refs/genome/GRCh38.primary.fa" "${BASE_DIR}/tools/beagle/beagle.27Feb25.75f.jar" \
         "${BASE_DIR}/tools/preseq-2.0.3/bin/bam2mr" "${BASE_DIR}/tools/preseq-2.0.3/bin/preseq" \
         "${BASE_DIR}/tools/sequoia/Sequoia" "${BASE_DIR}/tools/mpboot"; do
  [[ -e "$f" ]] && echo "[OK] $f" || echo "[MISS] $f"
done
log "Setup complete. Point wdl/inputs.json at the files under ${BASE_DIR}."
