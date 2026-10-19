#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User args
# -----------------------------
BASE_DIR="${1:-$PWD/resolveome_align}"
ENV_NAME="${2:-resolveome-align}"
THREADS="${3:-16}"
VEP_SPECIES="${4:-homo_sapiens}"
VEP_ASSEMBLY="${5:-GRCh38}"
ENSEMBL_RELEASE="${6:-112}"

# LSF controls for heavy work
QUEUE="${QUEUE:-standard}"
NCORES="${NCORES:-8}"
MEM_MB="${MEM_MB:-64000}"
WALL="${WALL:-08:00}"

# mode:
#   auto   -> if not already in LSF job, submit heavy part to LSF and exit
#   local  -> run heavy part directly in current shell
MODE="${MODE:-auto}"

# -----------------------------
# Helpers
# -----------------------------
have(){ command -v "$1" >/dev/null 2>&1; }
log(){ echo "[INFO] $*"; }
warn(){ echo "[WARN] $*"; }
err(){ echo "[ERROR] $*" >&2; }

download_if_missing() {
  local url="$1"; local out="$2"
  if [[ -s "$out" ]]; then
    echo "[OK] exists: $out"
    return 0
  fi
  log "downloading: $url"
  wget -O "$out" "$url"
}

# -----------------------------
# Base checks (safe on login)
# -----------------------------
log "BASE_DIR=${BASE_DIR}"
log "ENV_NAME=${ENV_NAME}"
log "THREADS=${THREADS}"
log "VEP_SPECIES=${VEP_SPECIES}"
log "VEP_ASSEMBLY=${VEP_ASSEMBLY}"
log "ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
log "MODE=${MODE}"

have conda || { err "conda not found"; exit 1; }
have bsub || { err "bsub not found"; exit 1; }

# singularity-only requirement
have singularity || { err "singularity not found (required)"; exit 1; }

mkdir -p "${BASE_DIR}"/{env,bin,logs,containers,refs/{genome,annotation,indexes/{bwa,bwa-mem2,samtools,picard},variants/{common,gatk/known_sites,manta,delly,lumpy,deepvariant},cnv,sv},fastq,results}

# -----------------------------
# Build heavy worker script
# -----------------------------
HEAVY_SCRIPT="${BASE_DIR}/bin/.setup_heavy_worker.sh"
cat > "${HEAVY_SCRIPT}" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$1"
ENV_NAME="$2"
THREADS="$3"
VEP_SPECIES="$4"
VEP_ASSEMBLY="$5"
ENSEMBL_RELEASE="$6"

have(){ command -v "$1" >/dev/null 2>&1; }
log(){ echo "[INFO] $*"; }

cd "${BASE_DIR}"

cat > env/environment.yml <<YAML
name: ${ENV_NAME}
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - bwa=0.7.18
  - bwa-mem2=2.2.1
  - samtools=1.20
  - bcftools=1.20
  - htslib=1.20
  - tabix
  - picard=3.2.0
  - fastqc=0.12.1
  - multiqc=1.21
  - seqkit=2.8.2
  - pigz=2.8
  - parallel
  - openjdk=17
  - curl
  - wget
  - jq
  - gatk4
  - delly
  - cnvkit
  - snpeff
YAML

log "Creating/updating conda env..."
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  conda env update -n "${ENV_NAME}" -f env/environment.yml --prune
else
  conda env create -f env/environment.yml
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# tool check
log "Tool check (conda env)"
for t in bwa bwa-mem2 samtools bcftools picard fastqc multiqc gatk delly cnvkit.py snpEff; do
  if have "$t"; then
    echo "  [OK] $t -> $(command -v "$t")"
  else
    echo "  [MISS] $t"
  fi
done

# refs
REF_DIR="${BASE_DIR}/refs"
GENOME_DIR="${REF_DIR}/genome"
ANN_DIR="${REF_DIR}/annotation"
IDX_DIR="${REF_DIR}/indexes"
mkdir -p "${GENOME_DIR}" "${ANN_DIR}" "${IDX_DIR}"/{bwa,bwa-mem2,samtools,picard}

download_if_missing() {
  local url="$1"; local out="$2"
  if [[ -s "$out" ]]; then
    echo "[OK] exists: $out"
    return 0
  fi
  log "downloading: $url"
  wget -O "$out" "$url"
}

FA_GZ="${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FA="${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_FA="${GENOME_DIR}/GRCh38.primary.fa"
DICT="${GENOME_DIR}/GRCh38.primary.dict"

download_if_missing \
  "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  "${FA_GZ}"

if [[ ! -s "${FA}" ]]; then
  log "Decompressing reference FASTA..."
  pigz -dc "${FA_GZ}" > "${FA}"
else
  echo "[OK] exists: ${FA}"
fi

ln -sf "${FA}" "${REF_FA}"

[[ -s "${REF_FA}.fai" ]] || samtools faidx "${REF_FA}"
[[ -s "${DICT}" ]] || picard CreateSequenceDictionary R="${REF_FA}" O="${DICT}"
[[ -s "${REF_FA}.bwt" ]] || bwa index "${REF_FA}"

if [[ ! -s "${REF_FA}.0123" && ! -s "${REF_FA}.bwt.2bit.64" ]]; then
  bwa-mem2 index "${REF_FA}"
fi

ln -sf "${REF_FA}" "${IDX_DIR}/bwa/GRCh38.primary.fa"
ln -sf "${REF_FA}" "${IDX_DIR}/bwa-mem2/GRCh38.primary.fa"
ln -sf "${REF_FA}.fai" "${IDX_DIR}/samtools/GRCh38.primary.fa.fai"
ln -sf "${DICT}" "${IDX_DIR}/picard/GRCh38.primary.dict"

# variant refs
ln -sf "${REF_FA}" "${REF_DIR}/variants/common/GRCh38.primary.fa"
ln -sf "${REF_FA}.fai" "${REF_DIR}/variants/common/GRCh38.primary.fa.fai"
ln -sf "${DICT}" "${REF_DIR}/variants/common/GRCh38.primary.dict"
ln -sf "${REF_FA}" "${REF_DIR}/variants/delly/reference.fa"
ln -sf "${REF_FA}.fai" "${REF_DIR}/variants/delly/reference.fa.fai"
ln -sf "${REF_FA}" "${REF_DIR}/variants/manta/reference.fa"
ln -sf "${REF_FA}.fai" "${REF_DIR}/variants/manta/reference.fa.fai"
ln -sf "${REF_FA}" "${REF_DIR}/variants/lumpy/reference.fa"
ln -sf "${REF_FA}.fai" "${REF_DIR}/variants/lumpy/reference.fa.fai"
ln -sf "${REF_FA}" "${REF_DIR}/variants/deepvariant/reference.fa"
ln -sf "${REF_FA}.fai" "${REF_DIR}/variants/deepvariant/reference.fa.fai"

# known sites (best effort only)
KNOWN="${REF_DIR}/variants/gatk/known_sites"
mkdir -p "${KNOWN}"
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

# containers (singularity only)
DV_SIF="${BASE_DIR}/containers/deepvariant.sif"
VEP_SIF="${BASE_DIR}/containers/vep.r${ENSEMBL_RELEASE}.sif"

if [[ -s "${DV_SIF}" ]]; then
  echo "[OK] DeepVariant container exists: ${DV_SIF}"
else
  log "Pulling DeepVariant container -> ${DV_SIF}"
  singularity pull "${DV_SIF}" docker://google/deepvariant:1.6.1 || echo "[WARN] DeepVariant pull failed"
fi

if [[ -s "${VEP_SIF}" ]]; then
  echo "[OK] VEP container exists: ${VEP_SIF}"
else
  log "Pulling VEP container -> ${VEP_SIF}"
  singularity pull "${VEP_SIF}" "docker://ensemblorg/ensembl-vep:release_${ENSEMBL_RELEASE}.0" || echo "[WARN] VEP pull failed"
fi

# snpEff DB
if have snpEff; then
  SNPEFF_DB_NAME="GRCh38.${ENSEMBL_RELEASE}"
  log "Attempting snpEff DB download: ${SNPEFF_DB_NAME}"
  set +e
  snpEff download -v "${SNPEFF_DB_NAME}"
  set -e
fi

mkdir -p "${ANN_DIR}/vep_cache"

cat > bin/check_variant_containers.sh <<'CHK'
#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
for i in manta.sif lumpy.sif deepvariant.sif; do
  [[ -s "$BASE_DIR/containers/$i" ]] && echo "[OK] $BASE_DIR/containers/$i" || echo "[MISS] $BASE_DIR/containers/$i"
done
ls -1 "$BASE_DIR/containers"/vep.r*.sif 2>/dev/null || echo "[MISS] VEP container (vep.r*.sif)"
CHK
chmod +x bin/check_variant_containers.sh

log "Final check"
for t in bwa bwa-mem2 samtools bcftools picard fastqc multiqc gatk delly cnvkit.py snpEff; do
  have "$t" && echo "  [OK] $t -> $(command -v "$t")" || echo "  [MISS] $t"
done
[[ -s "${DV_SIF}" ]] && echo "  [OK] deepvariant.sif -> ${DV_SIF}" || echo "  [MISS] deepvariant.sif"
[[ -s "${VEP_SIF}" ]] && echo "  [OK] vep.sif -> ${VEP_SIF}" || echo "  [MISS] vep.sif"

log "Heavy setup complete."
BASH
chmod +x "${HEAVY_SCRIPT}"

# -----------------------------
# Dispatch strategy
# -----------------------------
if [[ "${MODE}" == "local" ]]; then
  log "Running heavy setup locally in current shell."
  bash "${HEAVY_SCRIPT}" "${BASE_DIR}" "${ENV_NAME}" "${THREADS}" "${VEP_SPECIES}" "${VEP_ASSEMBLY}" "${ENSEMBL_RELEASE}"
  exit 0
fi

# MODE=auto: if already in LSF allocation, run directly; else submit to LSF
if [[ -n "${LSB_JOBID:-}" ]]; then
  log "Detected LSF job context (LSB_JOBID=${LSB_JOBID}); running heavy setup now."
  bash "${HEAVY_SCRIPT}" "${BASE_DIR}" "${ENV_NAME}" "${THREADS}" "${VEP_SPECIES}" "${VEP_ASSEMBLY}" "${ENSEMBL_RELEASE}"
else
  log "Submitting heavy setup to LSF compute node..."
  mkdir -p "${BASE_DIR}/logs"
  bsub -J "resolveome_setup" \
       -q "${QUEUE}" \
       -n "${NCORES}" \
       -M "${MEM_MB}" \
       -R "rusage[mem=${MEM_MB}]" \
       -W "${WALL}" \
       -oo "${BASE_DIR}/logs/setup.%J.out" \
       -eo "${BASE_DIR}/logs/setup.%J.err" \
       "bash ${HEAVY_SCRIPT} ${BASE_DIR} ${ENV_NAME} ${THREADS} ${VEP_SPECIES} ${VEP_ASSEMBLY} ${ENSEMBL_RELEASE}"

  log "Submitted. Monitor with:"
  echo "  bjobs -w -J resolveome_setup"
  echo "  tail -f ${BASE_DIR}/logs/setup.*.out"
fi
