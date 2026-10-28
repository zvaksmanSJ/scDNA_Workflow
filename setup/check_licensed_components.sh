#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# check_licensed_components.sh
#
# Validates whether each BioSkryb/BaseJumper-LICENSED component the
# production Python/LSF pipeline (python_lsf_pipeline/run_all_resolveome.py)
# needs is present at its configured location, and gives specific,
# actionable guidance when one is missing. See LICENSED_COMPONENTS.md
# (workflow_DNA/ root) for the full background on what each component is,
# why it's licensed rather than free/redistributable, and how to obtain
# your own legitimate copy from BioSkryb.
#
# This script does NOT:
#   - download, build, bundle, or redistribute any licensed component
#   - replicate, approximate, or otherwise work around any license
#   - modify anything -- it only checks and reports
#
# It honors the same 4 environment-variable overrides run_all_resolveome.py
# itself understands, so running this BEFORE a real pipeline run (with the
# same env vars set, if you're using them) tells you exactly what that run
# would fail on, before it fails:
#
#   SCDNA_BASEJ_PIPELINES_DIR         (default: <workflow_DNA>/resolveome_align/basej-public-pipelines)
#   SCDNA_BASEJ_CONTAINERS_DIR        (default: <workflow_DNA>/resolveome_align/containers/basej-somatic)
#   SCDNA_DEEPVARIANT_CHECKPOINT_DIR  (default: <project_root>/pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/DeepVariant/Illumina/bioskryb-af-20241102)
#   SCDNA_BASEJ_GENOMES_DIR           (default: <workflow_DNA>/resolveome_align/basej_genomes_nochr)
#
# Usage:
#   bash setup/check_licensed_components.sh
#   # or, pointing at components obtained under your own license, placed elsewhere:
#   SCDNA_BASEJ_PIPELINES_DIR=/path/to/your/basej-public-pipelines \
#   SCDNA_BASEJ_CONTAINERS_DIR=/path/to/your/basej-somatic-containers \
#   SCDNA_DEEPVARIANT_CHECKPOINT_DIR=/path/to/your/bioskryb-af-checkpoint \
#   SCDNA_BASEJ_GENOMES_DIR=/path/to/your/basej-genomes-bundle \
#     bash setup/check_licensed_components.sh
# ============================================================================

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_DNA=$(cd -- "$HERE/.." && pwd)
PROJECT_ROOT=$(cd -- "$WORKFLOW_DNA/.." && pwd)

BASEJ_PIPELINES_DIR="${SCDNA_BASEJ_PIPELINES_DIR:-$WORKFLOW_DNA/resolveome_align/basej-public-pipelines}"
BASEJ_CONTAINERS_DIR="${SCDNA_BASEJ_CONTAINERS_DIR:-$WORKFLOW_DNA/resolveome_align/containers/basej-somatic}"
DEEPVARIANT_CHECKPOINT_DIR="${SCDNA_DEEPVARIANT_CHECKPOINT_DIR:-$PROJECT_ROOT/pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/DeepVariant/Illumina/bioskryb-af-20241102}"
BASEJ_GENOMES_DIR="${SCDNA_BASEJ_GENOMES_DIR:-$WORKFLOW_DNA/resolveome_align/basej_genomes_nochr}"

FAILS=0

section() { echo; echo "==== $1 ===="; }

check_path() {
  local path="$1" label="$2"
  if [[ -e "$path" ]]; then
    echo "  [OK]   $label"
    echo "         -> $path"
  else
    echo "  [MISS] $label"
    echo "         -> $path (not found)"
    FAILS=$((FAILS + 1))
  fi
}

echo "ResolveOME DNA workflow -- licensed-component preflight"
echo "(This script only checks/reports. It never downloads, builds, or"
echo " bundles anything. See LICENSED_COMPONENTS.md for how to obtain each"
echo " component legitimately from BioSkryb/BaseJumper.)"

# ----------------------------------------------------------------------------
section "1. basej-public-pipelines (basej-somatic, basej-lineage, basej-dnaqc Nextflow pipelines)"
echo "Directory: \$SCDNA_BASEJ_PIPELINES_DIR (or its default location)"
check_path "$BASEJ_PIPELINES_DIR/basej-somatic/main.nf" "basej-somatic Nextflow entry point"
check_path "$BASEJ_PIPELINES_DIR/basej-lineage/main.nf" "basej-lineage Nextflow entry point"
check_path "$BASEJ_PIPELINES_DIR/basej-dnaqc/assets/plot_qc_config.json" "basej-dnaqc QC-threshold config"
if [[ ! -e "$BASEJ_PIPELINES_DIR/basej-somatic/main.nf" ]]; then
  cat <<'EOF'
       This is a BioSkryb/BaseJumper-licensed Nextflow pipeline bundle, NOT
       free/public software -- it is not something this project can supply,
       and no amount of local re-implementation substitutes for an actual
       license. To obtain it:
         1. Contact BioSkryb Genomics (https://www.bioskryb.com/) or your
            existing BioSkryb/BaseJumper account representative to establish
            or confirm your organization's license for the BaseJumper
            platform's underlying Nextflow pipeline source.
         2. Once licensed, BioSkryb will provide access to the
            "basej-public-pipelines" bundle (or equivalent delivery
            mechanism they use for licensees) containing basej-somatic/,
            basej-lineage/, and basej-dnaqc/ as Nextflow pipeline source.
         3. Place (or symlink) that bundle so its basej-somatic/,
            basej-lineage/, and basej-dnaqc/ subdirectories sit directly
            under the directory this check just reported above -- or set
            SCDNA_BASEJ_PIPELINES_DIR to point at wherever you placed it.
       Without a license, use main.wdl's free/open-source parity path
       instead (see workflow_DNA/README.md) -- it does not need this.
EOF
fi

# ----------------------------------------------------------------------------
section "2. basej-somatic Apptainer/Singularity containers (.sif)"
echo "Directory: \$SCDNA_BASEJ_CONTAINERS_DIR (or its default location)"
for image in utility-lite-0.2.0 variantannotation_0.4 r_ado_1.2 custom_snp_somatic_filter_sequoia_feb2026; do
  check_path "$BASEJ_CONTAINERS_DIR/${image}.sif" "$image.sif"
done
if [[ ! -e "$BASEJ_CONTAINERS_DIR/custom_snp_somatic_filter_sequoia_feb2026.sif" ]]; then
  cat <<'EOF'
       These are the pre-built Apptainer/Singularity container images (.sif)
       BioSkryb provides to BaseJumper licensees to actually execute the
       Nextflow pipelines above -- these are container IMAGES, not something
       you build yourself from public Dockerfiles (they bundle BioSkryb's
       own proprietary tool builds/configurations). Obtain these through the
       same BioSkryb/BaseJumper licensing channel as item 1 above; place the
       four .sif files directly under the directory reported above, or set
       SCDNA_BASEJ_CONTAINERS_DIR to point at wherever you placed them.
EOF
fi

# ----------------------------------------------------------------------------
section "3. DeepVariant single-cell checkpoint (BioSkryb-private model)"
echo "Directory: \$SCDNA_DEEPVARIANT_CHECKPOINT_DIR (or its default location)"
check_path "$DEEPVARIANT_CHECKPOINT_DIR/model/sc_af.model.ckpt.index" "sc_af.model.ckpt (index)"
check_path "$DEEPVARIANT_CHECKPOINT_DIR/model/sc_af.model.ckpt.data-00000-of-00001" "sc_af.model.ckpt (weights)"
check_path "$DEEPVARIANT_CHECKPOINT_DIR/population_vcfs_nochr" "population_vcfs_nochr/ (matching population-AF VCF bundle)"
if [[ ! -e "$DEEPVARIANT_CHECKPOINT_DIR/model/sc_af.model.ckpt.index" ]]; then
  cat <<'EOF'
       This is BioSkryb's own privately-trained single-cell DeepVariant
       checkpoint (and its matching population-allele-frequency VCF
       bundle), tuned specifically for BioSkryb's PTA/WGA chemistry -- it is
       not Google's public DeepVariant model and is not redistributable.
       Obtain it through your BioSkryb/BaseJumper license (same channel as
       above); place the checkpoint's "model/" and "population_vcfs_nochr/"
       directories directly under the directory reported above, or set
       SCDNA_DEEPVARIANT_CHECKPOINT_DIR to point at wherever you placed them.
       Without this checkpoint, main.wdl's free path (see workflow_DNA/
       README.md) falls back to DeepVariant's own free bundled WGS model
       automatically -- no checkpoint needed there at all.
EOF
fi

# ----------------------------------------------------------------------------
section "4. basej-genomes reference bundle (nf-core-style genome layout for the Nextflow stages)"
echo "Directory: \$SCDNA_BASEJ_GENOMES_DIR (or its default location)"
check_path "$BASEJ_GENOMES_DIR/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/GATK_bundle/Sequence/genome.fa.fai" \
  "genome.fa.fai (within the basej-genomes layout)"
if [[ ! -e "$BASEJ_GENOMES_DIR/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/GATK_bundle/Sequence/genome.fa.fai" ]]; then
  cat <<'EOF'
       This directory holds the same public GRCh38 reference genome used
       everywhere else in this project, but re-organized into the specific
       nf-core-style directory layout basej-somatic/basej-lineage's own
       Nextflow code expects as --genomes_base. BioSkryb typically provides
       this pre-arranged bundle alongside the pipeline bundle in item 1
       above (so you do not need to manually reorganize a reference genome
       yourself) -- obtain it through the same license channel. Place it
       so its genomes/ subdirectory sits directly under the directory
       reported above, or set SCDNA_BASEJ_GENOMES_DIR accordingly.
EOF
fi

# ----------------------------------------------------------------------------
echo
echo "==== Summary ===="
if [[ $FAILS -eq 0 ]]; then
  echo "ALL LICENSED COMPONENTS PRESENT -- the production Python/LSF pipeline"
  echo "(run_all_resolveome.py) should be able to find everything it needs."
  exit 0
else
  echo "$FAILS licensed component(s) missing -- see the guidance printed"
  echo "above each [MISS] line for how to legitimately obtain and place it."
  echo
  echo "If you do not have (or cannot obtain) a BioSkryb/BaseJumper license,"
  echo "use the free/open-source parity pipeline instead: main.wdl, run via"
  echo "Cromwell -- see workflow_DNA/README.md and workflow_DNA/wdl/README.md."
  echo "It needs NONE of the components checked above."
  exit 1
fi
