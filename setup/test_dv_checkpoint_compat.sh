#!/usr/bin/env bash
set -euo pipefail

# ---- edit these ----
REF="/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/workflow_DNA/resolveome_align/refs/genome/GRCh38.primary.fa"
BAM="/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/results/DNA_Ilumina_output_1/bam/SY5Y_Cell10_DNA_SJNBL046426_C855_DNA_C10/SY5Y_Cell10_DNA_SJNBL046426_C855_DNA_C10.dedup.bam"
CKPT_DIR="/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/workflow_DNA/resolveome_align/refs/genome/deepvariant_pta"
CKPT="${CKPT_DIR}/sc_af.model.ckpt"
REGION="1:1-2000000"
WORKROOT="$PWD/dv_compat_test"
THREADS=4

# Candidate SIFs to test (add/remove as needed)
CANDIDATES=(
  "/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/workflow_DNA/resolveome_align/containers/deepvariant.sif"
  # "/path/to/deepvariant_1.8.0.sif"
  # "/path/to/deepvariant_1.7.0.sif"
)

APP=$(command -v apptainer || command -v singularity || true)
[[ -n "$APP" ]] || { echo "ERROR: apptainer/singularity not found"; exit 2; }

mkdir -p "$WORKROOT"

for SIF in "${CANDIDATES[@]}"; do
  echo "=== TESTING: $SIF ==="
  [[ -s "$SIF" ]] || { echo "  -> missing SIF, skip"; continue; }

  TDIR="$WORKROOT/$(basename "$SIF" .sif)"
  rm -rf "$TDIR"
  mkdir -p "$TDIR"
  export TMPBASE="$TDIR/tmp"
  mkdir -p "$TMPBASE"

  set +e
  "$APP" exec -B /research_jude:/research_jude --env TMPBASE="$TMPBASE" "$SIF" bash -lc "
    set -eo pipefail
    /opt/deepvariant/bin/make_examples \
      --mode calling \
      --ref '$REF' \
      --reads '$BAM' \
      --regions '$REGION' \
      --examples '$TDIR/examples.tfrecord@${THREADS}.gz' \
      --checkpoint '$CKPT'

    /opt/deepvariant/bin/call_variants \
      --examples '$TDIR/examples.tfrecord@${THREADS}.gz' \
      --outfile '$TDIR/call_variants_output.tfrecord.gz' \
      --checkpoint '$CKPT'
  " >"$TDIR/test.out" 2>"$TDIR/test.err"
  RC=$?
  set -e

  if [[ $RC -eq 0 ]]; then
    echo "  -> PASS (compatible)"
  else
    echo "  -> FAIL (incompatible or runtime issue)"
    tail -n 40 "$TDIR/test.err" || true
  fi
done
