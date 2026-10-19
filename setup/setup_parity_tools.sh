#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate resolveome-align
for tool in fastp bedtools preseq; do
  command -v "$tool" >/dev/null || {
    echo "[ERROR] Missing $tool. Install the additions in parity_tools.environment.yml into resolveome-align." >&2
    exit 2
  }
done

# Keep the reference's bam2mr/gc_extrap pair separate from the installed preseq 3.
TARGET="$HERE/../resolveome_align/tools/preseq-2.0.3"  # resolveome_align/ lives at workflow_DNA root, one level above this script's own setup/ subfolder
if [[ ! -x "$TARGET/bin/bam2mr" || ! -x "$TARGET/bin/preseq" ]]; then
  mkdir -p "$TARGET"
  BUILD=$(mktemp -d "$TARGET/build.XXXXXX")
  curl --fail --location --retry 3 \
    https://github.com/smithlabcode/preseq/releases/download/v2.0.3/preseq_v2.0.3.tar.bz2 \
    --output "$BUILD/preseq_v2.0.3.tar.bz2"
  tar -xjf "$BUILD/preseq_v2.0.3.tar.bz2" -C "$BUILD"
  (
    cd "$BUILD/preseq"
    make all CXXFLAGS="-O2 -Wall -fno-common" CFLAGS="-O2 -fcommon"
  )
  mkdir -p "$TARGET/bin"
  install -m 755 "$BUILD/preseq/preseq" "$TARGET/bin/preseq"
  install -m 755 "$BUILD/preseq/bam2mr" "$TARGET/bin/bam2mr"
  sha256sum "$BUILD/preseq_v2.0.3.tar.bz2" > "$TARGET/source_archive.sha256"
  rm -rf "$BUILD"
fi
ln -sfn "$TARGET/bin/bam2mr" "$CONDA_PREFIX/bin/bam2mr"
ln -sfn "$TARGET/bin/preseq" "$CONDA_PREFIX/bin/preseq-bioskryb"
test -x "$CONDA_PREFIX/bin/bam2mr"
test -x "$CONDA_PREFIX/bin/preseq-bioskryb"
echo "[OK] Reference preseq 2.0.3 and bam2mr installed in resolveome-align; validate with the workflow preflight/smoke tests"
