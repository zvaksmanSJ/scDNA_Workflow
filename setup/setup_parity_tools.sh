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
#
# This builds preseq v2.0.3 (GPLv3, fully public -- https://github.com/
# smithlabcode/preseq, also on bioconda) directly from its official GitHub
# release tarball below. It is NOT a private/proprietary BioSkryb binary and
# has no license restriction of its own; it is pinned to exactly v2.0.3
# because that release's `bam2mr`/`gc_extrap` CLI is what this project's
# pipeline (and BioSkryb's own ResolveDNA protocol) was built against --
# preseq 3.x renamed/restructured this command (`bam2mr` -> `to-mr`) in an
# incompatible way, so the newer bioconda `preseq` package (also installed
# in this same environment for anything else that wants it) cannot simply
# be substituted here.
#
# TARGET_DIR defaults to this project's own existing resolveome_align/
# installation (preserving this script's original, already-in-production-
# use behavior for direct/standalone invocation); setup.sh instead passes
# an explicit BASE_DIR-relative path so a fresh clone's run installs into
# the portable wdl_runtime/ tree instead of creating a new resolveome_align/
# directory out of thin air.
TARGET="${TARGET_DIR:-$HERE/../resolveome_align/tools/preseq-2.0.3}"
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
