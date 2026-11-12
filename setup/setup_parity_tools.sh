#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wdl-resolveome-align
# NOTE: deliberately NOT checking for a plain "preseq" binary here -- that
# package used to be in parity_tools.environment.yml but was removed (every
# build requires zlib<1.3, conflicting with samtools=1.20's newer zlib, and
# nothing in this pipeline ever invokes a bare `preseq` call anyway -- see
# that file's own comment). This script builds its own preseq 2.0.3 from
# source below regardless, so requiring the (now absent) conda one first
# would only block that unnecessarily.
for tool in fastp bedtools; do
  command -v "$tool" >/dev/null || {
    echo "[ERROR] Missing $tool. Install the additions in parity_tools.environment.yml into wdl-resolveome-align." >&2
    exit 2
  }
done

# Build our own bam2mr/gc_extrap pair rather than relying on bioconda's
# newer preseq package.
#
# This builds preseq v2.0.3 (GPLv3, fully public -- https://github.com/
# smithlabcode/preseq, also on bioconda) directly from its official GitHub
# release tarball below. It is NOT a private/proprietary BioSkryb binary and
# has no license restriction of its own; it is pinned to exactly v2.0.3
# because that release's `bam2mr`/`gc_extrap` CLI is what this project's
# pipeline (and BioSkryb's own ResolveDNA protocol) was built against --
# preseq 3.x renamed/restructured this command (`bam2mr` -> `to-mr`) in an
# incompatible way, so the newer bioconda `preseq` package cannot simply be
# substituted here -- which is also why that package isn't installed into
# this environment at all (see parity_tools.environment.yml's own comment).
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

  # Build fully static (glibc + GSL + zlib all statically linked). This
  # script typically runs interactively on a login/head node, but the
  # resulting binaries are later executed via LSF on whatever compute node
  # the scheduler picks -- on a cluster with a heterogeneous node fleet
  # (mixed RHEL7/8/9 images etc.) a dynamically-linked build only works on
  # nodes with a glibc new enough to match the login node, and fails with
  # "version `GLIBC_2.3x' not found" everywhere older. Static linking makes
  # the binary's own node-of-origin irrelevant. CONDA_PREFIX/lib already
  # has libgsl.a/libgslcblas.a/libz.a (pulled in as transitive deps of other
  # packages), and glibc-static is required on the build host for libc.a.
  CONDA_LIBGSL_A="$CONDA_PREFIX/lib/libgsl.a"
  CONDA_LIBGSLCBLAS_A="$CONDA_PREFIX/lib/libgslcblas.a"
  CONDA_LIBZ_A="$CONDA_PREFIX/lib/libz.a"
  STATIC_OK=1
  for a in "$CONDA_LIBGSL_A" "$CONDA_LIBGSLCBLAS_A" "$CONDA_LIBZ_A"; do
    [[ -s "$a" ]] || STATIC_OK=0
  done
  (
    cd "$BUILD/preseq"
    if [[ "$STATIC_OK" == "1" ]] && \
       make all CXXFLAGS="-O2 -Wall -fno-common -static" CFLAGS="-O2 -fcommon -static" \
         LIBS="$CONDA_LIBGSL_A $CONDA_LIBGSLCBLAS_A $CONDA_LIBZ_A" 2>/tmp/preseq_static_build.$$.log; then
      echo "[OK] Built bam2mr/preseq fully statically (portable across heterogeneous cluster nodes)"
    else
      echo "[WARN] Static build unavailable (missing glibc-static or conda's static GSL/zlib archives)." >&2
      echo "[WARN] Falling back to a dynamically-linked build -- this WILL fail with" >&2
      echo "[WARN] 'GLIBC_2.3x not found' on any LSF compute node whose glibc is older" >&2
      echo "[WARN] than this build host's. See $TARGET/build log above/at /tmp/preseq_static_build.$$.log" >&2
      make clean >/dev/null 2>&1 || true
      make all CXXFLAGS="-O2 -Wall -fno-common" CFLAGS="-O2 -fcommon"
    fi
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
echo "[OK] Reference preseq 2.0.3 and bam2mr installed in wdl-resolveome-align; validate with the workflow preflight/smoke tests"
