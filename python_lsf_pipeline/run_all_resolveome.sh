#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate resolveome-align
if [[ "${1:-}" == "--submit" ]]; then
  shift
  # Submission is deliberately explicit; --preflight and --plan never start analyses.
  bsub -J resolveome_parity_driver -q standard -n 2 -M 16000 \
    -R 'rusage[mem=16000]' -W 168:00 \
    -oo "$HERE/logs/resolveome_parity_driver.%J.out" \
    -eo "$HERE/logs/resolveome_parity_driver.%J.err" \
    bash "$HERE/run_all_resolveome.sh" --run "$@"
else
  exec "$CONDA_PREFIX/bin/python" "$HERE/run_all_resolveome.py" "$@"
fi
