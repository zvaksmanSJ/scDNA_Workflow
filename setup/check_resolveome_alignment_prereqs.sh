#!/bin/bash
set -euo pipefail

# Usage:
#   bash check_resolveome_alignment_prereqs.sh \
#     /path/to/GRCh38.primary.fa \
#     /path/to/dna_cells.tsv \
#     resolveome-align

REF="${1:-}"
CELLS="${2:-}"
ENV_NAME="${3:-resolveome-align}"

if [[ -z "$REF" || -z "$CELLS" ]]; then
  echo "Usage: $0 <ref.fa> <dna_cells.tsv> [conda_env]"
  exit 2
fi

FAILS=()

check_file() {
  local f="$1"
  local label="$2"
  if [[ -s "$f" ]]; then
    echo "[PASS] $label : $f"
  else
    echo "[FAIL] $label : $f"
    FAILS+=("$label :: $f")
  fi
}

check_cmd_in_env() {
  local cmd="$1"
  local label="$2"
  if conda run -n "$ENV_NAME" bash -lc "command -v $cmd >/dev/null 2>&1"; then
    echo "[PASS] tool $label ($cmd) in env $ENV_NAME"
  else
    echo "[FAIL] tool $label ($cmd) missing in env $ENV_NAME"
    FAILS+=("tool $label ($cmd) missing in env $ENV_NAME")
  fi
}

echo "==== ResolveOME Illumina alignment preflight ===="
echo "REF      : $REF"
echo "CELLS TSV: $CELLS"
echo "ENV      : $ENV_NAME"
echo

# Core input files
check_file "$REF" "Reference FASTA"
check_file "${REF}.fai" "samtools FASTA index (.fai)"
check_file "$CELLS" "Cell manifest TSV"

# Reference dict (optional for pure alignment, but expected in your setup)
if [[ -s "${REF%.fa}.dict" || -s "${REF%.fasta}.dict" ]]; then
  echo "[PASS] Sequence dictionary (.dict) present"
else
  echo "[WARN] Sequence dictionary (.dict) missing (not required for bwa-mem2 align, but recommended)"
fi

# bwa-mem2 index files (at least one complete style should exist)
if [[ -s "${REF}.bwt.2bit.64" || -s "${REF}.0123" ]]; then
  echo "[PASS] bwa-mem2 index appears present"
else
  echo "[FAIL] bwa-mem2 index missing (expected e.g. ${REF}.bwt.2bit.64 or ${REF}.0123)"
  FAILS+=("bwa-mem2 index missing for $REF")
fi

# Conda/tool checks
if command -v conda >/dev/null 2>&1; then
  echo "[PASS] conda available"
else
  echo "[FAIL] conda not found in PATH"
  FAILS+=("conda not found in PATH")
fi

check_cmd_in_env bwa-mem2 "Aligner"
check_cmd_in_env samtools "BAM tools"
check_cmd_in_env fastqc "QC"
check_cmd_in_env multiqc "QC aggregation"

# Validate cells TSV columns + paired files
echo
echo "Checking cells TSV rows and FASTQ pairs..."
python - "$CELLS" <<'PY'
import csv, sys, pathlib
cells = pathlib.Path(sys.argv[1])
required = {"sample_id","r1","r2"}
fails = []
with cells.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    hdr = set(r.fieldnames or [])
    miss = required - hdr
    if miss:
        print(f"[FAIL] cells TSV missing columns: {sorted(miss)}")
        sys.exit(10)
    n = 0
    for i,row in enumerate(r, start=2):
        sid=row["sample_id"].strip(); r1=row["r1"].strip(); r2=row["r2"].strip()
        if not sid:
            fails.append(f"line {i}: empty sample_id")
            continue
        if not pathlib.Path(r1).is_file():
            fails.append(f"{sid}: missing R1 -> {r1}")
        if not pathlib.Path(r2).is_file():
            fails.append(f"{sid}: missing R2 -> {r2}")
        n += 1
print(f"[PASS] cells TSV header OK")
print(f"[INFO] samples parsed: {n}")
if fails:
    print("[FAIL] FASTQ issues:")
    for x in fails:
        print(f"  - {x}")
    sys.exit(11)
else:
    print("[PASS] all FASTQ paths exist")
PY
py_status=$?
if [[ $py_status -ne 0 ]]; then
  FAILS+=("cells TSV/FASTQ validation failed")
fi

echo
echo "==== Summary ===="
if [[ ${#FAILS[@]} -eq 0 ]]; then
  echo "ALL REQUIRED COMPONENTS PRESENT"
  exit 0
else
  echo "FAILED COMPONENTS:"
  for f in "${FAILS[@]}"; do
    echo " - $f"
  done
  exit 1
fi
