#!/bin/bash
set -euo pipefail

REF="/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/workflow_DNA/resolveome_align/refs/genome/GRCh38.primary.fa"
LOGDIR="$(dirname "$REF")/index_logs"
mkdir -p "$LOGDIR"

bsub <<EOF
#!/bin/bash
#BSUB -J build_ref_idx
#BSUB -q standard
#BSUB -n 8
#BSUB -M 128000
#BSUB -R "rusage[mem=128000]"
#BSUB -W 08:00
#BSUB -oo ${LOGDIR}/build_ref_idx.%J.out
#BSUB -eo ${LOGDIR}/build_ref_idx.%J.err

set -euo pipefail
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate resolveome-align

[[ -s "${REF}" ]] || { echo "Missing REF"; exit 1; }

# build only if missing
[[ -s "${REF}.fai" ]] || samtools faidx "${REF}"
[[ -s "${REF%.fa}.dict" ]] || picard CreateSequenceDictionary R="${REF}" O="${REF%.fa}.dict"

if [[ ! -s "${REF}.bwt.2bit.64" && ! -s "${REF}.0123" ]]; then
  bwa-mem2 index "${REF}"
fi

echo "Index build complete"
ls -lh "${REF}".*
EOF
