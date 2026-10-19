#!/bin/bash
set -euo pipefail

REF="/research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/workflow_DNA/resolveome_align/refs/genome/GRCh38.primary.fa"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate resolveome-align

# required side files
samtools faidx "$REF"
picard CreateSequenceDictionary R="$REF" O="${REF%.fa}.dict" || true

# critical for your failure
bwa-mem2 index "$REF"

# quick check
ls -lh "${REF}".*
