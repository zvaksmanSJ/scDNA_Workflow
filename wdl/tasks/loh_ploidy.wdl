## Per-cell CNA + loss-of-heterozygosity (and, optionally, best-fit ploidy
## selection) via Control-FREEC (Boeva et al. 2011/2012, GPL, bioconda).
## `chr_len_file` MUST be restricted to exactly the chromosomes present in
## `chrom_fasta_dir` (bug found and fixed this session: the full genome
## .fai also lists Y, MT, and ~170 decoy/alt/unplaced scaffolds that are not
## split into per-chromosome FASTAs, and FREEC requires a FASTA for every
## contig listed in chrLenFile -- it is NOT the full-genome .fai).
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: its `input {}`
## fields are filled in by main.wdl's TWO calls to this same task --
## `call loh_t.ControlFreec as freec_fixed { ... }` (ploidy = "2") and
## `call loh_t.ControlFreec as freec_ploidy { ... }` (ploidy =
## freec_candidate_ploidies, e.g. "2,3,4,5,6,7,8") -- using values from
## main.wdl's single inputs.json; there is no separate JSON per cell/call.
##
## KNOWN CAVEAT (found during this project's own real-cohort run, kept here
## rather than hidden): widening the candidate ploidy range can push MORE
## cells to the new ceiling rather than converging on a stable answer --
## finer copy-number discretization at higher ploidy mechanically improves
## FREEC's fit score without necessarily reflecting real biology. Both the
## fixed-ploidy=2 and candidate-selected results are kept and surfaced
## separately in the final dashboard specifically so this is visible, not
## papered over.
version 1.0

# Control-FREEC's own config file (freec_config.txt) is built by this
# command block from the WDL inputs below, then `freec` itself is invoked
# against it; it writes BAF/CNV/ratio/info output files named after the
# input BAM's basename, which are then renamed to this cell's own naming
# convention for WDL's `output {}` globbing to pick up.
task ControlFreec {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File chrom_fasta_dir_tar     # tar of per-chromosome FASTAs (chr1-22+X), each with .fai
    File chr_len_file            # samtools .fai restricted to chr1-22+X ONLY (see module docstring)
    File snp_positions_vcf       # site-only 1000 Genomes VCF, used for both BAF pileup and SNP positions
    File snp_positions_vcf_tbi
    String ploidy = "2"          # e.g. "2" (fixed) or "2,3,4,5,6,7,8" (candidate selection, FREEC's own feature)
    String sex = "XY"
    Int window = 50000
    Int min_coverage_per_position = 5
    Int threads = 4
    String freec_conda_env = "freec"
  }
  command <<<
    set -euo pipefail
    # untar the per-chromosome FASTAs that chrFiles below will point at
    CHROM_FASTA_DIR=$(mktemp -d)
    tar -xf ~{chrom_fasta_dir_tar} -C "$CHROM_FASTA_DIR"

    # Control-FREEC's own INI-style config file, built from the WDL inputs:
    #   [general] -- chromosome lengths/FASTAs, ploidy candidate(s), window size
    #   [sample]  -- the BAM to call CNAs/LOH from
    #   [BAF]     -- the site-only VCF used both to build a pileup and as the
    #                known-SNP-position list for B-allele-frequency (LOH) calling
    cat > freec_config.txt <<EOF
[general]
chrLenFile = ~{chr_len_file}
chrFiles = $CHROM_FASTA_DIR
ploidy = ~{ploidy}
sex = ~{sex}
window = ~{window}
breakPointType = 2
outputDir = $(pwd)
maxThreads = ~{threads}

[sample]
mateFile = ~{bam}
inputFormat = BAM
mateOrientation = FR

[BAF]
makePileup = ~{snp_positions_vcf}
fastaFile = ~{ref_fasta}
SNPfile = ~{snp_positions_vcf}
minimalCoveragePerPosition = ~{min_coverage_per_position}
shiftInQuality = 33
EOF

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{freec_conda_env}
    freec -conf freec_config.txt

    # FREEC names its own output files after the input BAM's basename;
    # rename them to this cell's own sample_id-based naming convention so
    # WDL's `output {}` paths below can find them.
    BAM_BASENAME=$(basename ~{bam})
    test -s "${BAM_BASENAME}_CNVs"
    mv "${BAM_BASENAME}_CNVs" "~{sample_id}.dedup.bam_CNVs"
    mv "${BAM_BASENAME}_ratio.txt" "~{sample_id}.dedup.bam_ratio.txt"
    mv "${BAM_BASENAME}_BAF.txt" "~{sample_id}.dedup.bam_BAF.txt" || true
    mv "${BAM_BASENAME}_info.txt" "~{sample_id}.dedup.bam_info.txt"
    rm -rf "$CHROM_FASTA_DIR"
  >>>
  output {
    File cnvs = "~{sample_id}.dedup.bam_CNVs"
    File ratio = "~{sample_id}.dedup.bam_ratio.txt"
    File info = "~{sample_id}.dedup.bam_info.txt"
  }
  runtime { cpu: threads, memory: "32 GB" }
}
