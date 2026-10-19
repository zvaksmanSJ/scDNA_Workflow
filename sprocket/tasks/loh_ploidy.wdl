## Per-cell CNA + loss-of-heterozygosity (and, optionally, best-fit ploidy
## selection) via Control-FREEC (Boeva et al. 2011/2012, GPL, bioconda).
## Sprocket port of ../../wdl/tasks/loh_ploidy.wdl -- see that file's header
## comment for the `chr_len_file`-must-be-filtered bug writeup and the
## candidate-ploidy model-selection caveat (both unchanged here).
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: no `conda activate` -- `freec`
## itself must already be on PATH inside the container named by
## `freec_container` (see ../containers/freec.Dockerfile).
version 1.2

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
    String freec_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-freec:latest"
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

    freec -conf freec_config.txt

    # FREEC names its own output files after the input BAM's basename;
    # rename them to this cell's own sample_id-based naming convention so
    # the `output {}` paths below can find them.
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
  requirements {
    container: freec_container
    cpu: threads
    memory: "32 GB"
  }
}
