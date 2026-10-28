## Per-cell, per-autosome statistical phasing via Beagle 5.5 (GPLv3, free)
## against the 1000 Genomes GRCh38 high-coverage phased reference panel.
## chrX is intentionally excluded: this is appropriate for male samples
## (hemizygous non-PAR chrX) where the panel's chrX reference VCF mixes
## haploid/diploid genotype representations that Beagle cannot parse in one
## pass -- confirmed this session to fail identically on all cells tested.
## If your cohort includes female samples needing chrX phasing, split chrX
## into PAR1/non-PAR/PAR2 regions and phase them separately (not implemented
## here, since this cohort is male and did not need it).
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: both tasks' inputs
## are filled in by main.wdl's `call phasing_t.PhaseChromosome as
## phase_chrom { ... }` (inside a nested per-cell x per-autosome scatter) and
## `call phasing_t.ConcatPhasedChromosomes as concat_phased { ... }` blocks,
## using values from main.wdl's single inputs.json; there is no separate
## JSON per cell/chromosome in normal cohort-wide use.
version 1.0

# One PhaseChromosome call per (cell, autosome) pair. Steps: subset the raw
# VCF to just this chromosome, rename it from "1" to "chr1" (matching the
# 1000G panel's own contig-naming convention), run Beagle (statistical
# phasing only -- `impute=false`, since we only want existing calls phased,
# not new sites imputed), then rename the phased result back to "1" so every
# downstream file in this pipeline stays on one consistent naming scheme.
task PhaseChromosome {
  input {
    String sample_id
    File raw_vcf             # un-annotated merged DeepVariant VCF (GT-only input)
    File raw_vcf_csi
    String chrom              # "1".."22" (no chrX; see module docstring)
    File ref_panel_vcf        # 1000G GRCh38 panel VCF for this chromosome
    File ref_panel_vcf_tbi
    File genetic_map          # plink-format genetic map for this chromosome
    File beagle_jar
    Int mem_gb = 48
    Int threads = 4
    String conda_env = "wdl-resolveome-align"
  }
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}

    echo -e "~{chrom}\tchr~{chrom}" > rename_to_chr.tsv
    echo -e "chr~{chrom}\t~{chrom}" > rename_from_chr.tsv

    bcftools view -r ~{chrom} ~{raw_vcf} -Oz -o input_chrom.vcf.gz
    bcftools annotate --rename-chrs rename_to_chr.tsv -Oz -o input_renamed.vcf.gz input_chrom.vcf.gz
    bcftools index -f input_renamed.vcf.gz

    java -Xmx~{mem_gb - 4}g -jar ~{beagle_jar} \
      gt=input_renamed.vcf.gz \
      ref=~{ref_panel_vcf} \
      map=~{genetic_map} \
      chrom=chr~{chrom} \
      out=phased_chrom \
      nthreads=~{threads} \
      impute=false

    test -s phased_chrom.vcf.gz
    bcftools index -f phased_chrom.vcf.gz
    bcftools annotate --rename-chrs rename_from_chr.tsv -Oz -o ~{sample_id}.phased.chr~{chrom}.vcf.gz phased_chrom.vcf.gz
    bcftools index -f ~{sample_id}.phased.chr~{chrom}.vcf.gz
  >>>
  output {
    File phased_chrom_vcf = "~{sample_id}.phased.chr~{chrom}.vcf.gz"
    File phased_chrom_vcf_csi = "~{sample_id}.phased.chr~{chrom}.vcf.gz.csi"
  }
  runtime { cpu: threads, memory: "~{mem_gb} GB" }
}

task ConcatPhasedChromosomes {
  input {
    String sample_id
    Array[File] phased_chrom_vcfs   # chr1..chr22, in any order
    String conda_env = "wdl-resolveome-align"
  }
  # Gathers this cell's 22 per-chromosome phased VCFs (from the scatter
  # above) back into one whole-genome phased VCF.
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}

    # bcftools concat requires its inputs sorted by contig/position; sort
    # the chromosome-numbered inputs numerically (1..22) regardless of the
    # scatter-completion order WDL hands them back in.
    python3 - <<'PYEOF' ~{sep=" " phased_chrom_vcfs} > sorted_list.txt
import re, sys
files = sys.argv[1:]
def key(f):
    m = re.search(r"\.chr(\d+)\.vcf\.gz$", f)
    return int(m.group(1)) if m else 999
for f in sorted(files, key=key):
    print(f)
PYEOF

    bcftools concat -Oz -o ~{sample_id}.phased.vcf.gz -f sorted_list.txt
    bcftools index -f ~{sample_id}.phased.vcf.gz
  >>>
  output {
    File phased_vcf = "~{sample_id}.phased.vcf.gz"
    File phased_vcf_csi = "~{sample_id}.phased.vcf.gz.csi"
  }
  runtime { cpu: 2, memory: "8 GB" }
}
