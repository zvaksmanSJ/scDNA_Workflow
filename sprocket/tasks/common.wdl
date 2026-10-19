## Small shared helper task(s), used across multiple other task files in
## this Sprocket port -- imported directly by ../main.wdl (WDL tasks cannot
## call other tasks directly, only workflows can, so this file is imported
## alongside the others rather than being imported BY cnv_sv.wdl/
## variant_calling.wdl themselves).
version 1.2

# Bcftools-indexes any VCF that didn't already come with an index from its
# own tool. Used as a fallback for Manta (see cnv_sv.wdl's MantaCall -- only
# called when Manta's own `.tbi` is missing) and unconditionally for every
# DeepVariant/VEP chunk output in variant_calling.wdl (those tools never
# produce their own index). Produces a `.csi` (bcftools index's own default
# format), matching the original Cromwell-version tasks' exact behavior.
task IndexVcf {
  input {
    File vcf
    String out_name    # the desired output basename, e.g. "~{sample_id}.manta.vcf.gz"
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    cp ~{vcf} ~{out_name}
    bcftools index -f ~{out_name}
  >>>
  output {
    File indexed_vcf = out_name
    File vcf_csi = "~{out_name}.csi"
  }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}
