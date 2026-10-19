## Cohort-level calculations, summaries, and dashboard generation -- all
## genuine Python/R logic (no external bioinformatics tool invocation
## beyond what each script's own docstring already documents). Sprocket
## port of ../../wdl/tasks/extra_qc_dashboard.wdl.
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: no `conda activate` -- every
## task below runs in exactly ONE container (named by
## `resolveome_align_container` for everything except MutationalSignatures,
## which needs `mutsig_container` instead -- see
## ../containers/resolveome-align.Dockerfile and
## ../containers/mutsig.Dockerfile).
##
## Common pattern across most tasks below: since each script expects a
## per-sample-subfolder directory layout (not a flat Array[File] list),
## every task first runs a small inline `python3 - <<'PYEOF'` block to stage
## its Array[File] inputs into that expected layout (named `*_staged/`),
## alongside a simple one-column `manifest.tsv` of sample IDs -- then calls
## the real calculation script against that staged directory.
version 1.2

# Median Absolute Pairwise Difference (MAPD): a per-cell CNV
# noise/uniformity QC metric, computed directly from CNVkit's .cnr bins.
task MapdQc {
  input {
    Array[File] cnr_files
    Array[String] sample_ids
    File mapd_qc_script
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p cnv_staged
    python3 - <<'PYEOF'
import shutil, os
files = """~{sep="\n" cnr_files}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, f in zip(samples, files):
    d = os.path.join("cnv_staged", sample)
    os.makedirs(d, exist_ok=True)
    shutil.copy(f, os.path.join(d, f"{sample}.dedup.cnr"))
PYEOF
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{mapd_qc_script} --cnv-dir cnv_staged --manifest manifest.tsv --output mapd_qc.tsv
  >>>
  output { File mapd_qc_tsv = "mapd_qc.tsv" }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}

# False-discovery-rate estimate for the small-variant calls, cross-checked
# against each cell's own Picard WGS coverage metrics.
task FdrEstimate {
  input {
    Array[File] final_annotated_vcfs
    Array[File] wgs_metrics_files
    Array[String] sample_ids
    File fdr_estimate_script
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p variants_staged picard_staged
    python3 - <<'PYEOF'
import shutil, os
vcfs = """~{sep="\n" final_annotated_vcfs}""".splitlines()
wgs = """~{sep="\n" wgs_metrics_files}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, vcf, w in zip(samples, vcfs, wgs):
    vd = os.path.join("variants_staged", sample); os.makedirs(vd, exist_ok=True)
    shutil.copy(vcf, os.path.join(vd, f"{sample}.deepvariant.raw.vep.snpeff.vcf.gz"))
    pd = os.path.join("picard_staged", sample); os.makedirs(pd, exist_ok=True)
    shutil.copy(w, os.path.join(pd, f"{sample}.wgs_metrics"))
PYEOF
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{fdr_estimate_script} --variants-dir variants_staged --picard-dir picard_staged \
      --manifest manifest.tsv --output fdr_estimate.tsv
  >>>
  output { File fdr_estimate_tsv = "fdr_estimate.tsv" }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}

# Recurrent-artifact panel: flags variant sites that recur suspiciously
# often across the whole cohort (likely technical artifacts rather than
# real biology). Its per-site "singleton" output (sites_tsv) also doubles
# as a proxy private-mutation set for MutationalSignatures below.
task ArtifactPanel {
  input {
    Array[File] final_annotated_vcfs
    Array[String] sample_ids
    File artifact_panel_script
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p variants_staged
    python3 - <<'PYEOF'
import shutil, os
vcfs = """~{sep="\n" final_annotated_vcfs}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, vcf in zip(samples, vcfs):
    d = os.path.join("variants_staged", sample); os.makedirs(d, exist_ok=True)
    shutil.copy(vcf, os.path.join(d, f"{sample}.deepvariant.raw.vep.snpeff.vcf.gz"))
PYEOF
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{artifact_panel_script} --variants-dir variants_staged --manifest manifest.tsv \
      --output-summary artifact_panel_summary.tsv --output-sites artifact_panel_sites.tsv
  >>>
  output {
    File summary_tsv = "artifact_panel_summary.tsv"
    File sites_tsv = "artifact_panel_sites.tsv"
  }
  requirements {
    container: resolveome_align_container
    cpu: 2
    memory: "8 GB"
  }
}

# Cohort-wide summary table over Stage 4 (phasing.wdl)'s per-cell phased VCFs.
task PhasingSummary {
  input {
    Array[File] phased_vcfs
    Array[String] sample_ids
    File phasing_summary_script
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p phasing_staged
    python3 - <<'PYEOF'
import shutil, os
vcfs = """~{sep="\n" phased_vcfs}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, vcf in zip(samples, vcfs):
    d = os.path.join("phasing_staged", sample); os.makedirs(d, exist_ok=True)
    shutil.copy(vcf, os.path.join(d, f"{sample}.phased.vcf.gz"))
PYEOF
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{phasing_summary_script} --phasing-dir phasing_staged --manifest manifest.tsv \
      --output phasing_summary.tsv
  >>>
  output { File phasing_summary_tsv = "phasing_summary.tsv" }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}

# Cohort-wide LOH + ploidy summary tables, PLUS the cross-cell LOH-by-
# genomic-region recurrence heatmap -- all computed from Stage 5
# (loh_ploidy.wdl)'s candidate-ploidy Control-FREEC outputs (freec_ploidy,
# not freec_fixed, in main.wdl's call).
task LohAndPloidySummary {
  input {
    Array[File] cnvs_files
    Array[File] info_files
    Array[String] sample_ids
    File loh_summary_script
    File ploidy_summary_script
    File loh_by_region_script
    File chrom_len_file
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p freec_staged
    python3 - <<'PYEOF'
import shutil, os
cnvs = """~{sep="\n" cnvs_files}""".splitlines()
infos = """~{sep="\n" info_files}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, c, i in zip(samples, cnvs, infos):
    d = os.path.join("freec_staged", sample); os.makedirs(d, exist_ok=True)
    shutil.copy(c, os.path.join(d, f"{sample}.dedup.bam_CNVs"))
    shutil.copy(i, os.path.join(d, f"{sample}.dedup.bam_info.txt"))
PYEOF
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{loh_summary_script} --freec-dir freec_staged --manifest manifest.tsv --output loh_summary.tsv
    python3 ~{ploidy_summary_script} --freec-dir freec_staged --manifest manifest.tsv --output ploidy_summary.tsv
    python3 ~{loh_by_region_script} --freec-dir freec_staged --manifest manifest.tsv \
      --chrom-len-file ~{chrom_len_file} --output-tsv loh_by_region.tsv --output-png loh_by_region_heatmap.png
  >>>
  output {
    File loh_summary_tsv = "loh_summary.tsv"
    File ploidy_summary_tsv = "ploidy_summary.tsv"
    File loh_by_region_tsv = "loh_by_region.tsv"
    File loh_by_region_png = "loh_by_region_heatmap.png"
  }
  requirements {
    container: resolveome_align_container
    cpu: 2
    memory: "8 GB"
  }
}

# Clean per-cell QC overview table + figure, parsed directly out of
# MultiQC's own general_stats data file -- NOT re-staged into a per-sample
# directory like the other tasks above, since this script reads one single
# already-aggregated file (multiqc.data[0] in main.wdl -- see
# qc_aggregate.wdl's MultiQC task for why that index is always
# general_stats).
task QcSummary {
  input {
    File multiqc_general_stats
    Array[String] sample_ids
    File qc_summary_script
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    python3 ~{qc_summary_script} --multiqc-general-stats ~{multiqc_general_stats} \
      --manifest manifest.tsv --output-tsv qc_summary.tsv --output-png qc_overview.png
  >>>
  output {
    File qc_summary_tsv = "qc_summary.tsv"
    File qc_overview_png = "qc_overview.png"
  }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}

# COSMIC mutational-signature refitting via MutationalPatterns (R/
# Bioconductor, free), using the artifact panel's private/singleton
# variant sites (from ArtifactPanel above) as a proxy mutation set -- see
# mutational_signatures.R's own docstring for the explicit methodological
# caveats this proxy approach carries. Runs in `mutsig_container`, NOT
# `resolveome_align_container`, since it needs MutationalPatterns +
# BSgenome.Hsapiens.
task MutationalSignatures {
  input {
    File artifact_panel_sites_tsv
    File mutational_signatures_r_script
    String mutsig_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-mutsig:latest"
  }
  command <<<
    set -euo pipefail
    Rscript ~{mutational_signatures_r_script} ~{artifact_panel_sites_tsv} mutational_signatures_out
  >>>
  output {
    Array[File] outputs = glob("mutational_signatures_out/*")
    File signature_contributions_tsv = "mutational_signatures_out/signature_contributions.tsv"
    File cosine_similarity_tsv = "mutational_signatures_out/reconstruction_cosine_similarity.tsv"
    File mutation_spectrum_png = "mutational_signatures_out/mutation_spectrum_96.png"
    File signature_barplot_png = "mutational_signatures_out/signature_contribution_barplot.png"
  }
  requirements {
    container: mutsig_container
    cpu: 2
    memory: "16 GB"
  }
}

# The final deliverable: builds one tabbed HTML dashboard (Overview/QC/
# Phasing/LOH/LOH-by-Region/Lineage/Mutational-Signatures) stitching
# together every table/figure produced by every cohort-wide task above.
# First re-stages everything into one `extra_methods_staged/` directory
# layout matching what build_dashboard_script expects (mirroring this
# project's own on-disk extra_methods output layout), THEN calls the real
# Python dashboard-building script against it.
task BuildDashboard {
  input {
    File mapd_qc_tsv
    File fdr_estimate_tsv
    File artifact_panel_summary_tsv
    File phasing_summary_tsv
    File loh_summary_tsv
    File ploidy_summary_tsv
    File loh_by_region_tsv
    File loh_by_region_png
    File qc_summary_tsv
    File qc_overview_png
    File signature_contributions_tsv
    File cosine_similarity_tsv
    File mutation_spectrum_png
    File signature_barplot_png
    Array[File] lineage_outputs
    File bootstrap_support_report
    File verification_summary_json
    File build_dashboard_script
    String donor_id = "cohort"
    String multiqc_relpath = "../components/multiqc/multiqc_report.html"
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail
    EM=extra_methods_staged
    mkdir -p "$EM/lineage" "$EM/lineage_verification" "$EM/mutational_signatures"
    cp ~{mapd_qc_tsv} "$EM/mapd_qc.tsv"
    cp ~{fdr_estimate_tsv} "$EM/fdr_estimate.tsv"
    cp ~{artifact_panel_summary_tsv} "$EM/artifact_panel_summary.tsv"
    cp ~{phasing_summary_tsv} "$EM/phasing_summary.tsv"
    cp ~{loh_summary_tsv} "$EM/loh_summary.tsv"
    cp ~{ploidy_summary_tsv} "$EM/ploidy_summary.tsv"
    cp ~{loh_by_region_tsv} "$EM/loh_by_region.tsv"
    cp ~{loh_by_region_png} "$EM/loh_by_region_heatmap.png"
    cp ~{qc_summary_tsv} "$EM/qc_summary.tsv"
    cp ~{qc_overview_png} "$EM/qc_overview.png"
    cp ~{signature_contributions_tsv} "$EM/mutational_signatures/signature_contributions.tsv"
    cp ~{cosine_similarity_tsv} "$EM/mutational_signatures/reconstruction_cosine_similarity.tsv"
    cp ~{mutation_spectrum_png} "$EM/mutational_signatures/mutation_spectrum_96.png"
    cp ~{signature_barplot_png} "$EM/mutational_signatures/signature_contribution_barplot.png"
    for f in ~{sep=" " lineage_outputs}; do cp "$f" "$EM/lineage/"; done
    cp ~{bootstrap_support_report} "$EM/lineage_verification/bootstrap_support_report.tsv"
    cp ~{verification_summary_json} "$EM/lineage_verification/verification_summary.json"

    python3 ~{build_dashboard_script} --extra-methods-dir "$EM" --donor-id ~{donor_id} \
      --multiqc-relpath ~{multiqc_relpath} --output dashboard.html
  >>>
  output { File dashboard_html = "dashboard.html" }
  requirements {
    container: resolveome_align_container
    cpu: 1
    memory: "4 GB"
  }
}
