## Cohort-wide MultiQC aggregation over all per-cell QC/variant/CNV outputs.
## Sprocket port of ../../wdl/tasks/qc_aggregate.wdl.
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: no `conda activate` -- MultiQC
## must already be on PATH inside the container named by
## `resolveome_align_container` (see
## ../containers/resolveome-align.Dockerfile).
version 1.2

task MultiQC {
  input {
    Array[File] input_files   # flattened list of every QC/log file to scan
    String resolveome_align_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest"
  }
  command <<<
    set -euo pipefail

    mkdir -p multiqc_inputs multiqc_out
    # MultiQC recognizes tool output by file name/content, not location, so
    # staging everything into one flat scan directory (via symlink) is safe.
    i=0
    for f in ~{sep=" " input_files}; do
      i=$((i+1))
      ln -s "$f" "multiqc_inputs/file_${i}_$(basename "$f")"
    done

    # Scan the flat staging directory and write one combined HTML report +
    # a multiqc_data/ directory of parsed per-tool TSVs (general_stats.txt
    # is the one qc_summary.py later reads -- see extra_qc_dashboard.wdl).
    multiqc -f -o multiqc_out multiqc_inputs
  >>>
  output {
    File report_html = "multiqc_out/multiqc_report.html"
    Array[File] data = glob("multiqc_out/multiqc_data/*")
  }
  requirements {
    container: resolveome_align_container
    cpu: 2
    memory: "8 GB"
  }
}
