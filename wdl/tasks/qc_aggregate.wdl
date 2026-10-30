## Cohort-wide MultiQC aggregation over all per-cell QC/variant/CNV outputs.
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: its `input {}`
## fields are filled in by main.wdl's `call qc_t.MultiQC as multiqc { ... }`
## block, using values from main.wdl's own single inputs.json -- there is no
## separate JSON file for this task in normal cohort-wide use (see main.wdl's
## header comment for the full explanation, and for how to run it standalone
## for debugging if you ever need to).
version 1.0

# Wraps MultiQC (free, MIT-licensed, bioconda). Runs ONCE for the whole
# cohort (not scattered per cell) -- main.wdl passes it every cell's
# fastp/fastqc/picard QC files plus SnpEff reports, flattened into one list.
task MultiQC {
  input {
    Array[File] input_files   # flattened list of every QC/log file to scan
    String conda_env = "wdl-resolveome-align"
  }
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}

    mkdir -p multiqc_inputs multiqc_out
    # MultiQC recognizes tool output by file name/content, not location, so
    # staging everything into one flat scan directory (via symlink) is safe
    # and keeps the WDL file-localization model simple.
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
  runtime {
    cpu: 2
    memory: "8 GB"
  }
}

