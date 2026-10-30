## Per-cell pass/Borderline/FAIL QC gate, ported from BaseJumper's own
## composite-score formula (see extra_methods/cell_qc_gate.py's own
## docstring for the full citation + the exact 5 metrics/cutoffs/scoring
## rule). Runs immediately after AlignTrimQC and BEFORE every expensive
## downstream stage (DeepVariant/VEP/SnpEff, Manta, Beagle phasing,
## Control-FREEC LOH/ploidy, Sequoia lineage) -- main.wdl wraps all of those
## in `if (gate.qc_pass)` so a FAILing cell's analysis stops right here,
## exactly as BaseJumper's own basej-dnaqc gate runs before basej-wgs/
## basej-somatic/basej-lineage.
##
## THIS FILE IS IMPORTED BY main.wdl, NOT RUN STANDALONE: its `input {}`
## fields are filled in by main.wdl's `call qc_gate_t.CellQcGate as gate
## { ... }` block (inside the same per-cell scatter as AlignTrimQC and
## CnvkitCall), using values from main.wdl's own single inputs.json.
version 1.0

# Runs a lightweight CNVkit pass purely to get MAPD/SKEW for the gate (the
# full-depth CNV/SV stage elsewhere in main.wdl runs its own, separate
# CnvkitCall for the real analysis -- see cnv_sv.wdl's own docstring and
# cell_qc_gate.py's module docstring for why this mirrors BaseJumper's own
# architecture rather than being wasted duplicate work).
task CellQcGate {
  input {
    String sample_id
    File bam
    File bam_bai
    File ref_fasta
    File ref_fasta_fai
    File? cnv_reference_cnn
    File idxstats
    File alignment_summary_metrics
    File preseq_result
    File cell_qc_gate_script
    # NOTE: must keep a decimal point (3500000000.0, not 3500000000) -- a
    # bare integer-looking literal this large fails to parse at all under
    # Cromwell 87's WDL 1.0 grammar even for a Float-typed declaration
    # (see main.wdl's matching qc_gate_cutoff_preseq for the full story).
    Float cutoff_preseq = 3500000000.0
    Float cutoff_chimeras = 0.20
    Float cutoff_chrm = 0.20
    Float cutoff_mapd = 0.25
    Float cutoff_skew = 0.25
    Int threads = 4
    String conda_env = "wdl-resolveome-align"
  }
  command <<<
    set -euo pipefail
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ~{conda_env}
    mkdir -p cnv_for_gate

    if [ -n "~{cnv_reference_cnn}" ]; then
      cnvkit.py batch ~{bam} \
        --method wgs \
        -r ~{cnv_reference_cnn} \
        --output-dir cnv_for_gate \
        --processes ~{threads}
    else
      cnvkit.py batch ~{bam} \
        --method wgs \
        --fasta ~{ref_fasta} \
        --output-dir cnv_for_gate \
        --processes ~{threads}
    fi

    sample_base=$(basename ~{bam} .bam)
    cnr="cnv_for_gate/${sample_base}.cnr"
    cns="cnv_for_gate/${sample_base}.cns"
    call_cns="cnv_for_gate/~{sample_id}.qc_gate.call.cns"
    cnvkit.py call "$cns" -o "$call_cns"

    python3 ~{cell_qc_gate_script} \
      --sample-id ~{sample_id} \
      --idxstats ~{idxstats} \
      --alignment-summary-metrics ~{alignment_summary_metrics} \
      --preseq-file ~{preseq_result} \
      --cnr "$cnr" \
      --call-cns "$call_cns" \
      --cutoff-preseq ~{cutoff_preseq} \
      --cutoff-chimeras ~{cutoff_chimeras} \
      --cutoff-chrm ~{cutoff_chrm} \
      --cutoff-mapd ~{cutoff_mapd} \
      --cutoff-skew ~{cutoff_skew} \
      --output-tsv "~{sample_id}.qc_gate.tsv" \
      --output-explanation "~{sample_id}.qc_gate_explanation.txt" \
      --output-pass-flag "~{sample_id}.qc_gate_status.txt"
  >>>
  output {
    File qc_gate_tsv = "~{sample_id}.qc_gate.tsv"
    File qc_gate_explanation = "~{sample_id}.qc_gate_explanation.txt"
    String qc_status = read_string("~{sample_id}.qc_gate_status.txt")
    Boolean qc_pass = qc_status != "FAIL"
  }
  runtime {
    cpu: threads
    memory: "16 GB"
  }
}
