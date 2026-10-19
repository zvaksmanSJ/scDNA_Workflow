## Cohort-wide lineage reconstruction (Sequoia, GPLv3, free, standalone from
## github.com/TimCoorens/Sequoia -- NOT the licensed BaseJumper/BioSkryb
## container) plus independent bootstrap/resample verification via
## ape::dist.topo (Robinson-Foulds distance). Sprocket port of
## ../../wdl/tasks/lineage.wdl.
##
## KEY DIFFERENCE FROM THE CROMWELL VERSION: the Cromwell version activates
## TWO conda environments in one task (`resolveome-align` for python3/
## bcftools, then prepends `mutsig`'s bin/ to PATH for R/Sequoia) -- Sprocket
## runs the whole command block inside exactly ONE container per task, so
## both toolsets are combined into a single `lineage_container` image (see
## ../containers/lineage.Dockerfile) instead of being two separate conda
## environments.
version 1.2

task BuildLineageTree {
  input {
    Array[File] final_annotated_vcfs   # one per cell, *.deepvariant.raw.vep.snpeff.vcf.gz
    Array[String] sample_ids
    File lineage_track_script
    File sequoia_build_phylogeny_script
    File mpboot_binary
    String donor_id = "cohort"
    String gender = "male"
    Int mpboot_bootstrap = 1000
    Int max_shared_sites = 20000
    Int random_seed = 42
    String lineage_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-lineage:latest"
  }
  command <<<
    set -euo pipefail
    mkdir -p variants_staged
    python3 - <<'PYEOF'
import shutil, os
vcfs = """~{sep="\n" final_annotated_vcfs}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, vcf in zip(samples, vcfs):
    d = os.path.join("variants_staged", sample)
    os.makedirs(d, exist_ok=True)
    shutil.copy(vcf, os.path.join(d, f"{sample}.deepvariant.raw.vep.snpeff.vcf.gz"))
PYEOF

    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    MPBOOT_DIR=$(mktemp -d)
    cp ~{mpboot_binary} "$MPBOOT_DIR/"
    chmod +x "$MPBOOT_DIR"/*

    python3 ~{lineage_track_script} \
      --variants-dir variants_staged \
      --manifest manifest.tsv \
      --output-dir lineage_out \
      --donor-id ~{donor_id} \
      --gender ~{gender} \
      --mpboot-bootstrap ~{mpboot_bootstrap} \
      --sequoia-script ~{sequoia_build_phylogeny_script} \
      --mpboot-path "$MPBOOT_DIR/" \
      --max-shared-sites ~{max_shared_sites} \
      --random-seed ~{random_seed}
  >>>
  output {
    Array[File] lineage_outputs = glob("lineage_out/*")
    File nv_matrix = glob("lineage_out/*_NV.tsv")[0]
    File nr_matrix = glob("lineage_out/*_NR.tsv")[0]
  }
  requirements {
    container: lineage_container
    cpu: 4
    memory: "32 GB"
  }
}

task VerifyLineage {
  input {
    Array[File] final_annotated_vcfs
    Array[String] sample_ids
    File lineage_track_script
    File verify_lineage_script
    File verify_lineage_r_script
    File sequoia_build_phylogeny_script
    File mpboot_binary
    Array[File] primary_lineage_outputs   # all files from BuildLineageTree, for the tree file
    String donor_id = "cohort"
    String gender = "male"
    Int mpboot_bootstrap = 1000
    Int max_shared_sites = 20000
    Int rerun_seed = 123
    Float bootstrap_threshold = 70.0
    String lineage_container = "docker://ghcr.io/zvaksmanSJ/scdna-workflow-lineage:latest"
  }
  # Independent verification of the primary tree built above: re-runs the
  # same lineage-tracking pipeline with a different random seed (rerun_seed)
  # and resampled sites, then compares the resulting tree topology against
  # the primary tree via ape::dist.topo (Robinson-Foulds distance) to report
  # per-branch bootstrap support -- a genuine sanity check that the tree
  # isn't an artifact of one particular site-sampling draw.
  command <<<
    set -euo pipefail
    mkdir -p variants_staged primary_lineage
    python3 - <<'PYEOF'
import shutil, os
vcfs = """~{sep="\n" final_annotated_vcfs}""".splitlines()
samples = """~{sep="\n" sample_ids}""".splitlines()
for sample, vcf in zip(samples, vcfs):
    d = os.path.join("variants_staged", sample)
    os.makedirs(d, exist_ok=True)
    shutil.copy(vcf, os.path.join(d, f"{sample}.deepvariant.raw.vep.snpeff.vcf.gz"))
primary = """~{sep="\n" primary_lineage_outputs}""".splitlines()
for f in primary:
    shutil.copy(f, os.path.join("primary_lineage", os.path.basename(f)))
PYEOF

    printf "sample_id\n" > manifest.tsv
    for s in ~{sep=" " sample_ids}; do printf "%s\n" "$s" >> manifest.tsv; done

    MPBOOT_DIR=$(mktemp -d)
    cp ~{mpboot_binary} "$MPBOOT_DIR/"
    chmod +x "$MPBOOT_DIR"/*

    python3 ~{verify_lineage_script} \
      --lineage-track-script ~{lineage_track_script} \
      --variants-dir variants_staged \
      --manifest manifest.tsv \
      --primary-lineage-dir primary_lineage \
      --donor-id ~{donor_id} \
      --gender ~{gender} \
      --mpboot-bootstrap ~{mpboot_bootstrap} \
      --sequoia-script ~{sequoia_build_phylogeny_script} \
      --mpboot-path "$MPBOOT_DIR/" \
      --max-shared-sites ~{max_shared_sites} \
      --rerun-seed ~{rerun_seed} \
      --bootstrap-threshold ~{bootstrap_threshold} \
      --verify-r-script ~{verify_lineage_r_script} \
      --output-dir verification_out
  >>>
  output {
    File bootstrap_support_report = "verification_out/bootstrap_support_report.tsv"
    File verification_summary_json = "verification_out/verification_summary.json"
  }
  requirements {
    container: lineage_container
    cpu: 4
    memory: "32 GB"
  }
}
