# workflow_DNA/

The ResolveOME single-cell DNA workflow, in two complementary, fully
cross-referenced forms. **For the full rationale, specification, and
history, see [README_RFC_DNA_Workflow.md](./README_RFC_DNA_Workflow.md)** --
this file is just a quick navigation map of the folder layout below.

## Folder layout

| Path | What's in it |
|---|---|
| **`main.wdl`** | The general-purpose WDL workflow's entry point. Run this (with Cromwell -- see `wdl/README.md`) for a free-tools-only, portable version of the pipeline. |
| **`README_RFC_DNA_Workflow.md`** | The full RFC: motivation, design rationale, dependencies, complete tool-by-tool command reference, known issues/fixes, and verification status for both pipelines below. |
| `wdl/` | `main.wdl`'s supporting task modules (`tasks/*.wdl`), WDL-support Python scripts (`scripts/`), the required Cromwell config (`cromwell_local.conf`), and its own usage `README.md`. |
| `sprocket/` | A Sprocket-native (https://github.com/stjude-rust-labs/sprocket, St. Jude's own next-gen WDL engine) port of the same workflow -- its own `main.wdl`, task files, 4 new Dockerfiles, `sprocket.toml`, and its own `README.md` explaining exactly how and why it differs from the Cromwell-oriented `main.wdl`/`wdl/` above. **Written but not run, built, or pushed anywhere** -- see that README's own "What this port does NOT do" section. |
| `setup/` | Every setup/provisioning script for *both* pipelines -- most notably `setup.sh`, the WDL pipeline's free-tools-only installer (conda envs, reference downloads, containers). |
| `wdl_runtime/` | Where `setup/setup.sh` installs everything it downloads by default (refs, conda envs, containers, tools) -- kept entirely separate from `resolveome_align/` below. Empty until you run `setup.sh`. |
| `extra_methods/` | Calculation scripts (Python/R) used by **both** pipelines: Beagle phasing, Control-FREEC LOH/ploidy, Sequoia lineage + verification, MAPD/FDR/artifact-panel QC, mutational signatures, and the final tabbed dashboard builder. |
| `python_lsf_pipeline/` | The original, **still-active production** Python/LSF pipeline (`step1_*.py` .. `step7_*.py`, orchestrated by `run_all_resolveome.py`/`.sh`), its tests, and its own `logs/`. |
| `docs/` | The workflow diagram (`workflow_diagram.{png,svg,dot}`). |
| `resolveome_align/` | **Left exactly as-is** (per explicit instruction): this project's real reference genome, conda environment, and container bundle for its own cohort run. Not touched or duplicated by anything above. |
| `Backup_code/` | Historical snapshots of older script versions -- not active workflow code. |
| `dv_compat_test/` | Scratch/results directory for `setup/test_dv_checkpoint_compat.sh`. |

## Which pipeline do I want?

- **Running this project's own real cohort, on this cluster, with its
  licensed BioSkryb DeepVariant model?** Use `python_lsf_pipeline/` (the
  original pipeline) -- see `README_RFC_DNA_Workflow.md`'s
  ["Workflow (Python/LSF production pipeline)"](./README_RFC_DNA_Workflow.md#workflow-pythonlsf-production-pipeline)
  section.
- **Deploying this workflow somewhere new, or want a free/portable path with
  no licensed dependencies?** Use `main.wdl` -- see `wdl/README.md` for
  exact run commands, and `setup/setup.sh` to provision everything it needs.
- **Running on St. Jude's own Sprocket engine instead of Cromwell?** Use
  `sprocket/main.wdl` -- see `sprocket/README.md` for exactly how and why it
  differs from the Cromwell version above (container-per-task model means
  4 new Docker images are needed; written but not built/run).

Both pipelines share `extra_methods/` for every phasing/LOH/lineage/QC/
mutational-signature/dashboard calculation; neither one duplicates that
logic.

## Where are my results?

Both pipelines produce a lot of per-cell, per-stage output; **the dashboard
HTML is always the right place to start** -- everything below is a map for
when you need a specific underlying file instead.

### WDL pipeline (`main.wdl` or `sprocket/main.wdl`)

Run `wdl/scripts/organize_outputs.py` against the engine's `outputs.json`
(see `wdl/README.md`) to get this tree (works unchanged for Sprocket's
output too):

```
<dest>/
  dashboards/dashboard.html              <- OPEN THIS FIRST
  cohort/
    qc_gate/qc_gate_summary.tsv          <- every cell's PASS/Borderline/FAIL, one row each
    lineage/                             <- final tree, NR/NV matrices, bootstrap support
    qc/, phasing/, loh_ploidy/, mutational_signatures/   <- cohort-wide summary tables/figures
  per_cell/<sample_id>/
    qc_gate/<sample_id>.qc_gate_explanation.txt   <- WHY this specific cell passed/failed
    variants/<sample_id>...vcf.gz          <- final DeepVariant+VEP+SnpEff calls (passed cells only)
    cnv/, sv/                              <- CNVkit + Manta calls (passed cells only)
    phasing/, loh_ploidy/                  <- Beagle phasing + Control-FREEC LOH/ploidy (passed cells only)
    alignment/, qc/                        <- dedup BAM, flagstat/idxstats, fastp/fastqc/Picard/preseq (ALL cells)
  logs/                                   <- safe to delete; not a result
  run_scripts/                            <- inputs.json/outputs.json used; safe to delete once audited
```

**The most important handful of files, if you only look at a few:**

| File | Why |
|---|---|
| `dashboards/dashboard.html` | The single-page overview: QC, phasing, LOH, lineage, mutational signatures, all in one tabbed report. |
| `cohort/qc_gate/qc_gate_summary.tsv` | Every cell's PASS/Borderline/FAIL status and whether its downstream analysis ran, in one glance. |
| `per_cell/<sample_id>/qc_gate/<sample_id>.qc_gate_explanation.txt` | For any specific FAILed (or Borderline) cell, exactly which of the 5 QC metrics tripped and by how much. |
| `per_cell/<sample_id>/variants/<sample_id>...vcf.gz` | That cell's final annotated small-variant calls -- only exists for cells that passed the QC gate. |
| `cohort/lineage/` | The final Sequoia lineage tree + its bootstrap-support verification report. |

### Python/LSF production pipeline

Run `extra_methods/organize_final_results.py` (symlink-based, non-
destructive, safe to re-run) to get:

```
DNA_results/final/
  dashboards/dashboard.html              <- OPEN THIS FIRST
  components/<name>/                     <- one symlinked folder per stage
                                             (variants_per_cell, cnv_cnvkit_per_cell,
                                             sv_manta_per_cell, bam_per_cell, qc_raw/,
                                             phasing_per_cell, loh_ploidy_per_cell,
                                             lineage, mutational_signatures, ...)
```

Step2's own `passed_cells_only_manifest.tsv` / `failed_cells_manifest.tsv` /
`cell_qc_metrics.tsv` (written directly under that run's Step2 output
directory, not under `final/`) are the production pipeline's equivalent of
the WDL gate's `qc_gate_summary.tsv` above -- see
["Cell selection and interpretation of QC"](./README_RFC_DNA_Workflow.md#cell-selection-and-interpretation-of-qc)
in the RFC for its hard-threshold gate (different formula from the WDL
pipeline's BaseJumper-derived gate -- see `wdl/README.md`'s "Per-cell QC
gate" section).

