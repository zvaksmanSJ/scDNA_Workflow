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

Both pipelines share `extra_methods/` for every phasing/LOH/lineage/QC/
mutational-signature/dashboard calculation; neither one duplicates that
logic.
