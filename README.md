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

## Getting started (WDL + Cromwell path)

This section covers `main.wdl` (the free, portable, deployable path) end to
end: what you need before you start, how to install everything this
workflow needs, and the fastest path to a running job. If you're running
the original production pipeline on this project's own cluster instead, see
`README_RFC_DNA_Workflow.md`'s
["Workflow (Python/LSF production pipeline)"](./README_RFC_DNA_Workflow.md#workflow-pythonlsf-production-pipeline)
section; if you want Sprocket instead of Cromwell, see `sprocket/README.md`
(same setup.sh prerequisites below, plus Sprocket itself and 4 Docker images
it needs that Cromwell doesn't).

### 1. Prerequisites (not installed by `setup/setup.sh`)

`setup/setup.sh` installs everything that's reasonably automatable --
conda environments, the reference genome, containers, Beagle, Sequoia, and
so on. A handful of things need to already be on your system first:

| Requirement | Why you need it | Get it |
|---|---|---|
| Linux (or WSL2 on Windows) | this workflow was built/tested on an HPC Linux cluster; Apptainer/Singularity require Linux | -- |
| conda, Miniconda, or Miniforge | `setup.sh` creates every conda environment with `conda env create` | https://docs.conda.io/en/latest/miniconda.html |
| Apptainer or Singularity | runs the DeepVariant/VEP/Manta containers | https://apptainer.org/docs/admin/main/installation.html |
| `git`, `curl`, `tar`, `gzip`/`gunzip`, `unzip` | `setup.sh` uses these to download/extract references and clone Sequoia | almost always preinstalled on Linux; otherwise `apt install git curl tar gzip unzip` (Debian/Ubuntu) or `yum install ...` (RHEL/CentOS) |
| Java 11 or newer (JRE or JDK) | required to run both Cromwell (the WDL engine) and Beagle (phasing -- itself a `.jar`) | https://adoptium.net/ (OpenJDK builds), or your distro's `openjdk-11-jre` / `openjdk-17-jre` package |
| Cromwell | the WDL execution engine -- **not installed by `setup.sh`** | download a release jar from https://github.com/broadinstitute/cromwell/releases (this workflow was built and validated against Cromwell 87) |
| Disk space, order of 100+ GB | reference genome + VEP cache + 1000 Genomes phasing panel, before any run outputs | varies mainly with the VEP cache and phasing-panel size; plan for more as you add cells/outputs |

### 2. Download and install

```bash
# 1. Clone this repo
git clone https://github.com/zvaksmanSJ/scDNA_Workflow.git
cd scDNA_Workflow

# 2. Run the setup script -- idempotent (every step is skip-if-already-
#    present) and safe to re-run. Installs into ./wdl_runtime/ by default;
#    see setup/setup.sh's own header comment for its BASE_DIR/THREADS/
#    ENSEMBL_RELEASE/VEP_CACHE_VERSION positional args and its
#    SKIP_<SECTION>=1 env-var overrides (SKIP_CONDA, SKIP_REFERENCE,
#    SKIP_KNOWN_SITES, SKIP_CONTAINERS, SKIP_PHASING, SKIP_FREEC_REFS,
#    SKIP_LINEAGE).
./setup/setup.sh
```

What this installs, all free/open (no licensed or private dependencies):
- 4 conda environments: `resolveome-align`, `freec`, `manta_env`, `mutsig`
- a GRCh38 reference FASTA (Ensembl) + bwa-mem2/samtools/Picard indices
- GATK's public legacy known-sites bundle
- stock `google/deepvariant` and `ensemblorg/ensembl-vep` containers, pulled
  with Apptainer/Singularity (not any licensed BioSkryb image)
- Beagle 5.5 + its public PLINK genetic maps + the 1000 Genomes GRCh38
  high-coverage phasing panel
- Sequoia + MPBoot (standalone, GPL-licensed lineage-tree building)

It does **not** install Cromwell, Java, conda, or Apptainer/Singularity
themselves (see Prerequisites above). It also does not provide the
BioSkryb-private DeepVariant checkpoint or the `preseq-bioskryb`/`bam2mr`
binaries this project's own real cohort run uses -- no public source exists
for either; `main.wdl` falls back to DeepVariant's free bundled model if you
don't supply a private checkpoint, and `wdl/README.md`/`sprocket/README.md`
document the preseq gap explicitly (you'd need to supply those binaries
yourself to use them).

> **If you override `ENSEMBL_RELEASE` away from its default (111) when
> running `setup.sh`**, also update `inputs.json`'s `snpeff_db` to match
> (`"GRCh38.<ENSEMBL_RELEASE>"`) -- `setup.sh` downloads the SnpEff database
> name-matched to whatever `ENSEMBL_RELEASE` you give it, and `main.wdl`'s
> `snpeff_db` input isn't derived from that automatically.

### 3. Quick start

```bash
# 1. Build your inputs.json. Start from the complete example and fill in
#    your own cohort (sample_ids/r1_fastqs/r2_fastqs/donor_id/gender) plus
#    every path setup.sh just produced (it prints each one as it runs) --
#    or use wdl/scripts/build_inputs_json.py to merge a FASTQ manifest TSV
#    into the example automatically (see wdl/README.md).
cp wdl/example_inputs.json inputs.json
# ... edit inputs.json ...

# 2. Run it with Cromwell (from the repo root -- main.wdl's own `import`
#    statements are written relative to that location)
java -Dconfig.file=wdl/cromwell_local.conf -jar cromwell-87.jar run \
  main.wdl -i inputs.json -m outputs.json

# 3. Reorganize the outputs into a clean, browsable results tree (see
#    "Where are my results?" below for what's in it)
python3 wdl/scripts/organize_outputs.py \
  --outputs-json outputs.json \
  --sample-ids cell1 cell2 ... \
  --dest /path/to/results \
  --inputs-json inputs.json

# 4. Open the dashboard
xdg-open /path/to/results/dashboards/dashboard.html   # or `open` on macOS
```

For a real multi-cell cohort at LSF scale (not just this quick-start local
run), see `wdl/README.md`'s "Running at LSF scale (production)" section.

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

