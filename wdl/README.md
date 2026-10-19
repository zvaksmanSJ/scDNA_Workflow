# ResolveOME single-cell DNA workflow -- WDL version

A general-purpose WDL 1.0 port of the Python/LSF pipeline in
`../python_lsf_pipeline/`, covering every stage: alignment+QC, DeepVariant/
VEP/SnpEff small variants, CNVkit/Manta CNV+SV, Beagle phasing, Control-FREEC
LOH+ploidy, Sequoia lineage+verification, and all `../extra_methods/`
calculations (MAPD/FDR/artifact-panel QC, phasing/LOH/ploidy/QC summaries,
mutational signatures, and the final tabbed dashboard).

Unlike the original pipeline, this version starts from already-prepared,
already-merged per-cell FASTQs (no multi-lane FASTQ-combining step), uses
only free/open tools and resources, and is meant to be runnable at any site
(not just this specific HPC cluster) -- see `../setup/setup.sh`.

## Layout

```
workflow_DNA/
  main.wdl                        -- top-level workflow (lives at workflow_DNA
                                      root, NOT inside wdl/ -- see note below)
  wdl/
    tasks/*.wdl                    -- one file per pipeline stage (8 files)
    scripts/make_region_chunks.py  -- DeepVariant region-chunking helper
    scripts/organize_outputs.py    -- post-run output reorganizer (see below)
    cromwell_local.conf            -- REQUIRED Cromwell config to actually
                                       run this (see "Why Cromwell" below)
  setup/
    setup.sh                       -- free-tools installer for this workflow
  wdl_runtime/                     -- setup.sh's default install location
                                      (refs/envs/containers/tools it
                                      downloads; separate from
                                      resolveome_align/, which stays untouched)
```

`main.wdl` deliberately sits at `workflow_DNA/` root (a sibling of this
`wdl/` folder, not inside it), so the folder structure reads as: top-level
entry point (`main.wdl`) + its supporting task library (`wdl/`) + its
installer (`setup/`) + its runtime data (`wdl_runtime/`).


## Validation performed this session

- `miniwdl check main.wdl` (run from `workflow_DNA/` root, miniwdl 1.16.0):
  **zero real errors**, only cosmetic WDL-1.1-style deprecation suggestions
  and 5 intentionally-unused bwa-mem2 index inputs (documented in
  `tasks/align_qc.wdl`; they exist only so the WDL engine stages those index
  files alongside the reference FASTA).
- A genuine local execution smoke test: `tasks/align_qc.wdl`'s `AlignTrimQC`
  task was run through Cromwell 87 + `cromwell_local.conf` against a tiny
  synthetic 50kb reference and ~2000 synthetic read pairs. It executed the
  *real* tools (bwa-mem2, fastp, fastqc, samtools, Picard, preseq) and
  produced a valid BAM (100% mapped, 99.95% properly paired) and real Picard
  `wgs_metrics`. This is a mechanics/plumbing test, not a full 15-cell
  production run -- see "What was NOT run" below.

## Why Cromwell, not miniwdl, to actually run this workflow

`miniwdl check` is excellent for *static* validation (type-checking, syntax,
unused-input detection) and was used throughout development. However,
**miniwdl cannot execute** this workflow: every miniwdl task runs inside its
own container (Docker/Singularity/Podman), with no "run directly on the
host" mode. This workflow's design deliberately does NOT containerize every
task -- only DeepVariant/VEP/Manta are containerized (via `apptainer exec`
calls *inside* their own command blocks); every other task assumes the named
conda environments from `../setup/setup.sh` (`resolveome-align`, `freec`,
`manta_env`, `mutsig`) are already on the host, exactly like the real
LSF/Python pipeline this ports from.

Cromwell's "Local" backend runs each task's command directly via `/bin/bash`
on the host, matching this design. **You must pass `cromwell_local.conf`**
(or an equivalent config) when running with Cromwell, and see the comments
in that file for two non-obvious gotchas that were hit and fixed during
validation (a `submit-docker` template quirk, and a `memory` runtime-attribute
type quirk).

## Running

```bash
# 1. One-time setup -- installs into ../wdl_runtime/ by default (see
#    ../setup/setup.sh for full options)
../setup/setup.sh

# 2. Build inputs.json pointing at the files setup.sh produced (see each
#    task's `input {}` block, or `../main.wdl`'s top-level `input {}` block,
#    for the full parameter list)

# 3. Run with Cromwell (download from https://github.com/broadinstitute/cromwell/releases)
#    -- run this from workflow_DNA/ root, since main.wdl's own `import`
#    statements are written relative to that location
cd ..   # i.e. workflow_DNA/
java -Dconfig.file=wdl/cromwell_local.conf -jar cromwell-87.jar run \
  main.wdl -i inputs.json -m outputs.json

# 4. Reorganize the outputs into a clean, browsable results tree
python3 wdl/scripts/organize_outputs.py \
  --outputs-json outputs.json \
  --sample-ids cell1 cell2 ... \
  --dest /path/to/results \
  --inputs-json inputs.json
```

`scripts/organize_outputs.py` lays every output out under `per_cell/<id>/`
(alignment/qc/variants/cnv/sv/phasing/loh_ploidy) and `cohort/`
(qc/phasing/loh_ploidy/mutational_signatures/lineage), with `dashboards/` and
`run_scripts/` kept separate -- see its own `--help` / docstring for details.

## Building inputs.json for a real cohort (many cells)

Hand-editing `sample_ids`/`r1_fastqs`/`r2_fastqs` JSON arrays for a dozen-plus
cells is tedious and error-prone. `scripts/build_inputs_json.py` instead
merges a simple two-column-plus-ID manifest TSV into a "base" inputs.json
(everything else -- reference, containers, Beagle panel, extra_methods
script paths -- filled in once, e.g. starting from `example_inputs.json`):

```bash
# manifest.tsv, tab-separated, one row per cell:
#   cell1	/path/to/cell1_R1.fastq.gz	/path/to/cell1_R2.fastq.gz
#   cell2	/path/to/cell2_R1.fastq.gz	/path/to/cell2_R2.fastq.gz
#   ...

python3 wdl/scripts/build_inputs_json.py \
  --manifest manifest.tsv \
  --base wdl/example_inputs.json \
  --donor-id my_cohort \
  --output inputs.json
```

This only replaces the 3 FASTQ-list keys; every other key in `--base` is
copied through unchanged, so you still need to fill in `example_inputs.json`
(or your own base file) once per deployment with real paths from
`setup/setup.sh`'s output and all 22 autosome entries in the 3 phasing Maps.

## Running at LSF scale (production)

Cromwell's "Local" backend (`cromwell_local.conf`, used above) runs every
task as a local subprocess on whichever single host runs Cromwell itself --
fine for a smoke test or 1-2 cells, but not how you'd want to run a real
multi-cell cohort (no distribution across the cluster). For real scale, use
`cromwell_lsf.conf` instead -- it submits every task via `bsub`, matching
this cluster's own conventions (`-q standard -n <cpu> -M <mem> -R
"span[hosts=1] rusage[mem=<mem>]" -W <walltime>`, the same flags used
throughout `../extra_methods/*_lsf.py`):

```bash
java -Dconfig.file=wdl/cromwell_lsf.conf -jar cromwell-87.jar run \
  main.wdl -i inputs.json -m outputs.json
```

No changes to any `.wdl` file are needed -- only the `-Dconfig.file=` flag
changes between `cromwell_local.conf` and `cromwell_lsf.conf`. Every task
gets a 48-hour LSF walltime by default (comfortably more than any single
cell/chunk/chromosome unit of work should need, since the WDL already
parallelizes each stage down to that granularity); override per task with
`runtime { lsf_walltime: "72:00" }` in the relevant `tasks/*.wdl` file if a
specific stage genuinely needs longer. **`cromwell_lsf.conf` was written
against Cromwell's own documented LSF backend pattern and this cluster's
proven `bsub` conventions, but has NOT been live-tested against a real
`bsub` submission** (doing so would itself submit a real LSF job) -- if you
want extra confidence before a full cohort run, consider first trying it
against a minimal 1-2-cell `inputs.json`.

## What was NOT run this session (and why)

A full 15-cell production run through this WDL workflow was **not**
attempted. The handoff constraints for this session explicitly prohibit
submitting/altering production jobs, and re-running the full cohort through
a second pipeline implementation would duplicate the multi-day real-data run
already completed (and reported on) through the original Python/LSF
pipeline earlier this session. What *was* validated: static correctness
(`miniwdl check`), and genuine local mechanics (a real task executing real
tools on synthetic data, per above). Anyone adopting this workflow for a new
cohort should expect normal first-run debugging at full scale (resource
sizing, LSF backend wiring, container availability), the same as any new
pipeline deployment.
