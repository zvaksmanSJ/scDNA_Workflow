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
    tasks/*.wdl                    -- one file per pipeline stage (9 files,
                                       including the per-cell QC gate --
                                       see "Per-cell QC gate" below)
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


## Per-cell QC gate (BaseJumper-derived composite score)

Right after Stage 1 (`AlignTrimQC`), `main.wdl` runs a new `CellQcGate` task
(`wdl/tasks/cell_qc_gate.wdl`) for every cell, **before** any of the
expensive downstream stages (DeepVariant/VEP/SnpEff, Manta, Beagle phasing,
Control-FREEC LOH/ploidy, Sequoia lineage). A cell that scores <=2/5 is
marked **FAIL**, and every one of those downstream stages is skipped
entirely for that cell -- not run, not silently stubbed out, actually
absent from the scatter for that cell's index.

This replaces the original `python_lsf_pipeline/step2` script's hard
thresholds (1,000,000 reads / 0.70 mapping rate / 0.60 proper-pair rate)
with the same composite-score formula BaseJumper's own `basej-dnaqc`
pipeline uses, so this workflow's pass/fail decision matches BaseJumper's
rather than inventing a third, different standard:

| # | Metric | Cutoff | Source |
|---|--------|--------|--------|
| 1 | `preseq_count` (library complexity) | `> 3,500,000,000` | AlignTrimQC's preseq output |
| 2 | `PCT_CHIMERAS` | `< 0.20` | AlignTrimQC's Picard `alignment_summary_metrics` |
| 3 | chrM read fraction | `< 0.20` | AlignTrimQC's `idxstats` |
| 4 | `MAPD_CNV_Log2` (CNV bin noise) | `< 0.25` | a lightweight CNVkit pass CellQcGate runs itself |
| 5 | `SKEW_CNV` (CNV directional bias) | `< 0.25` | the same lightweight CNVkit pass |

`CompositeScore` = how many of the 5 pass (0-5). `>=4` is **PASS**, `3` is
**Borderline** (still runs every downstream stage, just flagged for
review), `<=2` is **FAIL**. A cell with zero total reads is always FAIL
regardless of score. See
[`../extra_methods/cell_qc_gate.py`](../extra_methods/cell_qc_gate.py)'s own
module docstring for:
- the exact citation (BaseJumper's `basej-dnaqc` pipeline, vendored for
  reference under `../resolveome_align/basej-public-pipelines/basej-dnaqc/`)
- why `SKEW_CNV` here is an *adaptation* of BaseJumper's formula (same
  metric definition and cutoff, but computed from CNVkit's `.cnr`/`.call.cns`
  output instead of Ginkgo's, since this workflow uses CNVkit, not Ginkgo)
- why the CNVkit pass this gate runs is a genuinely separate, lightweight
  call from the full-depth CNV/SV stage later in the pipeline (mirrors
  BaseJumper's own architecture: `basej-dnaqc` runs its own dedicated CNV
  call purely for QC scoring, distinct from `basej-wgs`/`basej-somatic`'s
  real CNV/SV calling)

Every cell -- pass, Borderline, or FAIL -- gets a `qc_gate_explanation.txt`
(workflow output `qc_gate_explanations`) listing each of the 5 metrics'
actual value, its cutoff, and whether it passed, plus a final verdict line.
For a FAILed cell this is the first place to look to understand why.
Workflow outputs `qc_statuses` and `qc_passed_sample_ids` give the
cohort-wide pass/fail list and the filtered sample ID list, respectively;
cohort-wide outputs (small variants, phasing, LOH/ploidy, lineage) only
ever contain PASSED cells, while alignment/QC-stage outputs (`dedup_bams`,
`wgs_metrics`, `multiqc_report`, etc.) still cover every cell, pass or fail,
since QC visibility into *why* a cell failed is the whole point of keeping
those around.


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
conda environments from `../setup/setup.sh` (`wdl-resolveome-align`, `freec`,
`mutsig`) are already on the host, exactly like the real
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

This replaces the 3 FASTQ-list keys, AND regenerates the 3 phasing Maps
(`phasing_ref_panel_by_chrom`, `phasing_ref_panel_tbi_by_chrom`,
`genetic_map_by_chrom`) to a complete chr1-22 set by globbing the real
per-chromosome files `setup/setup.sh` downloads under `--wdl-runtime-dir`
(default: `wdl_runtime`, matching `setup.sh`'s own default `BASE_DIR`) --
**not** by copying whatever `--base` happens to contain for them.
`wdl/example_inputs.json` deliberately only illustrates those 3 keys with 2
example chromosomes (1, 2) to stay readable; copying that through unchanged
(the old behavior, still available via `--skip-autosome-regen` if you've
hand-built a complete `--base`) produces an `inputs.json` that runs cleanly
through alignment/variant-calling for every cell, then fails partway
through, at the phasing stage's chr3 scatter iteration -- after the
expensive stages already ran. Every other key in `--base` besides those 3
Maps and the FASTQ lists is still copied through unchanged, so you still
need to fill in `example_inputs.json` (or your own base file) once per
deployment with the rest of `setup/setup.sh`'s real output paths.

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
