# ResolveOME single-cell DNA workflow -- Sprocket port

A Sprocket-native (https://github.com/stjude-rust-labs/sprocket) adaptation
of the Cromwell-oriented general-purpose WDL workflow at
[`../main.wdl`](../main.wdl). Same pipeline, same stages, same bug fixes --
only the container-handling mechanics differ, because Sprocket's execution
model is genuinely different from Cromwell's.

**This port was written but not run, built, or pushed anywhere** -- see
"What this port does NOT do" at the bottom of this file before using it.

## Why a separate port, not just a new Cromwell config?

Everything else in this project's WDL rewrite (`../main.wdl`, `../wdl/`,
`../setup/setup.sh`) only ever needed a new *engine config* file to switch
from Cromwell's Local backend to its LSF backend (`../wdl/cromwell_lsf.conf`)
-- the `.wdl` files themselves never had to change. Sprocket is different:
it has a fundamentally different task-execution model that required real
changes to the WDL task definitions themselves, not just a new config file.
See the next section for exactly what and why.

## The one architectural difference that matters: one container per task

Cromwell's "Local"/"LSF" backends (used by `../main.wdl`) can run a task's
command block directly on the host, `conda activate`-ing a named
environment, and only specific tasks manually `apptainer exec` a second,
separate container mid-command for tools that need one (DeepVariant/VEP/
Manta). **Sprocket's execution model has no such option**: every task runs
inside exactly one container, chosen by that task's own
`requirements { container: ... }` value -- there is no "run bare on the
host" mode, and no way to invoke a second container from inside a task's
own command block. See https://sprocket.bio/concepts/containers for the
authoritative statement of this ("Sprocket runs each task inside a
container... When a task does not specify an image, Sprocket uses its own
default container").

This has two concrete consequences, both implemented in `./tasks/*.wdl`:

1. **Every tool stack this workflow used to reach via `conda activate` now
   needs an actual container image instead.** The Cromwell version's three
   conda environments (`resolveome-align`, `freec`, `mutsig`) become three
   new Dockerfiles under `./containers/` -- plus a fourth, `lineage.Dockerfile`,
   combining what used to be two separate conda environments activated in
   sequence within one Cromwell task (`BuildLineageTree`/`VerifyLineage`
   needed both `resolveome-align`'s python3/bcftools AND `mutsig`'s R/ape,
   reached via a PATH-prepend hack) into one image, since Sprocket only
   allows one container per task. **Building and pushing these four images
   is a required one-time setup step this port does NOT perform for you**
   -- see "What this port does NOT do" below.

2. **Tasks that used to run TWO toolsets in one command block had to be
   restructured.** Tracing exactly which of the 8 original task files
   actually did this (not just assumed) found only one real case:
   `cnv_sv.wdl`'s `MantaCall`, which manually `apptainer exec`s the Manta
   container, then falls back to a `conda activate`d `bcftools index` only
   if Manta didn't produce its own `.tbi`. This is now two tasks:
   `MantaCall` (runs entirely inside the Manta container, with a `File?`
   index output) and a new shared `IndexVcf` helper (`tasks/common.wdl`),
   called from `main.wdl` only via an `if (!defined(manta.manta_vcf_tbi))`
   block -- preserving the exact original fallback behavior rather than
   always blindly re-indexing. `variant_calling.wdl`'s `DeepVariantChunk`/
   `VepAnnotateChunk` also manually `apptainer exec`'d a container and then
   `conda activate`d for a trailing `bcftools index` -- but tracing
   `../main.wdl`'s actual data flow shows that index was **never consumed
   downstream** (`MergeVcfChunks` re-sorts and re-indexes every chunk
   itself regardless), so those two tasks (renamed `DeepVariantCall`/
   `VepAnnotateCall` here) simply dropped that unused step and run entirely
   inside their own tool container -- no splitting needed. `lineage.wdl`'s
   dual-conda-env tasks needed no splitting either, since WDL tasks can't
   call other tasks directly (only workflows can) -- building one combined
   `lineage` container was the correct fix, not a workflow-level restructure.

Every other task (`align_qc.wdl`, `cnv_sv.wdl`'s `CnvkitCall`, `phasing.wdl`,
`loh_ploidy.wdl`, `qc_aggregate.wdl`, `extra_qc_dashboard.wdl`,
`variant_calling.wdl`'s `PlanChunks`/`MergeVcfChunks`/`SnpEffAnnotate`) only
ever used ONE conda environment throughout its whole command block, so those
needed nothing more than a straightforward container swap (remove
`conda activate`, add `requirements { container: ... }`).

## WDL version: 1.0 -> 1.2

`requirements { container: ... }` is WDL 1.2 syntax (`runtime { docker: ... }`
is the pre-1.2 equivalent). Since this port needed every task file touched
anyway, every file here declares `version 1.2` -- a deliberate, one-way
upgrade from the Cromwell port's `version 1.0` (chosen there specifically
because Cromwell 87 cannot parse `version 1.1`/`1.2` declarations at all).
Sprocket has no such limitation.

## Layout

```
sprocket/
  main.wdl                 -- top-level workflow (this port's entry point)
  sprocket.toml             -- LSF+Apptainer backend config (see below)
  example_inputs.json       -- complete example inputs.json (see below)
  tasks/
    align_qc.wdl            -- container swap only
    variant_calling.wdl      -- PlanChunks/MergeVcfChunks/SnpEffAnnotate: swap only.
                                 DeepVariantCall/VepAnnotateCall: dropped their
                                 unused trailing bcftools-index step (see above)
    cnv_sv.wdl               -- CnvkitCall: swap only. MantaCall: split (see above)
    phasing.wdl              -- container swap only
    loh_ploidy.wdl            -- container swap only
    lineage.wdl               -- container swap only (one NEW combined image)
    qc_aggregate.wdl          -- container swap only
    extra_qc_dashboard.wdl    -- container swap only (2 different images)
    common.wdl                -- NEW: the shared IndexVcf helper task
  containers/
    resolveome-align.Dockerfile   -- NEW: replaces the `resolveome-align` conda env
    freec.Dockerfile              -- NEW: replaces the `freec` conda env
    mutsig.Dockerfile              -- NEW: replaces the `mutsig` conda env (MutationalSignatures only)
    lineage.Dockerfile             -- NEW: combines resolveome-align's python3/bcftools
                                       WITH mutsig's R/ape, PLUS Sequoia's own R
                                       dependencies pre-installed at build time
                                       (see that Dockerfile's own header comment)
```

## `sprocket.toml`

Configures Sprocket's LSF + Apptainer backend, based directly on
https://sprocket.bio/reference/backends/lsf and
https://sprocket.bio/guides/lsf (the official upstream reference), with
values tailored to match this project's own existing `bsub` conventions
(the same ones `../extra_methods/*_lsf.py` already use). **Several values
in this file are explicit placeholders you must edit before use** --
`default_lsf_queue.name`, the resource-limit values, and especially
`apptainer.image_cache_dir` (set to a real shared-filesystem path, or every
run will re-pull every container image from scratch) -- see the comments
inside the file itself for exactly which.

## `example_inputs.json`

A complete, working example covering every one of `main.wdl`'s 50 input
keys (annotated with `_comment_*` placeholder keys documenting where each
value comes from), structured exactly like `../wdl/example_inputs.json`
(the Cromwell version) with the differences described above: no more
`deepvariant_sif`/`vep_sif`/`manta_sif` `File` inputs or
`conda_env`/`freec_conda_env`/`r_conda_env` `String` inputs; new
`*_container` `String` inputs instead, defaulting to:

| Container | Default value | Status |
|---|---|---|
| `deepvariant_container` | `docker://google/deepvariant:1.6.1` | free, stock, public -- already used by `../setup/setup.sh` |
| `vep_container` | `docker://ensemblorg/ensembl-vep:release_111.0` | free, stock, public -- already used by `../setup/setup.sh` |
| `manta_container` | `docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6` | free, stock, public -- confirmed to exist on quay.io, pinned to the exact same build `../setup/setup.sh`'s `manta_env` conda recipe already uses |
| `resolveome_align_container` | `docker://ghcr.io/zvaksmanSJ/scdna-workflow-resolveome-align:latest` | **placeholder -- you must build + push this (see below)** |
| `freec_container` | `docker://ghcr.io/zvaksmanSJ/scdna-workflow-freec:latest` | **placeholder -- you must build + push this** |
| `lineage_container` | `docker://ghcr.io/zvaksmanSJ/scdna-workflow-lineage:latest` | **placeholder -- you must build + push this** |
| `mutsig_container` | `docker://ghcr.io/zvaksmanSJ/scdna-workflow-mutsig:latest` | **placeholder -- you must build + push this** |

(You can instead point any of these at a `file:///path/to/local.sif` if you
pre-pull/build and place a `.sif` on shared storage rather than using a
registry -- see https://sprocket.bio/concepts/containers for every
supported protocol.)

## How you would actually run this (not done as part of writing this port)

```bash
# 1. Build and push the 4 new container images this port needs (one-time;
#    NOT done by this port -- see containers/*.Dockerfile, especially
#    resolveome-align.Dockerfile's header comment about the bam2mr/
#    preseq-bioskryb gap)
docker build -f containers/resolveome-align.Dockerfile -t ghcr.io/<you>/scdna-workflow-resolveome-align:latest .
docker push ghcr.io/<you>/scdna-workflow-resolveome-align:latest
# ... repeat for freec.Dockerfile, mutsig.Dockerfile, lineage.Dockerfile

# 2. Install Sprocket (NOT done by this port)
#    https://github.com/stjude-rust-labs/sprocket/releases, or:
cargo install sprocket --locked

# 3. Run ../setup/setup.sh (same as the Cromwell path -- installs the
#    reference genome, Beagle, 1000G panel, Sequoia/MPBoot, etc. into
#    ../wdl_runtime/)
../setup/setup.sh

# 4. Build your real inputs.json from example_inputs.json (same manifest-
#    merging helper as the Cromwell path also works here, since both engines
#    accept the same JSON key convention):
python3 ../wdl/scripts/build_inputs_json.py \
  --manifest manifest.tsv --base example_inputs.json \
  --donor-id my_cohort --output inputs.json

# 5. Edit sprocket.toml's placeholder values (queue name, resource limits,
#    image_cache_dir -- see the file's own comments)

# 6. Run it (from workflow_DNA/, i.e. one directory up from this file)
cd ..
sprocket run --config sprocket/sprocket.toml sprocket/main.wdl @inputs.json -o /path/to/output

# 7. Reorganize the outputs (works unchanged against Sprocket's outputs.json
#    too -- see wdl/scripts/organize_outputs.py's own docstring)
python3 wdl/scripts/organize_outputs.py \
  --outputs-json /path/to/output/runs/ResolveOmeDnaWorkflow/<timestamp>/outputs.json \
  --sample-ids cell1 cell2 cell3 --dest /path/to/results
```

## Validation performed on this port

`miniwdl check main.wdl` (run from this `sprocket/` directory, miniwdl
1.16.0, which does support WDL 1.2): **zero real errors** -- only the same
cosmetic WDL-1.1-style placeholder-option deprecation notices as the
Cromwell port, plus the same 5 intentionally-unused bwa-mem2 index inputs.
The new conditional `if (!defined(manta.manta_vcf_tbi))` / `select_all` /
`flatten` wiring for the Manta index fallback type-checks cleanly.

**No execution smoke test was performed** for this port (unlike the
Cromwell port, which had a genuine local execution smoke test performed
against `tasks/align_qc.wdl`). Running anything through Sprocket requires
the 4 new container images to actually exist, which requires a `docker
build`/`docker push` -- real infrastructure actions this port deliberately
does not perform, consistent with "write the code, don't actually do it."

## What this port does NOT do

- **Does not install Sprocket.** (step 2 above)
- **Does not build or push any of the 4 new Docker images.** (step 1 above)
  `resolveome-align.Dockerfile` additionally has a known, documented gap:
  `bam2mr`/`preseq-bioskryb` have no public package source and must be
  supplied by you from this project's existing installation -- see that
  Dockerfile's own header comment. This gap already existed in the Cromwell
  port's `../setup/setup.sh` (which also never installs these two
  binaries); it is not something this Sprocket port introduces.
- **Does not run `../setup/setup.sh`.** (step 3 above -- same script as the
  Cromwell path; nothing Sprocket-specific needed here)
- **Does not run the workflow**, smoke-test it against real Sprocket/LSF
  execution, or submit any job of any kind.
- **Does not verify `sprocket.toml`'s placeholder values** (queue name,
  resource limits, `image_cache_dir`) against this cluster's real LSF
  configuration -- these are documented placeholders, not verified facts.
