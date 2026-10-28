# Licensed components -- what's licensed, why, and how to use your own

This project's **production** Python/LSF pipeline
(`python_lsf_pipeline/run_all_resolveome.py`) uses four components licensed
from BioSkryb Genomics (the company behind the BaseJumper platform and the
ResolveDNA single-cell WGA protocol this project's real cohort data comes
from). **None of these are bundled, redistributed, or downloadable through
this repository or its `setup/` scripts** -- they require your own valid
BioSkryb/BaseJumper license, obtained directly from BioSkryb.

**If you don't have (or can't get) such a license**, skip to
["No license? Use the free/open-source path instead"](#no-license-use-the-freeopen-source-path-instead)
below -- this project also has a complete, independently-built pipeline that
needs none of these four components and produces broadly comparable (not
identical) results using only free/open-source tools.

## The four licensed components

| # | Component | What it is | Why it's licensed |
|---|---|---|---|
| 1 | `basej-public-pipelines/` (`basej-somatic/`, `basej-lineage/`, `basej-dnaqc/`) | BaseJumper's own Nextflow pipeline source: somatic-variant statistical filtering, source-lineage/phylogeny reconstruction, and QC-threshold definitions | BioSkryb's own pipeline implementation (filtering algorithms, statistical models, QC cutoff rationale) -- proprietary source code, not public |
| 2 | `containers/basej-somatic/*.sif` | 4 pre-built Apptainer/Singularity container images that actually execute pipeline #1's Nextflow code | Pre-built images bundling BioSkryb's own proprietary tool builds/configurations -- not something you can rebuild from public Dockerfiles |
| 3 | DeepVariant single-cell checkpoint (`sc_af.model.ckpt` + `population_vcfs_nochr/`) | A DeepVariant model BioSkryb privately trained specifically for their PTA/WGA chemistry, plus its matching population-allele-frequency VCF bundle | Trained model weights are BioSkryb's own IP -- this is NOT Google's public DeepVariant model |
| 4 | `basej_genomes_nochr/` reference bundle | The public GRCh38 reference genome, pre-organized into the exact nf-core-style directory layout component #1's Nextflow code expects | The underlying genome data is public; this specific repackaging is typically provided by BioSkryb alongside component #1 for licensee convenience |

None of the rest of this project's production pipeline is licensed --
alignment (BWA-MEM2), the free DeepVariant/VEP/Manta containers, CNVkit,
Control-FREEC, Beagle phasing, Sequoia (lineage tool itself, as opposed to
BioSkryb's wrapper around it), and everything in `extra_methods/` are all
free/open-source, and already installed automatically by
[`setup/setup.sh`](setup/setup.sh) regardless of whether you have a
BioSkryb license at all.

## How to obtain your own copy of each component

There is no public download for any of the four components above -- they
are provided directly by BioSkryb to organizations with a valid license.
If your organization already has one:

1. Contact BioSkryb Genomics (https://www.bioskryb.com/) or your existing
   BaseJumper platform account representative to confirm your license
   covers the underlying Nextflow pipeline source and containers (as
   opposed to only the hosted BaseJumper platform UI/API) -- the exact
   delivery mechanism BioSkryb uses for this (e.g. a private Git
   repository, a direct file transfer, platform-specific export tooling)
   is between you and BioSkryb and is outside this project's knowledge or
   control.
2. Once you have each component in hand, place it according to the
   directory layout below (or point this project's environment-variable
   overrides at wherever you've placed it instead -- see the next
   section). Nothing in this repo needs to change to use your own copy.

## Directory layout and environment-variable overrides

The production pipeline (`run_all_resolveome.py`, plus
`step2_library_qc_dashboard.py`, `step3_submit_variants_cnv_sv_annot_lsf.py`,
and `step7_ginkgo_cnv.py`, which each independently reference one or more of
these components) defaults to this project's own existing real-cohort
locations for all four components, all under `resolveome_align/` (or, for
component 3, a sibling `pipeline_resources/` directory at the project
root) -- preserved exactly as-is for this project's own runs. Each default
can be overridden independently via an environment variable, so you can
place your own licensed copies anywhere (a separate licensed-software area
on shared storage, for example) without editing any script:

| Component | Environment variable | Default (this project's own layout) |
|---|---|---|
| 1. `basej-public-pipelines/` | `SCDNA_BASEJ_PIPELINES_DIR` | `workflow_DNA/resolveome_align/basej-public-pipelines` |
| 2. `containers/basej-somatic/*.sif` | `SCDNA_BASEJ_CONTAINERS_DIR` | `workflow_DNA/resolveome_align/containers/basej-somatic` |
| 3. DeepVariant checkpoint + population VCFs | `SCDNA_DEEPVARIANT_CHECKPOINT_DIR` | `pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/DeepVariant/Illumina/bioskryb-af-20241102` (project root) |
| 4. `basej_genomes_nochr/` | `SCDNA_BASEJ_GENOMES_DIR` | `workflow_DNA/resolveome_align/basej_genomes_nochr` |

Expected contents of each directory, once populated:

```
<SCDNA_BASEJ_PIPELINES_DIR>/
  basej-somatic/main.nf
  basej-lineage/main.nf
  basej-dnaqc/assets/plot_qc_config.json
  ... (the rest of each pipeline's own source tree)

<SCDNA_BASEJ_CONTAINERS_DIR>/
  utility-lite-0.2.0.sif
  variantannotation_0.4.sif
  r_ado_1.2.sif
  custom_snp_somatic_filter_sequoia_feb2026.sif

<SCDNA_DEEPVARIANT_CHECKPOINT_DIR>/
  model/sc_af.model.ckpt.index
  model/sc_af.model.ckpt.data-00000-of-00001
  population_vcfs_nochr/cohort-chr<1-22,X,Y>.release_missing2ref.no_calls.vcf.gz

<SCDNA_BASEJ_GENOMES_DIR>/
  genomes/Homo_sapiens/NCBI/GRCh38/Annotation/GATK_bundle/Sequence/genome.fa.fai
  ... (the rest of the nf-core-style genome bundle)
```

When an environment variable is **unset**, the corresponding default above
is used -- this project's own production runs have never set any of these
four variables, so unset is the common case and always preserves exactly
today's existing behavior.

## Validating your setup

Before attempting a real production run with your own licensed components,
run the preflight validator -- it checks whether each component is present
at its configured location (honoring the same 4 environment variables
above) and prints specific, actionable guidance for whichever one is
missing. **It never downloads, builds, or bundles anything itself** -- it
only checks and reports:

```bash
bash setup/check_licensed_components.sh

# or, with your own components placed elsewhere:
SCDNA_BASEJ_PIPELINES_DIR=/path/to/your/basej-public-pipelines \
SCDNA_BASEJ_CONTAINERS_DIR=/path/to/your/basej-somatic-containers \
SCDNA_DEEPVARIANT_CHECKPOINT_DIR=/path/to/your/bioskryb-af-checkpoint \
SCDNA_BASEJ_GENOMES_DIR=/path/to/your/basej-genomes-bundle \
  bash setup/check_licensed_components.sh
```

Exit code `0` means every licensed component was found; `1` means at least
one is missing (with guidance printed for each).

## No license? Use the free/open-source path instead

If you don't have, or can't obtain, a BioSkryb/BaseJumper license, use
**`main.wdl`** (run via Cromwell) instead of the production Python/LSF
pipeline -- see [`README.md`](README.md)'s "Getting started (WDL + Cromwell
path)" section and [`wdl/README.md`](wdl/README.md) for exact run commands.
`setup/setup.sh` installs everything that path needs automatically, with
**zero licensed dependencies of any kind**:

- Alignment/variant-calling: BWA-MEM2, DeepVariant's own free bundled model
  (not BioSkryb's checkpoint), stock VEP/Manta containers
- Phasing: Beagle 5.5 (GPLv3) + the public 1000 Genomes GRCh38 panel
- LOH/ploidy: Control-FREEC (GPL)
- Lineage: Sequoia + MPBoot (GPLv3), run standalone via
  `extra_methods/lineage_track.py` -- the same lineage-building tool
  BioSkryb's own `basej-lineage` wraps, just invoked directly instead of
  through their licensed Nextflow pipeline
- A cohort-internal somatic/artifact filter in place of `basej-somatic`'s
  licensed statistical pipeline

This free path is a genuinely different, independently-built
implementation -- not an attempt to clone BioSkryb's proprietary method --
so its results are **broadly comparable, not identical**, to a real
BioSkryb-licensed production run. Every place the two paths diverge
methodologically is documented in
[`README_RFC_DNA_Workflow.md`](README_RFC_DNA_Workflow.md)'s methods
comparison tables.
