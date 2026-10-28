# ResolveOME Single-Cell DNA Workflow RFC <!-- omit in toc -->

**Status:** Draft
**RFC number:** 0000 (placeholder; see [RFC process note](#rfc-process-note))
**Workflow version:** 3.0.0 (supersedes the v2.0.0 draft in `SingleCell_DNA_Workflow_python/docs/README_RFC_DNA_Illumina_Workflow_v2.0.0.md`)

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
  - [Full-depth alignment and QC](#full-depth-alignment-and-qc)
  - [Cell selection and interpretation of QC](#cell-selection-and-interpretation-of-qc)
  - [Variant and copy-number analyses](#variant-and-copy-number-analyses)
  - [Somatic analysis, lineage, and SV similarity](#somatic-analysis-lineage-and-sv-similarity)
  - [Extra-methods additions: phasing, LOH/ploidy, lineage verification, and the tabbed dashboard](#extra-methods-additions-phasing-lohploidy-lineage-verification-and-the-tabbed-dashboard)
  - [The general-purpose WDL rewrite](#the-general-purpose-wdl-rewrite)
  - [Reference compatibility and limitations](#reference-compatibility-and-limitations)
- [Specification](#specification)
  - [Dependencies](#dependencies)
  - [Reference files](#reference-files)
  - [Workflow (Python/LSF production pipeline)](#workflow-pythonlsf-production-pipeline)
  - [Tool-by-tool commands, in execution order](#tool-by-tool-commands-in-execution-order)
  - [Workflow (general-purpose WDL pipeline)](#workflow-general-purpose-wdl-pipeline)
  - [Output organization](#output-organization)
  - [Execution and failure handling](#execution-and-failure-handling)
- [Appendix](#appendix)
  - [Verification and current implementation status](#verification-and-current-implementation-status)
  - [Known issues and fixes baked into this workflow](#known-issues-and-fixes-baked-into-this-workflow)
  - [RFC process note](#rfc-process-note)

# Introduction

This RFC describes a reproducible, full-depth workflow for single-cell DNA
Illumina data, in two complementary forms:

1. **The production Python/LSF pipeline**
   (`workflow_DNA/python_lsf_pipeline/step1*.py` through `step7*.py`,
   orchestrated by `run_all_resolveome.py`/`.sh` in that same folder), which
   combines the ResolveOME analysis steps with selected St. Jude Cloud/BaseJ
   DNA methods, plus an extensive set of **extra-methods additions** (phasing,
   loss-of-heterozygosity + ploidy, lineage reconstruction + independent
   verification, additional QC metrics, mutational-signature analysis, and a
   consolidated tabbed dashboard) developed to extend the original pipeline's
   analytical depth.
2. **A general-purpose WDL 1.0 port** (`workflow_DNA/main.wdl` +
   `workflow_DNA/wdl/`) of the same pipeline (steps 1-3 plus every
   extra-methods addition), built so the workflow can be run at any site
   with free/open tools only, starting from already-prepared per-cell
   FASTQs, with no dependency on this project's specific HPC cluster or any
   licensed/proprietary resource.

`workflow_DNA/` is organized so its two top-level items are the WDL entry
point (`main.wdl`) and this RFC (`README_RFC_DNA_Workflow.md`); everything
else lives in a clearly-named subfolder -- see
[Output organization](#output-organization) for the full layout map.

This file is a local draft modeled on the
[St. Jude Cloud RFC process](https://github.com/stjudecloud/rfcs) and its
[RNA-Seq v2 pipeline RFC](https://github.com/stjudecloud/rfcs/blob/main/text/0001-rnaseq-workflow-v2.0.0.md).
It is not an accepted RFC or a submission to that repository. The repository
currently limits RFC initiation to St. Jude Cloud team members and asks
community members to contact the team before investing in a proposal.

This document consolidates and supersedes three earlier project documents,
now preserved for historical/operational reference in
`SingleCell_DNA_Workflow_python/docs/`:

- `README_RFC_DNA_Illumina_Workflow_v2.0.0.md` -- the prior RFC draft this
  document is built on and expands.
- `REVISED_WORKFLOW_STEPS.md` -- the exhaustive command-line reference for
  every launcher flag and tool invocation in the Python/LSF pipeline.
- `AGENT_HANDOFF_FULL_STEP3.md` -- a deep-dive handoff on step3 (variant/CNV/
  SV calling + annotation) internals, including troubleshooting notes for
  several bugs fixed during development (see
  [Known issues and fixes](#known-issues-and-fixes-baked-into-this-workflow)).

This RFC summarizes the workflow specification and its rationale; the three
documents above remain the authoritative exhaustive reference for every
launcher flag and historical troubleshooting note.

# Motivation

- **Use full-depth reads for production analysis.** Low-depth subsets are
  useful for QC and development, but are not a substitute for full-depth DNA
  variant calling. Production alignment therefore uses all available trimmed
  reads.
- **Make processing reproducible and auditable.** The workflow records
  inputs, tool checks, manifests, job identifiers, stage state, and logs.
  Failed or missing jobs must not be reported as successful stages.
- **Preserve useful project analyses while improving reference parity.**
  Existing small-variant, structural-variant, and copy-number analyses remain
  available alongside the source-derived somatic and lineage workflows.
- **Keep biological and technical decisions visible.** Cell eligibility,
  cohort-level caller gates, sex-chromosome handling, and unavailable metrics
  are reported explicitly instead of being hidden in scripts or converted to
  success-shaped defaults. The same principle extends to this session's
  extra-methods additions: known methodological caveats (e.g. the
  Control-FREEC candidate-ploidy boundary effect, or using artifact-panel
  singleton sites as a mutational-signature proxy) are documented in the
  dashboard itself, not hidden.
- **Compare complementary copy-number and lineage views.** CNVkit and Ginkgo
  are both retained, and small-variant lineage is kept distinct from
  exploratory structural-variant similarity.
- **Extend analytical depth without re-deriving the core pipeline.** Haplotype
  phasing, LOH/ploidy estimation, lineage tree construction with independent
  statistical verification, and mutational-signature refitting all build on
  the existing per-cell outputs rather than requiring a new primary pipeline.
- **Make the workflow runnable anywhere, by anyone, starting from free tools
  only.** The production pipeline above depends on this project's specific
  HPC environment and, for DeepVariant specifically, a licensed BioSkryb
  single-cell checkpoint and population-VCF bundle. The WDL port removes
  both dependencies: it defaults to DeepVariant's own free, bundled WGS
  model and stock Ensembl/GATK/1000-Genomes reference resources, uses
  Sequoia's standalone GPLv3 script instead of the licensed BioSkryb
  container for lineage reconstruction, and ships a `setup/setup.sh` that installs
  every tool, environment, and reference a new site would need.

# Discussion

## Full-depth alignment and QC

The workflow validates input FASTQ gzip streams, performs raw and trimmed
read QC, and trims paired reads with fastp. Adapter detection and forced
poly-G trimming are intentional for merged inputs whose headers may not
identify the instrument cycle. The production BAM is aligned from the
complete trimmed read pairs; no fixed two-million-read subsampling is used
for variant analysis.

Duplicates are removed with `samtools markdup -r`. Alignment and library
metrics are collected before and after duplicate removal where appropriate.
Undefined metrics, such as an insert-size statistic that cannot be estimated,
are reported as unavailable with a reason; they are not represented as zero.

## Cell selection and interpretation of QC

The existing cell gate is retained: at least 1,000,000 reads, mapping rate at
least 0.70, and proper-pair rate at least 0.60. Mapping and pairing rates are
computed from the trimmed alignment before duplicate removal, preserving the
workflow's established cell-selection semantics. Read trimming may change
which cells pass; membership is not forced to match an earlier cohort.

The library dashboard reports FastQC, fastp, Picard, samtools, preseq, and
Ginkgo metrics. Its reference comparison score is informational and does not
replace the explicit cell gate. A component that cannot be computed is marked
`UNKNOWN` with a reason.

**Note (WDL rewrite only):** the hard-threshold gate above describes the
original Python/LSF pipeline. The WDL rewrite (`main.wdl`/`wdl/` and its
Sprocket port) instead gates on a BaseJumper-derived composite score (5
metrics, 0-5 scale, PASS/Borderline/FAIL) -- a deliberate, requested change
from the simple-threshold approach above, matching BaseJumper's own
`basej-dnaqc` pipeline rather than this project's historical gate. See "The
general-purpose WDL rewrite" below and `wdl/README.md`'s "Per-cell QC gate"
section for the full formula.

## Variant and copy-number analyses

Per-cell small-variant calling uses the existing BioSkryb DeepVariant model,
allele-frequency inputs, HP channel, and chunking configuration. SnpEff and
VEP annotation, Manta structural-variant calling, and CNVkit are retained.
The candidate-preparation step supports the historical restricted analysis
but the revised workflow selects the reference-compatible, genome-wide mode;
coding, per-cell depth, and per-cell VAF filters are not applied before the
source cohort caller gate.

CNVkit and Ginkgo are complementary analyses, not interchangeable outputs.
CNVkit retains both 1-Mb and 200-kb cohort summaries. Ginkgo uses the
matching local GRCh38 binning and GC/boundary resources and preserves its
source segmentation and export behavior. Exon/intron coverage and the
library dashboard provide additional QC rather than changing variant calls.

## Somatic analysis, lineage, and SV similarity

The BaseJ somatic workflow retains its source-defined quality, artifact,
statistical, and shared-variant filters. The cohort-level caller gate
requires at least one passing cell (`N_PASS(MLV==0) >= 1`). Male
sex-chromosome correction is the selected configuration for the current
analysis; chromosomes 1-22, X, and Y are used for somatic processing.

Source lineage receives small-variant-only matrices and genotype tables.
Symbolic SVs are not represented as indels. SNV, indel, and combined lineage
modes use the source VAF presence/absence settings and 1,000 MPBoot
bootstrap replicates. If no small variants remain, lineage is reported as
`NOT_COMPUTABLE`.

The project's existing SV similarity analysis remains exploratory. Its
all-feature result is accompanied by a shared-feature sensitivity analysis
that requires detection in at least two cells. These similarity results are
not ancestry inference.

## Extra-methods additions: phasing, LOH/ploidy, lineage verification, and the tabbed dashboard

Beyond the original v2.0.0 scope, this workflow adds (under
`workflow_DNA/extra_methods/`):

- **Haplotype phasing** (`phase_variants_lsf.py`) via Beagle 5.5 against the
  1000 Genomes GRCh38 high-coverage phased reference panel, **autosomes
  only**. chrX is excluded because the panel's chrX reference VCF mixes
  haploid/diploid genotype representations that Beagle cannot parse in one
  pass for hemizygous male samples -- confirmed to fail identically on every
  cell tested. A cohort that includes female samples needing chrX phasing
  would need to split chrX into PAR1/non-PAR/PAR2 regions and phase them
  separately (not implemented, since this cohort did not need it).
- **Loss-of-heterozygosity and ploidy estimation** (`loh_freec_lsf.py`) via
  Control-FREEC, run **twice** per cell: once at a fixed ploidy of 2 (a
  simple, conservative baseline) and once letting Control-FREEC select the
  best-fitting ploidy from a candidate range (`--ploidy 2,3,4,5,6,7,8`).
  **Known caveat, surfaced in the dashboard rather than hidden:** widening
  the candidate range can push *more* cells to the new ceiling rather than
  converging on a stable answer -- finer copy-number discretization at
  higher ploidy mechanically improves Control-FREEC's own fit score without
  necessarily reflecting real biology. Both results are kept and shown
  separately for this reason.
- **Cohort-wide LOH-by-region recurrence** (`loh_by_region.py`): a heatmap
  and table of how often each genomic region shows LOH across the cohort,
  useful for spotting both real recurrent LOH and likely artifacts (e.g. a
  pericentromeric region flagged as LOH in essentially every cell is more
  likely a mapping/reference artifact than real shared biology).
- **Per-cell CNV burden + genome-wide scatter plots** (`cnv_overview.py`):
  summarizes CNVkit's own segment-level integer copy-number calls
  (`.dedup.call.cns`, already produced by the existing `variants_cnv_sv`
  stage) into a per-cell burden table (percent of genome gained/lost/
  altered, segment count, min/max copy number), and re-presents CNVkit's
  own genome-wide scatter plot per cell in a grid for quick visual
  comparison across the cohort. No new CNV caller is introduced; this is
  purely a summary/visualization layer over existing CNVkit output.
- **Ploidy-determination narrative** (`ploidy_determination.py`): parses
  Control-FREEC's own stdout (from the candidate-ploidy run) to report,
  per cell, the RSS (residual sum of squares) fit score and percent of
  genome explained for every candidate ploidy it tried (2 through 8), plus
  its own plain-text reasoning for the final choice -- including the cases
  where it deliberately overrides its own best-RSS candidate in favor of a
  lower, nearly-as-good-fitting one (FREEC's own documented tie-breaking
  behavior). This directly answers "how was this cell's non-diploid ploidy
  actually determined", not just "what number did it report."
- **Lineage reconstruction with independent verification**
  (`lineage_track.py`, `verify_lineage.py`): builds a Sequoia lineage tree
  from the cohort's annotated VCFs, then independently re-runs the same
  pipeline with a different random seed and resampled sites, comparing the
  two tree topologies via `ape::dist.topo` (Robinson-Foulds distance) to
  report real per-branch bootstrap support -- a genuine check that the
  primary tree is not an artifact of one particular site-sampling draw, not
  just a restatement of Sequoia's own internal bootstrap. Note:
  `ape::dist.topo`'s branch-score distance can legitimately come back as NA
  for rooted trees under the PH85 method (confirmed via `dist.topo`'s own
  "Some trees were rooted" warning) -- the dashboard reports this
  explicitly rather than crashing or silently showing a wrong number (a
  real bug found and fixed this session: `build_extra_dashboard.py` used to
  crash formatting a `None` branch-score distance as a float).
- **SV gene-fusion annotation** (`sv_gene_fusion_annotation.py`): intersects
  the SAME Manta structural-variant calls already produced by the
  `variants_cnv_sv` stage (`sv_features.tsv`/`sv_binary.tsv`) against a gene
  model -- no new SV caller, no re-calling. Each BND (translocation-type)
  event's two breakpoints are checked, gene-pair by gene-pair, against a
  WHITELIST of known cancer-associated gene fusions (default: the Mitelman
  Database of Chromosome Aberrations and Gene Fusions in Cancer, via
  Arriba's MIT-licensed redistribution) and a BLACKLIST of gene pairs
  recurrently reported as fusion transcripts in normal, non-cancer tissue
  (default: Babiceanu et al. 2016 + Greger et al. 2014) -- both pluggable
  via config for a custom additional list without replacing the defaults
  (see `extra_methods/refs/sv_gene_fusion/README.md` for full citations/
  licensing). DEL/DUP events are separately annotated with every gene their
  span overlaps; no whitelist/blacklist check applies to non-BND SV types
  (both reference lists are RNA-fusion-transcript evidence, which only a
  translocation/breakend can plausibly produce -- documented in the
  script's own docstring). Produces a full annotated table (every SV), a
  filtered "main" table (gene-overlapping, recurrent SVs), a ranked top-50
  table with a plain-English explanation per row, and an OncoScan-style
  figure (chromosome schematic + per-cell presence strip) for the
  highest-ranked subset of those.
- **Additional QC metrics**: MAPD (Median Absolute Pairwise Difference, a
  per-cell CNV noise/uniformity metric from CNVkit's own bins), an FDR
  estimate for the small-variant calls, and a recurrent-artifact panel that
  flags variant sites recurring suspiciously often across the cohort
  (likely technical artifacts) -- whose per-site "singleton" output doubles
  as a private-mutation proxy set for mutational-signature analysis.
- **Mutational-signature refitting** (`mutational_signatures.R`) via
  MutationalPatterns (R/Bioconductor, free) against COSMIC signatures, using
  the artifact panel's singleton sites as a proxy for private mutations.
  This is an explicit methodological approximation (private mutations are
  not independently called here), documented as a caveat alongside the
  results rather than presented as a definitive private-mutation call set.
- **A consolidated tabbed dashboard** (`build_extra_dashboard.py`) with tabs
  for Overview, QC (embedding the real MultiQC report), Library QC
  (embedding `step2_library_qc_dashboard.py`'s own standalone informational
  comparison against the licensed source pipeline's five pass/fail
  thresholds, with the rationale for each threshold and, per cell, why any
  FAIL/UNKNOWN call was made), Phasing, CNV, LOH, LOH-by-Region, Lineage, SV
  Gene Fusions, and Mutational Signatures -- the single entry point for
  reviewing every result above alongside the original pipeline's own
  outputs.
- **A browsable final-results folder** (`organize_final_results.py`),
  non-destructive and symlink-based, that lays out every component's output
  under `DNA_results/final/` for easy review without duplicating data.

## The general-purpose WDL rewrite

`workflow_DNA/wdl/` is a complete WDL 1.0 port of every stage above (steps
1-3 plus every extra-methods addition), designed to run at any site, not
just this one. Key design decisions:

- **Starts from already-prepared FASTQs.** No multi-lane FASTQ-merge step is
  included; that is a site-specific concern for whoever generates the input
  FASTQs, not part of this general-purpose workflow.
- **Mixed tool-invocation model, by design.** Only DeepVariant, VEP, and
  Manta are individually containerized (via `apptainer exec` calls made
  *inside* their own task command blocks); every other task assumes the
  named conda environments created by `setup/setup.sh` (`resolveome-align`,
  `freec`, `manta_env`, `mutsig`) are already on the execution host --
  mirroring exactly how the real LSF/Python pipeline already runs on this
  HPC cluster. **This has a real execution-engine consequence:** miniwdl
  cannot run this workflow, because every miniwdl task is fully
  containerized with no "run on the host directly" mode. Cromwell's "Local"
  backend, which runs each task's rendered command script directly via
  `/bin/bash` on the host, is required instead -- see
  `wdl/cromwell_local.conf` and `wdl/README.md` for the exact configuration
  (including two non-obvious Cromwell config gotchas found and fixed during
  validation this session).
- **DeepVariant defaults to the free, general-purpose path.** The real
  cohort run in this project uses a BioSkryb-private single-cell checkpoint
  (`sc_af.model.ckpt`) and a matching population-VCF bundle, neither of
  which are free/redistributable. The WDL task supports both: if
  `deepvariant_model_dir_tar` is left unset (the default), it runs the
  simple, standard `run_deepvariant --model_type=WGS` command using
  DeepVariant's own bundled free model, with zero extra downloads; only if
  that input is explicitly supplied does it fall back to the lower-level
  make_examples/call_variants/postprocess_variants pipeline against the
  custom checkpoint, exactly matching the real cohort's own invocation.
- **Sequoia runs standalone, not in a licensed container.** Sequoia itself
  (github.com/TimCoorens/Sequoia) is GPLv3 and fully free; only the
  BaseJumper/BioSkryb `basej-somatic` apptainer container this project's
  real run uses to invoke it is licensed. `lineage_track.py` and
  `verify_lineage.py` were extended with a `--sequoia-script`/
  `--mpboot-path` direct-invocation mode, alongside the existing
  `--sequoia-sif` apptainer mode (mutually exclusive, the real cohort run
  keeps using its proven container path unchanged), so the WDL workflow can
  use the free standalone tool.
- **`setup/setup.sh` is free-tools-only and idempotent.** It installs: the four
  conda environments above (pinned to the exact versions verified against
  this project's own real installs); a GRCh38 FASTA from Ensembl plus
  bwa-mem2/samtools/Picard indices; GATK's public legacy known-sites bundle;
  stock Docker Hub containers for DeepVariant/VEP/Manta (not BioSkryb's
  private images); Beagle 5.5 and its public PLINK genetic maps; the 1000
  Genomes GRCh38 high-coverage phased panel; and Sequoia + MPBoot from
  GitHub. Every step is skip-if-already-present and safe to re-run. It does
  not touch `resolveome_align/refs` or `resolveome_align/env` (this
  project's own existing installation), by design -- see
  [Reference compatibility and limitations](#reference-compatibility-and-limitations).
- **Validated, not merely syntax-checked.** `miniwdl check` passes with zero
  real errors (only cosmetic WDL-1.1-style deprecation notices and 5
  intentionally-unused bwa-mem2 index-staging inputs, both documented
  in-line). Beyond that static check, a real local execution smoke test was
  run: the `AlignTrimQC` task was executed through Cromwell + its Local
  backend against a tiny synthetic reference and synthetic read pairs, using
  the real tools (bwa-mem2, fastp, fastqc, samtools, Picard, preseq) --
  producing a genuinely valid BAM (100% mapped, 99.95% properly paired) and
  real Picard WGS metrics. See `wdl/README.md`'s "Validation performed this
  session" section for the full account, including what was *not* run (a
  full production-scale cohort through the WDL path, which this session's
  handoff constraints do not authorize and which would duplicate the
  multi-day real-cohort run already completed through the Python/LSF path).
- **No hardcoded paths.** Every file the workflow needs (reference, model,
  containers, scripts) is a WDL `File`/`File?` input supplied via a single
  `inputs.json`; string defaults are limited to generic values (conda env
  names, a model filename, a ploidy range, `"male"`) that carry no
  site-specific assumption.
- **Per-cell QC gate: BaseJumper-derived composite score, not the original
  pipeline's hard thresholds.** Right after alignment, a new `CellQcGate`
  task scores every cell on 5 metrics (library complexity, chimeric-read
  rate, mitochondrial-read fraction, CNV bin noise, CNV directional bias),
  each with BaseJumper's own published cutoff, and a FAILing cell (<=2/5)
  has every downstream stage (DeepVariant/VEP/SnpEff, Manta, phasing,
  LOH/ploidy, lineage) skipped entirely for that cell -- not just flagged.
  See `wdl/README.md`'s "Per-cell QC gate" section and
  `extra_methods/cell_qc_gate.py`'s own docstring for the full formula,
  citation, and the one adaptation needed (CNVkit's bins instead of
  Ginkgo's, for the CNV-based metrics -- this project uses CNVkit, not
  Ginkgo, for its real CNV/SV analysis).

### A parallel Sprocket port

`sprocket/` contains a second, Sprocket-native port of this same workflow
(https://github.com/stjude-rust-labs/sprocket, St. Jude's own next-generation
WDL engine, used in place of Cromwell at some St. Jude sites). It required
real changes beyond a new engine config, because Sprocket's execution model
has no "run bare on the host with a conda env activated" option -- every
task runs inside exactly one container. See `sprocket/README.md` for the
full explanation of what changed and why (container-per-task task
restructuring, four new Dockerfiles replacing the three conda environments
above, a `sprocket.toml` backend config). **This port was written and
statically validated (`miniwdl check` passes clean) but not built, pushed,
or run** -- it was explicitly requested as a "write the code, don't actually
do it" exercise, distinct in scope from the WDL rewrite described in the
rest of this section.

## Reference compatibility and limitations

The workflow is configured for GRCh38 and uses a coordinated set of
reference assets: FASTA and indexes, annotation, DeepVariant model and
population resources, CNVkit reference, Ginkgo bins/GC/boundaries, and
annotation caches. Changing only the FASTA path does not change or validate
the rest of that bundle. A different genome build requires matching assets
and configuration.

Some runtime dependencies differ from the original source environment. In
particular, the Ginkgo runtime uses R 4.4.0 rather than the source image's
R 3.6.3, and CNVkit is pinned to 0.9.12 to avoid an incompatibility with
newer pandas. Numerical and plot-rendering identity with the source
environment is not claimed.

For the WDL path specifically: `setup/setup.sh` installs its own, separate set of
conda environments, reference downloads, and containers into a user-chosen
`BASE_DIR`, and never modifies this project's existing
`resolveome_align/refs` or `resolveome_align/env` directories -- those
remain exactly as they are for this project's own real-cohort run. A site
adopting the general WDL workflow should expect the normal first-run
debugging of any new pipeline deployment (resource sizing, LSF-backend
wiring if running at production scale, container availability), since only
a local plumbing/mechanics smoke test -- not a full-scale production run --
has been performed against the WDL path to date.

# Specification

## Dependencies

**Python/LSF production pipeline.** Runs in the existing `resolveome-align`
environment. LSF is used for per-cell and cohort compute jobs; Apptainer
containers and the local Nextflow installation are used for source
pipelines. The workflow performs a preflight of executables, containers,
reference assets, model/cache files, and required R packages before
production submission.

Explicitly versioned dependencies include fastp 0.20.1, CNVkit 0.9.12, and
the reference-compatible preseq 2.0.3 build (exposed as `preseq-bioskryb`).
Other executable and container versions are verified from the installed
environment and recorded by preflight.

**Extra-methods additions** use two further conda environments: `freec`
(Control-FREEC 11.6, samtools 1.24, bedtools 2.31) and `mutsig`
(MutationalPatterns 3.20.1, BSgenome.Hsapiens.UCSC.hg38 1.4.5,
VariantAnnotation, R 4.5). Beagle 5.5 (`beagle.27Feb25.75f.jar`) and its
PLINK genetic maps are installed under `extra_methods/tools/beagle/`.

**General-purpose WDL pipeline.** `setup/setup.sh` creates four conda
environments with these exact pinned versions (verified against this
project's own installs):

| Environment | Key packages |
|---|---|
| `resolveome-align` | python 3.11, bwa 0.7.18, bwa-mem2 2.2.1, samtools 1.20, bcftools 1.20, htslib 1.20, picard 3.2.0, fastqc 0.12.1, multiqc 1.21, cnvkit 0.9.12, snpeff 5.1, gatk4 4.6.2.0, bedtools |
| `freec` | control-freec 11.6, samtools 1.24, bedtools 2.31, bioconductor-rsamtools |
| `manta_env` | python 2.7, manta 1.6.0 |
| `mutsig` | r-base 4.5, bioconductor-mutationalpatterns, bioconductor-bsgenome.hsapiens.ucsc.hg38, bioconductor-variantannotation |

It also pulls stock (not BioSkryb-private) containers --
`google/deepvariant:1.6.1` and `ensemblorg/ensembl-vep:release_111.0` --
plus an attempted `halllab/manta:latest` pull (CNVkit/Manta can also run
from the `manta_env` conda environment directly if the container pull is
unavailable at a given site). WDL execution itself requires either miniwdl
(for static `miniwdl check` validation only) or Cromwell 87+ with its Local
backend configured per `wdl/cromwell_local.conf` (for actual execution).

## Reference files

The specification assumes a GRCh38 analysis and matching local resources:

- A primary GRCh38 FASTA with matching FAI and sequence dictionary.
- GENCODE v46 annotation for coverage and gene-level reporting.
- The BioSkryb DeepVariant `sc_af.model.ckpt` checkpoint and
  chromosome-specific population VCFs (production pipeline only; the WDL
  pipeline defaults to DeepVariant's own free bundled model instead, see
  above).
- The matching CNVkit reference, Ginkgo variable-bin/GC/boundary resources,
  VEP caches, and SnpEff database.
- BaseJ genome resources and source pipeline checkouts used by somatic and
  lineage stages (production pipeline only).
- For phasing and LOH/ploidy (both pipelines): the 1000 Genomes GRCh38
  high-coverage phased panel (`CCDG_14151_B01_GRM_WGS_2020-08-05_chr*`,
  EBI FTP, public), Beagle's own PLINK-format genetic maps, and a
  23-contig-only (chr1-22+X) `.fai`/per-chromosome FASTA split for
  Control-FREEC (it requires a FASTA for every contig in its `chrLenFile`
  and does not tolerate the ~170 decoy/alt/unplaced scaffolds the full
  genome `.fai` lists).
- For lineage (WDL pipeline): Sequoia's `build_phylogeny.R` (cloned from
  GitHub) and the MPBoot binary (downloaded from its GitHub release).

These resources are a bundle. A different assembly or replacement
annotation must be validated as a compatible set before use.

## Workflow (Python/LSF production pipeline)

**This pipeline uses 4 BioSkryb/BaseJumper-licensed components** (the
DeepVariant single-cell checkpoint, the basej-somatic/basej-lineage/
basej-dnaqc Nextflow pipeline source, their container images, and a
genome-bundle layout they expect) -- see
**[`LICENSED_COMPONENTS.md`](./LICENSED_COMPONENTS.md)** for exactly what
each is, why it's licensed, how to point this pipeline at your own
separately-obtained copies via environment-variable overrides (no code
changes needed), and a preflight script
(`setup/check_licensed_components.sh`) that validates your setup before a
real run. If you don't have such a license, see that same document's final
section for the free/open-source `main.wdl` alternative instead.

The following stages define the intended processing order. Each stage
consumes the validated outputs of its predecessors and records its inputs,
outputs, commands, and status.

1. **Preflight and run planning.** Validate the input cell manifest, FASTQ
   readability, executables, containers, references, and caches.
2. **Full-depth preprocessing, alignment, and library metrics.** fastp ->
   FastQC -> BWA-MEM2 -> samtools fixmate/sort/markdup -> Picard metrics ->
   preseq.
3. **Cell calling.** Apply the three explicit cell thresholds (reads,
   mapping rate, proper-pair rate).
4. **Five mutually-independent per-cell tracks, submitted together and
   polled concurrently (2026-10-21 rewrite; see
   ["Known issues and fixes"](#known-issues-and-fixes-baked-into-this-workflow)).**
   A failure anywhere in one track is recorded and never blocks or crashes
   the others:
   - *Small variants:* DeepVariant (8-channel, HP + allele-frequency) ->
     VEP -> SnpEff.
   - *CNV/SV:* CNVkit + Manta.
   - *Ginkgo* (orthogonal CNV cross-check).
   - *Coverage QC* (breadth/depth + whole-genome coverage map).
   - *Control-FREEC* (LOH + ploidy, fixed=2 and candidate=2,3,4,5,6,7,8
     runs).
   Two further per-cell steps fire as soon as their single real prerequisite
   resolves, independent of the other four tracks: Beagle phasing (needs
   only the small-variants track's raw VCF) and the composite,
   BaseJumper-derived cell QC gate (needs only the CNV/SV track's CNVkit
   `.cnr`/`.call.cns` -- informational only, does not exclude any cell).
5. **Ginkgo, coverage, and library dashboard.** Ginkgo/coverage are part of
   step 4 above; the library QC dashboard runs immediately after step 3
   (it only needs step 2's QC files, so it does not wait on step 4 at all).
6. **Candidate preparation and somatic filtering.** BaseJ somatic workflow.
7. **Lineage and SV comparison.** Source lineage (SNV/indel/combined) plus
   SV similarity.
8. **Integrated reporting.**
9. **Extra methods (now wired directly into `run_all_resolveome.py`, not a
   separate manual step).** Phasing, LOH/ploidy (fixed + candidate), LOH-by-
   region, per-cell CNV burden + scatter plots, per-cell ploidy-determination
   narrative (parsed from Control-FREEC's own candidate-ploidy reasoning),
   MAPD/FDR/artifact-panel QC, mutational signatures, and the extra-methods
   tabbed dashboard -- see
   [the extra-methods discussion](#extra-methods-additions-phasing-lohploidy-lineage-verification-and-the-tabbed-dashboard)
   above. `lineage_track.py`/`verify_lineage.py` are deliberately NOT part of
   this automated orchestration: `source_lineage` (step 7, BaseJ lineage
   Nextflow) already builds and bootstraps a Sequoia/MPBoot tree from the
   same annotated VCFs when it runs, and running both would duplicate
   expensive compute for no new information. For a cohort whose
   `source_lineage` stage has NOT been run (e.g. `DNA_results`/Batch 1, where
   the somatic-filtering->lineage Nextflow chain was never executed),
   `lineage_track.py`/`verify_lineage.py` were run manually against that
   cohort instead, specifically so the extra-methods dashboard's Lineage tab
   has real data rather than showing "Pending" indefinitely.

## Tool-by-tool commands, in execution order

The following are command templates matching the generated worker code (not
every temporary-file move, lock, or manifest write). Independent per-cell/
chunk jobs run in parallel; the ordering reflects dependencies, not serial
execution.

**1. FASTQ integrity and raw FastQC**
```bash
gzip -t "$R1"; gzip -t "$R2"
fastqc -t 2 -o "$O/qc/fastqc_raw/$SID" "$R1" "$R2"
```

**2. fastp (0.20.1): adapter, poly-G and default quality filtering**
```bash
fastp --thread 8 --in1 "$R1" --in2 "$R2" \
  --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
  --json "$SID.fastp.json" --html "$SID.fastp.html" \
  --detect_adapter_for_pe --trim_poly_g
```

**3. BWA-MEM2: trimmed FASTQ to name-sorted BAM**
```bash
bwa-mem2 mem -M -Y -K 2500000000 -t 6 \
  -R "@RG\tID:$SID\tSM:$SID\tPL:ILLUMINA\tLB:DNA_lib\tPU:$SID.unit1" \
  "$REF" "$TRIM_R1" "$TRIM_R2" \
  | samtools sort -n -@ 2 -m 3G -o "$NAME_BAM" -
```

**4. samtools: fixmate, coordinate sort, dedup, index**
```bash
samtools fixmate -m -@ 2 "$NAME_BAM" "$FIXMATE_BAM"
samtools sort -@ 2 -m 3G -o "$SORTED_BAM" "$FIXMATE_BAM"
samtools markdup -r -s -f "$MARKDUP_STATS" -@ 2 "$SORTED_BAM" "$DEDUP_BAM"
samtools index -@ 2 "$DEDUP_BAM"
samtools quickcheck -v "$DEDUP_BAM"
```

**5. Picard: library and WGS metrics**
```bash
picard -Xmx8g CollectMultipleMetrics R="$REF" I="$DEDUP_BAM" O="$PREFIX" \
  PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectGcBiasMetrics \
  PROGRAM=MeanQualityByCycle PROGRAM=CollectQualityYieldMetrics
picard -Xmx8g CollectWgsMetrics R="$REF" I="$DEDUP_BAM" O="$PREFIX.wgs_metrics"
```

**6. samtools depth + preseq: library complexity QC**
```bash
samtools depth -b "$METRICS_INTERVALS" "$DEDUP_BAM" | awk '...' > "$COV_SUMMARY"
bam2mr -seg_len 100000 -o "$SID.mr" "$DEDUP_BAM"
preseq-bioskryb gc_extrap -o "$SID.curve" "$SID.mr"
```

**7. DeepVariant (per chunk; production pipeline's BioSkryb checkpoint path)**
```bash
apptainer exec --cleanenv --no-home --no-mount /programs -B "$(pwd):$(pwd)" \
  "$DV_SIF" bash -c '
    /opt/deepvariant/bin/make_examples --mode calling --ref REF --reads BAM \
      --regions REGION --examples "$T/me.tfrecord@8.gz" \
      --population_vcfs POP_VCF --use_allele_frequency
    /opt/deepvariant/bin/call_variants --examples "$T/me.tfrecord@8.gz" \
      --outfile "$T/cv.tfrecord.gz" --checkpoint "$MODELBASE/sc_af.model.ckpt"
    /opt/deepvariant/bin/postprocess_variants --ref REF \
      --infile "$T/cv.tfrecord.gz" --outfile CHUNK_VCF
  '
```
*(WDL pipeline's free default instead: a single `run_deepvariant
--model_type WGS --ref REF --reads BAM --regions REGION --output_vcf OUT`
call, no checkpoint/population-VCF needed -- see
[the WDL rewrite discussion](#the-general-purpose-wdl-rewrite).)*

**8. VEP (per chunk)**
```bash
apptainer exec ... "$VEP_SIF" vep --input_file CHUNK_VCF --output_file STDOUT \
  --format vcf --vcf --compress_output bgzip --offline --cache \
  --dir_cache CACHE --species homo_sapiens --assembly GRCh38 \
  --fasta REF --everything --fork 4 --cache_version 111
```

**9. bcftools merge + SnpEff**
```bash
bcftools concat -a -Oz -o "$MERGED_VCF" "${CHUNKS[@]}"
snpEff -Xmx8000m -v -stats "$SID.snpeff.html" -csvStats "$SID.snpeff.genes.csv" \
  GRCh38.105 "$MERGED_VEP_VCF" | bgzip -c > "$ANNOTATED_VCF"
```

**10. CNVkit**
```bash
cnvkit.py batch "$DEDUP_BAM" --method wgs --fasta "$REF" -r "$CNV_REF_CNN" \
  --output-dir cnv_out --diagram --scatter
cnvkit.py call "$CNS" -o "$CALL_CNS"
cnvkit.py genemetrics "$CNR" -s "$CNS" > "$GENEMETRICS_TSV"
```

**11. Manta**
```bash
apptainer exec ... "$MANTA_SIF" bash -c '
  configManta.py --bam BAM --referenceFasta REF --runDir RUNDIR
  python2 RUNDIR/runWorkflow.py -m local -j THREADS
'
```

**12. Beagle (phasing, per autosome; extra-methods addition)**
```bash
java -Xmx44g -jar beagle.27Feb25.75f.jar gt=input_renamed.vcf.gz \
  ref=REF_PANEL_VCF map=GENETIC_MAP chrom=chr$CHROM out=phased_chrom \
  nthreads=4 impute=false
```

**13. Control-FREEC (LOH + ploidy; extra-methods addition, run twice)**
```bash
freec -conf freec_config.txt   # once with ploidy=2, once with ploidy=2,3,4,5,6,7,8
```

**14. Sequoia lineage (extra-methods addition)**
```bash
Rscript rscript_sequoia_build_phylogeny_only.R \
  --mpboot_path "$MPBOOT_DIR/" --mpboot_bootstrap 1000 ...
```

The complete, unabridged command reference -- including Ginkgo, BaseJ
somatic, SV similarity, and every launcher flag -- remains in
`SingleCell_DNA_Workflow_python/docs/REVISED_WORKFLOW_STEPS.md`.

## Workflow (general-purpose WDL pipeline)

`main.wdl` wires together 8 task modules under `wdl/tasks/`, scattered
over every cell (and, for small variants, nested-scattered over genome
chunks within each cell):

1. `align_qc.wdl` (`AlignTrimQC`) -- fastp, FastQC, BWA-MEM2, samtools
   dedup, Picard metrics, preseq. One call per cell, in parallel.
2. `variant_calling.wdl` (`PlanChunks`, `DeepVariantChunk`,
   `VepAnnotateChunk`, `MergeVcfChunks` x2, `SnpEffAnnotate`) -- DeepVariant
   + VEP run per (cell, chunk) pair in parallel; results merged and
   SnpEff-annotated per cell.
3. `cnv_sv.wdl` (`CnvkitCall`, `MantaCall`) -- one call per cell, in
   parallel, reading the same dedup BAM from step 1.
4. `phasing.wdl` (`PhaseChromosome`, `ConcatPhasedChromosomes`) -- Beagle,
   per (cell, autosome) pair in parallel, then concatenated per cell.
5. `loh_ploidy.wdl` (`ControlFreec`, called twice per cell) -- fixed
   ploidy=2 and candidate-ploidy-selected runs, per cell, in parallel.
6. `qc_aggregate.wdl` (`MultiQC`) -- once for the whole cohort, after every
   cell's step 1/2 QC files are available.
7. `extra_qc_dashboard.wdl` (`MapdQc`, `FdrEstimate`, `ArtifactPanel`,
   `PhasingSummary`, `LohAndPloidySummary`, `QcSummary`,
   `MutationalSignatures`, `BuildDashboard`) -- cohort-wide summaries and the
   final dashboard, each run once.
8. `lineage.wdl` (`BuildLineageTree`, `VerifyLineage`) -- Sequoia tree build
   and independent bootstrap verification, each run once for the cohort.

Every input to every task is supplied through exactly one `inputs.json` for
the whole workflow (keys prefixed `ResolveOmeDnaWorkflow.<input_name>`); the
imported `wdl/tasks/*.wdl` files are library modules, not independently run
against their own JSON in normal use. See `main.wdl`'s own header
comment and `wdl/README.md` for the full explanation, example commands, and
the required `wdl/cromwell_local.conf` execution config.

## Output organization

### `workflow_DNA/` top-level folder layout

```
workflow_DNA/
  main.wdl                       <- WDL entry point (only file, besides the
  README_RFC_DNA_Workflow.md        RFC, that lives directly at top level)

  wdl/                            <- WDL task modules + WDL-support scripts
    tasks/*.wdl
    scripts/{make_region_chunks.py, organize_outputs.py}
    cromwell_local.conf
    README.md

  extra_methods/                 <- calculation scripts (Python/R), used by
                                     BOTH the WDL pipeline and the production
                                     pipeline below (phasing, LOH/ploidy,
                                     lineage, QC/signature summaries, the
                                     dashboard builder, plus their own
                                     LSF-launcher scripts)

  python_lsf_pipeline/            <- the original production Python/LSF
                                     pipeline: step1_*.py .. step7_*.py,
                                     run_all_resolveome.py/.sh,
                                     submit_resolveome_workflow_...py,
                                     every test_*.py, and this pipeline's
                                     own logs/

  setup/                         <- every setup/provisioning script, for
                                     BOTH pipelines: setup.sh (WDL's free-
                                     tools installer), setup_resolveome_
                                     alignment.sh, setup_parity_tools.sh,
                                     check_resolveome_alignment_prereqs.sh,
                                     fix_bwa-mem2_ref.sh,
                                     submit_build_ref_index_lsf.sh,
                                     test_dv_checkpoint_compat.sh, and their
                                     supporting config files

  wdl_runtime/                   <- setup/setup.sh's default install
                                     location (refs/envs/containers/tools it
                                     downloads for the WDL pipeline) --
                                     deliberately separate from
                                     resolveome_align/ below

  docs/                           <- workflow_diagram.{png,svg,dot}

  resolveome_align/               <- UNTOUCHED: this project's real
                                     reference genome, conda env, and
                                     container bundle for its own cohort run
  Backup_code/                    <- UNTOUCHED: historical script snapshots
  dv_compat_test/                 <- UNTOUCHED: test_dv_checkpoint_compat.sh's
                                     own scratch/results directory
```

**Production pipeline:** `extra_methods/organize_final_results.py` builds a
non-destructive, symlink-based `DNA_results/final/` tree: per-component
subfolders, plus `dashboards/dashboard.html` as the single entry point.

**WDL pipeline:** `wdl/scripts/organize_outputs.py` consumes either engine's
output JSON (miniwdl's flat format or Cromwell's `{"outputs": {...}}`
nested format) and lays every result out under:

```
<dest>/
  per_cell/<sample_id>/{alignment,qc,qc_gate,variants,cnv,sv,phasing,loh_ploidy}/
  cohort/{qc,qc_gate,phasing,loh_ploidy,mutational_signatures,lineage}/
  dashboards/dashboard.html
  logs/                 <- logs only; safe to delete without losing results
  run_scripts/           <- inputs.json/outputs.json/extra scripts; safe to
                            delete once the run no longer needs to be audited
```

`qc_gate` (both per-cell and cohort) holds the per-cell QC gate's own
outputs -- see "The general-purpose WDL rewrite" above and
`wdl/README.md`'s "Per-cell QC gate" section. `cohort/qc_gate/
qc_gate_summary.tsv` is the fastest way to see every cell's PASS/
Borderline/FAIL status and whether its downstream analysis ran, in one
file; `per_cell/<sample_id>/qc_gate/<sample_id>.qc_gate_explanation.txt` is
the per-cell breakdown of which of the 5 metrics passed/failed and why.

## Execution and failure handling

The shell launcher defaults to non-submitting preflight behavior when no
mode is selected. `--plan` and `--preflight` do not submit production jobs.
Production submission requires an explicit operator decision.

The driver records per-stage state and LSF job IDs and persists it after
every transition. Foundational stages (alignment, cell calling) remain
strictly fail-fast: a failure there is a real error, not a success-shaped
fallback, since literally nothing downstream can proceed without them.
Starting at the per-cell variant/CNV/SV/phasing/LOH stage (2026-10-21
rewrite), the five mutually-independent tracks described above are
submitted concurrently and polled together: each is resolved (`DONE` or, if
any of its jobs exited, `FAILED_PARTIAL`) the moment its *own* job set
reaches a terminal state, and that resolution immediately triggers its own
specific downstream dependents -- a DeepVariant/VEP/SnpEff failure can
never block recognition of the CNVkit/Manta/Ginkgo/coverage/Control-FREEC
tracks' success (or vice versa). Every stage from that point onward
(`combined_cnv`, the cohort summary analyses, somatic/lineage, SV
similarity, both dashboards) runs in a tolerant mode: a failure is recorded
in `workflow_state.json` with its error message and does not crash the
driver: strictly-downstream stages are skipped with an explicit
`NOT_COMPUTABLE` status rather than attempted against missing input.
Completed outputs and prior result roots must not be overwritten or
invalidated to force a rerun.

Review-only commands (run from `workflow_DNA/python_lsf_pipeline/`):
```bash
bash ./run_all_resolveome.sh --help
bash ./run_all_resolveome.sh --plan
bash ./run_all_resolveome.sh --preflight
```

For the WDL pipeline, `setup/setup.sh` is similarly idempotent (every step is
skip-if-already-present), and individual sections can be skipped with
`SKIP_<SECTION>=1` environment variables; see `setup/setup.sh`'s own header
comment for the full list.

# Appendix

## Verification and current implementation status

**Production pipeline:** reported preflight and synthetic validation
includes 52 focused regression tests, full tool/reference preflight,
synthetic paired-read integration through trimming/alignment/dedup/
Picard/preseq, synthetic source lineage, and a synthetic Ginkgo test with
local/LSF agreement. These checks validate workflow behavior on test data;
they do not establish the biological validity or completion of a
production cohort.

**Extra-methods additions:** now wired directly into `run_all_resolveome.py`
(2026-10-21) rather than run as separate manual scripts. Run end to end
against the real, completed 15-cell Batch 1 cohort (`DNA_results/`) to
validate the wiring: MAPD QC, FDR estimate, artifact panel, mutational
signatures, Beagle phasing (15/15 cells), both Control-FREEC runs (fixed
ploidy=2 and candidate ploidy=2-8, 15/15 cells each), phasing summary, LOH
summary, ploidy summary, LOH-by-region, and QC summary all completed
successfully, with real results on disk under `DNA_results/extra_methods/`.
The candidate-ploidy run reproduces the documented ceiling-pushing caveat
in practice (ploidy distribution across the cohort: 8=4 cells, 6=4, 7=3,
2=3, 3=1 -- most cells were pushed toward the upper end of the 2-8 range
rather than converging on a stable value, exactly the known FREEC behavior
described above). The extra-methods dashboard renders the complete result
set, and `organize_final_results.py` was run to refresh the symlinked
final-results tree under `DNA_results/final/`.

**WDL pipeline:** `miniwdl check main.wdl` passes with zero real errors.
A genuine local execution smoke test (not a full production run) was
performed: `wdl/tasks/align_qc.wdl`'s `AlignTrimQC` task, run through Cromwell
87 + `wdl/cromwell_local.conf` against a tiny synthetic 50kb reference and
~2000 synthetic read pairs, executed the real tools end to end and produced
a valid BAM (100% mapped, 99.95% properly paired) and real Picard
`wgs_metrics`. A full multi-cell production run through the WDL path has
not been performed; see
[the WDL rewrite discussion](#the-general-purpose-wdl-rewrite) for why.

A run is complete only when its stage manifests, job outcomes, required
output checks, and downstream reports have been audited. This draft does
not claim that every stage has completed for a real cohort beyond what is
explicitly stated above.

## Known issues and fixes baked into this workflow

These were found and fixed during development; both the Python/LSF and WDL
code already incorporate the fixes (noted here so they are not
rediscovered/reverted):

- **Single fail-fast "variants" stage crashed the whole driver on one
  flaky chunk (2026-10-20/21).** The original orchestrator submitted
  DeepVariant+VEP+SnpEff+CNVkit+Manta as one combined stage and waited on
  it as a single blocking unit; one transient chunk failure anywhere in
  that ~10,000-job population raised immediately, aborting the driver
  before even fully independent stages (Ginkgo, coverage) got a chance to
  run. Fix: split into `variants_small` (DeepVariant/VEP/SnpEff) and
  `variants_cnv_sv` (CNVkit/Manta), submitted concurrently with Ginkgo,
  coverage, and Control-FREEC via a shared tolerant poller
  (`poll_multi_stage`) that resolves and reacts to each stage independently
  -- see ["Execution and failure handling"](#execution-and-failure-handling).
- **DeepVariant memory headroom for the HP/haplotype channel.** No chunk
  has actually hit `TERM_MEMLIMIT` at the prior 24 GB floor across either
  real cohort checked (2,962 + 10,381 DeepVariant chunk logs), but
  `--dv-add-hp-channel` is now always on (required by the
  `sc_af.model.ckpt` checkpoint's 8-channel shape -- see the channel-count
  fix below), and an extra input channel means a larger per-example tensor
  for `call_variants`' CNN to hold than the 7-channel (no-HP) configuration
  briefly used needed. Fix: raised the floor to 32 GB proactively
  (`max(32000, --mem-mb // 2)`), rather than waiting for a real OOM kill.
- **DeepVariant 8-channel/HP-channel misdiagnosis.** `--dv-add-hp-channel`
  had been flipped to `default=False` based on a misdiagnosis blaming MT
  chunk failures; the real cause was MT's missing population-VCF entry
  (already handled separately). With HP disabled, every non-MT chunk
  produced only 7 of the 8 channels the `sc_af.model.ckpt` checkpoint was
  trained on, failing 100% of the time with a tensor-shape mismatch. Fixed
  by restoring `default=True`.
- **SnpEff OOM under LSF.** LSF reserved memory but Java was not explicitly
  given enough heap, so SnpEff crashed loading the GRCh38.105 database with
  default heap. Fix: explicit `-Xmx`, with the LSF memory request set to at
  least that heap plus 8000 MB of overhead.
- **CNVkit `--fasta`/`-r` mutual exclusion.** CNVkit's `batch` subcommand
  rejects `--fasta` together with `-r/--reference`: `-r` implies reusing an
  existing flat reference, so `--fasta` (only needed to build a new
  reference on the fly) must be omitted whenever a reference `.cnn` is
  supplied.
- **Apptainer `--no-mount /programs` hang.** Prevents a hang when the host
  has a `/programs` mount that collides with the container's own
  filesystem (affects DeepVariant, VEP, and Manta container invocations).
- **DeepVariant MT exclusion.** No MT entry exists in the allele-frequency
  population VCF bundle, causing non-deterministic tensor shapes in
  `make_examples`/`call_variants`; MT is excluded from the "all
  chromosomes" sweep for the BioSkryb-checkpoint path.
- **Serial VCF-chunk merge bottleneck.** The merge step now normalizes and
  indexes each chunk in parallel before `bcftools concat`, rather than
  serially.
- **Control-FREEC `chrLenFile` must be filtered.** The full-genome `.fai`
  lists ~170 decoy/alt/unplaced scaffolds plus Y/MT that are not split into
  per-chromosome FASTAs; Control-FREEC requires a FASTA for every contig
  listed in `chrLenFile`, so a filtered 23-contig (chr1-22+X) `.fai` must be
  used.
- **Beagle chrX phasing failure.** The 1000G panel's chrX reference VCF
  mixes haploid/diploid genotype representations that Beagle cannot parse
  in one pass for hemizygous male samples; chrX is excluded from phasing.
- **Alignment/cell-call jobs scattered across hosts.** `step1`/`step2`'s
  `bsub -R` requests were missing `span[hosts=1]` (step3 already had it),
  so LSF could satisfy an `-n 8` request by scattering single cores across
  up to 8 different hosts once the cluster saturated. Since fastp/bwa-mem2/
  samtools are single-host multithreaded, only the one host actually
  running the process did real work -- measured 6-13% CPU efficiency on
  scattered jobs vs. near-full utilization once fixed.
- **MultiQC dependency string exceeded Linux's per-argument limit.** At
  this project's scale (~10,000+ jobs in the combined small-variant/CNV/SV
  population), an LSF `-w` dependency expression built from one job ID per
  prerequisite reached ~194 KB, over the 128 KB `MAX_ARG_STRLEN` limit on a
  single argument (a different, smaller limit than the overall `ARG_MAX`),
  crashing the orchestrator's own subprocess call. Fix: MultiQC now runs as
  its own standalone, dependency-free stage once `wait_jobs()`/
  `poll_multi_stage()` has already confirmed the prerequisite stage's jobs
  are terminal -- no LSF `-w` expression is ever built for it.
- **NFS write-corruption pattern.** Some tools intermittently corrupted
  output when writing directly to the shared NFS mount; the fix is a
  local-then-copy pattern (write to local scratch, then copy the completed
  file to its final NFS destination).
- **`verify_lineage.py` VCF suffix mismatch.** An early version looked for
  the wrong final-VCF filename suffix; fixed to match
  `*.deepvariant.raw.vep.snpeff.vcf.gz` exactly.
- **`build_extra_dashboard.py` crash on a null branch-score distance.**
  `ape::dist.topo()`'s branch-score distance can legitimately return NA for
  rooted trees under the PH85 method (confirmed via `dist.topo`'s own
  "Some trees were rooted" warning, not a scripting bug) -- the dashboard
  builder formatted this value as a float unconditionally and crashed the
  first time a real run produced NA. Fixed to detect a non-numeric value
  and report "not computable for this tree pair (NA from ape::dist.topo)"
  instead of crashing.
- **Candidate-ploidy Control-FREEC jobs submitted with the wrong ploidy
  range.** A first attempt at the candidate-ploidy run used `--ploidy
  2,3,4` instead of the documented `2,3,4,5,6,7,8` range; caught before
  any downstream result depended on it, the affected jobs were killed and
  resubmitted with the correct range.
- **Incorrect "no public source" claim for `bam2mr`/`preseq-bioskryb`.**
  Despite the pipeline-internal "-bioskryb" naming (used only to
  distinguish this exact build from the newer, CLI-incompatible bioconda
  `preseq` 3.x, which renamed `bam2mr` to `to-mr`), these are preseq
  v2.0.3's own public, GPLv3-licensed `bam2mr`/`gc_extrap` commands
  (github.com/smithlabcode/preseq, also on bioconda) -- not a private
  BioSkryb binary. `setup/setup_parity_tools.sh` already built these from
  the official public release tarball, but it was never wired into
  `setup/setup.sh`'s own flow, and three docs (this file's predecessor
  claims in `README.md`, `sprocket/README.md`, and
  `sprocket/containers/resolveome-align.Dockerfile`) incorrectly stated no
  public source existed. Fixed: `setup.sh` now has its own
  `SKIP_PARITY_TOOLS`-gated section that updates the `resolveome-align`
  conda env with `parity_tools.environment.yml`'s packages and calls
  `setup_parity_tools.sh` (parameterized via a new `TARGET_DIR` override so
  a fresh clone installs under the portable `wdl_runtime/` tree instead of
  this project's own `resolveome_align/`), the Sprocket port's Dockerfile
  now builds the same public source directly instead of requiring a manual
  `COPY`, and all three docs are corrected.
- **WDL `version 1.1` incompatibility with Cromwell.** Cromwell 87 only
  recognizes `version 1.0` (and `draft-2`); a `version 1.1` declaration
  fails to parse at all. Every file under `wdl/` declares `version 1.0` for
  this reason (the workflow never actually used any 1.1-exclusive syntax --
  only the older, still-fully-valid 1.0-style placeholder-option syntax).
- **Cromwell Local-backend config gotchas.** (1) A custom
  `runtime-attributes` override must still declare `String? docker` /
  `String? docker_user`, because Cromwell's built-in Local/SFS backend
  config always includes a `submit-docker` template referencing `${docker}`
  even when unused. (2) Do not redeclare `cpu`/`memory` in that same
  override -- Cromwell's shared-filesystem backend provides built-in
  parsing/validation for those two attributes, and overriding them with a
  plain `String?` breaks `runtime { memory: "~{mem_gb} GB" }` validation.
  Both are documented and fixed in `wdl/cromwell_local.conf`.

## RFC process note

The St. Jude Cloud RFC repository describes a draft, stakeholder-discussion,
and broader-community-discussion process. Before submission, this local
document should be reviewed with the repository maintainers. If accepted
for discussion, the RFC repository's maintainers assign the RFC number; the
placeholder `0000` should then be replaced in the filename, title, and any
resource paths as appropriate.
