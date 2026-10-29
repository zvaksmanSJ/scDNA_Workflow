# Quick start: download, install, and run this project from scratch

This is a condensed, linear checklist for a brand-new user/machine with
nothing installed yet. For full rationale and every option, see
[`README.md`](./README.md) and [`README_RFC_DNA_Workflow.md`](./README_RFC_DNA_Workflow.md).

## Which path do I want?

This repo has **two** ways to run the same pipeline:

| | Free/portable path (`main.wdl`) | Production path (`python_lsf_pipeline/`) |
|---|---|---|
| Needs a BioSkryb license? | **No** -- 100% free/open tools | **Yes** -- 4 BioSkryb/BaseJumper components (see [`LICENSED_COMPONENTS.md`](./LICENSED_COMPONENTS.md)) |
| Runs anywhere? | Yes (any Linux host with conda+Apptainer+Java) | Built for this project's own LSF cluster |
| What this guide covers | **Everything below** | Steps 1-2, then see [`LICENSED_COMPONENTS.md`](./LICENSED_COMPONENTS.md) + `README_RFC_DNA_Workflow.md`'s ["Workflow (Python/LSF production pipeline)"](./README_RFC_DNA_Workflow.md#workflow-pythonlsf-production-pipeline) section instead of steps 3+ |

If you don't have a BioSkryb license (or just want to try the workflow),
follow this guide straight through -- it uses only free/open tools and
reference data.

## Prerequisites (install these yourself first)

| Need | Get it |
|---|---|
| Linux (or WSL2) | -- |
| conda / Miniconda / Miniforge | https://docs.conda.io/en/latest/miniconda.html |
| Apptainer or Singularity | https://apptainer.org/docs/admin/main/installation.html |
| `git`, `curl`, `tar`, `gzip`, `unzip` | usually preinstalled; otherwise `apt install`/`yum install` |
| Java 11+ | https://adoptium.net/ |
| Cromwell (the WDL engine) | download a release jar: https://github.com/broadinstitute/cromwell/releases (built/validated against Cromwell 87) |
| ~100+ GB free disk | reference genome + VEP cache + 1000 Genomes phasing panel, before any run output |

## 1. Clone

```bash
git clone https://github.com/zvaksmanSJ/scDNA_Workflow.git
cd scDNA_Workflow
```

## 2. Run setup.sh (installs everything else, free tools only)

```bash
./setup/setup.sh
```

Idempotent and safe to re-run (every step is skip-if-already-present).
Installs into `./wdl_runtime/` by default: 4 conda envs, a GRCh38 reference
+ indices, stock (non-licensed) DeepVariant/VEP/Manta containers, Beagle 5.5 +
its genetic maps + the 1000 Genomes phasing panel, and Sequoia+MPBoot for
lineage trees. Takes a while (large downloads) -- you can skip sections
you don't need yet with `SKIP_<SECTION>=1` env vars (`SKIP_CONDA`,
`SKIP_REFERENCE`, `SKIP_CONTAINERS`, `SKIP_PHASING`,
`SKIP_FREEC_REFS`, `SKIP_PARITY_TOOLS`, `SKIP_LINEAGE`) -- see the script's
own header comment for details. It does **not** install Cromwell/Java/
conda/Apptainer themselves (see Prerequisites above).

## 3. Build your inputs.json

Make a manifest TSV (tab-separated, no header needed): one row per cell --
`sample_id`, path to R1 FASTQ, path to R2 FASTQ.

```bash
cat > my_cells.tsv << 'EOF'
cell1	/path/to/cell1_R1.fastq.gz	/path/to/cell1_R2.fastq.gz
cell2	/path/to/cell2_R1.fastq.gz	/path/to/cell2_R2.fastq.gz
EOF

python3 wdl/scripts/build_inputs_json.py \
  --manifest my_cells.tsv \
  --base wdl/example_inputs.json \
  --donor-id my_cohort \
  --output inputs.json
```

This merges your FASTQ list into `wdl/example_inputs.json`'s fixed
infrastructure paths (reference, containers, Beagle panel, script paths --
all pointing at what `setup.sh` just installed), **and** regenerates the
full chr1-22 phasing/genetic-map entries from `setup.sh`'s real output (the
bundled example only illustrates those with 2 chromosomes, so this
regeneration step is required, not optional -- see the script's own
`--help` for why).

Then open `inputs.json` and check/edit the handful of things that are
still per-deployment choices, not auto-filled:
- `gender` (`male`/`female` -- the example defaults to `male`)
- `metrics_intervals_bed` (defaults to a placeholder `.fai`; supply your
  own BED of regions of interest for real coverage QC)
- `dv_use_allele_frequency` / DeepVariant checkpoint overrides, only if you
  have your own private checkpoint (otherwise leave as-is: DeepVariant's
  free bundled model is used)

## 4. Run it with Cromwell

```bash
java -Dconfig.file=wdl/cromwell_local.conf -jar cromwell-87.jar run \
  main.wdl -i inputs.json -m outputs.json
```

Run this from the repo root (`main.wdl`'s own `import` statements are
relative to that location). `cromwell_local.conf` runs every task as a
local subprocess -- fine for a smoke test or a couple of cells. For a real
multi-cell cohort on this cluster's LSF scheduler instead, swap in
`wdl/cromwell_lsf.conf` -- see `wdl/README.md`'s "Running at LSF scale
(production)" section (only the `-Dconfig.file=` flag changes).

## 5. Organize outputs and open the dashboard

```bash
python3 wdl/scripts/organize_outputs.py \
  --outputs-json outputs.json \
  --sample-ids cell1 cell2 \
  --dest /path/to/results \
  --inputs-json inputs.json

xdg-open /path/to/results/dashboards/dashboard.html   # `open` on macOS
```

See `README.md`'s "Where are my results?" section for what every folder
under `<dest>/` contains.

## Worked example: a quick smoke test with small subsampled FASTQs

If you just want to confirm the whole chain works end to end before
committing to a full-size run, use a handful of small, subsampled FASTQs
(a few million reads each run through alignment/variant-calling/etc in
minutes instead of hours). Build `my_cells.tsv` by pointing at your own
subsampled files, e.g.:

```bash
cat > my_cells.tsv << 'EOF'
test_Cells_C1	/path/to/test_Cells_C1_R1.fastq.gz	/path/to/test_Cells_C1_R2.fastq.gz
test_Cells_C2	/path/to/test_Cells_C2_R1.fastq.gz	/path/to/test_Cells_C2_R2.fastq.gz
EOF
```

then follow steps 3-5 exactly as above. This is a good first run to try
before scaling up to a real cohort.

## If you have a BioSkryb license instead

Stop after step 2 (or skip it entirely) and follow
[`LICENSED_COMPONENTS.md`](./LICENSED_COMPONENTS.md) to point
`python_lsf_pipeline/` at your own licensed components, then
`README_RFC_DNA_Workflow.md`'s production-pipeline section for how to run
it -- that path produces this project's own reference results and is what
the real cohort in `DNA_results_Expended/` was generated with.
