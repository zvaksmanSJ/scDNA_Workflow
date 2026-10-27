# Required fixes for `basej-lineage` (licensed BioSkryb component)

`run_all_resolveome.py`'s `source_lineage` stage runs the BaseJumper/BioSkryb
`basej-lineage` Nextflow pipeline, expected at
`workflow_DNA/resolveome_align/basej-public-pipelines/basej-lineage/`. This
directory (like its sibling `basej-somatic/`) is a **licensed BioSkryb
component, not covered by `setup/setup.sh` and not tracked in this repo** --
you obtain it separately under your own BioSkryb license and place it at that
path yourself.

As first exercised end-to-end against this project's real 31-cell cohort
(2026-10-27), a freshly-obtained copy of `basej-lineage` has **three** real
bugs when run under a modern Nextflow release (confirmed against 26.04.6;
`basej-somatic`, which is otherwise structured very similarly, does not hit
any of these). None of these are specific to this project's data --
all three are either Nextflow-version compatibility issues or missing
edge-case handling in the vendored pipeline's own code, and will recur for
anyone else who obtains a fresh copy and runs it as-is. Apply all three
patches below before using `basej-lineage` in production.

## 1. Workflow-block-local closures lose scope (compilation failure)

**Symptom:** `nextflow run .../basej-lineage/main.nf` fails to even compile:
```
Error main.nf:327:34: `resolveParam` is not defined
Error subworkflows/lineage_from_matrices.nf:217:26: `findPileupTree` is not defined
Error subworkflows/lineage_from_matrices.nf:217:66: `findPhyloPlacedVariantsAll` is not defined
```

**Root cause:** `main.nf` and `subworkflows/lineage_from_matrices.nf` each
define a small helper as a **local closure variable** inside the `workflow {
}` block (`def resolveParam = { ... }`, `def findPileupTree = { ... }`, `def
findPhyloPlacedVariantsAll = { ... }`). Under this Nextflow version, such a
closure loses scope either (a) after a subworkflow/process invocation
statement further down the same workflow body (`resolveParam`, called again
after `LINEAGE_FROM_MATRICES(...)`), or (b) when referenced from inside a
nested `.map { ... }` operator closure (`findPileupTree` /
`findPhyloPlacedVariantsAll`). `basej-somatic` never uses this
closure-as-local-variable pattern anywhere, which is why it doesn't hit this.

**Fix:** convert all three to plain top-level script functions (`def
resolveParam(String paramName, Object paramValue, Map manifestParams) { ...
}`, `def findPileupTree(outputs, subdir, mutKind) { ... }`, `def
findPhyloPlacedVariantsAll(outputs, subdir) { ... }`), defined outside any
`workflow { }` block, alongside the file's `include { ... }` statements.
Top-level functions are plain script methods with no such scope limitation.
`resolveParam` additionally needs its previously-closed-over `manifest_params`
Map passed in explicitly as a third parameter at every call site, since a
top-level function can't close over a workflow-block-local variable the way
the original closure did.

## 2. `file('/dev/null')` as an input sentinel breaks every container that receives it

**Symptom:** a process that received this sentinel as one of its declared
`path` inputs fails on container startup with:
```
open /dev/null: permission denied
.command.run: line 158: /dev/null: Permission denied
Command 'ps' required by nextflow to collect task metrics cannot be found
```
100% reproducible (not transient) on every attempt, regardless of which
compute node it lands on.

**Root cause:** `genotype_bin_file`, `mandatory_qc_file`, and both
`findPileupTree`/`findPhyloPlacedVariantsAll` helpers (see fix 1) use
`file('/dev/null')` as the literal path value for "this optional input was
not provided." Declaring the real `/dev/null` device node as a Nextflow
`path` input makes Nextflow auto-bind-mount its parent directory (the
**host's own live `/dev`**) into the container so the file is visible there
-- confirmed via `singularity exec -B /dev ... echo test > /dev/null` vs. the
same command without `-B /dev` on the same node/container/host. Bind-mounting
`/dev` wholesale into an unprivileged apptainer/singularity container breaks
`/dev/null`'s own device semantics for the rest of that task.

**Fix:** replace every `file('/dev/null')` sentinel in `main.nf` and
`subworkflows/lineage_from_matrices.nf` with a real, empty, regular file
committed under `assets/EMPTY_SENTINEL` (0 bytes) instead, e.g. `def
EMPTY_SENTINEL = file("${projectDir}/assets/EMPTY_SENTINEL")` (top level, or
local, as appropriate) and use that value everywhere `file('/dev/null')` was
used before. A plain regular file needs no special container device
handling, so Nextflow's auto-bind only ever touches ordinary project paths.

## 3. `MULTIQC_LINEAGE` requires a report even when there is nothing to show

**Symptom:** `source_lineage` fails at its very last step, after the real
phylogeny tree(s) have already been built successfully:
```
ERROR ~ Error executing process > 'MULTIQC_LINEAGE'
Caused by:
  Missing output file(s) `multiqc_report.html` expected by process `MULTIQC_LINEAGE`
```

**Root cause:** `MULTIQC_LINEAGE`'s script writes up to three
`*_mqc.html` custom-content sections (digital genotype heatmaps, VAF
heatmaps, signature bar graphs), each only if its corresponding PNG(s)
exist, then runs `multiqc . --config ... (run_modules: [custom_content])`.
Whenever **all three** are legitimately empty -- e.g. variant placement
correctly found zero variants left to place beyond what the phylogeny step
already used (common/expected, not an error), and/or
`--run_phylogeny_mutsig false` -- `multiqc` has zero inputs, logs "No
analysis results found", and correctly (from its own perspective) writes no
report at all. But the process declares `multiqc_report.html` as a
mandatory output, so Nextflow then fails the whole task over what is
actually a successful run with nothing left to visualize.

**Fix:** after the three conditional `_write_png_section(...)` calls in
`MULTIQC_LINEAGE`'s script, add a fallback: if no `*_mqc.html` file was
written by any of them, write one minimal placeholder `*_mqc.html` section
explaining why (e.g. `lineage_no_results_mqc.html`), so `multiqc` always has
at least one custom-content input and therefore always produces
`multiqc_report.html`. This does not change behavior for any cohort that
does have real heatmap/signature data -- it only prevents a cosmetic
reporting step from failing an otherwise-successful lineage run.

---

All three fixes were verified end-to-end against the real 31-cell
`DNA_results_Expended` cohort on 2026-10-27: `source_lineage` reached `DONE`
with a genuine Sequoia/MPBoot phylogeny tree (882 variants) and a real
`multiqc_report.html`.
