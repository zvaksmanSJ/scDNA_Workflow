# SV gene-fusion whitelist / blacklist reference data

These two files are the **default** gene-pair reference lists used by
`extra_methods/sv_gene_fusion_annotation.py` to flag BND (translocation-type)
structural variants whose two breakpoint genes form a known cancer-associated
fusion (whitelist) or a known recurrent fusion seen in normal, non-cancer
tissue (blacklist, usually a transcriptional read-through/splicing artifact,
not a somatic driver event). See that script's own docstring for exactly how
they're applied, and `run_all_resolveome.py`'s `sv_gene_annotation` stage for
how custom additional lists can be supplied (config JSON `custom_whitelists`/
`custom_blacklists`) without replacing these defaults.

## `whitelist_mitelman_cancer_fusions.tsv` -- 542 gene pairs

**Source:** Mitelman Database of Chromosome Aberrations and Gene Fusions in
Cancer (Mitelman F, Johansson B, Mertens F, Eds., 2020;
https://mitelmandatabase.isb-cgc.org) -- the standard, long-maintained
reference for cytogenetically/molecularly confirmed cancer gene fusions.

**Obtained via:** the Arriba RNA-fusion caller's own redistributed copy
(`database/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz`,
https://github.com/suhrig/arriba, release v2.5.0). Arriba's `LICENSE` file
explicitly states its MIT/Expat license covers "arriba (code & executable),
documentation, **and database files**" -- so this redistribution is
unambiguously free to use and re-redistribute with attribution, unlike
querying the Mitelman Database's own web interface directly (which has no
documented bulk-download/reuse terms). Extracted and reduced to unique
`gene1 / gene2` pairs only (Arriba's own copy also carries genomic
coordinates per breakpoint pair, used by Arriba for RNA-seq read evidence
matching -- not needed here, since we're annotating DNA-level Manta SV
breakpoints against gene *identity*, not exact transcript junctions).

## `blacklist_normal_tissue_fusions.tsv` -- 325 gene pairs

**Sources** (two independent, peer-reviewed RNA-seq surveys of fusion
transcripts recurrently detected in **normal, non-cancer** tissue --
i.e. "these show up even in healthy samples, treat recurrence of these
specific pairs as a reason for *lower*, not higher, suspicion of a somatic
driver event"):

- Babiceanu M, Qin F, Xie Z, Jia Y, Lopez K, Janus N, Facemire L, Kumar S,
  Pang Y, Qi Y, Lazar IM, Li H. "Recurrent chimeric fusion RNAs in
  non-cancer tissues and cells." *Nucleic Acids Research* 2016;
  44(6):2859-72. PMID: 26837576. (239 of the 325 pairs)
- Greger L, Su J, Rung J, Ferreira PG, GEUVADIS Consortium, Lappalainen T,
  Dermitzakis ET, Brazma A. "Tandem RNA chimeras contribute to transcriptome
  diversity in human population and are associated with intronic genetic
  variants." *PLoS ONE* 2014; 9(8):e104567. PMID: 25133550. (86 of the 325
  pairs)

**Obtained via:** the Trinity CTAT project's `FusionAnnotator/CTAT_HumanFusionLib`
data repo (https://github.com/FusionAnnotator/CTAT_HumanFusionLib,
`fusion_lib_data/data/Babiceanu_Normal.dat` and `Greger_Normal.txt`), which is
the standard community aggregator these two specific published gene-pair
lists are distributed through (consumed by STAR-Fusion/FusionInspector's own
official reference-building pipeline). That repo itself carries no explicit
license file; the underlying content being redistributed here is factual
data (gene-pair findings) extracted from the two peer-reviewed papers above,
not the aggregator's own original creative work, and is used here for
internal, non-commercial research QC purposes with full citation to the
primary literature. If stricter license provenance is ever required,
replace this file (or supply an additional one via `custom_blacklists`) with
an independently-sourced equivalent -- see the sv_gene_fusion_annotation.py
docstring for the exact expected format.

## Important caveat (documented here and repeated in the dashboard)

Both lists were built for **RNA-level fusion transcripts** (chimeric
mRNAs detected by RNA-seq), not DNA-level structural-variant breakpoints.
Using them to annotate Manta-called BND events is a reasonable, standard
adaptation (a genomic BND SV whose two breakpoints land in two different
genes is exactly the kind of DNA event that *could* produce such a fusion
transcript if both ends are transcribed/spliced together), but it is an
adaptation, not an exact match to either list's original validation context.
Treat whitelist/blacklist hits as a prioritization aid for manual review,
not a standalone pathogenicity call.

## Custom lists

Both are plain TSVs with a `gene1<TAB>gene2` header (blacklist has an
optional third `source` column, ignored by the loader beyond that). Any
additional file in the same two-column-minimum format works as a custom
whitelist/blacklist; see `run_all_resolveome.py`'s `sv_gene_annotation` stage
and the JSON config's `custom_whitelists`/`custom_blacklists` arrays.
