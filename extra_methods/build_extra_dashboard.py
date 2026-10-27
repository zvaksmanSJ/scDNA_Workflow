#!/usr/bin/env python3
"""Build a self-contained, tabbed HTML dashboard presenting every result
from the ResolveOME DNA extra-methods workflow:
  - Overview: BaseJumper-vs-reference-paper-vs-our-workflow methods/QC
    comparison.
  - QC: MultiQC-style per-cell overview (figure + table) plus MAPD,
    calibrated FDR, and cross-cell artifact-panel metrics, plus the real
    native MultiQC report embedded for the complete interactive view.
  - Phasing: per-cell Beagle phasing results.
  - LOH: Control-FREEC loss-of-heterozygosity calls, both with ploidy fixed
    at 2 and with Control-FREEC selecting the best-fit ploidy per cell from
    a candidate list (its own built-in feature).
  - LOH by region: cohort-wide LOH recurrence heatmap + table.
  - Lineage: Sequoia tree + independent bootstrap/resample verification.
  - Mutational signatures: COSMIC signature refitting on cell-private SNVs.

No external network/CDN dependency: figures are embedded as base64 PNG data
URIs, the Newick-tree renderer is inline vanilla-JS + SVG, and tab
switching is plain JS -- the dashboard works fully offline as one file
(the embedded MultiQC report is the one exception, included via <iframe>
pointing at a relative path since inlining a multi-MB report is
impractical; everything else is self-contained).

Uses string.Template ($-style placeholders) rather than f-strings for the
HTML/CSS/JS template, since the embedded JS/CSS use literal `{ }`
extensively and f-string brace-escaping for a block this size is
error-prone.
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from pathlib import Path
from string import Template


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table_html(rows: list[dict], columns: list[str] | None = None, max_rows: int | None = None) -> str:
    if not rows:
        return "<p class='missing'>Not yet available.</p>"
    cols = columns or list(rows[0].keys())
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    shown = rows[:max_rows] if max_rows else rows
    body = ""
    for row in shown:
        body += "<tr>" + "".join(f"<td>{html.escape(str(row.get(c, '')))}</td>" for c in cols) + "</tr>\n"
    note = ""
    if max_rows and len(rows) > max_rows:
        note = f"<p class='section-note'>Showing top {max_rows} of {len(rows)} rows.</p>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>{note}"


def img_data_uri(path: Path) -> str | None:
    if not path.is_file():
        return None
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def img_tag(path: Path, alt: str, max_width: str = "100%") -> str:
    uri = img_data_uri(path)
    if uri is None:
        return "<p class='missing'>Figure not yet available.</p>"
    return f"<img src='{uri}' alt='{html.escape(alt)}' style='max-width:{max_width};border:1px solid #ccc;'>"


def cnv_scatter_grid(cnv_overview_rows: list[dict]) -> str:
    """Per-cell CNVkit genome-wide scatter plots (already produced by the
    existing variants_cnv_sv stage -- CNVkit's own `.dedup-scatter.png`),
    laid out as a labeled grid so every cell's CNV profile is visible at a
    glance alongside its burden metrics."""
    if not cnv_overview_rows:
        return "<p class='missing'>Not yet available.</p>"
    cards = []
    for row in cnv_overview_rows:
        png_path = Path(row.get("scatter_png", "")) if row.get("scatter_png") else None
        fig = img_tag(png_path, f"{row['sample_id']} CNV scatter", max_width="100%") if png_path else \
            "<p class='missing'>Figure not available.</p>"
        altered = row.get("pct_genome_altered", "")
        cards.append(
            "<div style='display:inline-block;width:300px;vertical-align:top;margin:6px;"
            "border:1px solid #ddd;border-radius:6px;padding:6px;background:white;'>"
            f"<div style='font-size:0.8rem;font-weight:bold;margin-bottom:4px;'>{html.escape(row['sample_id'])}</div>"
            f"{fig}"
            f"<div style='font-size:0.78rem;color:#555;margin-top:4px;'>"
            f"{altered}% genome altered &middot; {row.get('n_segments', '')} segments</div>"
            "</div>"
        )
    return "".join(cards)


def ploidy_narrative_html(narrative_rows: list[dict], candidate_rows: list[dict]) -> str:
    """Per-cell table of Control-FREEC's own candidate-ploidy reasoning
    (parsed verbatim from its stdout by ploidy_determination.py -- no new
    selection logic is added here), plus the raw per-candidate RSS/percent-
    explained fit values underneath so the quantitative basis for each
    cell's final ploidy is visible, not just the summary narrative."""
    if not narrative_rows:
        return "<p class='missing'>Not yet available.</p>"
    by_sample: dict[str, list[dict]] = {}
    for row in candidate_rows:
        by_sample.setdefault(row["sample_id"], []).append(row)

    parts = [table_html(narrative_rows, ["sample_id", "final_ploidy", "narrative", "status"])]
    parts.append("<h4>Per-candidate fit values (Control-FREEC's own RSS score and percent of genome "
                  "explained for every ploidy it tried)</h4>")
    for row in narrative_rows:
        sample = row["sample_id"]
        cands = sorted(by_sample.get(sample, []), key=lambda r: int(r["candidate_ploidy"]))
        if not cands:
            continue
        parts.append(f"<p class='section-note' style='margin-bottom:2px;'><b>{html.escape(sample)}</b> "
                      f"(final ploidy: {html.escape(str(row.get('final_ploidy', '')))})</p>")
        parts.append(table_html(cands, ["candidate_ploidy", "rss_score", "pct_genome_explained"]))
    return "".join(parts)


def find_tree_file(lineage_dir: Path, donor_id: str) -> Path | None:
    candidate = lineage_dir / f"{donor_id}_snv_tree_with_branch_length_selectedscheme.tree"
    if candidate.is_file():
        return candidate
    matches = list(lineage_dir.glob(f"{donor_id}_*tree_with_branch_length_selectedscheme.tree"))
    return matches[0] if matches else None


def find_official_basej_lineage_tree(project_root: Path) -> Path | None:
    """Fallback source for the Lineage tab's tree, used when the standalone
    lineage_track.py output (find_tree_file(), above) isn't present. When
    the real, official BaseJumper `source_lineage` Nextflow stage succeeds
    (run_all_resolveome.py's normal/preferred path -- lineage_track.py is a
    deliberate standalone approximation used only as a fallback when the
    official stage can't be run), its real tree lives under
    basej_lineage/.../secondary_analyses/lineage_artifacts/
    lineage_artifacts_<cohort>/trees/, named by that pipeline's own
    --cohort_id (e.g. "somatic_cells"), not this project's --donor-id --
    so find_tree_file()'s donor-id-based glob never matches it. Prefer the
    'both' (SNV+indel combined) tree if present, then fall back to snv/indel
    alone, matching --sequoia_phylogeny_mode snv,indel,both's own priority.
    """
    trees_glob = "basej_lineage/workflow_outputs/*/*/secondary_analyses/lineage_artifacts/lineage_artifacts_*/trees"
    for trees_dir in project_root.glob(trees_glob):
        for suffix in ("both", "snv", "indel"):
            matches = list(trees_dir.glob(f"*_{suffix}_tree_with_branch_length_selectedscheme.tree"))
            if matches:
                return matches[0]
    return None


def sv_similarity_section(project_root: Path) -> str:
    """Secondary, orthogonal lineage signal: structural-variant (Manta)
    sharing clustering, already computed by the existing sv_similarity /
    sv_similarity_shared stages (step6_sv_similarity.py) -- never previously
    surfaced in this dashboard at all. This is ADDITIONAL to, and does NOT
    replace, the primary Sequoia/MPBoot SNV phylogeny tree rendered above:
    it is average-linkage Jaccard-distance clustering of shared Manta SV
    breakpoints across cells, not a rooted phylogenetic inference. Two
    variants are shown, exactly as the pipeline produces them:
      - "all features": every detected SV breakpoint, including private
        (single-cell) ones -- more features, more sensitive to per-cell
        calling noise.
      - "shared features only" (>=2 supporting cells): restricted to SVs
        seen in more than one cell -- fewer features, but each one is by
        definition not a single-cell-only artifact.
    IMPORTANT CAVEAT (reproduced verbatim from step6_sv_similarity.py's own
    documented limitation, not new editorializing added here): SV detection
    sensitivity correlates with per-cell sequencing depth/coverage, so two
    cells can cluster together because they share similar coverage/calling
    characteristics rather than because they share true common ancestry.
    This view should be read as a complementary, exploratory cross-check
    alongside the primary SNV tree -- concordant clustering across both
    methods is a stronger signal than either alone; discordance is not
    automatically evidence against the SNV tree, since this method has its
    own distinct failure modes.
    """
    variants = [
        ("all SV features (including private, single-cell-only calls)", project_root / "sv_similarity"),
        ("shared SV features only (>=2 supporting cells)", project_root / "sv_similarity_shared"),
    ]
    blocks = []
    any_present = False
    for label, directory in variants:
        dendro = directory / "sv_similarity.dendrogram.png"
        heatmap = directory / "sv_similarity.distance_heatmap.png"
        summary_rows = read_tsv(directory / "sample_sv_summary.tsv")
        report_html = directory / "sv_similarity.report.html"
        if not (dendro.is_file() or summary_rows):
            blocks.append(f"<h4>{html.escape(label)}</h4><p class='missing'>Not available "
                           "(no SV features met the detection/prevalence threshold for this cohort).</p>")
            continue
        any_present = True
        report_link = (
            f"<p class='section-note'>Full interactive report: "
            f"<a href='../{directory.name}/sv_similarity.report.html' target='_blank'>"
            f"{html.escape(directory.name)}/sv_similarity.report.html</a></p>"
            if report_html.is_file() else ""
        )
        blocks.append(
            f"<h4>{html.escape(label)}</h4>"
            "<div style='display:flex;gap:16px;flex-wrap:wrap;align-items:flex-start;'>"
            f"<div style='flex:1;min-width:320px;'>{img_tag(dendro, f'{label} -- dendrogram')}</div>"
            f"<div style='flex:1;min-width:320px;'>{img_tag(heatmap, f'{label} -- pairwise distance heatmap')}</div>"
            "</div>"
            + report_link
            + table_html(summary_rows, ["sample_id", "detected_features", "private_features",
                                        "shared_features", "nearest_cell", "nearest_jaccard_distance"],
                         max_rows=50)
        )
    if not any_present:
        return "<p class='missing'>Not available (sv_similarity / sv_similarity_shared stages did not " \
               "detect enough SV features in this cohort to cluster).</p>"
    return "".join(blocks)


def sv_prevalence_histogram(project_root: Path) -> str:
    """Pure HTML/CSS horizontal bar chart (no new plotting dependency) of
    how many distinct SV breakpoints were detected in exactly how many
    cells, bucketed from sv_similarity's own sv_feature_prevalence.tsv
    (every detected SV site x its SVTYPE and detected_cells count).
    Dominated by private (single-cell) calls is normal/expected -- most
    SV calls from single-cell WGA data are private artifacts or
    cell-specific somatic events; features shared across many cells are
    the ones worth prioritizing as likely-clonal/lineage-informative."""
    path = project_root / "sv_similarity" / "sv_feature_prevalence.tsv"
    if not path.is_file():
        return "<p class='missing'>Not available.</p>"
    buckets = [(1, 1, "1 (private)"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"),
               (5, 9, "5-9"), (10, 19, "10-19"), (20, 10**9, "20+")]
    counts = {label: 0 for _, _, label in buckets}
    svtype_counts: dict[str, int] = {}
    total = 0
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            svtype_counts[row["SVTYPE"]] = svtype_counts.get(row["SVTYPE"], 0) + 1
            n = int(row["detected_cells"])
            for low, high, label in buckets:
                if low <= n <= high:
                    counts[label] += 1
                    break
    if total == 0:
        return "<p class='missing'>Not available.</p>"
    max_count = max(counts.values()) or 1
    bars = []
    for _, _, label in buckets:
        count = counts[label]
        width_pct = 100 * count / max_count
        bars.append(
            "<div style='display:flex;align-items:center;gap:8px;margin:2px 0;'>"
            f"<div style='width:70px;font-size:0.82rem;text-align:right;'>{html.escape(label)}</div>"
            f"<div style='background:#2c7fb8;height:16px;width:{width_pct:.1f}%;min-width:2px;'></div>"
            f"<div style='font-size:0.8rem;color:#555;'>{count:,} ({100*count/total:.1f}%)</div>"
            "</div>"
        )
    svtype_table = "".join(
        f"<tr><td>{html.escape(t)}</td><td>{c:,}</td><td>{100*c/total:.1f}%</td></tr>"
        for t, c in sorted(svtype_counts.items(), key=lambda kv: -kv[1])
    )
    return (
        f"<p class='section-note'>{total:,} distinct SV breakpoints detected across the cohort "
        "(Manta, per-cell, before any similarity filtering). \"Detected in N cells\" buckets how "
        "many of those are private vs. recurrently shared.</p>"
        "<div style='display:flex;gap:24px;flex-wrap:wrap;'>"
        f"<div style='flex:2;min-width:360px;'><h4>Detected in how many cells</h4>{''.join(bars)}</div>"
        "<div style='flex:1;min-width:220px;'><h4>By SV type</h4><table><tr><th>SVTYPE</th>"
        f"<th>Count</th><th>%</th></tr>{svtype_table}</table></div>"
        "</div>"
    )


def sv_fusion_provenance_boxes(provenance: dict) -> str:
    """Render sv_gene_fusion_annotation.py's own provenance.json as labeled
    boxes (same visual pattern as the library-QC tab's provenance boxes),
    so it's clear at a glance exactly which whitelist/blacklist file(s)
    this run actually used -- bundled default(s), any user-supplied
    custom additional file(s), and how many pairs each contributed."""
    if not provenance:
        return "<p class='missing'>Not yet available.</p>"
    boxes = []
    for key in ("whitelist", "blacklist"):
        value = provenance.get(key)
        if not isinstance(value, dict):
            continue
        body = "<br>".join(f"<b>{html.escape(k)}:</b> {html.escape(str(v))}" for k, v in value.items())
        boxes.append(f"<div class='prov-box'><div class='prov-box-title'>{html.escape(key)}</div>"
                     f"<div class='prov-box-body'>{body}</div></div>")
    for key in ("n_svs_total", "n_bnd_events", "n_bnd_fusion_pair_checked", "n_whitelist_hits",
                "n_blacklist_hits", "n_svs_in_main_table", "min_supporting_cells_for_main_table",
                "figures_generated"):
        if key in provenance:
            boxes.append(f"<div class='prov-box'><div class='prov-box-title'>{html.escape(key)}</div>"
                         f"<div class='prov-box-body'>{html.escape(str(provenance[key]))}</div></div>")
    return "<div class='prov-grid'>" + "".join(boxes) + "</div>"


def sv_fusion_top_table_html(top_rows: list[dict], sv_fusion_dir: Path) -> str:
    """Top-N ranked SV table: unlike table_html()'s generic column-dump
    renderer, this embeds each row's OncoScan-style figure (when one was
    generated for it) directly alongside its text explanation."""
    if not top_rows:
        return "<p class='missing'>Not yet available.</p>"
    parts = ["<table><thead><tr><th>Rank</th><th>Type</th><th>Region</th><th>Cells sharing</th>"
             "<th>Whitelist</th><th>Blacklist</th><th>Explanation</th><th>Figure</th></tr></thead><tbody>"]
    for row in top_rows:
        figure_relpath = row.get("figure_png", "")
        figure_html = img_tag(sv_fusion_dir / figure_relpath, f"SV rank {row.get('rank')}", max_width="360px") \
            if figure_relpath else "<span class='missing'>(not generated for this rank)</span>"
        whitelist_mark = "&#10003; WHITELIST" if row.get("whitelist_hit") == "True" else ""
        blacklist_mark = "&#10003; BLACKLIST" if row.get("blacklist_hit") == "True" else ""
        parts.append(
            "<tr>"
            f"<td>{html.escape(str(row.get('rank', '')))}</td>"
            f"<td>{html.escape(row.get('SVTYPE', ''))}</td>"
            f"<td>{html.escape(row.get('region', ''))}</td>"
            f"<td>{html.escape(row.get('supporting_cells', ''))}/{html.escape(row.get('n_cohort_cells', ''))} "
            f"({html.escape(row.get('pct_cells', ''))})</td>"
            f"<td>{whitelist_mark}</td><td>{blacklist_mark}</td>"
            f"<td style='max-width:420px'>{html.escape(row.get('explanation', ''))}</td>"
            f"<td style='min-width:360px'>{figure_html}</td>"
            "</tr>"
        )
    parts.append("</tbody></table>")
    return "".join(parts)


METHODS_COMPARISON_ROWS = [
    ("Alignment", "Sentieon bwa mem (licensed) or BWA-MEM2 open-source fallback", "BWA-MEM2 (open source; no Sentieon license available)", "BWA-MEM2 (unchanged)", "Sentieon (licensed BWA-MEM reimplementation)"),
    ("Post-alignment processing", "Sentieon LocusCollector/Dedup (or samtools markdup fallback for BWA-MEM2 path)", "samtools fixmate/markdup only", "samtools fixmate/markdup only (unchanged)", "Sentieon Dedup + Realigner (indel realignment) + QualCal (BQSR)"),
    ("Phasing", "Not performed", "Not performed", "Beagle 5.5 (GPLv3, free) + 1000 Genomes GRCh38 panel &mdash; <b>new in this workflow</b>", "Eagle v2.4.1 + 1000 Genomes panel"),
    ("Small-variant calling", "DeepVariant (basej-google-deepvariant), no HP channel, same sc_af.model.ckpt checkpoint", "DeepVariant (BioSkryb checkpoint)", "DeepVariant (unchanged)", "SCAN2 (bespoke PTA-artifact statistical model wrapping GATK HaplotypeCaller)"),
    ("Somatic calling", "Sentieon DNAscope (basej-somatic primary caller)", "Population AF prior + cohort-internal BaseJ-style filter (no matched bulk)", "Unchanged", "Matched bulk tissue from the same donor, jointly genotyped"),
    ("False-positive control", "Population AF database only", "Population AF database only", "+ Cross-cell recurrence panel from our own cohort &mdash; <b>new</b>", "Cross-sample panel: 52 cells + 18 bulks from 17 other individuals"),
    ("CNV/aneuploidy", "Ginkgo only (no CNVkit in BaseJumper source)", "CNVkit + Ginkgo (two independent callers) &mdash; <b>our own addition vs. BaseJumper</b>", "Unchanged", "SCAN2 depth-binning (QC-oriented, no dedicated segmentation caller described)"),
    ("Ploidy estimation", "Not described", "Not performed (CNVkit assumed fixed diploid baseline)", "Control-FREEC multi-candidate ploidy selection (its own built-in BIC-like best-fit feature) &mdash; <b>new</b>", "Not described"),
    ("Loss of heterozygosity (LOH)", "Not present in BaseJumper source", "Not performed", "Control-FREEC native BAF-based LOH genotyping &mdash; <b>new</b>", "SCAN2 het-site tracking (QC use, not a dedicated LOH caller)"),
    ("Structural variants", "Not present in BaseJumper source", "Manta &mdash; <b>our own addition vs. BaseJumper</b>", "Unchanged", "Not described in available Methods"),
    ("Annotation", "Not present (basej-somatic has optional VEP for germline filtering only)", "VEP + SnpEff", "Unchanged", "None described"),
    ("Lineage reconstruction", "Sequoia (basej-lineage), from Sentieon DNAscope calls", "Not previously performed directly (official basej-lineage stage still pending in main pipeline)", "Sequoia, directly from per-cell DeepVariant calls, plus independent bootstrap/resample verification via <code>ape::dist.topo</code> &mdash; <b>new, lightweight, same-session track</b>", "Sequoia, from SCAN2 calls + matched-bulk comparison"),
    ("SV similarity/clustering", "Not present in BaseJumper", "Our own exploratory SV-similarity analysis &mdash; <b>our own addition vs. BaseJumper</b>", "Unchanged", "Not described"),
    ("Mutational signatures", "Not present", "Not performed", "MutationalPatterns (Bioconductor, free) COSMIC-signature NNLS refitting on cell-private SNVs &mdash; <b>new</b>", "HDP (hierarchical Dirichlet process)"),
    ("Amplification-evenness QC", "Not computed", "Not computed", "MAPD (median adjacent-bin log2 diff, from existing CNVkit bins) &mdash; <b>new</b>", "MAPD (from SCAN2 summary objects)"),
    ("Calibrated FDR reporting", "Not computed", "PASS/FAIL gate only, no error-rate estimate", "Per-cell calibrated FDR from literature PTA false-positive rate &mdash; <b>new</b>", "Per-cell FDR from synthetic-data-derived false-positive rate"),
    ("Orthogonal validation", "None described", "Internal software/parity QC only", "Unchanged", "6 independent duplex-sequencing technologies + bulk-lineage cross-check"),
]


PAGE_TEMPLATE = Template(r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ResolveOME DNA: Extra Methods QC &amp; Lineage Dashboard</title>
<style>
  body { font-family: -apple-system, Segoe UI, Arial, sans-serif; margin: 0; color: #1a1a1a; background: #fafafa; }
  header { padding: 1rem 2rem 0.3rem 2rem; }
  h1 { font-size: 1.5rem; margin-bottom: 0.2rem; }
  h2 { margin-top: 1.8rem; border-bottom: 2px solid #ddd; padding-bottom: 0.3rem; }
  h3 { margin-top: 1.4rem; }
  table { border-collapse: collapse; width: 100%; margin: 1rem 0; background: white; }
  th, td { border: 1px solid #ccc; padding: 6px 10px; font-size: 0.82rem; text-align: left; }
  th { background: #2c3e50; color: white; position: sticky; top: 0; }
  tr:nth-child(even) { background: #f4f4f4; }
  .missing { color: #999; font-style: italic; }
  .new-tag { background: #2ecc71; color: white; padding: 1px 6px; border-radius: 3px; font-size: 0.75rem; }
  .caveat { background: #fff8e1; border-left: 4px solid #f39c12; padding: 0.8rem 1rem; margin: 1rem 0; font-size: 0.9rem; }
  .section-note { color: #555; font-size: 0.9rem; margin-bottom: 0.5rem; }
  #tree-svg-container { background: white; border: 1px solid #ccc; padding: 1rem; overflow-x: auto; }
  .status-badge { display: inline-block; padding: 2px 10px; border-radius: 10px; font-size: 0.8rem; color: white; }
  .status-complete { background: #2ecc71; }
  .status-pending { background: #e67e22; }

  .tabbar { display: flex; flex-wrap: wrap; gap: 2px; padding: 0 2rem; background: #fafafa;
            position: sticky; top: 0; z-index: 10; border-bottom: 2px solid #2c3e50; }
  .tablink { background: #dfe6e9; border: none; padding: 10px 18px; cursor: pointer; font-size: 0.95rem;
             border-radius: 6px 6px 0 0; color: #2c3e50; }
  .tablink:hover { background: #b2bec3; }
  .tablink.active { background: #2c3e50; color: white; font-weight: bold; }
  .tabcontent { display: none; padding: 0.5rem 2rem 2rem 2rem; }
  .tabcontent.active { display: block; }
  iframe.multiqc-embed { width: 100%; height: 1400px; border: 1px solid #ccc; }
  .fig-wrap { margin: 1rem 0; }
  .prov-grid { display: flex; flex-wrap: wrap; gap: 10px; margin: 0.5rem 0 1.5rem 0; }
  .prov-box { background: #f4f6f8; border: 1px solid #ccd3d9; border-radius: 6px; padding: 8px 12px;
              min-width: 220px; max-width: 420px; font-size: 0.82rem; }
  .prov-box-title { font-weight: bold; color: #2c3e50; margin-bottom: 4px; word-break: break-word; }
  .prov-box-body { color: #333; word-break: break-word; }
</style>
</head>
<body>
<header>
<h1>ResolveOME DNA &mdash; Extra Methods QC &amp; Cell-Lineage Dashboard</h1>
<p class="section-note">Cohort: 15 cells (SY5Y_Cell1-16, excluding Cell4 per explicit instruction: corrupt R2 FASTQ).
Generated alongside, not replacing, the official BaseJumper-parity workflow dashboard.</p>
</header>

<div class="tabbar">
  <button class="tablink active" onclick="showTab(event,'tab-overview')">Overview</button>
  <button class="tablink" onclick="showTab(event,'tab-qc')">QC</button>
  <button class="tablink" onclick="showTab(event,'tab-libraryqc')">Library QC</button>
  <button class="tablink" onclick="showTab(event,'tab-phasing')">Phasing</button>
  <button class="tablink" onclick="showTab(event,'tab-cnv')">CNV</button>
  <button class="tablink" onclick="showTab(event,'tab-loh')">LOH</button>
  <button class="tablink" onclick="showTab(event,'tab-lohregion')">LOH by Region</button>
  <button class="tablink" onclick="showTab(event,'tab-lineage')">Lineage</button>
  <button class="tablink" onclick="showTab(event,'tab-sv')">Structural Variants</button>
  <button class="tablink" onclick="showTab(event,'tab-svfusion')">SV Gene Fusions</button>
  <button class="tablink" onclick="showTab(event,'tab-mutsig')">Mutational Signatures</button>
</div>

<div id="tab-overview" class="tabcontent active">
<h2>Methods comparison: our prior workflow vs. new additions vs. reference paper</h2>
<p class="section-note">Reference paper: Park lab/consortium PTA single-cell lung &amp; colon study
(bioRxiv 2025.10.31.685648). Rows marked <span class="new-tag">new</span> in the methods table were added to our
workflow as a direct result of this comparison.</p>
<table>
<thead><tr><th>Stage</th><th>BaseJumper (official)</th><th>Our workflow (prior to today)</th><th>Our workflow (current, with new additions)</th><th>Reference paper</th></tr></thead>
<tbody>
$methods_rows_html
</tbody>
</table>

<div class="caveat">
<b>Scope notes on the new additions:</b>
<ul>
<li><b>Phasing</b>: run per-autosome (chr1-22) with Beagle 5.5 against the official 1000 Genomes GRCh38
high-coverage phased reference panel (IGSR/EBI, public/open, no license restriction). Beagle is GPLv3 (fully free
for any use); Eagle2 was not used due to its more restrictive non-commercial academic license terms. chrX is
intentionally excluded: this cohort is male, so chrX is hemizygous outside the small pseudoautosomal regions, and
the panel's chrX reference VCF mixes haploid/diploid genotype representations that Beagle cannot parse in one pass.</li>
<li><b>Ploidy / LOH</b>: Control-FREEC (Boeva et al. 2011/2012, GPL) computes BAF directly from BAM + the same
1000 Genomes SNP positions used for phasing; LOH is read from FREEC's own native genotype field (a single
repeated parental-allele letter), not a custom classifier. Ploidy was estimated two ways: a fixed diploid run
(matching common practice / BaseJumper/CNVkit's implicit assumption) and a multi-candidate run letting FREEC's own
built-in selection choose the best-fit value per cell &mdash; both are shown in the LOH tab so the assumption's
impact is visible rather than hidden.</li>
<li><b>Cross-cell artifact panel</b>: built from our own 15-cell clonal cohort (no matched bulk or unrelated
individuals available), so unlike the reference paper's cross-individual panel, this cannot fully distinguish
systematic technical artifacts from true clonal germline variants shared by allele-dropout noise. Sites are
reported in three explicit buckets (UNIVERSAL / SINGLETON / RECURRENT_SUBSET) rather than silently merged.</li>
<li><b>Calibrated FDR</b>: applies the published Luquette et al. synthetic-data PTA false-positive rate
(0.013 SNV/Mb, 0.00073 indel/Mb) to each cell's own callable genome size (Picard WGS metrics) and its full
PASS call set (germline + somatic combined) &mdash; not a curated private-somatic-only call set as in the
reference paper, so these FDR values describe the overall call set's expected artifact burden, not a
somatic-mutation-specific error rate.</li>
<li><b>Cell lineage</b>: built directly from each cell's own per-cell DeepVariant calls (NV/NR read-count
matrices), not from a full per-cell pileup across all BAMs as in the reference paper, and restricted to a
reproducibly-subsampled set of lineage-informative (non-invariant, non-near-universal) SNV sites for
MPBoot bootstrap tractability.</li>
<li><b>Mutational signatures</b>: there is no matched-bulk sample, so there is no curated/validated
private-somatic call set available. The SINGLETON bucket from the artifact panel (sites PASS in exactly one of
the 15 cells) is used as the best available proxy for cell-private mutations; this will still contain residual
allele-dropout/amplification noise a matched bulk would otherwise let us subtract. Refitting (not de novo
extraction) was used against the published COSMIC catalog because per-cell mutation counts are too low to support
reliable de novo signature discovery.</li>
</ul>
</div>
</div>

<div id="tab-qc" class="tabcontent">
<h2>QC overview</h2>
<p class="section-note">Re-presents existing MultiQC aggregate data (not a re-run of MultiQC) as a quick per-cell
overview figure and table, plus the three extra_methods QC additions. The full native interactive MultiQC report
is embedded below for the complete per-tool view (FastQC, fastp, Picard, samtools, SnpEff, etc.).</p>
<div class="fig-wrap">$qc_overview_fig</div>
<h3>Per-cell QC summary table</h3>
$qc_summary_table

<h3>MAPD (amplification-evenness QC)</h3>
<p class="section-note">Median Absolute Pairwise Difference between adjacent CNVkit bins' log2 ratios &mdash;
lower is more even amplification. Computed directly from existing per-cell CNVkit .cnr output (no new upstream
tool required).</p>
$mapd_table

<h3>Calibrated FDR (SNV/indel)</h3>
$fdr_table

<h3>Cross-cell artifact panel summary</h3>
$artifact_table

<h3>Full native MultiQC report</h3>
<p class="section-note">Embedded via iframe from <code>$multiqc_relpath</code> (relative to this dashboard's final
location). If figures below appear blank, open that file directly.</p>
<iframe class="multiqc-embed" src="$multiqc_relpath"></iframe>
</div>

<div id="tab-libraryqc" class="tabcontent">
<h2>CUSTOM library QC comparison</h2>
<p class="section-note">Embeds the existing <code>library_dashboard</code> stage's own standalone page
(<code>qc/library_dashboard/source_threshold_dashboard.html</code>) -- an independent, per-cell comparison of
this workflow's own measured library-prep/sequencing metrics against the licensed source pipeline's five
configured pass/fail thresholds (preseq library complexity, chimeric-read rate, mitochondrial depth fraction,
and CNV MAPD/skew noise). <b>Informational only</b>: this comparison never replaces or overrides this
workflow's own cell-QC gate (<a href="#" onclick="showTab(event,'tab-qc');return false;">see the QC tab</a>
for the gate actually used). The embedded page itself documents, per cell and per metric, exactly why any
FAIL or UNKNOWN call was made, the rationale behind each of the five configured cutoffs, and full run
provenance.</p>
<iframe class="multiqc-embed" src="$library_qc_relpath"></iframe>
</div>

<div id="tab-phasing" class="tabcontent">
<h2>Phasing (Beagle 5.5 + 1000 Genomes GRCh38 panel)</h2>
<p class="section-note">Per-cell statistical phasing of heterozygous SNVs, run per-autosome (chr1-22) against the
official 1000 Genomes GRCh38 high-coverage phased reference panel. chrX is intentionally excluded (see Overview
tab caveats). <code>het_phased_fraction</code> is the fraction of heterozygous genotype calls that received a
real (non-"/") phase &mdash; this is the het-SNP backbone the LOH tab's genotype calls build on.</p>
$phasing_table
</div>

<div id="tab-cnv" class="tabcontent">
<h2>CNV burden &amp; per-cell profiles (CNVkit)</h2>
<p class="section-note">Genome-wide copy-number burden computed directly from CNVkit's own segment-level integer
copy-number calls (<code>.dedup.call.cns</code>, already produced by the existing <code>variants_cnv_sv</code>
stage -- no new upstream caller). <code>pct_genome_altered</code> is the fraction of the analyzed genome in a
segment with a called copy number other than 2 (diploid); <code>pct_genome_gain</code>/<code>pct_genome_loss</code>
split that into amplified vs. deleted. A handful of very short, very-high-copy-number segments (large
<code>max_cn</code> values) are common CNVkit single-cell noise/artifact and should be read alongside the
genome-wide fraction, not in isolation.</p>
$cnv_overview_table
<h3>Per-cell CNV scatter plots</h3>
<p class="section-note">CNVkit's own genome-wide log2-ratio scatter plot per cell (already generated during
variant calling; re-presented here for a quick side-by-side view across the cohort).</p>
<div>$cnv_scatter_grid</div>
</div>

<div id="tab-loh" class="tabcontent">
<h2>Loss of heterozygosity &amp; ploidy (Control-FREEC)</h2>
<p class="section-note">Control-FREEC (Boeva et al., <i>Bioinformatics</i> 2011/2012; GPL license, bioconda) computes
a B-allele-frequency profile directly from each cell's BAM plus the 1000 Genomes SNP positions (no matched bulk
required for WGS) and reports LOH natively in its own genotype field (a single repeated parental-allele letter,
e.g. "A" or "B", rather than a mixed "AB"-type genotype). No custom LOH classifier was written &mdash; this table
only aggregates FREEC's own native output. <code>loh_genome_fraction</code> is the fraction of the analyzed genome
(chr1-22+X) falling in a region FREEC assigned a single-parental-allele genotype.</p>

<h3>LOH calls (ploidy fixed at 2)</h3>
$loh_table

<h3>Ploidy estimation (Control-FREEC candidate selection)</h3>
<p class="section-note">Control-FREEC was re-run letting it choose the best-fit ploidy per cell from candidates
{2,3,4,5,6,7,8} (its own built-in feature: "In case of doubt, you can set different values and Control-FREEC will
select the one that explains most observed CNAs" &mdash; FREEC manual). <code>sample_purity</code> is also FREEC's
own estimate (1.0 = no normal-cell contamination modeled, appropriate for single-cell input with no control).</p>
<div class="caveat">
<b>Important methodological caveat:</b> widening the candidate range from {2,3,4,5} to {2,3,4,5,6,7,8} did not
converge on a stable answer &mdash; <i>more</i> cells moved to the new upper boundary (5 of 15 selected the
maximum candidate, ploidy=8) rather than settling on an interior value. This is a known failure mode of BAF/
read-depth ploidy-selection criteria: higher candidate ploidies offer finer copy-number discretization, which can
mechanically improve the fit statistic without reflecting genuine biology, especially in noisy low-coverage
single-cell data. SY5Y is a neuroblastoma line and plausibly is aneuploid/near-triploid or more complex, but these
specific per-cell point estimates should be treated as suggestive, not definitive, without an orthogonal
validation (karyotyping, DNA-content flow cytometry, or a matched-bulk sample) &mdash; exactly the kind of
independent validation the reference paper has access to (6 technologies) and we do not. Both the fixed-ploidy=2
and the candidate-selected tables are shown so this assumption's impact is visible rather than hidden either way.</div>
$ploidy_table

<h3>How each cell's ploidy was determined</h3>
<p class="section-note">For cells where <code>output_ploidy</code> above is not 2 (not diploid), this shows
Control-FREEC's own stated reasoning -- the per-candidate RSS (residual sum of squares) fit score and percent of
observed copy-number changes explained for every ploidy it tried (2 through 8), followed by its own decision
narrative (including the cases where it deliberately overrode its own best-RSS pick in favor of a lower, nearly-
as-good-fitting ploidy -- FREEC's own documented tie-breaking behavior, not a custom rule added here).</p>
$ploidy_narrative
</div>


<div id="tab-lohregion" class="tabcontent">
<h2>LOH by region (cohort-wide recurrence)</h2>
<p class="section-note">Genome binned into 5Mb windows; a cell is marked LOH-positive in a bin if any
Control-FREEC segment overlapping that bin carries a single-parental-allele genotype. Rows near chromosome
centromeres can show artificially high recurrence because informative heterozygous SNPs are naturally sparse
there (less confident BAF signal), not necessarily shared biology &mdash; interpret very high-recurrence
pericentromeric bins with that in mind.</p>
<div class="fig-wrap">$loh_region_fig</div>
<h3>Most recurrent LOH regions (top 50 by fraction of cohort)</h3>
$loh_region_table
</div>

<div id="tab-lineage" class="tabcontent">
<h2>Cell lineage <span class="status-badge status-$lineage_status_class">$lineage_status</span></h2>
<p class="section-note">$lineage_method_note</p>
<div id="tree-svg-container">
  <p class="missing" id="tree-missing-msg">$tree_missing_msg</p>
  <svg id="tree-svg" width="900" height="500" style="display:$tree_display"></svg>
</div>

<h3>Lineage verification</h3>
<p class="section-note">Two independent checks, neither of which required writing a custom tree-comparison
algorithm: (1) per-branch bootstrap support from the primary tree's own MPBoot replicates, flagged against a
$bootstrap_threshold% threshold; (2) the primary tree compared to a second, independently-built tree from a
different random subsample of lineage-informative sites (different <code>--random-seed</code>), using
<code>ape::dist.topo()</code> (Robinson-Foulds distance) &mdash; <code>ape</code> is an existing dependency already
present in the Sequoia container, so no new package was installed.</p>
$verification_summary_html
<h4>Per-branch bootstrap support (primary tree)</h4>
$bootstrap_table

<h3>Structural-variant similarity clustering (orthogonal cross-check)</h3>
<p class="section-note">Additional, independent lineage-adjacent signal -- does <b>not</b> replace the
primary Sequoia/MPBoot SNV tree above. Average-linkage hierarchical clustering on the Jaccard distance
between cells' shared Manta structural-variant (SV) breakpoints (already computed by the
<code>sv_similarity</code>/<code>sv_similarity_shared</code> stages). <b>Important caveat:</b> SV detection
sensitivity correlates with per-cell sequencing depth/coverage, so two cells can cluster together because
they share similar coverage/calling characteristics rather than true common ancestry -- this is exploratory
similarity clustering, <b>not</b> a rooted phylogenetic inference. Read it as a complementary cross-check:
concordant clustering with the SNV tree above is a stronger signal than either method alone; discordance is
not automatically evidence against the SNV tree, since this method has its own distinct failure modes.</p>
$sv_similarity_section
</div>

<div id="tab-sv" class="tabcontent">
<h2>Structural variants (Manta, per-cell)</h2>
<p class="section-note">Raw structural-variant call data from the existing per-cell Manta calling
(<code>variants_cnv_sv</code> stage) and the <code>sv_similarity</code>/<code>sv_similarity_shared</code>
stages' own feature tables -- no new upstream caller. This tab is about the SV calls themselves (counts,
types, how widely shared each is); for the similarity-clustering <i>use</i> of this same data as an
orthogonal lineage cross-check, see the Lineage tab's "Structural-variant similarity clustering" section
above (same underlying calls, different analysis/purpose -- not duplicated here); for gene/fusion-identity
annotation of these same calls, see the SV Gene Fusions tab next to this one.</p>
<h3>SV prevalence across the cohort</h3>
$sv_prevalence_histogram
<h3>Per-cell SV summary</h3>
<p class="section-note">Per-cell counts of detected SV breakpoints (<code>detected_features</code>), how
many are private to that cell vs. shared with at least one other cell, the per-type breakdown (Manta's own
BND/DEL/DUP/INS classes), and each cell's nearest neighbor by SV-sharing Jaccard distance (same metric used
for the Lineage tab's clustering, surfaced per-cell here for direct inspection).</p>
$sv_percell_table
</div>

<div id="tab-svfusion" class="tabcontent">
<h2>SV gene-fusion annotation <span class="status-badge status-$svfusion_status_class">$svfusion_status</span></h2>
<p class="section-note">Computed directly from the SAME Manta structural-variant calls already produced by the
existing <code>variants_cnv_sv</code> stage (<code>sv_features.tsv</code> / <code>sv_binary.tsv</code>) &mdash;
<b>no new SV caller, no re-calling</b>. Each BND (translocation-type) event's two breakpoints are checked,
gene-pair by gene-pair, against a WHITELIST of known cancer-associated gene fusions (default: the Mitelman
Database of Chromosome Aberrations and Gene Fusions in Cancer, via Arriba's MIT-licensed redistribution) and a
BLACKLIST of gene pairs recurrently reported as fusion transcripts in normal, non-cancer tissue (default:
Babiceanu et al. 2016 + Greger et al. 2014). DEL/DUP spans are separately annotated with every gene they
overlap; no whitelist/blacklist check applies to DEL/DUP/INS (both reference lists are built from RNA-level
fusion-transcript evidence, which only a translocation/breakend can plausibly produce &mdash; see
<code>sv_gene_fusion_annotation.py</code>'s own docstring for the full reasoning). See the Lineage tab for the
separate, exploratory SV-sharing <i>similarity-clustering</i> analysis &mdash; this tab is a different kind of
analysis (gene/fusion identity, not cell-to-cell similarity).</p>
<div class="caveat"><b>Caveat:</b> both reference lists were built from RNA-level fusion-transcript evidence,
not DNA-level SV breakpoints. Applying them to Manta's DNA-level breakpoints here is a standard, reasonable
adaptation in the field, but it is an adaptation &mdash; treat a whitelist/blacklist hit as a prioritization
aid for manual review, not a standalone pathogenicity or benignity call.</div>

<h3>Whitelist / blacklist used this run</h3>
$svfusion_provenance_boxes

<h3>Top-ranked SVs (with explanations and figures)</h3>
<p class="section-note">Ranking: (1) BND events whose two breakpoint genes match the cancer-fusion whitelist,
most-recurrent first; then (2) every other non-blacklisted SV, also most-recurrent first. Blacklisted SVs
(known normal-tissue fusions) are excluded from this ranked list entirely &mdash; shown, flagged, in the full
table below instead. An OncoScan-style figure (chromosome schematic with breakpoint/gene labels, plus a
per-cell presence strip across the whole cohort) is generated for the highest-ranked subset of these rows
(bounds dashboard size/runtime); every ranked row gets a text explanation regardless of whether it has a
figure.</p>
$svfusion_top_table

<h3>All gene-overlapping, recurrent SVs</h3>
<p class="section-note">$svfusion_main_note</p>
$svfusion_main_table
</div>

<div id="tab-mutsig" class="tabcontent">
<h2>Mutational signatures</h2>
<p class="section-note">COSMIC reference-signature refitting (MutationalPatterns, Bioconductor, GPL-3, free) via
NNLS backward-elimination (<code>fit_to_signatures_strict</code>) on cell-private SNVs (see Overview tab caveat on
what "cell-private" means without a matched bulk). Reconstruction cosine similarity &ge;0.9 for all cells
indicates the selected signature combination explains the observed 96-trinucleotide-context spectrum well.</p>
<div class="fig-wrap">$mutsig_spectrum_fig</div>
<div class="fig-wrap">$mutsig_barplot_fig</div>
<h3>Reconstruction quality</h3>
$mutsig_cosine_table
<h3>Signature contributions (absolute mutation counts per signature)</h3>
<p class="section-note">Columns are COSMIC SBS signature IDs; cells omitted below a signature's column had zero
estimated contribution after backward-elimination refitting.</p>
$mutsig_contrib_table
</div>

<script>
function showTab(evt, tabId) {
  document.querySelectorAll('.tabcontent').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tablink').forEach(el => el.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  evt.currentTarget.classList.add('active');
}

const newickString = $newick_json;

function parseNewick(s) {
  const tokens = s.replace(/\s/g, '').match(/([^():,;]+|[():,;])/g);
  let idx = 0;
  function parseNode() {
    const node = { name: '', length: 0, children: [] };
    if (tokens[idx] === '(') {
      idx++;
      node.children.push(parseNode());
      while (tokens[idx] === ',') { idx++; node.children.push(parseNode()); }
      idx++; // ')'
    }
    if (tokens[idx] && tokens[idx] !== ':' && tokens[idx] !== ',' && tokens[idx] !== ';' && tokens[idx] !== ')') {
      node.name = tokens[idx]; idx++;
    }
    if (tokens[idx] === ':') { idx++; node.length = parseFloat(tokens[idx]) || 0; idx++; }
    return node;
  }
  return parseNode();
}

function layoutTree(root) {
  const leaves = [];
  (function collectLeaves(n) { if (!n.children.length) leaves.push(n); else n.children.forEach(collectLeaves); })(root);
  leaves.forEach((leaf, i) => { leaf.y = i; });
  function assignY(n) {
    if (!n.children.length) return n.y;
    const ys = n.children.map(assignY);
    n.y = (Math.min(...ys) + Math.max(...ys)) / 2;
    return n.y;
  }
  assignY(root);
  function assignX(n, depth) {
    n.x = depth;
    n.children.forEach(c => assignX(c, depth + (c.length || 0.01) * 8));
  }
  assignX(root, 0);
  return [root, leaves];
}

if (newickString) {
  const tree = parseNewick(newickString);
  const [root, leaves] = layoutTree(tree);
  const svg = document.getElementById('tree-svg');
  const marginLeft = 40, marginRight = 220, marginTop = 20;
  const rowHeight = 28;
  const height = marginTop * 2 + leaves.length * rowHeight;
  svg.setAttribute('height', height);
  const maxX = Math.max(...leaves.map(l => l.x), 1);
  const plotWidth = 900 - marginLeft - marginRight;
  const ns = 'http://www.w3.org/2000/svg';
  function sx(x) { return marginLeft + (x / maxX) * plotWidth; }
  function sy(y) { return marginTop + y * rowHeight; }
  function drawNode(n) {
    n.children.forEach(c => {
      const hline = document.createElementNS(ns, 'line');
      hline.setAttribute('x1', sx(n.x)); hline.setAttribute('x2', sx(c.x));
      hline.setAttribute('y1', sy(c.y)); hline.setAttribute('y2', sy(c.y));
      hline.setAttribute('stroke', '#333'); hline.setAttribute('stroke-width', '1.5');
      svg.appendChild(hline);
      drawNode(c);
    });
    if (n.children.length > 1) {
      const ys = n.children.map(c => c.y);
      const vline = document.createElementNS(ns, 'line');
      vline.setAttribute('x1', sx(n.x)); vline.setAttribute('x2', sx(n.x));
      vline.setAttribute('y1', sy(Math.min(...ys))); vline.setAttribute('y2', sy(Math.max(...ys)));
      vline.setAttribute('stroke', '#333'); vline.setAttribute('stroke-width', '1.5');
      svg.appendChild(vline);
    }
    if (!n.children.length) {
      const text = document.createElementNS(ns, 'text');
      text.setAttribute('x', sx(n.x) + 6); text.setAttribute('y', sy(n.y) + 4);
      text.setAttribute('font-size', '11'); text.textContent = n.name;
      svg.appendChild(text);
      const dot = document.createElementNS(ns, 'circle');
      dot.setAttribute('cx', sx(n.x)); dot.setAttribute('cy', sy(n.y)); dot.setAttribute('r', 3);
      dot.setAttribute('fill', '#2980b9');
      svg.appendChild(dot);
    }
  }
  drawNode(root);
}
</script>

</body>
</html>
""")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--extra-methods-dir", required=True, type=Path)
    p.add_argument("--donor-id", default="SY5Y_cohort")
    p.add_argument("--multiqc-relpath", default="../components/multiqc/multiqc_report.html",
                    help="Path to the native MultiQC report, relative to this dashboard's own "
                         "final location (used as an <iframe> src, not copied/inlined).")
    p.add_argument("--library-qc-relpath", default="../qc/library_dashboard/source_threshold_dashboard.html",
                    help="Path to step2_library_qc_dashboard.py's own standalone page, relative to this "
                         "dashboard's own final location (used as an <iframe> src, not copied/inlined).")
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()

    EM = args.extra_methods_dir
    mapd_rows = read_tsv(EM / "mapd_qc.tsv")
    fdr_rows = read_tsv(EM / "fdr_estimate.tsv")
    artifact_rows = read_tsv(EM / "artifact_panel_summary.tsv")
    phasing_rows = read_tsv(EM / "phasing_summary.tsv")
    loh_rows = read_tsv(EM / "loh_summary.tsv")
    ploidy_rows = read_tsv(EM / "ploidy_summary.tsv")
    cnv_overview_rows = read_tsv(EM / "cnv_overview.tsv")
    ploidy_candidate_rows = read_tsv(EM / "ploidy_candidates.tsv")
    ploidy_narrative_rows = read_tsv(EM / "ploidy_narrative.tsv")
    qc_summary_rows = read_tsv(EM / "qc_summary.tsv")
    bootstrap_rows = read_tsv(EM / "lineage_verification" / "bootstrap_support_report.tsv")
    mutsig_contrib_rows = read_tsv(EM / "mutational_signatures" / "signature_contributions.tsv")
    mutsig_cosine_rows = read_tsv(EM / "mutational_signatures" / "reconstruction_cosine_similarity.tsv")

    loh_region_rows_all = read_tsv(EM / "loh_by_region.tsv")
    loh_region_rows_all.sort(key=lambda r: float(r.get("fraction_cells_loh", 0) or 0), reverse=True)
    loh_region_rows = [r for r in loh_region_rows_all if float(r.get("fraction_cells_loh", 0) or 0) > 0]

    sv_fusion_dir = EM / "sv_gene_fusion"
    sv_fusion_provenance_path = sv_fusion_dir / "provenance.json"
    sv_fusion_provenance = (
        json.loads(sv_fusion_provenance_path.read_text()) if sv_fusion_provenance_path.is_file() else {}
    )
    top_n = sv_fusion_provenance.get("top_n", 50)
    sv_fusion_top_rows = read_tsv(sv_fusion_dir / f"sv_gene_fusion_top{top_n}.tsv")
    sv_fusion_main_rows_all = read_tsv(sv_fusion_dir / "sv_gene_annotation_main.tsv")

    verification_json_path = EM / "lineage_verification" / "verification_summary.json"
    verification_summary = None
    if verification_json_path.is_file():
        verification_summary = json.loads(verification_json_path.read_text())

    lineage_dir = EM / "lineage"
    tree_file = find_tree_file(lineage_dir, args.donor_id)
    tree_is_official = False
    if tree_file is None:
        tree_file = find_official_basej_lineage_tree(EM.parent)
        tree_is_official = tree_file is not None
    newick = tree_file.read_text().strip() if tree_file else None
    lineage_method_note = (
        "Reconstructed by the official BaseJumper <code>source_lineage</code> Nextflow "
        "pipeline (Sequoia + MPBoot), run on the full somatic-filtering-&gt;"
        "integrated-matrices-&gt;lineage chain -- the same production path BaseJumper "
        "itself uses, not an approximation."
        if tree_is_official else
        "Reconstructed via Sequoia (same tool used by the official BaseJumper lineage "
        "pipeline), run directly on this cohort's own merged DeepVariant calls rather "
        "than waiting for the full somatic-filtering-&gt;integrated-matrices-&gt;lineage "
        "Nextflow chain."
    )

    methods_rows_html = "".join(
        f"<tr><td>{stage}</td><td>{basejumper}</td><td>{our_prev}</td><td>{our_new}</td><td>{paper}</td></tr>\n"
        for stage, basejumper, our_prev, our_new, paper in METHODS_COMPARISON_ROWS
    )

    if verification_summary is not None:
        rf = verification_summary.get("robinson_foulds_distance")
        bsd = verification_summary.get("branch_score_distance")
        # ape::dist.topo()'s branch-score distance can legitimately come back
        # as NA/null (observed with rooted trees under the PH85 method --
        # see dist.topo()'s own "Some trees were rooted" warning), so this
        # must not assume a numeric value is always present.
        bsd_text = f"{bsd:.4f}" if isinstance(bsd, (int, float)) else "not computable for this tree pair (NA from ape::dist.topo)"
        same_leaves = verification_summary.get("same_leaf_set")
        n_low = verification_summary.get("n_low_confidence_branches")
        n_branches = verification_summary.get("n_internal_branches")
        verification_summary_html = (
            "<ul>"
            f"<li>Internal branches: {n_branches}, of which {n_low} fall below the "
            f"{verification_summary.get('bootstrap_threshold')}% bootstrap-support threshold</li>"
            f"<li>Independent resample (different subsampled sites): "
            f"{'same 15-cell leaf set recovered' if same_leaves else 'LEAF SET MISMATCH -- see raw output'}</li>"
            f"<li>Robinson-Foulds topology distance between primary and resampled trees: {rf}</li>"
            f"<li>Branch-score distance: {bsd_text}</li>"
            "</ul>"
        )
    else:
        verification_summary_html = (
            "<p class='missing'>Not applicable: this tree is the official BaseJumper "
            "<code>source_lineage</code> pipeline's own output. That pipeline does not "
            "perform an independent resample/bootstrap-concordance check the way this "
            "project's standalone <code>lineage_track.py</code>/<code>verify_lineage.py</code> "
            "scripts do (used only when the official pipeline itself can't be run); its own "
            "per-branch bootstrap support (from MPBoot) is reported directly in the tree "
            "file instead.</p>"
            if tree_is_official else
            "<p class='missing'>Not yet available.</p>"
        )

    bootstrap_threshold = (
        verification_summary.get("bootstrap_threshold", 70) if verification_summary else 70
    )

    mutsig_contrib_cols = None
    if mutsig_contrib_rows:
        # Drop signature columns that are all-zero across every cell to keep the table readable.
        all_cols = list(mutsig_contrib_rows[0].keys())
        kept = ["sample_id"]
        for c in all_cols:
            if c == "sample_id":
                continue
            if any(float(r.get(c, 0) or 0) > 0 for r in mutsig_contrib_rows):
                kept.append(c)
        mutsig_contrib_cols = kept

    html_doc = PAGE_TEMPLATE.substitute(
        methods_rows_html=methods_rows_html,
        qc_overview_fig=img_tag(EM / "qc_overview.png", "Per-cell QC overview"),
        qc_summary_table=table_html(qc_summary_rows),
        mapd_table=table_html(mapd_rows, ["sample_id", "mapd", "n_bins", "n_adjacent_pairs", "status"]),
        fdr_table=table_html(fdr_rows, ["sample_id", "callable_mb", "n_snv_pass", "n_indel_pass",
                                        "expected_fp_snv", "expected_fp_indel", "fdr_snv", "fdr_indel", "status"]),
        artifact_table=table_html(artifact_rows, ["sample_id", "n_pass_total", "n_universal", "n_singleton",
                                                   "n_recurrent_subset_flagged_artifact"]),
        multiqc_relpath=args.multiqc_relpath,
        library_qc_relpath=args.library_qc_relpath,
        phasing_table=table_html(phasing_rows, ["sample_id", "n_records", "n_het", "n_het_phased",
                                                 "het_phased_fraction", "status"]),
        cnv_overview_table=table_html(cnv_overview_rows, ["sample_id", "n_segments", "pct_genome_gain",
                                                          "pct_genome_loss", "pct_genome_altered",
                                                          "max_cn", "min_cn", "status"]),
        cnv_scatter_grid=cnv_scatter_grid(cnv_overview_rows),
        loh_table=table_html(loh_rows, ["sample_id", "n_regions", "n_loh_regions", "loh_region_fraction",
                                        "loh_bp", "total_bp", "loh_genome_fraction", "status"]),
        ploidy_table=table_html(ploidy_rows, ["sample_id", "output_ploidy", "sample_purity",
                                              "good_polynomial_fit", "status"]),
        ploidy_narrative=ploidy_narrative_html(ploidy_narrative_rows, ploidy_candidate_rows),
        loh_region_fig=img_tag(EM / "loh_by_region_heatmap.png", "Cohort-wide LOH recurrence heatmap"),
        loh_region_table=table_html(loh_region_rows, ["chrom", "start", "end", "n_cells_loh",
                                                       "n_cells_total", "fraction_cells_loh",
                                                       "cells_with_loh"], max_rows=50),
        bootstrap_threshold=bootstrap_threshold,
        verification_summary_html=verification_summary_html,
        bootstrap_table=table_html(bootstrap_rows, ["internal_node_index", "bootstrap_support", "low_confidence"]),
        lineage_status="Complete" if newick else "Pending",
        lineage_status_class="complete" if newick else "pending",
        lineage_method_note=lineage_method_note,
        tree_missing_msg="" if newick else "Tree not yet available \u2014 lineage reconstruction is still running.",
        tree_display="block" if newick else "none",
        newick_json=json.dumps(newick),
        sv_similarity_section=sv_similarity_section(EM.parent),
        sv_prevalence_histogram=sv_prevalence_histogram(EM.parent),
        sv_percell_table=table_html(
            read_tsv(EM.parent / "sv_similarity" / "sample_sv_summary.tsv"),
            ["sample_id", "detected_features", "private_features", "shared_features",
             "BND", "DEL", "DUP", "INS", "nearest_cell", "nearest_jaccard_distance"],
            max_rows=100,
        ),
        svfusion_status="Complete" if sv_fusion_top_rows else "Pending",
        svfusion_status_class="complete" if sv_fusion_top_rows else "pending",
        svfusion_provenance_boxes=sv_fusion_provenance_boxes(sv_fusion_provenance),
        svfusion_top_table=sv_fusion_top_table_html(sv_fusion_top_rows, sv_fusion_dir),
        svfusion_main_note=(
            f"{len(sv_fusion_main_rows_all)} SVs hit at least one annotated gene and are shared by "
            f"&ge;{sv_fusion_provenance.get('min_supporting_cells_for_main_table', 2)} cells (showing the "
            "first 300, sorted by recurrence); the full, unfiltered annotation (every SV, including "
            "non-gene-overlapping and singleton ones) is at "
            "<code>extra_methods/sv_gene_fusion/sv_gene_annotation_full.tsv</code>."
            if sv_fusion_main_rows_all else "Not yet available."
        ),
        svfusion_main_table=table_html(
            sv_fusion_main_rows_all,
            ["VariantId", "SVTYPE", "region", "supporting_cells", "pct_cells",
             "genes_breakpoint1", "genes_breakpoint2", "genes_in_span", "whitelist_hit", "blacklist_hit"],
            max_rows=300,
        ),
        mutsig_spectrum_fig=img_tag(EM / "mutational_signatures" / "mutation_spectrum_96.png",
                                     "96-trinucleotide-context mutation spectrum"),
        mutsig_barplot_fig=img_tag(EM / "mutational_signatures" / "signature_contribution_barplot.png",
                                    "Relative signature contribution per cell"),
        mutsig_cosine_table=table_html(mutsig_cosine_rows),
        mutsig_contrib_table=table_html(mutsig_contrib_rows, mutsig_contrib_cols),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html_doc)
    print(f"[OK] Dashboard written -> {args.output}")


if __name__ == "__main__":
    main()
