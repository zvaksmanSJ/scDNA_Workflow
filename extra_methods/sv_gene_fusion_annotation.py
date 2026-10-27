#!/usr/bin/env python3
"""SV gene-fusion annotation, computed directly from Manta's already-called,
cohort-merged structural variants (`sv_features.tsv` / `sv_binary.tsv`, the
existing step5_filter_variant_inputs.py output) -- NO re-calling, no new SV
caller. This script only intersects already-called breakpoints against a
gene model and two curated gene-pair reference lists.

What this does, per SV, depends on its type:

  * BND (breakend / translocation-type, same- or cross-chromosome): looks up
    the gene(s), if any, overlapping each of its two breakpoints
    (CHROM1:POS1 and CHROM2:POS2). If both breakpoints land in >=1 gene,
    every (gene_at_bp1, gene_at_bp2) combination is checked, order-
    independent, against:
      - a WHITELIST of known cancer-associated gene-pair fusions (default:
        the Mitelman Database of Chromosome Aberrations and Gene Fusions in
        Cancer, via Arriba's MIT-licensed redistribution -- see
        refs/sv_gene_fusion/README.md for full provenance/citations), and
      - a BLACKLIST of gene pairs recurrently reported as fusion transcripts
        in NORMAL, non-cancer tissue (default: Babiceanu et al. 2016 +
        Greger et al. 2014, via the Trinity CTAT project's public data repo
        -- same README.md), i.e. likely transcriptional read-through/
        splicing artifacts rather than somatic driver events.
    This whitelist/blacklist check is intentionally restricted to BND
    events: both reference lists are built from RNA-level FUSION
    TRANSCRIPTS (two genes joined into one chimeric mRNA), which is the
    phenomenon a genomic translocation/breakend can plausibly produce if
    both breakpoints land in transcribed genes -- it is not a meaningful
    check for DEL/DUP/INS (a deletion or duplication does not fuse two
    genes into one transcript the way a translocation can).

  * DEL / DUP (same-chromosome, defines a genomic span): reports every gene
    whose interval overlaps the deleted/duplicated span [POS1, POS2] (not
    just genes exactly at the two breakpoints) -- the biologically relevant
    question for a copy-number-changing SV is "what genes does this
    region's copy number change affect", not "what single gene sits at
    each edge". No whitelist/blacklist check (see above).

  * INS (insertion): reports genes overlapping the (typically very short)
    insertion site only; no span/whitelist/blacklist logic applies.

CAVEAT (also repeated in the dashboard): both the whitelist and the
blacklist were built from RNA-seq/transcript-level evidence, not DNA-level
SV breakpoints. Using them to annotate Manta-called DNA breakpoints is a
standard, reasonable adaptation in the field, but it is an adaptation --
treat a whitelist/blacklist hit as a prioritization aid for manual review,
not a standalone pathogenicity or benignity call.

Outputs (all under --output-dir):
  sv_gene_annotation_full.tsv   -- every SV in sv_features.tsv, annotated
                                   (empty gene/fusion columns where nothing
                                   overlaps/applies); full transparency.
  sv_gene_annotation_main.tsv   -- the same, filtered to SVs that (a) hit at
                                   least one gene somewhere and (b) have
                                   supporting_cells >= --min-supporting-cells
                                   (default 2 -- singleton, single-cell-only
                                   SVs are, as elsewhere in this PTA/WGA
                                   pipeline, far more likely amplification
                                   artifacts than real recurrent events; the
                                   full TSV above still has them for anyone
                                   who wants to look).
  sv_gene_fusion_top<N>.tsv     -- the --top-n (default 50) highest-priority
                                   SVs, each with a plain-English explanation
                                   (see rank_svs()/build_explanation() for the
                                   exact ranking rule and wording), plus a
                                   relative path to an OncoScan-style figure
                                   for the --figures-top-n (default 20)
                                   highest-ranked of those (capped to bound
                                   dashboard size/runtime; all --top-n rows
                                   still get a text explanation regardless).
  figures/rankNN_<variant>.png  -- the OncoScan-style figures themselves.
  provenance.json               -- exactly which whitelist/blacklist file(s)
                                   were used (default and/or any
                                   user-supplied custom ones), counts, and
                                   summary statistics.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
DEFAULT_WHITELIST = HERE / "refs/sv_gene_fusion/whitelist_mitelman_cancer_fusions.tsv"
DEFAULT_BLACKLIST = HERE / "refs/sv_gene_fusion/blacklist_normal_tissue_fusions.tsv"
DEFAULT_GENES_BED = HERE.parent / "resolveome_align/refs/annotation/genes.bed"


def normalize_chrom(chrom: str) -> str:
    """genes.bed uses 'chr1'/'chrX'; sv_features.tsv and the reference .fai
    use '1'/'X' (no prefix) -- normalize everything to the latter so lookups
    always match regardless of which file a chromosome name came from."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def safe_filename(variant_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", variant_id)[:120]


@dataclass
class GeneIndex:
    """Per-chromosome sorted-by-start gene intervals plus a running max-end
    array, enabling O(log n + k) point/range overlap queries (k = number of
    genes actually returned) instead of an O(n) scan per query -- needed
    since this script runs tens of thousands of such queries (2 per SV
    breakpoint x ~100k SVs in a real cohort)."""

    by_chrom: dict

    @classmethod
    def from_bed(cls, path: Path) -> "GeneIndex":
        raw: dict[str, list[tuple[int, int, str]]] = {}
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            raw.setdefault(normalize_chrom(chrom), []).append((start, end, name))
        by_chrom = {}
        for chrom, entries in raw.items():
            entries.sort()
            starts = [e[0] for e in entries]
            running_max = []
            current_max = -1
            for start, end, _name in entries:
                current_max = max(current_max, end)
                running_max.append(current_max)
            by_chrom[chrom] = (entries, starts, running_max)
        return cls(by_chrom)

    def _candidates(self, chrom: str, query_end: int):
        chrom = normalize_chrom(chrom)
        if chrom not in self.by_chrom:
            return None, None, -1
        entries, starts, running_max = self.by_chrom[chrom]
        i = bisect.bisect_right(starts, query_end) - 1
        return entries, running_max, i

    def point(self, chrom: str, pos: int) -> list[str]:
        """Genes overlapping a single coordinate (an SV breakpoint)."""
        entries, running_max, i = self._candidates(chrom, pos)
        if i < 0 or running_max[i] < pos:
            return []
        hits = set()
        j = i
        while j >= 0:
            if running_max[j] < pos:
                break  # no earlier entry's end can reach pos (prefix-max is non-increasing as j decreases)
            start, end, name = entries[j]
            if start <= pos <= end:
                hits.add(name)
            j -= 1
        return sorted(hits)

    def span(self, chrom: str, start: int, end: int) -> list[str]:
        """Genes whose interval overlaps [start, end] (an SV's deleted/
        duplicated region), not just the two breakpoint coordinates."""
        entries, running_max, i = self._candidates(chrom, end)
        hits = set()
        j = i
        while j >= 0:
            if running_max[j] < start:
                break
            gene_start, gene_end, name = entries[j]
            if gene_start <= end and gene_end >= start:
                hits.add(name)
            j -= 1
        return sorted(hits)


@dataclass
class SVRow:
    variant_id: str
    svtype: str
    orientation: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    supporting_records: int
    supporting_cells: int
    breakpoint_tolerance: int


def load_sv_features(path: Path) -> list[SVRow]:
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for r in reader:
            rows.append(SVRow(
                variant_id=r["VariantId"], svtype=r["SVTYPE"], orientation=r["orientation"],
                chrom1=r["CHROM1"], pos1=int(r["POS1"]), chrom2=r["CHROM2"], pos2=int(r["POS2"]),
                supporting_records=int(r["supporting_records"]), supporting_cells=int(r["supporting_cells"]),
                breakpoint_tolerance=int(r["breakpoint_tolerance"]),
            ))
    return rows


def load_sv_binary(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    """Returns (cell_id_column_order, {VariantId: [cell_ids supporting it]})
    from the existing per-cell presence/absence matrix (0/1 per cell)."""
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        cells = header[1:]
        presence = {}
        for row in reader:
            variant_id = row[0]
            presence[variant_id] = [cells[i] for i, value in enumerate(row[1:]) if value == "1"]
    return cells, presence


def load_pair_set(paths: list[Path]) -> tuple[set[frozenset], dict[frozenset, str]]:
    """Union of gene-pair sets from one or more TSVs (header: gene1, gene2,
    optional source). Pairs are stored order-independent (frozenset) since
    a fusion A--B is the same biological event as B--A."""
    pairs: set[frozenset] = set()
    source_of: dict[frozenset, str] = {}
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Gene-pair list not found: {path}")
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            if not {"gene1", "gene2"}.issubset(fields):
                raise ValueError(f"{path}: expected at least columns gene1, gene2 (got {sorted(fields)})")
            for row in reader:
                key = frozenset((row["gene1"], row["gene2"]))
                pairs.add(key)
                source_of.setdefault(key, row.get("source") or path.name)
    return pairs, source_of


def annotate_sv(sv: SVRow, gene_index: GeneIndex,
                 whitelist_pairs: set[frozenset], blacklist_pairs: set[frozenset]) -> dict:
    genes_bp1 = gene_index.point(sv.chrom1, sv.pos1)
    genes_bp2 = gene_index.point(sv.chrom2, sv.pos2)
    same_chrom = normalize_chrom(sv.chrom1) == normalize_chrom(sv.chrom2)

    genes_in_span: list[str] = []
    if same_chrom and sv.svtype in ("DEL", "DUP"):
        low, high = sorted((sv.pos1, sv.pos2))
        genes_in_span = gene_index.span(sv.chrom1, low, high)

    fusion_pair_checked = sv.svtype == "BND" and bool(genes_bp1) and bool(genes_bp2)
    fusion_pair = None
    whitelist_hit = False
    whitelist_matched_pair = None
    blacklist_hit = False
    blacklist_matched_pair = None
    if fusion_pair_checked:
        for gene1 in genes_bp1:
            for gene2 in genes_bp2:
                if gene1 == gene2:
                    continue  # same gene at both ends: an intragenic rearrangement, not a fusion
                key = frozenset((gene1, gene2))
                if fusion_pair is None:
                    fusion_pair = (gene1, gene2)  # first candidate pair, reported even without a list match
                if key in whitelist_pairs and not whitelist_hit:
                    whitelist_hit = True
                    whitelist_matched_pair = (gene1, gene2)
                    fusion_pair = (gene1, gene2)
                if key in blacklist_pairs and not blacklist_hit:
                    blacklist_hit = True
                    blacklist_matched_pair = (gene1, gene2)

    return {
        "genes_breakpoint1": genes_bp1, "genes_breakpoint2": genes_bp2, "genes_in_span": genes_in_span,
        "fusion_pair_checked": fusion_pair_checked, "fusion_pair": fusion_pair,
        "whitelist_hit": whitelist_hit, "whitelist_matched_pair": whitelist_matched_pair,
        "blacklist_hit": blacklist_hit, "blacklist_matched_pair": blacklist_matched_pair,
    }


def format_region(sv: SVRow) -> str:
    same_chrom = normalize_chrom(sv.chrom1) == normalize_chrom(sv.chrom2)
    if same_chrom:
        low, high = sorted((sv.pos1, sv.pos2))
        return f"chr{normalize_chrom(sv.chrom1)}:{low:,}-{high:,}"
    return f"chr{normalize_chrom(sv.chrom1)}:{sv.pos1:,} :: chr{normalize_chrom(sv.chrom2)}:{sv.pos2:,}"


def build_explanation(sv: SVRow, ann: dict, n_total_cells: int) -> str:
    """Plain-English, per-SV justification for the top-N table -- every
    sentence here is derived directly from already-computed fields (no
    hidden logic), so the explanation is always consistent with the rest of
    that row's columns."""
    parts = []
    if sv.svtype == "BND":
        if ann["whitelist_hit"]:
            gene1, gene2 = ann["whitelist_matched_pair"]
            parts.append(f"Translocation breakpoint joins {gene1} and {gene2}, a gene pair listed in the "
                          "Mitelman Database of cancer-associated gene fusions (WHITELIST match).")
        elif ann["fusion_pair"]:
            gene1, gene2 = ann["fusion_pair"]
            parts.append(f"Translocation breakpoint joins {gene1} and {gene2}; this specific gene pair is in "
                          "neither the cancer-fusion whitelist nor the normal-tissue blacklist.")
        else:
            parts.append("Translocation breakpoint; at least one side does not overlap an annotated gene, "
                          "so no gene-pair fusion check could be performed.")
    elif sv.svtype in ("DEL", "DUP"):
        verb = "Deletion" if sv.svtype == "DEL" else "Duplication"
        if ann["genes_in_span"]:
            shown = ann["genes_in_span"][:8]
            more = f", and {len(ann['genes_in_span']) - 8} more" if len(ann["genes_in_span"]) > 8 else ""
            parts.append(f"{verb} spans {len(ann['genes_in_span'])} annotated gene(s): {', '.join(shown)}{more}.")
        else:
            parts.append(f"{verb} does not overlap any annotated gene.")
    else:
        parts.append(f"{sv.svtype} event; no whitelist/blacklist check applies to this SV type (see script "
                      "docstring for why only BND events get a gene-pair fusion check).")
    if ann["blacklist_hit"]:
        gene1, gene2 = ann["blacklist_matched_pair"]
        parts.append(f"CAUTION: {gene1}/{gene2} also matches a gene pair reported as a recurrent fusion "
                      "transcript in normal, non-cancer tissue (BLACKLIST) -- more likely a transcriptional "
                      "read-through/splicing artifact than a somatic driver event.")
    fraction = sv.supporting_cells / n_total_cells if n_total_cells else 0.0
    parts.append(f"Detected in {sv.supporting_cells} of {n_total_cells} cells ({fraction:.0%}), "
                  f"{sv.supporting_records} total supporting read/record pileup(s) cohort-wide.")
    return " ".join(parts)


def rank_svs(annotated: list[dict], top_n: int) -> list[dict]:
    """Priority order for the top-N table: (1) BND events whose breakpoint
    genes match the cancer-fusion whitelist, ranked by recurrence
    (supporting_cells, most cells first); then (2) every other non-
    blacklisted SV, also ranked by recurrence. Blacklisted SVs (known
    normal-tissue fusions) are EXCLUDED from this ranked/explained table
    entirely -- by definition they are not a prioritization candidate --
    but remain fully visible, flagged, in the full and main TSVs."""
    eligible = [a for a in annotated if not a["ann"]["blacklist_hit"]]
    whitelisted = sorted((a for a in eligible if a["ann"]["whitelist_hit"]),
                         key=lambda a: -a["sv"].supporting_cells)
    rest = sorted((a for a in eligible if not a["ann"]["whitelist_hit"]),
                 key=lambda a: -a["sv"].supporting_cells)
    return (whitelisted + rest)[:top_n]


def draw_sv_figure(sv: SVRow, ann: dict, supporting_cells: list[str], cell_order: list[str],
                    chrom_lengths: dict[str, int], rank: int, out_path: Path) -> None:
    """OncoScan-style per-SV figure: a linear chromosome schematic (NOT a
    true cytogenetically-banded ideogram -- a simplified proportional bar)
    showing the breakpoint(s)/span and any overlapping gene labels, plus a
    per-cell presence strip across the whole cohort (dark = this cell
    supports the SV) underneath, matching OncoScan's convention of showing
    sample-level support alongside the genomic event."""
    same_chrom = normalize_chrom(sv.chrom1) == normalize_chrom(sv.chrom2)
    chrom1 = normalize_chrom(sv.chrom1)
    chrom2 = normalize_chrom(sv.chrom2)
    len1 = chrom_lengths.get(chrom1, max(sv.pos1, sv.pos2) * 1.1)
    len2 = chrom_lengths.get(chrom2, max(sv.pos1, sv.pos2) * 1.1)

    fig, (ax_main, ax_cells) = plt.subplots(2, 1, figsize=(10, 4.2), gridspec_kw={"height_ratios": [3, 1]})

    if same_chrom and sv.svtype != "BND":
        box_bottom, box_top = 0.45, 0.65
        ax_main.add_patch(mpatches.FancyBboxPatch((0, box_bottom), 1, box_top - box_bottom,
                                                    boxstyle="round,pad=0,rounding_size=0.02",
                                                    linewidth=1, edgecolor="#333", facecolor="#dfe6e9"))
        x1, x2 = sv.pos1 / len1, sv.pos2 / len1
        low, high = sorted((x1, x2))
        if sv.svtype in ("DEL", "DUP"):
            color = "#e74c3c" if sv.svtype == "DEL" else "#2980b9"
            ax_main.add_patch(mpatches.Rectangle((low, box_bottom), max(high - low, 0.002), box_top - box_bottom,
                                                  facecolor=color, alpha=0.65))
        for x in (x1, x2):
            ax_main.plot([x, x], [box_bottom - 0.03, box_top + 0.03], color="black", lw=1.2)
        # DEL/DUP: the span label (below) already covers every gene the
        # breakpoint-only labels would show, and is more informative (the
        # breakpoints often sit within ~1-2kb of each other on a much
        # larger chromosome, so their two point-gene labels would otherwise
        # render on top of one another). Only fall back to breakpoint-level
        # labels when there's no span result to show (INS has no span
        # concept at all; a DEL/DUP could in principle also hit zero genes).
        if not (sv.svtype in ("DEL", "DUP") and ann["genes_in_span"]):
            point_genes = sorted(set(ann["genes_breakpoint1"]) | set(ann["genes_breakpoint2"]))
            if point_genes:
                ax_main.text((low + high) / 2, 0.97, ",".join(point_genes[:3]), ha="center", va="center",
                             fontsize=7)
        if ann["genes_in_span"]:
            label = ",".join(ann["genes_in_span"][:6])
            if len(ann["genes_in_span"]) > 6:
                label += f" (+{len(ann['genes_in_span']) - 6} more)"
            ax_main.text((low + high) / 2, 0.30, label, ha="center", va="center", fontsize=6.5, color="#555")
        ax_main.text(0.5, 0.14, f"chr{chrom1}", ha="center", fontsize=9)
        ax_main.set_xlim(-0.02, 1.02)
    else:
        # Two stacked chromosome bars (breakpoint1 on top, breakpoint2 on
        # bottom). Gene labels are placed just INSIDE the axes' top/bottom
        # edges (y=0.97 / y=0.03) rather than above/below the boxes by a
        # fixed offset -- two breakpoints that happen to be close together
        # (e.g. a same-chromosome BND, common in this data) map to nearly
        # identical x-positions, so a fixed-offset label can land outside
        # the axes' own [0,1] y-range entirely; matplotlib text is not
        # clipped there by default, so it silently overlaps the axes title
        # instead of raising any error. Anchoring to the axes' own edges
        # keeps every label inside the visible plot regardless of x.
        row0_box, row1_box = (0.80, 0.95), (0.05, 0.20)  # (bottom, top) of each chromosome bar
        row0_label_y, row1_label_y = 0.97, 0.03
        row0_chrom_label_y, row1_chrom_label_y = 0.875, 0.125
        # facing_edge: the box edge pointing INTO the gap between the two
        # bars (bottom of the top bar, top of the bottom bar) -- the arc
        # connects these two points so it visibly runs through that gap.
        rows = (
            (chrom1, sv.pos1, len1, ann["genes_breakpoint1"], row0_box, row0_label_y, row0_chrom_label_y, row0_box[0]),
            (chrom2, sv.pos2, len2, ann["genes_breakpoint2"], row1_box, row1_label_y, row1_chrom_label_y, row1_box[1]),
        )
        anchor_points = []
        for (chrom, pos, length, genes, (box_bottom, box_top), label_y, chrom_label_y, facing_edge) in rows:
            ax_main.add_patch(mpatches.FancyBboxPatch((0, box_bottom), 1, box_top - box_bottom,
                                                        boxstyle="round,pad=0,rounding_size=0.02",
                                                        linewidth=1, edgecolor="#333", facecolor="#dfe6e9"))
            x = pos / length
            ax_main.plot([x, x], [box_bottom - 0.03, box_top + 0.03], color="black", lw=1.2)
            if genes:
                ax_main.text(x, label_y, ",".join(genes[:3]), ha="center", va="center", fontsize=7)
            ax_main.text(1.03, chrom_label_y, f"chr{chrom}", fontsize=9, va="center")
            anchor_points.append((x, facing_edge))
        arc_color = "#c0392b" if ann["whitelist_hit"] else ("#7f8c8d" if ann["blacklist_hit"] else "#2980b9")
        arc = mpatches.FancyArrowPatch(anchor_points[0], anchor_points[1], connectionstyle="arc3,rad=0.3",
                                        arrowstyle="-", color=arc_color, lw=2)
        ax_main.add_patch(arc)
        ax_main.set_xlim(-0.02, 1.15)

    ax_main.set_ylim(0, 1)
    ax_main.axis("off")

    badges = []
    if ann["whitelist_hit"]:
        badges.append("WHITELIST: known cancer fusion")
    if ann["blacklist_hit"]:
        badges.append("BLACKLIST: normal-tissue fusion")
    badge_text = "  |  ".join(badges) if badges else "no whitelist/blacklist match"
    fraction = sv.supporting_cells / len(cell_order) if cell_order else 0.0
    ax_main.set_title(f"#{rank}  {sv.svtype}  {sv.variant_id}\n{badge_text}  --  "
                      f"{sv.supporting_cells}/{len(cell_order)} cells ({fraction:.0%})", fontsize=9)

    present = set(supporting_cells)
    colors = ["#2980b9" if cell in present else "#ecf0f1" for cell in cell_order]
    ax_cells.bar(range(len(cell_order)), [1] * len(cell_order), color=colors, width=1.0,
                 edgecolor="white", linewidth=0.3)
    ax_cells.set_xlim(-0.5, max(len(cell_order) - 0.5, 0.5))
    ax_cells.set_yticks([])
    ax_cells.set_xticks([])
    ax_cells.set_xlabel(f"{len(cell_order)} cells in cohort (dark blue = supports this SV)", fontsize=8)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_annotation_tsv(path: Path, rows: list[dict], n_total_cells: int) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["VariantId", "SVTYPE", "orientation", "region", "CHROM1", "POS1", "CHROM2", "POS2",
                          "supporting_records", "supporting_cells", "pct_cells",
                          "genes_breakpoint1", "genes_breakpoint2", "genes_in_span",
                          "fusion_pair_checked", "fusion_pair",
                          "whitelist_hit", "whitelist_matched_pair", "blacklist_hit", "blacklist_matched_pair"])
        for item in rows:
            sv, ann = item["sv"], item["ann"]
            fraction = sv.supporting_cells / n_total_cells if n_total_cells else 0.0
            writer.writerow([
                sv.variant_id, sv.svtype, sv.orientation, format_region(sv), sv.chrom1, sv.pos1, sv.chrom2, sv.pos2,
                sv.supporting_records, sv.supporting_cells, f"{fraction:.4f}",
                ",".join(ann["genes_breakpoint1"]), ",".join(ann["genes_breakpoint2"]),
                ",".join(ann["genes_in_span"]),
                ann["fusion_pair_checked"], "-".join(ann["fusion_pair"]) if ann["fusion_pair"] else "",
                ann["whitelist_hit"], "-".join(ann["whitelist_matched_pair"]) if ann["whitelist_matched_pair"] else "",
                ann["blacklist_hit"], "-".join(ann["blacklist_matched_pair"]) if ann["blacklist_matched_pair"] else "",
            ])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sv-matrix-dir", required=True, type=Path,
                    help="Directory containing sv_features.tsv and sv_binary.tsv "
                         "(step5_filter_variant_inputs.py's structural_variants output); no new SV calling.")
    p.add_argument("--manifest", required=True, type=Path,
                    help="passed_cells_only_manifest.tsv -- fixes the per-cell presence strip's cell order; "
                         "falls back to sv_binary.tsv's own column order if this can't be read.")
    p.add_argument("--genes-bed", type=Path, default=DEFAULT_GENES_BED)
    p.add_argument("--chrom-sizes", type=Path, required=True,
                    help="Reference .fai (or any chrom<TAB>length file) -- used only to scale the "
                         "schematic figures proportionally; not involved in any gene-overlap logic.")
    p.add_argument("--whitelist", action="append", default=None,
                    help="Gene-pair whitelist TSV (header: gene1, gene2[, source]); repeatable. "
                         "Defaults to the bundled Mitelman-Database-derived list (see refs/sv_gene_fusion/"
                         "README.md) when omitted.")
    p.add_argument("--blacklist", action="append", default=None,
                    help="Gene-pair blacklist TSV; repeatable. Defaults to the bundled Babiceanu/Greger "
                         "normal-tissue-fusion list when omitted.")
    p.add_argument("--custom-whitelist", action="append", default=[],
                    help="ADDITIONAL whitelist file(s), merged with (never replacing) the default(s) above.")
    p.add_argument("--custom-blacklist", action="append", default=[],
                    help="ADDITIONAL blacklist file(s), merged with (never replacing) the default(s) above.")
    p.add_argument("--min-supporting-cells", type=int, default=2,
                    help="Minimum supporting_cells for an SV to appear in the MAIN (filtered) table. "
                         "Singleton, single-cell-only SVs are, as with small variants elsewhere in this "
                         "PTA/WGA pipeline, far more likely to be amplification artifacts than real "
                         "recurrent events; they remain visible in the FULL, unfiltered TSV. Default: 2.")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--figures-top-n", type=int, default=20,
                    help="How many of the --top-n ranked SVs additionally get an OncoScan-style figure "
                         "(bounds dashboard size/runtime; every one of --top-n still gets a text "
                         "explanation regardless of this cap). Default: 20.")
    p.add_argument("--output-dir", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest_cells = []
    if args.manifest.is_file():
        with args.manifest.open() as handle:
            manifest_cells = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]

    gene_index = GeneIndex.from_bed(args.genes_bed)

    chrom_lengths = {}
    for line in args.chrom_sizes.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            chrom_lengths[normalize_chrom(parts[0])] = int(parts[1])

    whitelist_paths = [Path(p) for p in (args.whitelist or [DEFAULT_WHITELIST])] + \
        [Path(p) for p in args.custom_whitelist]
    blacklist_paths = [Path(p) for p in (args.blacklist or [DEFAULT_BLACKLIST])] + \
        [Path(p) for p in args.custom_blacklist]
    whitelist_pairs, _whitelist_sources = load_pair_set(whitelist_paths)
    blacklist_pairs, _blacklist_sources = load_pair_set(blacklist_paths)

    svs = load_sv_features(args.sv_matrix_dir / "sv_features.tsv")
    binary_cell_order, presence = load_sv_binary(args.sv_matrix_dir / "sv_binary.tsv")
    cell_order = manifest_cells or binary_cell_order
    n_total_cells = len(cell_order)

    annotated = [{"sv": sv, "ann": annotate_sv(sv, gene_index, whitelist_pairs, blacklist_pairs)} for sv in svs]

    write_annotation_tsv(args.output_dir / "sv_gene_annotation_full.tsv", annotated, n_total_cells)

    main_rows = [
        item for item in annotated
        if item["sv"].supporting_cells >= args.min_supporting_cells
        and (item["ann"]["genes_breakpoint1"] or item["ann"]["genes_breakpoint2"] or item["ann"]["genes_in_span"])
    ]
    main_rows.sort(key=lambda item: -item["sv"].supporting_cells)
    write_annotation_tsv(args.output_dir / "sv_gene_annotation_main.tsv", main_rows, n_total_cells)

    ranked = rank_svs(annotated, args.top_n)
    figures_dir = args.output_dir / "figures"
    top_path = args.output_dir / f"sv_gene_fusion_top{args.top_n}.tsv"
    with top_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "VariantId", "SVTYPE", "region", "supporting_cells", "n_cohort_cells",
                          "pct_cells", "whitelist_hit", "blacklist_hit", "explanation", "figure_png"])
        n_figures = 0
        for rank, item in enumerate(ranked, start=1):
            sv, ann = item["sv"], item["ann"]
            explanation = build_explanation(sv, ann, n_total_cells)
            figure_relpath = ""
            if rank <= args.figures_top_n:
                figure_relpath = f"figures/rank{rank:02d}_{safe_filename(sv.variant_id)}.png"
                draw_sv_figure(sv, ann, presence.get(sv.variant_id, []), cell_order, chrom_lengths, rank,
                                args.output_dir / figure_relpath)
                n_figures += 1
            fraction = sv.supporting_cells / n_total_cells if n_total_cells else 0.0
            writer.writerow([rank, sv.variant_id, sv.svtype, format_region(sv), sv.supporting_cells,
                              n_total_cells, f"{fraction:.1%}", ann["whitelist_hit"], ann["blacklist_hit"],
                              explanation, figure_relpath])

    provenance = {
        "stage": "SV gene-fusion annotation (Manta BND/DEL/DUP/INS breakpoints vs. gene models + "
                 "cancer-fusion whitelist / normal-tissue-fusion blacklist); no new SV calling performed.",
        "gene_model": str(args.genes_bed),
        "n_svs_total": len(svs),
        "n_cells": n_total_cells,
        "min_supporting_cells_for_main_table": args.min_supporting_cells,
        "n_svs_in_main_table": len(main_rows),
        "top_n": args.top_n,
        "figures_top_n": args.figures_top_n,
        "figures_generated": n_figures,
        "whitelist": {
            "default_files_used": [] if args.whitelist else [str(DEFAULT_WHITELIST)],
            "explicit_files_used": [str(p) for p in args.whitelist] if args.whitelist else [],
            "custom_files_supplied": [str(p) for p in args.custom_whitelist],
            "total_unique_pairs_loaded": len(whitelist_pairs),
        },
        "blacklist": {
            "default_files_used": [] if args.blacklist else [str(DEFAULT_BLACKLIST)],
            "explicit_files_used": [str(p) for p in args.blacklist] if args.blacklist else [],
            "custom_files_supplied": [str(p) for p in args.custom_blacklist],
            "total_unique_pairs_loaded": len(blacklist_pairs),
        },
        "n_bnd_events": sum(1 for item in annotated if item["sv"].svtype == "BND"),
        "n_bnd_fusion_pair_checked": sum(1 for item in annotated if item["ann"]["fusion_pair_checked"]),
        "n_whitelist_hits": sum(1 for item in annotated if item["ann"]["whitelist_hit"]),
        "n_blacklist_hits": sum(1 for item in annotated if item["ann"]["blacklist_hit"]),
        "caveat": "Whitelist and blacklist were both built from RNA-level fusion-transcript evidence, not "
                  "DNA-level SV breakpoints; applying them here is a standard but imperfect adaptation -- "
                  "treat hits as a manual-review prioritization aid, not a standalone pathogenicity call.",
    }
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"[OK] Annotated {len(svs)} SVs ({provenance['n_bnd_events']} BND) -> "
          f"{len(main_rows)} in main table, top {len(ranked)} ranked, {n_figures} figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
