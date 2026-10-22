#!/usr/bin/env python3
"""CUSTOM library QC comparison; informational DNAQC cutoffs never replace cell calling.

Run after per-cell fastp, samtools markdup and Picard CollectMultipleMetrics.
The source's five-component consensus additionally needs preseq, interval-depth
chrM fraction and Ginkgo metrics; unavailable components remain UNKNOWN, not zero.

Source-equivalent producers (PREFIX = ROOT/qc/picard/SID/SID):
  picard CollectMultipleMetrics R=REF I=dedup.bam O=PREFIX PROGRAM=null
    PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics
    PROGRAM=CollectGcBiasMetrics PROGRAM=MeanQualityByCycle
    PROGRAM=CollectQualityYieldMetrics
  fastp --thread N --in1 R1 --in2 R2 --out1 TRIM_R1 --out2 TRIM_R2
    --json ROOT/qc/fastp/SID/SID.fastp.json --html FASTP_HTML --detect_adapter_for_pe
  samtools markdup -r -s -f PREFIX.markdup_stats.txt FIXMATE_SORTED_BAM dedup.bam
  bam2mr -seg_len 100000 -o SID.mr dedup.bam
  preseq gc_extrap -o SID.curve SID.mr
  tail -n 1 SID.curve | cut -f 2 | awk '{printf "%.0f\\n", $1}'
  samtools depth -b SOURCE_BASE_METRICS_INTERVALS dedup.bam
  Rscript SOURCE/container/custom_r_qcplots/scripts/cnvSummarizer.R
    --rds_file GINKGO_RDS --out_file CNV_SUMMARY

The optional --source-metrics TSV combines the final preseq estimate, depth-sum
chrM fraction and original CNV summary MAPD_CNV_Log2/SKEW_CNV columns, keyed by
sample_name. Missing/NA values require explanatory *_reason fields when supplied.
Missing required files/headers are errors, except a missing/empty insert-size file
with a matching PREFIX.insert_size.status.tsv containing exactly
"not_computable<TAB>no_eligible_pairs". Nonempty insert metrics must agree with
that marker (no positive READ_PAIRS). A failed preseq model fit needs a two-column
SID.preseq.status.tsv ("status", "reason") with one NOT_COMPUTABLE row and reason
"preseq model fit failed", plus a nonempty SID.preseq.log or SID.gc_extrap.log.
Successful numeric estimates conflict with a failure marker. Standard
FastQC reports under ROOT/qc/fastqc_raw and ROOT/qc/fastqc_trimmed are included in
MultiQC if present (qc/fastqc also supported). DNAQC's
PICARD_METRICS does not invoke CollectWgsMetrics; none is assumed here.
"""

import argparse
import base64
import csv
import hashlib
import html
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess
import sys


SOURCE_ROOT = Path(__file__).resolve().parent.parent / "resolveome_align/basej-public-pipelines/basej-dnaqc"  # resolveome_align/ lives at workflow_DNA root
DEFAULT_CONFIG = SOURCE_ROOT / "assets/plot_qc_config.json"
RULES = {
    "preseq_count": ("cutoff_preseq", ">"),
    "PCT_CHIMERAS": ("cutoff_chimeras", "<"),
    "chrM": ("cutoff_chrM", "<"),
    "MAPD_CNV_Log2": ("cutoff_cnv_mapd", "<"),
    "SKEW_CNV": ("cutoff_cnv_sk", "<"),
}
GATE_FIELDS = (
    "sample_id", "total_reads", "mapped_reads", "properly_paired",
    "map_rate", "proper_pair_rate", "pass_qc",
)
EXTRA_SOURCE_FIELDS = ("preseq_count", "chrM", "MAPD_CNV_Log2", "SKEW_CNV")


def number(value, context):
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{context}: expected a number, got {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{context}: non-finite numeric value")
    return result


def require_keys(mapping, names, context):
    missing = set(names) - mapping.keys()
    if missing:
        raise ValueError(f"{context}: missing fields {sorted(missing)}")


def read_tsv(path, required, id_field):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_keys(dict.fromkeys(reader.fieldnames or []), required, str(path))
        records = {}
        for row in reader:
            sid = row[id_field]
            if not sid or Path(sid).name != sid or sid in (".", "..") or sid in records:
                raise ValueError(f"{path}: invalid or duplicate sample identifier {sid!r}")
            if None in row or any(row[key] is None for key in required):
                raise ValueError(f"{path}: malformed TSV row for {sid}")
            records[sid] = row
    if not records:
        raise ValueError(f"{path}: no sample rows")
    return records


def picard_rows(path, header_key, required):
    """Read only the metrics block, preserving trailing empty metadata columns."""
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.split("\t")[0] != header_key:
            continue
        headers = line.split("\t")
        require_keys(dict.fromkeys(headers), required, str(path))
        if len(set(headers)) != len(headers):
            raise ValueError(f"{path}: duplicate metric columns")
        rows = []
        for data in lines[index + 1:]:
            if not data or data.startswith("#"):
                break
            values = data.split("\t")
            if len(values) != len(headers):
                raise ValueError(f"{path}: header/data field count mismatch")
            rows.append(dict(zip(headers, values)))
        return rows
    raise ValueError(f"{path}: missing expected Picard header {header_key}")


def select_row(rows, path, category=None):
    if category is not None:
        rows = [row for row in rows if row["CATEGORY"] == category]
    if len(rows) != 1:
        raise ValueError(f"{path}: expected one aggregate metrics row, got {len(rows)}")
    return rows[0]


def numeric_metrics(row, prefix):
    result = {}
    for key, value in row.items():
        if key not in {"CATEGORY", "SAMPLE", "LIBRARY", "READ_GROUP", "PAIR_ORIENTATION",
                       "ACCUMULATION_LEVEL", "READS_USED"}:
            if value in ("", "NaN", "?"):
                # Picard's CollectWgsMetrics uses the literal "?" sentinel
                # (distinct from "" or "NaN") specifically for
                # FOLD_80/90/95_BASE_PENALTY when the coverage distribution
                # doesn't reach that percentile threshold -- expected and
                # common for uneven/low-coverage single-cell PTA data, not
                # a tool or parsing error.
                result[prefix + key] = None
                result[prefix + key + "_reason"] = "Picard reported an undefined value; insufficient eligible observations"
            else:
                result[prefix + key] = number(value, prefix + key)
        else:
            result[prefix + key] = value
    return result


def ratio(numerator, denominator):
    return numerator / denominator if numerator is not None and denominator else None


def read_markdup(path):
    stats = {}
    for line in path.read_text().splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if key in stats:
            raise ValueError(f"{path}: duplicate statistic {key}")
        stats[key] = value.strip()
    required = ("EXAMINED", "DUPLICATE PAIR", "DUPLICATE SINGLE")
    require_keys(stats, required, str(path))
    result = {}
    for key in required:
        value = number(stats[key], f"{path}: {key}")
        if value < 0 or not value.is_integer():
            raise ValueError(f"{path}: {key} must be a nonnegative integer")
        result["markdup_" + key.lower().replace(" ", "_")] = value
    duplicates = result["markdup_duplicate_pair"] + result["markdup_duplicate_single"]
    examined = result["markdup_examined"]
    if duplicates > examined:
        raise ValueError(f"{path}: duplicate primary reads exceed examined reads")
    result["markdup_primary_duplicate_fraction"] = ratio(duplicates, examined)
    result["markdup_primary_duplicate_fraction_reason"] = (
        "" if examined else "No primary reads examined"
    )
    result["picard_estimated_library_size_reason"] = (
        "samtools markdup statistics are not Picard DuplicationMetrics; no interchangeable library-size estimate"
    )
    result["picard_estimated_library_size"] = None
    return result


def read_pre_dedup_flagstat(path):
    patterns = {
        "total_reads": r"^(\d+)\s+\+\s+(\d+)\s+in total\b",
        "mapped_reads": r"^(\d+)\s+\+\s+(\d+)\s+mapped\s+\(",
        "properly_paired": r"^(\d+)\s+\+\s+(\d+)\s+properly paired\s+\(",
    }
    counts = {}
    for line in path.read_text().splitlines():
        for key, pattern in patterns.items():
            match = re.match(pattern, line)
            if match:
                # Original gate parser uses the QC-passing column, not QC-fail sum.
                counts[key] = int(match.group(1))
    require_keys(counts, patterns, str(path))
    counts["map_rate"] = ratio(counts["mapped_reads"], counts["total_reads"])
    counts["proper_pair_rate"] = ratio(counts["properly_paired"], counts["total_reads"])
    return {"pre_dedup_" + key: value for key, value in counts.items()}


def read_fastp(path):
    data = json.loads(path.read_text())
    require_keys(data, ("command", "summary", "adapter_cutting", "duplication"), str(path))
    if not isinstance(data["command"], str) or not data["command"]:
        raise ValueError(f"{path}: missing fastp producer command needed by MultiQC")
    summary = data["summary"]
    require_keys(summary, ("before_filtering", "after_filtering"), str(path))
    before = summary["before_filtering"]
    before_keys = ("total_reads", "total_bases", "q20_rate", "q30_rate", "gc_content",
                   "read1_mean_length", "read2_mean_length")
    require_keys(before, before_keys, str(path))
    require_keys(summary["after_filtering"], ("total_reads", "q30_rate"), str(path))
    require_keys(data["adapter_cutting"], ("adapter_trimmed_reads", "adapter_trimmed_bases"), str(path))
    require_keys(data["duplication"], ("rate",), str(path))
    result = {"fastp_" + key: number(before[key], f"{path}: {key}") for key in before_keys}
    result["fastp_reads_after_filter"] = number(summary["after_filtering"]["total_reads"], str(path))
    result["fastp_q30_rate_after"] = number(summary["after_filtering"]["q30_rate"], str(path))
    for key in ("adapter_trimmed_reads", "adapter_trimmed_bases"):
        result["fastp_" + key] = number(data["adapter_cutting"][key], str(path))
    result["fastp_duplication_rate"] = number(data["duplication"]["rate"], str(path))
    return result


def source_scores(metrics, config):
    states = {}
    for metric, (cutoff, operator) in RULES.items():
        value = metrics[metric]
        if value is None:
            states[metric] = "UNKNOWN"
        else:
            passed = value > config[cutoff] if operator == ">" else value < config[cutoff]
            states[metric] = "PASS" if passed else "FAIL"
    known = sum(state != "UNKNOWN" for state in states.values())
    passed = sum(state == "PASS" for state in states.values())
    score = passed if known == len(RULES) else None
    status = "UNKNOWN" if score is None else "PASS" if score >= 4 else "Borderline" if score == 3 else "FAIL"
    return states, score, status, known, passed


def preseq_input(root, sid, supplied):
    directory = root / "qc/preseq" / sid
    estimate_path = directory / f"{sid}_preseq.txt"
    status_path = directory / f"{sid}.preseq.status.tsv"
    logs = [directory / f"{sid}.preseq.log", directory / f"{sid}.gc_extrap.log"]
    text = estimate_path.read_text().strip() if estimate_path.is_file() else ""
    numeric = bool(re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", text))
    estimate = number(text, str(estimate_path)) if numeric else None
    if estimate is not None and estimate < 0:
        estimate = None
    if status_path.is_file():
        with status_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["status", "reason"]:
                raise ValueError(f"{status_path}: expected exactly status and reason columns")
            rows = list(reader)
        if rows != [{"status": "NOT_COMPUTABLE", "reason": "preseq model fit failed"}]:
            raise ValueError(f"{status_path}: expected one NOT_COMPUTABLE / preseq model fit failed row")
        if estimate is not None or supplied.get("preseq_count") not in (None, "", "NA"):
            raise ValueError(f"{status_path}: failure marker conflicts with successful preseq estimate")
        if not any(path.is_file() and path.stat().st_size > 0 for path in logs):
            raise ValueError(f"{status_path}: NOT_COMPUTABLE requires nonempty raw preseq model log")
        return {"preseq_count": None, "preseq_count_reason": "preseq model fit failed",
                "preseq_status": "NOT_COMPUTABLE"}
    if "preseq_count" in supplied:
        if estimate_path.is_file() and estimate is None:
            raise ValueError(f"{estimate_path}: invalid estimate without NOT_COMPUTABLE status")
        return {}
    if estimate is not None:
        return {"preseq_count": estimate, "preseq_status": "COMPUTABLE"}
    if directory.is_dir():
        raise ValueError(f"{directory}: missing/invalid estimate requires NOT_COMPUTABLE status and raw logs")
    return {}


def optional_source_inputs(root, sid, supplied):
    """Use named source-compatible products, never infer coverage from read counts."""
    result = dict(supplied)
    result.update(preseq_input(root, sid, supplied))
    if "chrM" not in result:
        path = root / "qc/picard" / sid / f"{sid}.cov_sentieonmetrics.sample_interval_summary"
        if path.is_file():
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                require_keys(dict.fromkeys(reader.fieldnames or []), ("Target", "total_coverage"), str(path))
                total = mitochondrial = 0
                targets = set()
                for row in reader:
                    target = row["Target"]
                    if target in targets:
                        raise ValueError(f"{path}: duplicate interval target {target}")
                    targets.add(target)
                    depth = number(row["total_coverage"], f"{path}: {target}")
                    if depth < 0:
                        raise ValueError(f"{path}: negative depth sum")
                    total += depth
                    if target.split(":")[0] in {"chrM", "MT"}:
                        mitochondrial += depth
            if not targets:
                raise ValueError(f"{path}: no interval-depth rows")
            result["chrM"] = ratio(mitochondrial, total)
            if total == 0:
                result["chrM_reason"] = "No aligned bases within source coverage intervals"
    return result


def collect_cell(root, sid, gate, config, extra):
    extra = optional_source_inputs(root, sid, extra)
    row = {"sample_id": sid}
    row["preseq_status"] = extra.get("preseq_status", "UNKNOWN")
    for key in GATE_FIELDS[1:-1]:
        row[key] = number(gate[key], f"{sid}: flagstat {key}")
    if gate["pass_qc"].lower() not in {"true", "false"}:
        raise ValueError(f"{sid}: pass_qc must be true or false")
    row["flagstat_gate_pass_qc"] = gate["pass_qc"].lower()
    prefix = root / "qc/picard" / sid / sid
    pre_dedup_flagstat = Path(str(prefix) + ".pre_dedup.flagstat.txt")
    if pre_dedup_flagstat.is_file():
        row.update(read_pre_dedup_flagstat(pre_dedup_flagstat))
        row["pre_dedup_flagstat_reason"] = ""
    else:
        row["pre_dedup_flagstat_reason"] = "Optional preserved pre-dedup flagstat not supplied"
    alignment_path = Path(str(prefix) + ".alignment_summary_metrics")
    aligned_rows = picard_rows(alignment_path, "CATEGORY",
                              ("CATEGORY", "TOTAL_READS", "PF_READS", "PF_READS_ALIGNED",
                               "READS_ALIGNED_IN_PAIRS", "PCT_CHIMERAS"))
    category = "PAIR" if any(r["CATEGORY"] == "PAIR" for r in aligned_rows) else "UNPAIRED"
    alignment = select_row(aligned_rows, alignment_path, category)
    row.update(numeric_metrics(alignment, "alignment_"))
    total = number(alignment["TOTAL_READS"], str(alignment_path))
    aligned = number(alignment["PF_READS_ALIGNED"], str(alignment_path))
    paired = number(alignment["READS_ALIGNED_IN_PAIRS"], str(alignment_path))
    row["picard_pct_aligned"] = ratio(aligned, total)
    row["picard_pct_pf"] = ratio(number(alignment["PF_READS"], str(alignment_path)), total)
    row["PCT_CHIMERAS"] = (
        number(alignment["PCT_CHIMERAS"], str(alignment_path))
        if category == "PAIR" and paired > 0 and alignment["PCT_CHIMERAS"] not in ("", "NaN") else None
    )
    row["PCT_CHIMERAS_reason"] = (
        "" if row["PCT_CHIMERAS"] is not None else
        "No eligible paired alignments or Picard reported an undefined chimera fraction"
    )
    suites = (
        (".insert_size_metrics", "MEDIAN_INSERT_SIZE",
         ("MEDIAN_INSERT_SIZE", "MEAN_INSERT_SIZE", "STANDARD_DEVIATION", "READ_PAIRS"), "insert_", paired),
        (".gc_bias.summary_metrics", "ACCUMULATION_LEVEL",
         ("ACCUMULATION_LEVEL", "TOTAL_CLUSTERS", "ALIGNED_READS", "AT_DROPOUT", "GC_DROPOUT"), "gc_", aligned),
        (".quality_yield_metrics", "TOTAL_READS",
         ("TOTAL_READS", "TOTAL_BASES", "Q20_BASES", "Q30_BASES"), "yield_", total),
    )
    insert_status_path = Path(str(prefix) + ".insert_size.status.tsv")
    insert_no_pairs = False
    if insert_status_path.is_file():
        status = insert_status_path.read_text().strip()
        if status != "not_computable\tno_eligible_pairs":
            raise ValueError(f"{insert_status_path}: expected not_computable TAB no_eligible_pairs")
        insert_no_pairs = True
    for suffix, header, required, namespace, eligible in suites:
        path = Path(str(prefix) + suffix)
        if namespace == "insert_" and (not path.is_file() or path.stat().st_size == 0):
            if not insert_no_pairs:
                raise ValueError(f"{path}: missing/empty collector output requires explicit no_eligible_pairs status")
            for key in required:
                row[namespace + key] = None
            row[namespace + "reason"] = "no_eligible_pairs (explicit matching collector status)"
            continue
        metric_rows = picard_rows(path, header, required)
        if namespace == "insert_" and insert_no_pairs:
            if any(number(item["READ_PAIRS"], str(path)) > 0 for item in metric_rows):
                raise ValueError(f"{path}: positive READ_PAIRS conflicts with no_eligible_pairs status")
        if not metric_rows and (eligible == 0 or namespace == "insert_"):
            for key in required:
                row[namespace + key] = None
            row[namespace + "reason"] = "No eligible reads/pairs for this collector; Picard emitted a header-only metrics block"
            continue
        if namespace == "insert_" and len(metric_rows) > 1:
            metric_rows = [r for r in metric_rows if r.get("PAIR_ORIENTATION") == "FR"]
        if namespace == "gc_":
            metric_rows = [r for r in metric_rows if r["ACCUMULATION_LEVEL"] == "All Reads"]
        selected = select_row(metric_rows, path)
        row.update(numeric_metrics(selected, namespace))
        row[namespace + "reason"] = (
            "Picard reported blank values; insufficient eligible observations"
            if any(row.get(namespace + key) is None for key in required) else ""
        )
    row["picard_q30_rate"] = (
        ratio(row["yield_Q30_BASES"], row["yield_TOTAL_BASES"])
        if row["yield_TOTAL_BASES"] is not None else None
    )
    row["picard_q30_rate_reason"] = "" if row["picard_q30_rate"] is not None else "No sequenced bases in QualityYield"
    wgs_path = Path(str(prefix) + ".wgs_metrics")
    if wgs_path.is_file():
        wgs = picard_rows(wgs_path, "GENOME_TERRITORY", ("GENOME_TERRITORY", "MEAN_COVERAGE"))
        row.update(numeric_metrics(select_row(wgs, wgs_path), "wgs_"))
    row.update(read_markdup(Path(str(prefix) + ".markdup_stats.txt")))
    row.update(read_fastp(root / "qc/fastp" / sid / f"{sid}.fastp.json"))
    for key in EXTRA_SOURCE_FIELDS:
        value = extra.get(key)
        if value in (None, "", "NA"):
            row[key] = None
            if key in extra and not extra.get(key + "_reason"):
                raise ValueError(f"{sid}: explicitly unavailable {key} needs {key}_reason")
            row[key + "_reason"] = extra.get(key + "_reason") or f"Upstream {key} not supplied; not derivable from library QC"
        else:
            row[key] = number(value, f"{sid}: {key}")
            row[key + "_reason"] = ""
    states, score, status, known, passed = source_scores(row, config)
    row.update({key + "_source_check": value for key, value in states.items()})
    row.update(source_composite_score=score, source_qc_informational_status=status,
               source_components_known=known, source_components_passed=passed)
    return row


def write_tsv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row[key] for key in columns})


def plot_embedded(rows, metric, cutoff=None, operator=None):
    known = [(i, row[metric]) for i, row in enumerate(rows) if row.get(metric) is not None]
    values = [value for _, value in known]
    if cutoff is not None:
        values.append(cutoff)
    lower = min(0, min(values)) if values else 0
    upper = max(values) if values else 1
    if upper <= lower:
        upper = lower + 1
    upper += (upper - lower) * .05

    def y_position(value):
        return 240 - 180 * (value - lower) / (upper - lower)

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="700" height="300" viewBox="0 0 700 300">',
        '<rect width="700" height="300" fill="white"/>',
        f'<text x="350" y="22" text-anchor="middle" font-family="sans-serif">{html.escape(metric)}</text>',
        '<path d="M90 60V240H680" fill="none" stroke="#333"/>',
        '<text x="350" y="286" text-anchor="middle" font-family="sans-serif" font-size="12">Cell index (table order, zero-based)</text>',
    ]
    for fraction in (0, .5, 1):
        value = lower + fraction * (upper - lower)
        position = y_position(value)
        parts.append(f'<text x="84" y="{position + 4:.2f}" text-anchor="end" font-size="11">{value:.3g}</text>')
    for index, value in known:
        x = 105 + 555 * index / max(1, len(rows) - 1)
        parts.append(f'<circle cx="{x:.2f}" cy="{y_position(value):.2f}" r="3" fill="#2475aa">'
                     f'<title>{html.escape(rows[index]["sample_id"])}: {value:g}</title></circle>')
    for index in sorted({0, len(rows) - 1}):
        x = 105 + 555 * index / max(1, len(rows) - 1)
        parts.append(f'<text x="{x:.2f}" y="260" text-anchor="middle" font-size="11">{index}</text>')
    if cutoff is not None:
        y = y_position(cutoff)
        parts.append(f'<path d="M90 {y:.2f}H680" stroke="red" stroke-dasharray="5 4"/>')
        label = html.escape(f"Source pass: {operator} {cutoff:g}")
        parts.append(f'<text x="100" y="45" font-size="12" fill="red">{label}</text>')
    if not known:
        parts.append('<text x="350" y="150" text-anchor="middle" font-size="13">'
                     'UNKNOWN — upstream metric not available</text>')
    parts.append("</svg>")
    encoded = base64.b64encode("".join(parts).encode()).decode("ascii")
    return f'<img alt="{html.escape(metric)}" src="data:image/svg+xml;base64,{encoded}">'


def render_dashboard(path, rows, config, provenance):
    metrics = list(RULES) + ["total_reads", "map_rate", "proper_pair_rate",
                            "fastp_q30_rate", "picard_q30_rate", "insert_MEDIAN_INSERT_SIZE",
                            "gc_AT_DROPOUT", "gc_GC_DROPOUT", "markdup_primary_duplicate_fraction"]
    figures = []
    for metric in metrics:
        if metric in RULES:
            cutoff, operator = RULES[metric]
            figures.append(plot_embedded(rows, metric, config[cutoff], operator))
        else:
            figures.append(plot_embedded(rows, metric))
    columns = ["sample_id", "flagstat_gate_pass_qc", "source_qc_informational_status",
               "source_composite_score", "source_components_known"]
    for metric in RULES:
        columns.extend([metric, metric + "_source_check", metric + "_reason"])
    table = "<tr>" + "".join(f"<th>{html.escape(key)}</th>" for key in columns) + "</tr>"
    for row in rows:
        table += "<tr>" + "".join(
            f"<td>{html.escape(str(row.get(key) if row.get(key) is not None else 'NA'))}</td>"
            for key in columns) + "</tr>"
    path.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>CUSTOM library QC</title>"
        "<style>body{font:14px sans-serif;margin:2em}img{max-width:49%}"
        "td,th{border:1px solid #ccc;padding:4px}table{border-collapse:collapse}"
        ".scroll{overflow:auto}</style></head><body><h1>CUSTOM library QC comparison</h1>"
        "<p>The flagstat read/mapping/proper-pair gate is retained verbatim from the cell-QC "
        "manifest. If its producer now uses a deduplicated BAM, read counts differ from the "
        "old pre-dedup gate. Preserved pre-dedup flagstat metrics, when available, are exported "
        "separately in the metrics TSV, not used to recalculate cell calls. Source checks are "
        "<b>informational</b>, not a replacement cell call. Source consensus sums five strict "
        "boolean comparisons (score ≥4 PASS, 3 Borderline, ≤2 FAIL). If any component is "
        "unavailable the score/status is NA/UNKNOWN, not an extrapolated partial score.</p>"
        "<p>chrM means mitochondrial summed interval depth / total interval depth (chrM or MT), "
        "not a mapped-read fraction. preseq_count is the final extrapolation-curve estimate, "
        "and exact source conversion uses preseq 2.0.3 bam2mr (-seg_len 100000); preseq 3.x "
        "to-mr is not invocation-compatible. A missing source-compatible estimate stays UNKNOWN, "
        "not samtools library size. Ginkgo MAPD/skew must use source CNV definitions. "
        "samtools duplicate fraction uses (DUPLICATE PAIR + DUPLICATE SINGLE)/EXAMINED "
        "(primary reads); it is not Picard PERCENT_DUPLICATION. Fastp is input-read QC; "
        "Picard here is post-dedup BAM QC. Denominators/stages differ; do not interpret "
        "their difference as sample failure. This is a custom comparison, not the original "
        "source CNV/clustering dashboard. FastQC is aggregated by MultiQC, not inferred from fastp.</p>"
        + "".join(figures) + "<h2>Cell checks and unavailable-metric reasons</h2><div class='scroll'><table>"
        + table + "</table></div><h2>Provenance</h2><pre>"
        + html.escape(json.dumps(provenance, indent=2)) + "</pre></body></html>"
    )


def run(outdir, ref, config_path=DEFAULT_CONFIG, threads=1, source_metrics=None, multiqc=None):
    root = Path(outdir).resolve()
    ref = Path(ref).resolve()
    if not ref.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {ref}")
    if threads < 1:
        raise ValueError("--threads must be positive")
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    require_keys(config, [rule[0] for rule in RULES.values()], str(config_path))
    for cutoff, _ in RULES.values():
        config[cutoff] = number(config[cutoff], f"{config_path}: {cutoff}")
    gate_path = root / "manifests/cell_qc_metrics.tsv"
    gates = read_tsv(gate_path, GATE_FIELDS, "sample_id")
    if source_metrics is None:
        for directory in ("ginkgo", "ginkgo_cnv"):
            default_metrics = root / directory / "source_metrics.tsv"
            if default_metrics.is_file():
                source_metrics = default_metrics
                break
    extras = read_tsv(Path(source_metrics), ("sample_name",), "sample_name") if source_metrics else {}
    if set(extras) - set(gates):
        raise ValueError(f"Source metrics contain unknown cells: {sorted(set(extras) - set(gates))}")
    rows = [collect_cell(root, sid, gates[sid], config, extras.get(sid, {})) for sid in sorted(gates)]
    if multiqc:
        multiqc_command = [multiqc]
    elif importlib.util.find_spec("multiqc"):
        multiqc_command = [sys.executable, "-m", "multiqc"]
    else:
        sibling = Path(sys.executable).parent / "multiqc"
        if not sibling.is_file():
            raise FileNotFoundError("MultiQC is missing from the selected Python environment; use --multiqc explicitly")
        multiqc_command = [str(sibling)]
    report_dir = root / "qc/library_dashboard"
    report_dir.mkdir(parents=True, exist_ok=True)
    provenance = {
        "stage": "CUSTOM source-threshold annotated library QC comparison",
        "reference": str(ref), "source_config": str(config_path),
        "source_config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "source_definitions": "basej-dnaqc/main.nf PICARD_METRICS and QC_PLOTS; "
                              "container/custom_r_qcplots/scripts/function_plot_qc_dna.R",
        "source_metrics": str(source_metrics) if source_metrics else None,
        "coverage_interval_scope": "Uses supplied depth-sum summary; producer canonical-FAI BED "
                                   "may differ from source pipeline's unavailable/empty interval BED",
        "figure_equivalence": "Custom embedded SVG plots with source thresholds, not original R/image figures",
        "threads": threads, "cutoffs": {key: config[key] for key, _ in RULES.values()},
    }
    write_tsv(report_dir / "library_qc_metrics.tsv", rows)
    checks = [
        {"check": "cells_total", "count": len(rows)},
        {"check": "flagstat_gate_pass", "count": sum(r["flagstat_gate_pass_qc"] == "true" for r in rows)},
    ]
    checks.extend({"check": "source_informational_" + status, "count": sum(
        r["source_qc_informational_status"] == status for r in rows
    )} for status in ("PASS", "Borderline", "FAIL", "UNKNOWN"))
    checks.extend({"check": key + "_" + state, "count": sum(
        r[key + "_source_check"] == state for r in rows
    )} for key in RULES for state in ("PASS", "FAIL", "UNKNOWN"))
    write_tsv(report_dir / "cohort_checks.tsv", checks)
    scan_dirs = [root / "qc/picard", root / "qc/fastp"]
    fastqc_dirs = [root / "qc" / name for name in ("fastqc_raw", "fastqc_trimmed", "fastqc")]
    present_fastqc = [path for path in fastqc_dirs if path.is_dir()]
    scan_dirs.extend(present_fastqc)
    provenance["fastqc"] = (
        [str(path) for path in present_fastqc] if present_fastqc
        else "Not present at qc/fastqc_raw, qc/fastqc_trimmed or qc/fastqc; MultiQC cannot supply FastQC results"
    )
    command = [*multiqc_command, *(str(p) for p in scan_dirs), "--force",
               "--outdir", str(report_dir), "--filename", "multiqc_library_qc.html"]
    provenance["multiqc_command"] = command
    render_dashboard(report_dir / "source_threshold_dashboard.html", rows, config, provenance)
    with (report_dir / "multiqc.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not (report_dir / "multiqc_library_qc.html").is_file():
        raise RuntimeError("MultiQC returned success without producing its combined HTML report")
    (report_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return report_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--source-asset", type=Path, default=DEFAULT_CONFIG,
                        help="DNAQC plot_qc_config.json (verified local asset by default)")
    parser.add_argument("--source-metrics", type=Path,
                        help="Optional TSV: sample_name, preseq_count, chrM, MAPD_CNV_Log2, SKEW_CNV; NA allowed")
    parser.add_argument("--threads", type=int, default=1,
                        help="Recorded producer resource setting; MultiQC/report rendering are serial")
    parser.add_argument("--multiqc", help="Existing MultiQC executable (default: selected Python's MultiQC module)")
    args = parser.parse_args()
    report = run(args.outdir, args.ref, args.source_asset, args.threads, args.source_metrics, args.multiqc)
    print(f"Library QC reports: {report}")


if __name__ == "__main__":
    main()
