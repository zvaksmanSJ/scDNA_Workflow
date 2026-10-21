#!/usr/bin/env python3
"""Explain, per cell, exactly HOW Control-FREEC determined its final ploidy
call -- not just the resulting number (already in ploidy_summary.tsv).

Control-FREEC's own stdout (captured in the candidate-ploidy run's LSF job
logs) already reports a per-candidate-ploidy RSS (residual sum of squares)
fit score and the percentage of observed copy-number changes explained by
each candidate, followed by its own stated reasoning for the final choice
(occasionally overriding its own best-RSS pick when a lower, similarly-
fitting ploidy explains the data almost as well with fewer unexplained
regions -- FREEC's own documented tie-breaking heuristic, not something
this script invents). This script only parses and reports that existing
output; no new ploidy-selection logic is implemented here.

Two outputs:
  - ploidy_candidates.tsv (long format): one row per sample x candidate
    ploidy, with that candidate's RSS score and percent-genome-explained.
  - ploidy_narrative.tsv: one row per sample with the plain-English
    sequence of FREEC's own decision lines and the resulting final ploidy
    (should match ploidy_summary.tsv's output_ploidy for the same run).
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

TABLE_HEADER_RE = re.compile(r"^Ploidy\s+RSS score\s+Percentage of Genome Explained\s*$")
TABLE_ROW_RE = re.compile(r"^(\d+)\t([\d.eE+-]+)\t([\d.eE+-]+)\s*$")
NARRATIVE_LINE_RES = [
    re.compile(r"^\.\.Best ploidy"),
    re.compile(r"^\.\.Changed ploidy"),
    re.compile(r"^unexplained regions for ploidy"),
    re.compile(r"^\.\.Running FREEC with ploidy set to\s+(\d+)"),
]
FINAL_PLOIDY_RE = re.compile(r"\.\.Running FREEC with ploidy set to\s+(\d+)")
BEST_PLOIDY_RE = re.compile(r"\.\.Best ploidy set to\s+(\d+)")


def pick_log(logs_dir: Path, sample: str) -> Path | None:
    """Multiple log files can exist per sample (e.g. a killed/resubmitted
    job retains its earlier, truncated log alongside the real completed
    one) -- the completed run's log is always larger, since a killed job's
    log is cut off mid-run. Picking the largest file is a simple, robust
    way to always select the real, complete log without depending on which
    specific job ID happened to succeed."""
    candidates = list(logs_dir.glob(f"freec_{sample}.*.out"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_size)


def parse_log(log_path: Path) -> tuple[list[dict], list[str], int | None]:
    lines = log_path.read_text(errors="replace").splitlines()
    candidates = []
    narrative = []
    final_ploidy = None
    best_ploidy = None
    in_table = False
    narrative_started = False
    for line in lines:
        if TABLE_HEADER_RE.match(line):
            in_table = True
            continue
        if in_table:
            m = TABLE_ROW_RE.match(line)
            if m:
                candidates.append({
                    "candidate_ploidy": int(m.group(1)),
                    "rss_score": float(m.group(2)),
                    "pct_genome_explained": float(m.group(3)),
                })
                continue
            in_table = False
            narrative_started = True
        if narrative_started:
            if any(r.match(line) for r in NARRATIVE_LINE_RES):
                narrative.append(line.lstrip("."))
                m = FINAL_PLOIDY_RE.search(line)
                if m:
                    final_ploidy = int(m.group(1))
                m = BEST_PLOIDY_RE.search(line)
                if m and best_ploidy is None:
                    best_ploidy = int(m.group(1))
                continue
            if candidates:
                # First non-matching line after the table ends the narrative block.
                break
    if final_ploidy is None:
        final_ploidy = best_ploidy
    return candidates, narrative, final_ploidy


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freec-logs-dir", required=True, type=Path,
                    help="--logs-dir passed to loh_freec_lsf.py for the CANDIDATE-ploidy run "
                         "(contains freec_<sample>.<jobid>.out per cell)")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--output-candidates", required=True, type=Path)
    p.add_argument("--output-narrative", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]
    if not samples:
        raise SystemExit("No samples found in manifest")

    candidate_rows = []
    narrative_rows = []
    n_ok = 0
    for sample in samples:
        log_path = pick_log(args.freec_logs_dir, sample)
        if log_path is None:
            narrative_rows.append({"sample_id": sample, "final_ploidy": "", "narrative": "",
                                   "status": "MISSING_LOG"})
            continue
        candidates, narrative, final_ploidy = parse_log(log_path)
        if not candidates:
            narrative_rows.append({"sample_id": sample, "final_ploidy": "", "narrative": "",
                                   "status": "UNPARSEABLE_LOG"})
            continue
        for row in candidates:
            candidate_rows.append({"sample_id": sample, **row})
        narrative_rows.append({
            "sample_id": sample,
            "final_ploidy": final_ploidy if final_ploidy is not None else "",
            "narrative": " ".join(narrative),
            "status": "OK",
        })
        n_ok += 1

    args.output_candidates.parent.mkdir(parents=True, exist_ok=True)
    with args.output_candidates.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "candidate_ploidy", "rss_score",
                                                "pct_genome_explained"], delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    args.output_narrative.parent.mkdir(parents=True, exist_ok=True)
    with args.output_narrative.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "final_ploidy", "narrative", "status"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(narrative_rows)

    print(f"[OK] Ploidy-determination narrative parsed for {n_ok}/{len(samples)} samples")
    print(f"[OK] Candidate-ploidy fit table -> {args.output_candidates}")
    print(f"[OK] Narrative -> {args.output_narrative}")


if __name__ == "__main__":
    main()
