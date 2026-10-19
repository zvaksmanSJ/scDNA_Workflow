#!/usr/bin/env python3
"""Build a complete main.wdl inputs.json from a simple FASTQ manifest.

Takes (1) a manifest TSV with one row per cell -- sample_id, r1_fastq,
r2_fastq -- and (2) a "base" JSON file already containing every fixed
infrastructure input (reference, containers, Beagle panel, extra_methods
script paths, etc. -- see wdl/example_inputs.json for a starting point you
customize ONCE per deployment), and merges the manifest's FASTQ list into
the base, producing a ready-to-run inputs.json. This avoids hand-editing
JSON arrays for every cell in a real cohort.

Manifest TSV format (tab-separated, no header required but recommended):
    sample_id	r1_fastq	r2_fastq
    cell1	/path/to/cell1_R1.fastq.gz	/path/to/cell1_R2.fastq.gz
    cell2	/path/to/cell2_R1.fastq.gz	/path/to/cell2_R2.fastq.gz
    ...

Usage:
    python3 build_inputs_json.py \\
      --manifest cells.tsv \\
      --base wdl/example_inputs.json \\
      --donor-id my_cohort \\
      --output inputs.json

Every non-FASTQ key already present in --base is copied through unchanged;
any "_comment*" keys (used for documentation in example_inputs.json) are
stripped. --donor-id overrides the base file's donor_id if given.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

WORKFLOW_NAME = "ResolveOmeDnaWorkflow"


def read_manifest(path: Path) -> tuple[list[str], list[str], list[str]]:
    sample_ids, r1s, r2s = [], [], []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 3:
                raise SystemExit(f"[ERROR] manifest row has fewer than 3 columns: {row!r}")
            sample_id, r1, r2 = row[0].strip(), row[1].strip(), row[2].strip()
            if sample_id.lower() in ("sample_id", "sample", "id") and r1.lower() in ("r1_fastq", "r1", "fastq_r1"):
                continue  # skip an optional header row
            for p, label in ((r1, "R1"), (r2, "R2")):
                if not Path(p).is_file():
                    print(f"[WARN] {label} for {sample_id} does not exist (yet?): {p}", file=sys.stderr)
            sample_ids.append(sample_id)
            r1s.append(r1)
            r2s.append(r2)
    if not sample_ids:
        raise SystemExit(f"[ERROR] no cells found in manifest: {path}")
    dupes = {s for s in sample_ids if sample_ids.count(s) > 1}
    if dupes:
        raise SystemExit(f"[ERROR] duplicate sample_id(s) in manifest: {sorted(dupes)}")
    return sample_ids, r1s, r2s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, type=Path, help="TSV: sample_id, r1_fastq, r2_fastq")
    ap.add_argument("--base", required=True, type=Path, help="base inputs.json with every fixed/infra input filled in")
    ap.add_argument("--donor-id", default=None, help="override the base file's donor_id")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    base = json.loads(args.base.read_text())
    base = {k: v for k, v in base.items() if not k.startswith("_comment")}

    sample_ids, r1s, r2s = read_manifest(args.manifest)
    base[f"{WORKFLOW_NAME}.sample_ids"] = sample_ids
    base[f"{WORKFLOW_NAME}.r1_fastqs"] = r1s
    base[f"{WORKFLOW_NAME}.r2_fastqs"] = r2s
    if args.donor_id:
        base[f"{WORKFLOW_NAME}.donor_id"] = args.donor_id

    args.output.write_text(json.dumps(base, indent=2) + "\n")
    print(f"[DONE] wrote {args.output} with {len(sample_ids)} cell(s): {', '.join(sample_ids)}")


if __name__ == "__main__":
    main()
