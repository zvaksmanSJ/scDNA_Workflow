#!/usr/bin/env python3
"""Plan and validate per-cell merges of DNA lane/flowcell FASTQ pairs.

DNA sibling of ../../workflow_RNA/step0_prepare_raw_fastqs.py -- same
rigor (true FASTQ read-pairing validation, physical run/lane dedup via the
sequencer's own read-header run ID rather than trusting filename lane
labels, idempotent provenance-tracked merges, never modifies/moves any
source file), adapted for:
  - DNA's filename convention (hyphens, not underscores, and a "P<n>"
    plate/pool token): SJNBL046426_C894-ResolOME-SY5Y_P3-DNA_C17_S1_L001_R1_001.fastq.gz
  - multiple --raw-dir inputs (a cohort's DNA is frequently split across an
    original sequencing run PLUS one or more later "TopOff" re-sequencing
    runs for extra depth; both commonly reuse the same L001/L002 filename
    labels for a *different* physical flowcell/lane, which is exactly why
    this script dedups by the read header's own run+lane ID, not the
    filename's lane label)
  - writing a manifest directly in the 3-column (sample_id, r1, r2) format
    ../scripts/build_inputs_json.py expects, instead of RNA's 4-column
    nf-core-style biosampleName/read1/read2/groups CSV
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import json
from pathlib import Path
import re
import shutil

# ---- generic FASTQ-format helpers (duplicated, not imported, from
# workflow_RNA/step1_prepare_samples.py -- stable, workflow-agnostic, and
# keeping workflow_DNA self-contained avoids a cross-workflow-folder import) ----

def fastq_record(stream):
    header = stream.readline()
    if not header:
        return None
    sequence, plus, quality = (stream.readline().rstrip("\r\n") for _ in range(3))
    if not header.startswith("@") or not plus.startswith("+") or not sequence \
            or len(sequence) != len(quality):
        raise ValueError(f"Malformed FASTQ record: {header.strip()}")
    name = header.split()[0][1:]
    return re.sub(r"/[12]$", "", name)


def validate_pairs(r1, r2):
    count = 0
    with gzip.open(r1, "rt") as left, gzip.open(r2, "rt") as right:
        while True:
            a, b = fastq_record(left), fastq_record(right)
            if a is None and b is None:
                break
            if a is None or b is None or a != b:
                raise ValueError(f"Unpaired FASTQ at pair {count + 1}: {a!r}, {b!r}")
            count += 1
    if not count:
        raise ValueError("Empty FASTQ pair")
    return count


# ---- DNA-specific planning ----

OUTER = re.compile(r"(.+)_S(\d+)_(L\d{3})_(R[12])_(\d+)\.fastq\.gz")
# label = "SJNBL046426_C894-ResolOME-SY5Y_P3-DNA_C17" -> patient, library,
# plate/pool, local cell number
INNER = re.compile(r"(SJNBL\d+)_(C\d+)-ResolOME-SY5Y_P(\d+)-DNA_C(\d+)")


def plan(raw_dirs, outdir):
    raw_dirs = [Path(d).resolve() for d in raw_dirs]
    outdir = Path(outdir).resolve()
    groups = {}
    for raw_dir in raw_dirs:
        for path in sorted(raw_dir.glob("*DNA*.fastq.gz")):
            match = OUTER.fullmatch(path.name)
            if not match:
                raise ValueError(f"Unsupported DNA FASTQ name: {path.name}")
            label, sample, lane, mate, part = match.groups()
            identity = INNER.fullmatch(label)
            if not identity:
                raise ValueError(f"Unsupported DNA library identity: {label}")
            patient, library, plate, cell = identity.groups()
            sid = f"SY5Y_Cell{int(cell)}_DNA_{patient}_{library}_DNA_C{int(cell)}"
            key = (sid, raw_dir, int(sample), lane, part)
            entry = groups.setdefault(key, {})
            if mate in entry:
                raise ValueError(f"Duplicate mate for {key}: {path}")
            entry[mate] = path
    if not groups:
        raise ValueError(f"No DNA FASTQs in {raw_dirs}")

    cells, physical = {}, set()
    for (sid, raw_dir, sample, lane, part), mates in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
        if set(mates) != {"R1", "R2"}:
            raise ValueError(f"Missing mate: {sid}, {raw_dir}, S{sample}, {lane}, {part}")
        with gzip.open(mates["R1"], "rt") as left, gzip.open(mates["R2"], "rt") as right:
            name = fastq_record(left)
            if name is None or name != fastq_record(right):
                raise ValueError(f"First read is not paired: {mates}")
        fields = name.split(":")
        if len(fields) < 7 or not fields[3].isdigit() or int(fields[3]) != int(lane[1:]):
            raise ValueError(f"Illumina read header and filename lane disagree: {mates['R1']}")
        # Dedup key uses the read header's own run+lane ID (fields[:4],
        # e.g. "INSTRUMENT:RUN:FLOWCELL:LANE"), NOT the filename's lane
        # label -- a TopOff re-sequencing run reuses "L001"/"L002" for a
        # genuinely different physical flowcell/lane, which must NOT be
        # treated as a duplicate of the original run's L001/L002.
        run_lane = ":".join(fields[:4])
        key = (sid, run_lane, part)
        if key in physical:
            raise ValueError(f"Duplicate physical run/lane input: {key}")
        physical.add(key)
        cell = cells.setdefault(sid, {"sample_id": sid, "pairs": [],
                                     "r1": str(outdir / "fastqs" / f"{sid}_R1.fastq.gz"),
                                     "r2": str(outdir / "fastqs" / f"{sid}_R2.fastq.gz")})
        cell["pairs"].append({
            "source_dir": str(raw_dir), "sample_number": sample, "lane": lane,
            "run_lane": run_lane, "part": part,
            **{mate: {"path": str(path), "size": path.stat().st_size,
                      "mtime_ns": path.stat().st_mtime_ns} for mate, path in mates.items()}})
    return {"raw_directories": [str(d) for d in raw_dirs], "cells": list(cells.values())}


def write_plan(manifest, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(manifest, indent=2) + "\n"
    path = outdir / "raw_inputs.json"
    if path.exists():
        if path.read_text() != encoded:
            raise ValueError("Raw input plan changed; use a new output directory")
    else:
        with path.open("x") as stream:
            stream.write(encoded)
    # 3-column, header-optional TSV, exactly the format
    # ../scripts/build_inputs_json.py's --manifest expects.
    table = outdir / "dna_cells.tsv"
    if not table.exists():
        temporary = table.with_suffix(".tsv.tmp")
        with temporary.open("x") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["sample_id", "r1_fastq", "r2_fastq"])
            for cell in manifest["cells"]:
                writer.writerow([cell["sample_id"], cell["r1"], cell["r2"]])
        temporary.replace(table)
    else:
        with table.open() as stream:
            existing = list(csv.reader(stream, delimiter="\t"))[1:]
        expected = [[cell["sample_id"], cell["r1"], cell["r2"]] for cell in manifest["cells"]]
        if existing != expected:
            raise ValueError("Existing DNA cell manifest differs from raw plan")


def merge_cell(manifest_path, index):
    manifest_path = Path(manifest_path)
    cells = json.loads(manifest_path.read_text())["cells"]
    if index < 1 or index > len(cells):
        raise ValueError(f"Cell index outside 1..{len(cells)}: {index}")
    cell = cells[index - 1]
    root = Path(cell["r1"]).parent
    root.mkdir(parents=True, exist_ok=True)
    marker = root / f'{cell["sample_id"]}.merge.json'
    with (root / f'{cell["sample_id"]}.lock').open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for pair in cell["pairs"]:
            for mate in ("R1", "R2"):
                item = pair[mate]
                stat = Path(item["path"]).stat()
                if (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"]):
                    raise ValueError(f'Raw input changed: {item["path"]}')
        if marker.exists():
            saved = json.loads(marker.read_text())
            if saved["input"] != cell:
                raise ValueError(f"Existing merge provenance changed: {marker}")
            for mate in ("r1", "r2"):
                stat = Path(cell[mate]).stat()
                if [stat.st_size, stat.st_mtime_ns] != saved["outputs"][mate]:
                    raise ValueError(f"Existing merged input changed: {cell[mate]}")
            print(f'[OK] Validated merge retained: {cell["sample_id"]}', flush=True)
            return
        targets = [Path(cell[m]) for m in ("r1", "r2")]
        temporary = [p.with_suffix(p.suffix + ".partial") for p in targets]
        if any(p.exists() for p in targets + temporary):
            raise FileExistsError(f'Uncommitted merge exists: {cell["sample_id"]}; inspect before retry')
        counts = []
        try:
            with temporary[0].open("xb") as left, temporary[1].open("xb") as right:
                for pair in cell["pairs"]:
                    count = validate_pairs(pair["R1"]["path"], pair["R2"]["path"])
                    counts.append(count)
                    for mate, output in (("R1", left), ("R2", right)):
                        with Path(pair[mate]["path"]).open("rb") as stream:
                            shutil.copyfileobj(stream, output, length=8 * 1024 * 1024)
                    print(f'[OK] {cell["sample_id"]} {pair["source_dir"]} {pair["run_lane"]}: {count} pairs', flush=True)
            merged_count = validate_pairs(*temporary)
            if merged_count != sum(counts):
                raise ValueError("Merged FASTQ count does not equal source lane counts")
            for source, target in zip(temporary, targets):
                source.replace(target)
            saved = {"input": cell, "source_pair_counts": counts, "merged_pairs": merged_count,
                     "outputs": {m: [Path(cell[m]).stat().st_size, Path(cell[m]).stat().st_mtime_ns]
                                 for m in ("r1", "r2")}}
            with marker.open("x") as stream:
                stream.write(json.dumps(saved, indent=2) + "\n")
        finally:
            for path in temporary:
                if path.exists():
                    path.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--raw-dir", type=Path, action="append",
                        help="repeatable: one or more source directories of *DNA*.fastq.gz "
                             "(e.g. the original run AND a TopOff re-sequencing run)")
    modes.add_argument("--merge-index", type=int)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    if args.raw_dir:
        if not args.outdir:
            parser.error("--raw-dir requires --outdir")
        manifest = plan(args.raw_dir, args.outdir)
        write_plan(manifest, args.outdir)
        print(f'[OK] Planned {len(manifest["cells"])} DNA cells; no raw input modified')
    else:
        if not args.manifest:
            parser.error("--merge-index requires --manifest")
        merge_cell(args.manifest, args.merge_index)


if __name__ == "__main__":
    main()
