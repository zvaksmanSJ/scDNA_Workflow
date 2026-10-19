#!/usr/bin/env python3
"""Plan fixed-size genomic chunks for per-region DeepVariant scatter.

Ported 1:1 from the chunking algorithm in
step3_submit_variants_cnv_sv_annot_lsf.py (make_chunks / read_fai_contigs):
each requested chromosome is split into consecutive windows of
--chunk-size-bp (default 20,000,000), covering the whole chromosome length
exactly once (last window may be shorter). MT is intentionally excluded
from "all" (see note below) -- same reasoning as the original: the
DeepVariant allele-frequency population VCF bundle has no MT/chrM entry,
and requesting the AF channel for a contig with no population reference
produced non-deterministic make_examples/call_variants tensor shapes.

Output: a tab-separated file (no header) with columns:
  fai_contig_name  canonical_chrom  start  end  region  chunk_id
suitable for Cromwell/miniwdl's read_tsv() to scatter over directly.
"""
from __future__ import annotations

import argparse
from pathlib import Path

CANONICAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def normalize_chrom_token(token: str) -> str:
    value = token.strip().upper()
    if value.startswith("CHR"):
        value = value[3:]
    if value == "M":
        value = "MT"
    return value


def chrom_matches_fai(contig_name: str, wanted: str) -> bool:
    contig = contig_name.upper()
    target = wanted.upper()
    if contig == target:
        return True
    if contig.startswith("CHR") and contig[3:] == target:
        return True
    if target == "MT" and contig in {"M", "CHRM", "CHRMT"}:
        return True
    return False


def parse_chroms(value: str) -> list[str]:
    if value.strip().lower() == "all":
        return [c for c in CANONICAL_CHROMS if c != "MT"]
    chromosomes = [normalize_chrom_token(item) for item in value.split(",") if item.strip()]
    invalid = [c for c in chromosomes if c not in CANONICAL_CHROMS]
    if invalid:
        raise SystemExit(f"[ERROR] Invalid chromosomes: {invalid}. Allowed: 1-22,X,Y,MT or all.")
    return chromosomes


def read_fai_contigs(fai_path: Path, chroms_arg: str):
    all_contigs = []
    with fai_path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                all_contigs.append((fields[0], int(fields[1])))
    requested = parse_chroms(chroms_arg)
    selected = []
    for chromosome in requested:
        match = None
        for contig_name, length in all_contigs:
            if chrom_matches_fai(contig_name, chromosome):
                match = (contig_name, length, chromosome)
                break
        if match is None:
            raise SystemExit(f"[ERROR] Chromosome {chromosome} was not found in {fai_path}")
        selected.append(match)
    return selected


def make_chunks(contigs, chunk_size: int):
    chunks = []
    for fai_name, length, canonical in contigs:
        start = 1
        while start <= length:
            end = min(length, start + chunk_size - 1)
            region = f"{fai_name}:{start}-{end}"
            chunk_id = f"{fai_name}_{start}_{end}"
            chunks.append((fai_name, canonical, start, end, region, chunk_id))
            start = end + 1
    return chunks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fai", required=True, type=Path)
    p.add_argument("--chroms", default="all")
    p.add_argument("--chunk-size-bp", type=int, default=20_000_000)
    p.add_argument("--output-tsv", required=True, type=Path)
    args = p.parse_args()

    contigs = read_fai_contigs(args.fai, args.chroms)
    chunks = make_chunks(contigs, args.chunk_size_bp)

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as fh:
        for fai_name, canonical, start, end, region, chunk_id in chunks:
            fh.write(f"{fai_name}\t{canonical}\t{start}\t{end}\t{region}\t{chunk_id}\n")
    print(f"[OK] {len(chunks)} chunks -> {args.output_tsv}")


if __name__ == "__main__":
    main()
