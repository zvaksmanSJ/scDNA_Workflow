#!/usr/bin/env python3
"""Submit exon/intron/gene coverage QC jobs for ResolveOME BAMs.

Coverage is computed with samtools depth -a, so zero-depth bases inside an
annotated interval are retained.  The report distinguishes exact zero coverage
from intervals with no base above the requested threshold (default: 3X).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shlex
import subprocess
import sys
import textwrap
from collections import defaultdict
from pathlib import Path


def q(value: str) -> str:
    return shlex.quote(str(value))


def normalize_contig(chrom: str, strip_chr_prefix: bool) -> str:
    if not strip_chr_prefix or not chrom.startswith("chr"):
        return chrom
    chrom = chrom[3:]
    return "MT" if chrom == "M" else chrom


def read_gtf(path: str, strip_chr_prefix: bool = False):
    opener = gzip.open if path.endswith(".gz") else open
    transcripts = defaultdict(list)
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = {}
            for item in fields[8].rstrip(";").split(";"):
                item = item.strip()
                if not item:
                    continue
                key, _, value = item.partition(" ")
                attrs[key] = value.strip().strip('"')
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id or not transcript_id:
                continue
            start, end = int(fields[3]), int(fields[4])
            chrom = normalize_contig(fields[0], strip_chr_prefix)
            transcripts[(chrom, gene_id, transcript_id, fields[6])].append(
                (start, end, attrs.get("gene_name", gene_id), attrs.get("exon_number", ""))
            )
    exons = []
    introns = []
    gene_exons = defaultdict(list)
    for (chrom, gene_id, transcript_id, strand), values in transcripts.items():
        values.sort(key=lambda item: item[0], reverse=strand == "-")
        for number, (start, end, gene_name, exon_number) in enumerate(values, 1):
            exon_id = f"exon_{gene_id}_{transcript_id}_{number}"
            exons.append((exon_id, "exon", chrom, start, end, gene_id, gene_name,
                          transcript_id, exon_number or str(number), strand))
            gene_exons[(chrom, gene_id, gene_name, strand)].append((start, end))
        genomic = sorted(values)
        for number, left in enumerate(genomic[:-1], 1):
            intron_start, intron_end = left[1] + 1, genomic[number][0] - 1
            if intron_start <= intron_end:
                intron_id = f"intron_{gene_id}_{transcript_id}_{number}"
                introns.append((intron_id, "intron", chrom, intron_start, intron_end,
                                gene_id, left[2], transcript_id, str(number), strand))
    gene_exon_unions = []
    for (chrom, gene_id, gene_name, strand), intervals in gene_exons.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for number, (start, end) in enumerate(merged, 1):
            feature_id = f"gene_exon_union_{gene_id}_{number}"
            gene_exon_unions.append(
                (feature_id, "gene_exon_union", chrom, start, end, gene_id,
                 gene_name, ".", str(number), strand)
            )
    return exons, introns, gene_exon_unions


def write_features(path: str, features):
    with open(path, "w", encoding="utf-8") as handle:
        for feature in features:
            _, _, chrom, start, end, *_ = feature
            handle.write(f"{chrom}\t{start - 1}\t{end}\n")


def submit(args, sample, bam, feature_bed, feature_tsv, output_dir, dependency=None):
    log_dir = os.path.join(args.logs_dir, sample)
    os.makedirs(log_dir, exist_ok=True)
    body = textwrap.dedent(
        f"""\
        set -euo pipefail
        mkdir -p {q(output_dir)}
        samtools depth -aa -q {args.min_mapping_quality} -Q {args.min_base_quality} \
          -b {q(feature_bed)} {q(bam)} | \
        python3 {q(args.worker)} \
          --depth - \
          --features {q(feature_tsv)} \
          --output-dir {q(output_dir)} \
          --sample {q(sample)} \
          --threshold {args.threshold}
        echo "[OK] Coverage QC completed for {sample}"
        """
    )
    command = [
        "bsub", "-J", f"step4_cov_{sample[:40]}",
        "-q", args.queue, "-n", str(args.threads),
        "-M", str(args.mem_mb), "-R", f"rusage[mem={args.mem_mb}]",
        "-W", args.wall,
        "-o", os.path.join(log_dir, "coverage.%J.out"),
        "-e", os.path.join(log_dir, "coverage.%J.err"),
    ]
    if dependency:
        command.extend(["-w", dependency])
    command.extend(["bash", "-lc", body])
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    print(result.stdout.strip())
    match = re.search(r"Job <(\d+)>", result.stdout)
    if not match:
        raise RuntimeError(f"Could not parse LSF job ID from: {result.stdout}")
    return match.group(1)


def submit_plot(args, job_ids):
    if not args.plot_script:
        return
    dependency = " && ".join(f"done({job_id})" for job_id in job_ids)
    body = textwrap.dedent(
        f"""\
        set -euo pipefail
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate {q(args.plot_env)}
        python {q(args.plot_script)} \
          --input-dir {q(args.output_dir)} \
          --output-prefix {q(args.plot_output_prefix)} \
          --chroms {q(args.chroms)} \
          --bin-size {args.plot_bin_size} \
          --max-gene-labels {args.plot_max_gene_labels} \
          --threshold {args.threshold}
        """
    )
    command = [
        "bsub", "-J", "step4_coverage_plots",
        "-q", args.queue, "-n", "1",
        "-M", str(args.mem_mb), "-R", f"rusage[mem={args.mem_mb}]",
        "-W", args.wall, "-w", dependency,
        "-o", os.path.join(args.logs_dir, "coverage_plots.%J.out"),
        "-e", os.path.join(args.logs_dir, "coverage_plots.%J.err"),
        "bash", "-lc", body,
    ]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    print(result.stdout.strip())


def write_worker(path: str):
    # The worker is kept in this file so the submitter is self-contained.
    if os.path.exists(path):
        return
    raise SystemExit(f"[ERROR] Worker script does not exist: {path}")


def run_worker(args):
    stats = {}
    features = []
    with open(args.features, encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            features.append(row)
            stats[row[0]] = [row[1], row[2], int(row[3]), int(row[4]), *row[5:], 0, 0, 0, 0, 0]
    by_chrom = defaultdict(list)
    for row in features:
        by_chrom[row[2]].append((int(row[3]), int(row[4]), row[0]))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    current_chrom = None
    intervals = []
    interval_index = 0
    active = []
    handle = sys.stdin if args.depth == "-" else open(args.depth, encoding="utf-8")
    try:
        for line in handle:
            chrom, pos_text, depth_text = line.split()[:3]
            pos, depth = int(pos_text), int(depth_text)
            if chrom != current_chrom:
                current_chrom = chrom
                intervals = by_chrom.get(chrom, [])
                interval_index = 0
                active = []
            while interval_index < len(intervals) and intervals[interval_index][0] <= pos:
                active.append(intervals[interval_index])
                interval_index += 1
            active = [interval for interval in active if interval[1] >= pos]
            for _, _, feature_id in active:
                item = stats[feature_id]
                item[-5] += 1
                item[-4] += depth
                if depth > args.threshold:
                    item[-3] += 1
                    item[-2] += depth
                item[-1] += depth == 0
    finally:
        if handle is not sys.stdin:
            handle.close()
    os.makedirs(args.output_dir, exist_ok=True)
    detail = os.path.join(args.output_dir, f"{args.sample}.feature_coverage.tsv")
    gene = defaultdict(lambda: [0, 0, 0, 0, 0, set()])
    with gzip.open(detail + ".gz", "wt", encoding="utf-8") as out:
        out.write("sample\tfeature_type\tgene_id\tgene_name\ttranscript_id\tfeature_id\t"
                  "chrom\tstart\tend\tlength\tbases_gt3x\tpercent_gt3x\t"
                  "mean_depth_gt3x\tmean_depth_all\tzero_bases\tstatus\n")
        for feature_id, item in stats.items():
            feature_type, chrom, start, end, gene_id, gene_name, transcript_id, number, strand, covered, total_depth, gt, gt_depth, zero = item
            length = end - start + 1
            # samtools -aa should provide every base; retain a defensive length bound.
            missing = max(0, length - covered)
            zero += missing
            mean_gt = (gt_depth / gt) if gt else 0.0
            mean_all = total_depth / length if length else 0.0
            status = "zero_coverage" if zero == length else ("no_gt3x" if gt == 0 else "covered")
            out.write(f"{args.sample}\t{feature_type}\t{gene_id}\t{gene_name}\t{transcript_id}\t"
                      f"{feature_id}\t{chrom}\t{start}\t{end}\t{length}\t{gt}\t"
                      f"{100.0 * gt / length:.4f}\t{mean_gt:.4f}\t{mean_all:.4f}\t{zero}\t{status}\n")
            if feature_type == "gene_exon_union":
                key = (gene_id, gene_name)
                g = gene[key]
                g[0] += length
                g[1] += gt
                g[2] += total_depth
                g[3] += gt_depth
                g[4] += zero
                g[5].add(feature_id)
    gene_path = os.path.join(args.output_dir, f"{args.sample}.gene_exon_coverage.tsv")
    with gzip.open(gene_path + ".gz", "wt", encoding="utf-8") as out:
        out.write("sample\tgene_id\tgene_name\tmerged_exon_interval_count\texon_bases\t"
                  "bases_gt3x\tpercent_exon_bases_gt3x\tmean_exon_depth_gt3x\t"
                  "mean_exon_depth_all\tzero_exon_bases\n")
        for (gene_id, gene_name), (length, gt, total, gt_depth, zero, exon_ids) in sorted(gene.items()):
            out.write(f"{args.sample}\t{gene_id}\t{gene_name}\t{len(exon_ids)}\t{length}\t{gt}\t"
                      f"{100.0 * gt / length if length else 0:.4f}\t"
                      f"{gt_depth / gt if gt else 0:.4f}\t"
                      f"{total / length if length else 0:.4f}\t{zero}\n")
    for suffix, predicate in (("no_coverage", "zero_coverage"), ("no_gt3x", "no_gt3x")):
        source = detail + ".gz"
        target = os.path.join(args.output_dir, f"{args.sample}.{suffix}.tsv")
        with gzip.open(source, "rt", encoding="utf-8") as src, gzip.open(target + ".gz", "wt", encoding="utf-8") as dst:
            header = next(src)
            dst.write(header)
            for line in src:
                if f"\t{predicate}\n" in line:
                    dst.write(line)
    summary = os.path.join(args.output_dir, f"{args.sample}.coverage_summary.tsv")
    with open(summary, "w", encoding="utf-8") as out:
        out.write("sample\tfeature_type\tfeatures\tzero_coverage_features\tno_gt3x_features\n")
        for feature_type in ("exon", "intron"):
            rows = [x for x in stats.values() if x[0] == feature_type]
            out.write(f"{args.sample}\t{feature_type}\t{len(rows)}\t"
                      f"{sum(x[-1] == x[3] - x[2] + 1 for x in rows)}\t"
                      f"{sum(x[-2] == 0 for x in rows)}\n")
    if args.delete_depth and args.depth != "-":
        os.unlink(args.depth)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="")
    parser.add_argument(
        "--manifest-base-dir",
        default="",
        help="Base directory for relative BAM paths in the manifest (default: current directory)",
    )
    parser.add_argument("--gtf")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--logs-dir")
    parser.add_argument("--queue", default="compute")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--mem-mb", type=int, default=8000)
    parser.add_argument("--wall", default="145:00")
    parser.add_argument("--threshold", type=int, default=3)
    parser.add_argument("--chroms", default="")
    parser.add_argument("--strip-chr-prefix", action="store_true")
    parser.add_argument("--min-mapping-quality", type=int, default=0)
    parser.add_argument("--min-base-quality", type=int, default=0)
    parser.add_argument("--worker", default=os.path.abspath(__file__))
    parser.add_argument("--depth")
    parser.add_argument("--features")
    parser.add_argument("--sample")
    parser.add_argument("--delete-depth", action="store_true")
    parser.add_argument("--max-concurrent", type=int, default=0)
    parser.add_argument("--plot-script", default="")
    parser.add_argument("--plot-output-prefix", default="")
    parser.add_argument("--plot-bin-size", type=int, default=50_000)
    parser.add_argument("--plot-max-gene-labels", type=int, default=28)
    parser.add_argument("--plot-env", default="resolveome-align")
    args = parser.parse_args()
    if args.depth:
        run_worker(args)
        return
    if not args.gtf or not args.logs_dir:
        parser.error("--gtf and --logs-dir are required when submitting jobs")
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.logs_dir, exist_ok=True)
    exons, introns, gene_exon_unions = read_gtf(args.gtf, args.strip_chr_prefix)
    if args.chroms:
        selected = {
            normalize_contig(item.strip(), args.strip_chr_prefix)
            for item in args.chroms.split(",")
            if item.strip()
        }
        exons = [feature for feature in exons if feature[2] in selected]
        introns = [feature for feature in introns if feature[2] in selected]
        gene_exon_unions = [
            feature for feature in gene_exon_unions if feature[2] in selected
        ]
    features = exons + introns + gene_exon_unions
    feature_tsv = os.path.join(args.output_dir, "annotation_features.tsv")
    feature_bed = os.path.join(args.output_dir, "annotation_features.bed")
    with open(feature_tsv, "w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(features)
    write_features(feature_bed, features)
    manifest = args.manifest
    if not manifest:
        raise SystemExit("--manifest is required when submitting jobs")
    job_ids = []
    with open(manifest, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("pass_qc", "true").lower() not in ("true", "1", "yes"):
                continue
            bam = row["bam"]
            if not os.path.isabs(bam):
                base_dir = args.manifest_base_dir or os.getcwd()
                bam = os.path.abspath(os.path.join(base_dir, bam))
            if not os.path.exists(bam):
                raise SystemExit(f"Missing BAM for {row['sample_id']}: {bam}")
            dependency = None
            if args.max_concurrent and len(job_ids) >= args.max_concurrent:
                dependency = f"ended({job_ids[-args.max_concurrent]})"
            job_ids.append(
                submit(
                    args,
                    row["sample_id"],
                    bam,
                    feature_bed,
                    feature_tsv,
                    os.path.join(args.output_dir, row["sample_id"]),
                    dependency,
                )
            )
    submit_plot(args, job_ids)


if __name__ == "__main__":
    main()
