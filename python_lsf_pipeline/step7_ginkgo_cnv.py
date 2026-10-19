#!/usr/bin/env python3
"""Standalone CUSTOM Ginkgo stage; does not replace or modify CNVkit."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent.parent  # workflow_DNA root (this script now lives one level deeper, in python_lsf_pipeline/)
SOURCE = HERE / "resolveome_align/basej-public-pipelines/basej-dnaqc"
REFERENCE = (HERE / "../pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/"
             "Annotation/Ginkgo_reference/GRCH38_RL_50_Bin1MB").resolve()
SOURCE_FILES = ("container/ginkgo/binUnsorted", "container/ginkgo/cnv_ginkgo.R",
                "container/ginkgo/CNVcaller", "container/ginkgo/rds_to_flat.R",
                "container/ginkgo/rscript_parse_rds_cnv_metrics.R",
                "container/ginko_parser/scripts/divide_clouds_by_sample.py",
                "container/custom_r_qcplots/scripts/cnvSummarizer.R")
PACKAGES = ("ctc", "DNAcopy", "inline", "gplots", "scales", "plyr", "ggplot2",
            "gridExtra", "dplyr", "magrittr", "optparse")


def command(argv, **kwargs):
    return subprocess.run([str(x) for x in argv], check=True, text=True, **kwargs)


def write_tsv(path, fields, rows):
    with Path(path).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def canonical(name):
    return name[3:] if name.startswith("chr") else name


TRUE = {"true", "t", "yes", "y", "1", "pass", "passed"}
FALSE = {"false", "f", "no", "n", "0", "fail", "failed"}


def resolve_input(value, manifest, root, label, line):
    path = Path(value).expanduser()
    if path.is_absolute():
        candidates = [path]
    elif root:
        candidates = [Path(root) / path]
    else:
        # Workflow manifests store paths relative to the results root, above manifests/.
        candidates = [base / path for base in (manifest.parent, *manifest.parents)]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    raise ValueError(f"Manifest line {line}: {label} not found: {value}")


def load_manifest(path, path_root=None):
    """Return (cells, excluded) from the workflow sample_id/bam/bai/pass_qc manifest."""
    path = Path(path).resolve()
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        sample_field = next((f for f in ("sample_id", "sample_name", "sample", "cell_id")
                             if f in fields), None)
        bam_field = next((f for f in ("bam", "bam_path", "aligned_bam") if f in fields), None)
        if not sample_field or not bam_field:
            raise ValueError(f"Manifest header {fields} requires sample_id and bam "
                             "(aliases sample_name/sample/cell_id and bam_path/aligned_bam)")
        cells, excluded, names = [], [], set()
        for line, row in enumerate(reader, 2):
            if None in row or any(row.get(f) is None for f in fields):
                raise ValueError(f"Manifest line {line}: column count differs from header")
            name, bam = row[sample_field].strip(), row[bam_field].strip()
            if not name or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
                raise ValueError(f"Manifest line {line}: unsafe or empty sample_id {name!r}")
            if name in names or not bam:
                raise ValueError(f"Manifest line {line}: duplicate sample or empty bam")
            names.add(name)
            if "pass_qc" in fields:
                flag = row["pass_qc"].strip().lower()
                if flag in FALSE:
                    excluded.append({"sample_name": name, "reason": f"pass_qc={row['pass_qc']}"})
                    continue
                if flag not in TRUE:
                    raise ValueError(f"Manifest line {line}: invalid pass_qc {row['pass_qc']!r}")
            cell = {"sample_name": name,
                    "bam_path": resolve_input(bam, path, path_root, "bam", line), "bai_path": ""}
            if "bai" in fields:
                cell["bai_path"] = resolve_input(row["bai"].strip(), path, path_root, "bai", line)
            cells.append(cell)
    if not cells:
        raise ValueError("Manifest has no passing cells")
    return cells, excluded


def load_references(directory, bin_size, read_length):
    suffix = f"variable_{bin_size}_{read_length}_bwa"
    paths = {key: Path(directory) / f"{prefix}{suffix}"
             for key, prefix in (("bins", ""), ("gc", "GC_"), ("bounds", "bounds_"))}
    for path in paths.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing/empty Ginkgo resource: {path}; no read-length fallback")
    with paths["bins"].open() as stream:
        table = list(csv.reader(stream, delimiter="\t"))
    if table[0] != ["CHR", "END"]:
        raise ValueError("Ginkgo bin header must be CHR<TAB>END")
    bins, last, chroms = [], {}, []
    for row in table[1:]:
        if len(row) != 2:
            raise ValueError("Malformed Ginkgo bin row")
        chrom, end = row[0], int(row[1])
        if end <= last.get(chrom, 0):
            raise ValueError(f"Non-increasing Ginkgo bin end on {chrom}")
        if chrom not in last:
            if canonical(chrom) in {canonical(c) for c in chroms}:
                raise ValueError("Ambiguous chromosome aliases in Ginkgo bins")
            chroms.append(chrom)
        elif chroms[-1] != chrom:
            raise ValueError("Ginkgo chromosome blocks must be contiguous")
        bins.append((chrom, end))
        last[chrom] = end
    gc = [float(x.strip()) for x in paths["gc"].read_text().splitlines()]
    if len(gc) != len(bins) or any(not 0 <= g <= 1 for g in gc):
        raise ValueError("GC resource must have one finite [0,1] value per bin")
    bounds = [line.split("\t") for line in paths["bounds"].read_text().splitlines()]
    expected = [(bins[i - 1][0], i + 1) for i in range(1, len(bins))
                if bins[i][0] != bins[i - 1][0]]
    if [(c, int(n)) for c, n in bounds] != expected:
        raise ValueError("Ginkgo bounds do not match bin chromosome transitions")
    return {k: str(v.resolve()) for k, v in paths.items()}, bins


def chromosome_map(header, bins, fai):
    """Require every modeled chromosome and check BAM/FASTA lengths; audit extras."""
    fasta = {}
    for line in Path(fai).read_text().splitlines():
        parts = line.split("\t")
        fasta[parts[0]] = int(parts[1])
    sequences = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            tags = dict(p.split(":", 1) for p in line.split("\t")[1:] if ":" in p)
            if "SN" not in tags or "LN" not in tags or tags["SN"] in sequences:
                raise ValueError("BAM header has missing/duplicate SN or missing LN")
            sequences[tags["SN"]] = int(tags["LN"])
    if not sequences:
        raise ValueError("BAM header has no @SQ records")
    for chrom, length in sequences.items():
        if fasta.get(chrom) != length:
            raise ValueError(f"BAM @SQ {chrom}:{length} disagrees with supplied FASTA .fai")
    aliases = {}
    for chrom in sequences:
        key = canonical(chrom)
        if key in aliases:
            raise ValueError(f"Ambiguous BAM chromosome aliases: {chrom}, {aliases[key]}")
        aliases[key] = chrom
    mapping = {}
    for chrom, end in bins:
        bam_chrom = aliases.get(canonical(chrom))
        if bam_chrom is None:
            raise ValueError(f"BAM header missing modeled chromosome {chrom}")
        if end > sequences[bam_chrom]:
            raise ValueError(f"Ginkgo bin {chrom}:{end} exceeds BAM/FASTA chromosome length")
        mapping[bam_chrom] = chrom
    return mapping, sequences


def tool_paths(prefix):
    return {name: str(Path(prefix) / name)
            for name in ("samtools", "bedtools", "Rscript", "python")}


def environment_report(tools, source):
    report, errors = {}, []
    for name, path in tools.items():
        try:
            argv = [path, "-c", "import pandas; print('pandas', pandas.__version__)"] \
                if name == "python" else [path, "--version"]
            run = command(argv, capture_output=True)
            report[name] = (run.stdout + run.stderr).splitlines()[0]
        except (OSError, subprocess.CalledProcessError) as exc:
            errors.append(f"{name}: {path}: {exc}")
    for name, usage_text in (("binUnsorted", "USAGE: binFile #bins readsFile nameReadsFile outfile"),
                             ("CNVcaller", "USAGE: ./temp CN_file outfile1 outfile2")):
        try:
            run = subprocess.run([str(Path(source) / "container/ginkgo" / name)], text=True,
                                 capture_output=True, timeout=10)
            usage = run.stdout + run.stderr
            if usage_text not in usage:
                raise ValueError(f"{name} executable smoke failed (exit {run.returncode}): {usage}")
            report[name] = usage.strip()
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            errors.append(str(exc))
    expr = ('p<-c(' + ",".join(json.dumps(p) for p in PACKAGES) + '); '
            'for(x in p) cat(x,if(requireNamespace(x,quietly=TRUE)) '
            'as.character(packageVersion(x)) else "MISSING","\\n"); '
            'cat("jpeg",capabilities("jpeg"),"\\n"); '
            'quit(status=as.integer(any(!vapply(p,requireNamespace,logical(1),quietly=TRUE))'
            ' || !capabilities("jpeg")))')
    try:
        run = subprocess.run([tools["Rscript"], "-e", expr], text=True, capture_output=True)
        report["R_packages"] = run.stdout.strip()
        if run.returncode:
            errors.append("R package/jpeg probe failed: " + (run.stdout + run.stderr).strip())
    except OSError as exc:
        errors.append(str(exc))
    return report, errors


def bin_cell(config, sample):
    stage = Path(config["stage"])
    cell = next(c for c in config["cells"] if c["sample_name"] == sample)
    work = stage / "cells" / sample
    status = {"sample_name": sample, "status": "ERROR"}
    tools = config["tools"]
    try:
        _, bins = load_references(config["reference_dir"], config["bin_size"], config["read_length"])
        command([tools["samtools"], "quickcheck", "-v", cell["bam_path"]], capture_output=True)
        header = command([tools["samtools"], "view", "-H", cell["bam_path"]], capture_output=True).stdout
        mapping, sequences = chromosome_map(header, bins, config["ref"] + ".fai")
        (work / "bam_header.sam").write_text(header)
        counts = {chrom: 0 for chrom in sequences}
        bed_path = work / "reads.bed"
        # binUnsorted uses chromosome/start; retain all alignments bamtobed emits.
        with bed_path.open("w") as bed, (work / "bamtobed.stderr").open("w") as err:
            proc = subprocess.Popen([tools["bedtools"], "bamtobed", "-i", cell["bam_path"]],
                                    stdout=subprocess.PIPE, stderr=err, text=True)
            try:
                for line in proc.stdout:
                    row = line.rstrip("\n").split("\t")
                    if len(row) < 3 or row[0] not in counts:
                        raise ValueError("Malformed BED or chromosome absent from BAM header")
                    counts[row[0]] += 1
                    if row[0] in mapping:
                        row[0] = mapping[row[0]]
                        bed.write("\t".join(row) + "\n")
            except BaseException:
                proc.terminate()
                proc.wait()
                raise
            finally:
                proc.stdout.close()
            if proc.wait():
                raise ValueError("bedtools bamtobed failed; see bamtobed.stderr")
        write_tsv(work / "chromosome_audit.tsv",
                  ["bam_chromosome", "ginkgo_chromosome", "bed_records", "disposition"],
                  [{"bam_chromosome": c, "ginkgo_chromosome": mapping.get(c, ""),
                    "bed_records": counts[c],
                    "disposition": "modeled" if c in mapping else "not_in_source_reference"}
                   for c in counts])
        modeled = sum(counts[c] for c in mapping)
        status.update(bed_records=sum(counts.values()), modeled_records=modeled,
                      unmodeled_records=sum(counts.values()) - modeled)
        if modeled == 0:
            raise ValueError("Zero modeled reads; cohort retained, CNV inference blocked")
        output = work / "binUnsorted.outdata"
        with (work / "binUnsorted.log").open("w") as log:
            # This executable counts the header in its size argument; otherwise chrY's last bin is lost.
            command([Path(config["source"]) / "container/ginkgo/binUnsorted",
                     config["resources"]["bins"], len(bins) + 1, bed_path, sample, output],
                    stdout=log, stderr=subprocess.STDOUT)
        values = read_counts(output, sample, len(bins))
        if sum(values) != modeled:
            raise ValueError(f"binUnsorted counted {sum(values)} of {modeled} modeled BED records")
        status.update(status="OK", reason="")
        bed_path.unlink()
    except Exception as exc:
        status["reason"] = str(exc)
        if isinstance(exc, subprocess.CalledProcessError):
            status["reason"] += ": " + (exc.stderr or "")
    (work / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    return status["status"] == "OK"


def read_counts(path, sample, number):
    rows = Path(path).read_text().splitlines()
    if not rows or rows[0].strip() != sample:
        raise ValueError(f"{path}: binUnsorted header differs from sample {sample}")
    values = [int(line.strip()) for line in rows[1:]]
    if len(values) != number or any(n < 0 for n in values):
        raise ValueError(f"{path}: expected {number} nonnegative bin counts, got {len(values)}")
    return values


def metric_values(summary):
    result = {}
    for key in ("MAPD_CNV_Log2", "SKEW_CNV"):
        value = summary[key]
        undefined = value in ("", "NA", "NaN", "Inf", "-Inf") or not math.isfinite(float(value))
        result[key] = "NA" if undefined else value
        result[f"{key}_reason"] = ("Original Ginkgo summarizer returned an undefined metric"
                                  if undefined else "")
    return result


def merge_cohort(config):
    stage = Path(config["stage"])
    statuses = []
    for cell in config["cells"]:
        status_path = stage / "cells" / cell["sample_name"] / "status.json"
        statuses.append(json.loads(status_path.read_text()) if status_path.exists() else
                        {"sample_name": cell["sample_name"], "status": "ERROR",
                         "reason": "Binning job did not produce status.json (see LSF logs)"})
    write_tsv(stage / "cohort_status.tsv",
              ["sample_name", "status", "reason", "bed_records", "modeled_records", "unmodeled_records"],
              [{k: s.get(k, "") for k in ("sample_name", "status", "reason", "bed_records",
                                         "modeled_records", "unmodeled_records")} for s in statuses])
    if any(s["status"] != "OK" for s in statuses):
        raise ValueError("Cohort inference blocked: see cohort_status.tsv; no cells silently removed")
    _, bins = load_references(config["reference_dir"], config["bin_size"], config["read_length"])
    names = [c["sample_name"] for c in config["cells"]]
    columns = [read_counts(stage / "cells" / name / "binUnsorted.outdata", name, len(bins))
               for name in names]
    collisions = [(a, b) for a in names for b in names if a != b and a in b]
    if collisions:
        # divide_clouds_by_sample.py selects CNV columns by substring.
        raise ValueError(f"Sample IDs collide in source per-sample parser: {collisions[:3]}")
    work = stage / "cohort"
    size = config["bin_size"]
    merged = f"binUnsorted.{size}.outdata"
    with (work / merged).open("w") as stream:
        stream.write("\t".join(names) + "\n")
        for row in zip(*columns):
            stream.write("\t".join(map(str, row)) + "\n")
    source = Path(config["source"])
    script_dir = source / "container/ginkgo"
    resources = config["resources"]
    rscript = config["tools"]["Rscript"]
    facs = 0
    if config["is_haplotype"] != 0:
        # GINKGO_SEGMENTATION_R writes one sample/ploidy prior per cell and uses FACS mode 1.
        (work / "ploidy.txt").write_text(
            "".join(f"{name}\t{config['is_haplotype']}\n" for name in names))
        facs = 1
    rds = f"ginkgo_res.binsize_{size}.RDS"
    segcopy = f"SegCopy.binsize_{size}.tsv"
    with (work / "ginkgo.log").open("w") as log:
        def run(argv, cwd=work):
            command(argv, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)

        run([rscript, script_dir / "cnv_ginkgo.R", merged, resources["bins"], resources["gc"],
             resources["bounds"], config["min_ploidy"], config["max_ploidy"],
             config["min_width"], facs])
        (work / "ginkgo_res.RDS").rename(work / rds)
        (work / "SegCopy").rename(work / segcopy)
        plots = sorted(work.glob("*.jpeg"))
        if len(plots) != len(names):
            raise ValueError(f"Expected {len(names)} Ginkgo JPEG plots, found {len(plots)}")
        renamed = []
        for plot in plots:
            renamed.append(plot.rename(work / f"cnv_binsize_{size}_{plot.name}").name)
        run(["tar", "-czf", f"cnv_plots_binsize_{size}.tar.gz", *renamed])
        run([script_dir / "CNVcaller", segcopy, f"CNV1_binsize_{size}.tsv",
             f"CNV2_binsize_{size}.tsv"])
        run([rscript, script_dir / "rds_to_flat.R", rds])
        per_sample = work / "per_sample"
        per_sample.mkdir(exist_ok=True)
        run([config["tools"]["python"],
             source / "container/ginko_parser/scripts/divide_clouds_by_sample.py",
             "-c", work / "clouds.tsv", "-v1", work / f"CNV1_binsize_{size}.tsv",
             "-v2", work / f"CNV2_binsize_{size}.tsv", "-b", resources["bins"]], cwd=per_sample)
        run([rscript, script_dir / "rscript_parse_rds_cnv_metrics.R", rds])
        (work / "AllSample-GinkgoSegmentSummary.txt").rename(
            work / "Legacy-GinkgoSegmentSummary.txt")
        # main.nf QC step overwrites the legacy summary with cnvSummarizer MAPD/SKEW.
        run([rscript, source / "container/custom_r_qcplots/scripts/cnvSummarizer.R",
             "--rds_file", rds, "--out_file", "AllSample-GinkgoSegmentSummary.txt"])
    expected = [rds, segcopy, merged, "SegNorm", "SegFixed", "SegBreaks", "SegStats",
                f"cnv_plots_binsize_{size}.tar.gz", f"CNV1_binsize_{size}.tsv",
                f"CNV2_binsize_{size}.tsv", "clouds.tsv", "rectangles.tsv",
                "sampleploidy.tsv", "Legacy-GinkgoSegmentSummary.txt",
                "AllSample-GinkgoSegmentSummary.txt"]
    for name in expected:
        if not (work / name).is_file() or (work / name).stat().st_size == 0:
            raise ValueError(f"Missing/empty cohort result: {name}")
    missing_split = [n for n in names if not (per_sample / f"{n}_clouds.tsv").is_file()]
    if missing_split:
        raise ValueError(f"Per-sample parser missed cells: {missing_split}")
    with (work / segcopy).open() as stream:
        header = stream.readline().rstrip("\n").split("\t")
    if header != ["CHR", "START", "END"] + names:
        raise ValueError("SegCopy lost or renamed cohort cells")
    with (work / "AllSample-GinkgoSegmentSummary.txt").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    summary = {row["SampleId"]: row for row in rows}
    if set(summary) != set(names) or len(rows) != len(names):
        raise ValueError("Ginkgo summary lost or renamed cohort cells")
    fields = ["sample_name"] + [f for f in rows[0] if f != "SampleId"]
    write_tsv(stage / "ginkgo_metrics.tsv", fields,
              [{"sample_name": n, **{f: summary[n][f] for f in fields[1:]}} for n in names])
    write_tsv(stage / "source_metrics.tsv",
              ["sample_name", "MAPD_CNV_Log2", "SKEW_CNV",
               "MAPD_CNV_Log2_reason", "SKEW_CNV_reason"],
              [{"sample_name": n, **metric_values(summary[n])} for n in names])
    (stage / "COMPLETE").write_text("CUSTOM Ginkgo; all manifest cells retained\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def submit(config, path, queue, memory):
    stage = Path(config["stage"])
    jobs = []
    env = dict(os.environ, PATH=config["tool_bin"] + os.pathsep + os.environ.get("PATH", ""))

    def launch(name, argv, dependency=None):
        script = stage / "jobs" / f"{name}.sh"
        script.write_text("#!/bin/bash\nset -euo pipefail\n"
                          f"export PATH={shlex.quote(env['PATH'])}\n"
                          + shlex.join([str(x) for x in argv]) + "\n")
        script.chmod(0o750)
        job = {"job_name": name, "sample_name": name[4:] if name.startswith("bin_") else "COHORT",
               "job_id": "", "status": "SUBMITTING", "dependency": dependency or "",
               "script": str(script)}
        jobs.append(job)
        write_tsv(stage / "job_manifest.tsv", list(job), jobs)
        args = ["bsub", "-q", queue, "-J", name, "-n", "1", "-M", str(memory),
                "-R", f"rusage[mem={memory}]", "-oo", str(stage / "logs" / f"{name}.out"),
                "-eo", str(stage / "logs" / f"{name}.err")]
        if dependency:
            args += ["-w", dependency]
        try:
            run = command(args + ["bash", str(script)], capture_output=True, env=env)
            found = re.search(r"Job <(\d+)>", run.stdout)
            if not found:
                raise ValueError("bsub returned no parseable job ID: " + run.stdout)
            job.update(job_id=found.group(1), status="SUBMITTED")
        except Exception:
            job["status"] = "SUBMIT_FAILED"
            raise
        finally:
            write_tsv(stage / "job_manifest.tsv", list(job), jobs)
        return job["job_id"]

    ids = [launch("bin_" + c["sample_name"],
                  [sys.executable, __file__, "--worker-config", path, "--worker-sample", c["sample_name"]])
           for c in config["cells"]]
    # ended() lets the cohort report every failed header/bin job, not hang on done().
    launch("ginkgo_cohort", [sys.executable, __file__, "--worker-config", path],
           " && ".join(f"ended({job})" for job in ids))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest")
    parser.add_argument("--outdir", help="New results root; creates ginkgo/ without overwriting")
    parser.add_argument("--path-root", help="Base for relative manifest BAM/BAI paths "
                        "(default: manifest directory, then its parents)")
    parser.add_argument("--ref", help="FASTA with existing .fai; never modifies the reference")
    parser.add_argument("--queue", default="standard")
    parser.add_argument("--mode", choices=("lsf", "local"), default="lsf")
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--reference-dir", type=Path, default=REFERENCE)
    parser.add_argument("--tool-bin", default=str(Path(sys.executable).parent),
                        help="Current resolveome-align bin directory; no PATH fallback")
    parser.add_argument("--bin-size", type=int, choices=(1000000,), default=1000000)
    parser.add_argument("--read-length", type=int, default=50)
    parser.add_argument("--min-ploidy", type=float, default=1.5)
    parser.add_argument("--max-ploidy", type=float, default=6.0)
    parser.add_argument("--min-width", type=int, default=5)
    parser.add_argument("--is-haplotype", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--memory-mb", type=int, default=8000)
    parser.add_argument("--check-environment", action="store_true")
    parser.add_argument("--worker-config", help=argparse.SUPPRESS)
    parser.add_argument("--worker-sample", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.worker_config:
        config = json.loads(Path(args.worker_config).read_text())
        os.environ["TMPDIR"] = str(Path(config["stage"]) / "tmp")
        if args.worker_sample:
            return 0 if bin_cell(config, args.worker_sample) else 1
        merge_cohort(config)
        return 0
    source = args.source_dir.resolve()
    resources, bins = load_references(args.reference_dir, args.bin_size, args.read_length)
    for name in SOURCE_FILES:
        if not (source / name).is_file():
            raise ValueError(f"Missing source Ginkgo file: {source / name}")
    tools = tool_paths(args.tool_bin)
    report, errors = environment_report(tools, source)
    if args.check_environment:
        print(json.dumps({"tools": report, "resources": resources, "bin_count": len(bins),
                          "errors": errors}, indent=2))
        return 1 if errors else 0
    if errors:
        raise ValueError("Ginkgo prerequisites unavailable:\n" + "\n".join(errors))
    if not args.manifest or not args.outdir or not args.ref:
        parser.error("--manifest, --outdir and --ref are required")
    if not 0 < args.min_ploidy <= args.max_ploidy or args.min_width < 2:
        parser.error("Invalid ploidy range or min-width")
    cells, excluded = load_manifest(args.manifest, args.path_root)
    ref = Path(args.ref).resolve()
    if not ref.is_file() or not Path(str(ref) + ".fai").is_file():
        raise ValueError(f"FASTA and existing .fai required: {ref}")
    outdir = Path(args.outdir).resolve()
    if "DNA_Ilumina_output_1" in outdir.parts:
        raise ValueError("Protected old output root: DNA_Ilumina_output_1")
    stage = outdir / "ginkgo"
    stage.mkdir(parents=True, exist_ok=False)
    for directory in ("jobs", "logs", "cells", "cohort", "tmp"):
        (stage / directory).mkdir()
    # R and sort scratch stay inside the run directory.
    os.environ["TMPDIR"] = str(stage / "tmp")
    for cell in cells:
        (stage / "cells" / cell["sample_name"]).mkdir()
    config = {"stage": str(stage), "cells": cells, "ref": str(ref), "source": str(source),
              "resources": resources, "reference_dir": str(args.reference_dir.resolve()),
              "tools": tools, "tool_bin": str(Path(args.tool_bin).resolve()),
              "bin_size": args.bin_size, "read_length": args.read_length,
              "min_ploidy": args.min_ploidy, "max_ploidy": args.max_ploidy,
              "min_width": args.min_width, "is_haplotype": args.is_haplotype,
              "method": "CUSTOM Ginkgo; CNVkit unchanged",
              "ploidy_prior": "Source GINKGO_SEGMENTATION_R: ploidy.txt sample<TAB>is_haplotype, FACS mode 1",
              "excluded_by_pass_qc": excluded,
              "environment": report, "manifest": str(Path(args.manifest).resolve()),
              "resource_sha256": {k: sha256(v) for k, v in resources.items()},
              "source_sha256": {n: sha256(source / n) for n in SOURCE_FILES}}
    path = stage / "config.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    write_tsv(stage / "input_manifest.tsv", ["sample_name", "bam_path", "bai_path"], cells)
    write_tsv(stage / "excluded_cells.tsv", ["sample_name", "reason"], excluded)
    if args.mode == "local":
        for cell in cells:
            bin_cell(config, cell["sample_name"])
        merge_cohort(config)
    else:
        submit(config, path, args.queue, args.memory_mb)
    print(f"CUSTOM Ginkgo {args.mode}: {stage}; cohort size={len(cells)}; "
          f"pass_qc excluded={len(excluded)}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
