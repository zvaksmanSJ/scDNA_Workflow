#!/usr/bin/env python3
"""Run the revised Illumina stages in order, preserving the previous results."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent  # ResolveOME_illumina project root (this script now lives one level deeper, in workflow_DNA/python_lsf_pipeline/)
TOOLS = HERE.parent / "resolveome_align"  # resolveome_align/ lives at workflow_DNA root, one level above this script's own python_lsf_pipeline/ subfolder
SOURCE = TOOLS / "basej-public-pipelines"
ANNOTATION = PROJECT / "pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Annotation"
CHROMS = ",".join([*(str(i) for i in range(1, 23)), "X", "Y"])
LINEAGE_SCRIPTS = (
    "rscript_sequoia_build_phylogeny_only.R", "rscript_minimal_variant_placement.R",
    "rscript_postprocess_vaf_tree_egg.R", "theme_ohchibi_pubr.R",
)


def run(argv, **kwargs):
    print("[CMD]", shlex.join(map(str, argv)), flush=True)
    try:
        return subprocess.run(list(map(str, argv)), check=True, text=True, **kwargs)
    except subprocess.CalledProcessError as exc:
        for output in (exc.stdout, exc.stderr):
            if output:
                print(output, file=sys.stderr, flush=True)
        raise


def require_file(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing/empty prerequisite: {path}")


def read_rows(path, delimiter="\t"):
    require_file(path)
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def preflight(args):
    """Probe executables/images, not merely filenames, before any submission."""
    from step1_submit_resolveome_alignment_IluminaOnly_lsf import read_cells_tsv, validate_inputs
    validate_inputs(args)
    rows = read_cells_tsv(args.cells)
    if not rows:
        raise ValueError("No input cells")
    required_tools = ("bsub", "bjobs", "apptainer", "fastp", "fastqc", "bwa-mem2",
                      "samtools", "bcftools", "picard", "multiqc", "cnvkit.py",
                      "snpEff", "bedtools", "Rscript", "bam2mr", "preseq-bioskryb")
    missing = [name for name in required_tools if not shutil.which(name)]
    if missing:
        raise RuntimeError(f"Missing tools in the current project environment: {missing}")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "python": sys.executable,
              "tools": {name: shutil.which(name) for name in required_tools}, "cells": len(rows)}
    probes = {
        "fastp": ["fastp", "--version"], "samtools": ["samtools", "--version"],
        "bcftools": ["bcftools", "--version"], "bedtools": ["bedtools", "--version"],
        "fastqc": ["fastqc", "--version"], "multiqc": ["multiqc", "--version"],
        "cnvkit": ["cnvkit.py", "version"], "snpEff": ["snpEff", "-version"],
        "nextflow": [TOOLS / "tools/nextflow", "-version"],
    }
    report["versions"] = {}
    for name, command in probes.items():
        result = run(command, capture_output=True)
        report["versions"][name] = (result.stdout + result.stderr).strip()
    result = subprocess.run(["picard", "CollectMultipleMetrics", "--version"],
                            capture_output=True, text=True)
    picard_version = (result.stdout + result.stderr).strip()
    if result.returncode not in (0, 1) or not re.search(r"Version:3\.\d+\.", picard_version):
        raise RuntimeError(f"Picard version probe failed: {picard_version}")
    report["versions"]["picard"] = picard_version
    if "0.20.1" not in report["versions"]["fastp"]:
        raise RuntimeError("Reference parity requires fastp 0.20.1")
    for name in ("deepvariant.sif", "manta.sif", "vep.r112.sif"):
        require_file(TOOLS / "containers" / name)
    model = ANNOTATION / "DeepVariant/Illumina/bioskryb-af-20241102/model"
    for suffix in (".index", ".data-00000-of-00001"):
        require_file(model / f"sc_af.model.ckpt{suffix}")
    populations = model.parent / "population_vcfs_nochr"
    for chromosome in [*(str(i) for i in range(1, 23)), "X", "Y"]:
        require_file(populations / f"cohort-chr{chromosome}.release_missing2ref.no_calls.vcf.gz")
    for cache in ("111_GRCh38", "115_GRCh38"):
        if not (ANNOTATION / "VEP/homo_sapiens" / cache).is_dir():
            raise FileNotFoundError(f"Required VEP cache {cache}")
    require_file(TOOLS / "refs/annotation/genes.bed")
    require_file(TOOLS / "refs/annotation/gencode.v46.annotation.gtf.gz")
    image_dir = TOOLS / "containers/basej-somatic"
    for image in ("utility-lite-0.2.0", "variantannotation_0.4", "r_ado_1.2",
                  "custom_snp_somatic_filter_sequoia_feb2026"):
        require_file(image_dir / f"{image}.sif")
    cache = TOOLS / "apptainer_cache"
    for name in ("bcftools-1.21--h8b25389_0", "bcftools-1.14--h88f3f91_0",
                 "samtools-1.21--h50ea8bc_0", "multiqc-1.33--pyhdfd78af_0",
                 "ensembl-vep-115.2--pl5321h2a3209d_1"):
        require_file(cache / f"quay.io-biocontainers-{name}.img")
    lineage_image = image_dir / "custom_snp_somatic_filter_sequoia_feb2026.sif"
    packages = ("optparse", "ape", "seqinr", "data.table", "treemut", "dplyr", "ggplot2",
                "ggtree", "reshape2", "tidytree", "egg", "scales")
    code = "for(p in c(" + ",".join(json.dumps(p) for p in packages) + (
        ")) { if(!requireNamespace(p,quietly=TRUE)) stop(paste('Missing',p));"
        "cat(p,as.character(packageVersion(p)),'\\n') }"
    )
    result = run(["apptainer", "exec", lineage_image, "Rscript", "-e", code], capture_output=True)
    report["lineage_packages"] = result.stdout
    for script in LINEAGE_SCRIPTS:
        run(["apptainer", "exec", lineage_image, "test", "-s", f"/usr/local/bin/scripts/{script}"])
    run(["apptainer", "exec", lineage_image, "/usr/local/bin/scripts/mpboot-avx", "-h"],
        capture_output=True)
    result = run([sys.executable, HERE / "step7_ginkgo_cnv.py", "--check-environment"],
                 capture_output=True)
    report["ginkgo"] = json.loads(result.stdout)
    run(["apptainer", "exec", TOOLS / "containers/manta.sif",
         "/usr/local/bin/configManta.py", "--help"],
        capture_output=True)
    require_file(TOOLS / "step3_bundle/resources/cnvkit/reference.cnn")
    run(["apptainer", "exec", TOOLS / "containers/vep.r112.sif", "vep", "--help"],
        capture_output=True)
    run(["apptainer", "exec",
         cache / "quay.io-biocontainers-ensembl-vep-115.2--pl5321h2a3209d_1.img",
         "vep", "--help"], capture_output=True)
    # The model/channel compatibility is unchanged from the validated workflow.
    for executable in ("make_examples", "call_variants", "postprocess_variants"):
        run(["apptainer", "exec", TOOLS / "containers/deepvariant.sif",
             "/bin/test", "-x", f"/opt/deepvariant/bin/{executable}"])
    require_file(TOOLS / "basej_genomes_nochr/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/GATK_bundle/Sequence/genome.fa.fai")
    report["reference_fai_sha256"] = hashlib.sha256(Path(args.ref + ".fai").read_bytes()).hexdigest()
    return report


def wait_jobs(job_ids, poll_seconds=60):
    remaining = set(job_ids)
    if not remaining or any(not str(jid).isdigit() for jid in remaining):
        raise ValueError("No valid submitted LSF job IDs")
    while remaining:
        pending = set()
        for offset in range(0, len(remaining), 100):
            batch = sorted(remaining)[offset:offset + 100]
            result = run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *batch],
                         capture_output=True)
            seen = {}
            for line in result.stdout.splitlines():
                fields = line.split()
                if len(fields) != 2 or not fields[0].isdigit():
                    raise RuntimeError(f"Unrecognized LSF status output: {line}")
                seen[fields[0]] = fields[1]
            if set(seen) != set(batch):
                raise RuntimeError(f"Cannot recover job status for {set(batch) - set(seen)}; inspect LSF history before resuming")
            for jid, status in seen.items():
                if status == "DONE":
                    continue
                if status not in {"PEND", "RUN", "WAIT"}:
                    raise RuntimeError(f"LSF job {jid}: {status}; dependent stages not started")
                pending.add(jid)
        remaining = pending
        if remaining:
            print(f"[WAIT] {len(remaining)} jobs outstanding", flush=True)
            time.sleep(poll_seconds)


def write_nf_config(root, kind, queue):
    config = root / "scripts" / f"{kind}_lsf_apptainer.config"
    image = TOOLS / "containers/basej-somatic/custom_snp_somatic_filter_sequoia_feb2026.sif"
    runtime = root / "scripts/lineage_runtime"
    bindings = ""
    if kind == "lineage":
        runtime.mkdir(parents=True, exist_ok=True)
        for script in LINEAGE_SCRIPTS:
            run(["apptainer", "exec", "-B", "/research_jude:/research_jude", image,
                 "cp", f"/usr/local/bin/scripts/{script}", runtime / script])
        bindings = " ".join(f"-B {runtime / script}:/usr/local/bin/{script}:ro" for script in LINEAGE_SCRIPTS)
    include = SOURCE / "basej-somatic/conf/local_lsf_apptainer.config"
    content = (
        (f"includeConfig '{include}'\n" if kind == "somatic" else "")
        + f"""
docker.enabled = false
singularity.enabled = true
singularity.autoMounts = true
singularity.cacheDir = '{TOOLS / "apptainer_cache"}'
singularity.runOptions = '--cleanenv --no-mount /programs -B /research_jude:/research_jude {bindings}'
process {{
    executor = 'lsf'
    queue = '{queue}'
    maxForks = 8
    withName: /.*/ {{
        errorStrategy = 'terminate'
        maxRetries = 0
    }}
"""
    )
    if kind == "somatic":
        content += f"""
    withName: 'FILTER_VEP_GERMLINE_CHR' {{
        container = '{TOOLS / "containers/basej-somatic/variantannotation_0.4.sif"}'
    }}
"""
    else:
        content += f"""
    withName: /SEQUOIA_.+|POSTPROCESS_.+|GATHER_LINEAGE_.+|CHECK_MATRIX_NON_EMPTY_.+/ {{
        container = '{image}'
    }}
"""
    content += "}\n"
    config.write_text(content)
    return config


def python_stage(name, *arguments):
    return [sys.executable, HERE / name, *map(str, arguments)]


def commands(args):
    root = args.outdir
    passed = root / "manifests/passed_cells_only_manifest.tsv"
    gtf = TOOLS / "refs/annotation/gencode.v46.annotation.gtf.gz"
    model = ANNOTATION / "DeepVariant/Illumina/bioskryb-af-20241102"
    return {
        "alignment": python_stage("step1_submit_resolveome_alignment_IluminaOnly_lsf.py",
                                 "--cells", args.cells, "--outdir", root, "--ref", args.ref,
                                 "--queue", args.queue, "--threads", args.threads, "--env", "resolveome-align"),
        "variants": python_stage("step3_submit_variants_cnv_sv_annot_lsf.py",
                                "--outdir", root, "--ref", args.ref, "--env", "resolveome-align",
                                "--queue", args.queue, "--threads", args.threads,
                                "--snv-tools", "deepvariant", "--cnv-tools", "cnvkit",
                                "--sv-tools", "manta", "--annotators", "snpeff,vep", "--chroms", "all",
                                "--deepvariant-sif", TOOLS / "containers/deepvariant.sif",
                                "--deepvariant-model-dir", model / "model",
                                "--deepvariant-checkpoint-filename", "sc_af.model.ckpt",
                                "--deepvariant-population-vcfs", model / "population_vcfs_nochr",
                                "--dv-use-allele-frequency",
                                "--manta-sif", TOOLS / "containers/manta.sif",
                                "--cnv-gene-bed", TOOLS / "refs/annotation/genes.bed",
                                "--cnv-reference-cnn", TOOLS / "step3_bundle/resources/cnvkit/reference.cnn",
                                "--snpeff-db", "GRCh38.105", "--vep-sif", TOOLS / "containers/vep.r112.sif",
                                "--vep-cache-dir", ANNOTATION / "VEP", "--vep-fasta", args.ref,
                                "--vep-cache-version", "111", "--run-multiqc-if-missing"),
        "ginkgo": python_stage("step7_ginkgo_cnv.py", "--manifest", passed,
                              "--outdir", root, "--ref", args.ref, "--queue", args.queue),
        "coverage": python_stage("step4_submit_coverage_qc_lsf.py", "--manifest", passed,
                                "--manifest-base-dir", root, "--gtf", gtf,
                                "--output-dir", root / "coverage_qc", "--logs-dir", root / "logs_step4",
                                "--chroms", CHROMS + ",MT", "--strip-chr-prefix", "--threshold", "3",
                                "--queue", args.queue, "--threads", "1", "--mem-mb", "8000",
                                "--wall", "145:00", "--max-concurrent", "8",
                                "--plot-script", HERE / "step4_plot_coverage_regions.py",
                                "--plot-output-prefix", root / "coverage_qc/whole_genome_coverage_map",
                                "--plot-env", "resolveome-align"),
        "somatic_inputs": python_stage("step5_filter_variant_inputs.py", "--manifest", passed,
                                      "--manifest-base-dir", root, "--variants-dir", root / "variants",
                                      "--sv-dir", root / "sv", "--output-dir", root / "somatic_inputs",
                                      "--output-csv", root / "somatic_inputs/basej_somatic_inputs.csv",
                                      "--small-variant-method", "reference", "--chroms", "all",
                                      "--group", "somatic_cells"),
    }


def perform(args):
    root = args.outdir
    root.mkdir(parents=True, exist_ok=True)
    (root / "scripts").mkdir(exist_ok=True)
    (root / "logs_driver").mkdir(exist_ok=True)
    state_path = root / "workflow_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    fingerprint = {
        "cells_sha256": hashlib.sha256(args.cells.read_bytes()).hexdigest(),
        "ref": args.ref, "gender": args.gender, "protocol": "parity_v2",
        "code_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(HERE.glob("step*.py"))
        },
    }
    if state.get("fingerprint", fingerprint) != fingerprint:
        raise ValueError("Existing run has different inputs/protocol; use a new output root")
    if not state and any((root / "bam").glob("*/*.bam")):
        raise ValueError("Existing untracked BAMs: use a new output root")
    state["fingerprint"] = fingerprint

    def persist():
        temp = state_path.with_suffix(".tmp")
        temp.write_text(json.dumps(state, indent=2) + "\n")
        temp.replace(state_path)

    def stage(name, action):
        status = state.get(name, {}).get("status")
        if status == "DONE":
            print(f"[RESUME] {name}", flush=True)
            return
        if status == "STARTED":
            raise RuntimeError(f"Stage {name} was interrupted. Inspect its job manifest/logs; do not resubmit blindly")
        state[name] = {"status": "STARTED", "started_utc": datetime.now(timezone.utc).isoformat()}
        persist()
        action()
        state[name]["status"] = "DONE"
        state[name]["finished_utc"] = datetime.now(timezone.utc).isoformat()
        persist()

    def submitted(name, command, manifest=None):
        with (root / "logs_driver" / f"{name}.submission.log").open("w") as log:
            result = run(command, stdout=log, stderr=subprocess.STDOUT)
        del result
        if manifest:
            ids = [row["job_id"] for row in read_rows(manifest)]
        else:
            text = (root / "logs_driver" / f"{name}.submission.log").read_text()
            ids = re.findall(r"Job <(\d+)>", text)
        state[name]["job_ids"] = ids
        persist()
        wait_jobs(ids)

    cmd = commands(args)
    stage("alignment", lambda: submitted("alignment", cmd["alignment"]))
    pointer = root / "manifests/aligned_bams_latest.tsv"
    require_file(pointer)
    aligned = Path(pointer.read_text().strip())
    stage("cell_qc", lambda: run(python_stage(
        "step2_submit_resolveome_barcode_cellcall_IluminaOnly_lsf.py",
        "--aligned-manifest", aligned, "--outdir", root, "--local")))
    cohort = read_rows(aligned)
    expected = {row["sample_id"] for row in read_rows(args.cells)}
    if len(cohort) != len(expected) or {row["sample_id"] for row in cohort} != expected:
        raise RuntimeError("Alignment manifest does not contain exactly the input cohort")
    for row in cohort:
        require_file(root / "bam" / row["sample_id"] / "alignment.complete")
        run(["samtools", "quickcheck", row["bam"]], capture_output=True)
    passing = read_rows(root / "manifests/passed_cells_only_manifest.tsv")
    if not passing:
        raise RuntimeError("No passing cells; no variants/CNV/lineage submitted")
    stage("variants", lambda: submitted("variants", cmd["variants"],
                                       root / "logs_step3/step3_job_manifest.tsv"))
    stage("ginkgo", lambda: submitted("ginkgo", cmd["ginkgo"], root / "ginkgo/job_manifest.tsv"))
    stage("coverage", lambda: submitted("coverage", cmd["coverage"]))

    def combined_cnv():
        for size, directory in ((1000000, "combined"), (200000, "combined_200kb")):
            prefix = root / "cnv" / directory / "combined_single_cell_cnv"
            run(python_stage("step4_combine_single_cell_cnv.py", "--cnv-dir", root / "cnv",
                             "--reference-fai", args.ref + ".fai", "--output-dir", prefix.parent,
                             "--bin-size", size))
            run(python_stage("step4_annotate_combined_cnv_genes.py", "--input-prefix", prefix,
                             "--gtf", TOOLS / "refs/annotation/gencode.v46.annotation.gtf.gz",
                             "--bin-size", size, "--recurrent-min-fraction", "0.20"))
    stage("combined_cnv", combined_cnv)
    stage("library_dashboard", lambda: run(python_stage(
        "step2_library_qc_dashboard.py", "--outdir", root, "--ref", args.ref)))
    stage("somatic_inputs", lambda: run(cmd["somatic_inputs"]))
    output = root / "basej_somatic"
    index = output / "workflow_outputs/demo/TEST_WORKFLOW_ID/index"

    def somatic():
        config = write_nf_config(root, "somatic", args.queue)
        run([TOOLS / "tools/nextflow", "-c", config, "run", SOURCE / "basej-somatic/main.nf",
             "-profile", "singularity", "--platform", "Illumina",
             "--input_csv", root / "somatic_inputs/basej_somatic_inputs.csv",
             "--genomes_base", TOOLS / "basej_genomes_nochr", "--genome", "GRCh38",
             "--chrs", CHROMS, "--gender", args.gender, "--architecture", "x86",
             "--model_vcf", "deepvariant", "--seq_mode", "wgs",
             "--vep_cache_dir", ANNOTATION / "VEP", "--outputDir", output,
             "--run_phylogeny_mutsig", "false", "-work-dir", root / "work_somatic",
             "-with-trace", root / "logs_driver/somatic_trace.tsv", "-resume"])
    stage("somatic", somatic)
    stage("integrated_matrices", lambda: run(python_stage(
        "step5_integrate_sv_lineage.py", "--lineage-inputs", index / "lineage_inputs.csv",
        "--sv-matrix-dir", root / "somatic_inputs/structural_variants",
        "--output-dir", index / "sv_integrated", "--group", "somatic_cells")))

    def lineage():
        # Source lineage consumes small variants only. Symbolic SVs are not indels.
        connector = read_rows(index / "lineage_inputs.csv", ",")
        nr = Path(next(row["path"] for row in connector if row["param"] == "nr_matrix"))
        if len(nr.read_text().splitlines()) <= 1:
            (root / "basej_lineage").mkdir(exist_ok=True)
            (root / "basej_lineage/status.json").write_text(json.dumps({
                "status": "NOT_COMPUTABLE", "reason": "No small variants passed; cannot infer a source lineage",
            }, indent=2) + "\n")
            print("[NOTE] No passing small variants: source lineage is explicitly not computable")
            return
        config = write_nf_config(root, "lineage", args.queue)
        reduced = root / "scripts/lineage_small_variants.csv"
        with reduced.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=["group", "param", "path"])
            writer.writeheader()
            writer.writerows(row for row in connector if row["param"] != "binary_matrix")
        run([TOOLS / "tools/nextflow", "-c", config, "run", SOURCE / "basej-lineage/main.nf",
             "-profile", "singularity", "--lineage_csv", reduced,
             "--cohort_id", "somatic_cells", "--gender", args.gender,
             "--sequoia_phylogeny_mode", "snv,indel,both",
             "--sequoia_mpboot_path", "/usr/local/bin/scripts/", "--sequoia_mpboot_bootstrap", "1000",
             "--outputDir", root / "basej_lineage", "-work-dir", root / "work_lineage",
             "-with-trace", root / "logs_driver/lineage_trace.tsv", "-resume"])
    stage("source_lineage", lineage)
    for prevalence, name in ((1, "sv_similarity"), (2, "sv_similarity_shared")):
        def sv_analysis(p=prevalence, n=name):
            features = read_rows(root / "somatic_inputs/structural_variants/sv_features.tsv")
            if not any(int(row["supporting_cells"]) >= p for row in features):
                directory = root / n
                directory.mkdir(exist_ok=True)
                (directory / "status.json").write_text(json.dumps({
                    "status": "NOT_COMPUTABLE", "reason": f"No SV features detected in at least {p} cells",
                }, indent=2) + "\n")
                return
            run(python_stage("step6_sv_similarity.py",
                            "--sv-matrix-dir", root / "somatic_inputs/structural_variants",
                            "--manifest", root / "manifests/passed_cells_only_manifest.tsv",
                            "--output-dir", root / n, "--min-prevalence", p))
        stage(name, sv_analysis)
    stage("final_dashboard", lambda: final_dashboard(root, state))


def final_dashboard(root, state):
    links = [
        ("Library QC and source thresholds", "qc/library_dashboard/source_threshold_dashboard.html"),
        ("FastQC / fastp / Picard MultiQC", "qc/library_dashboard/multiqc_library_qc.html"),
        ("SV similarity — all features", "sv_similarity/sv_similarity.report.html"),
        ("SV sensitivity — shared features only", "sv_similarity_shared/sv_similarity.report.html"),
    ]
    for _, relative in links:
        if not (root / relative).is_file():
            status = (root / relative).parent / "status.json"
            require_file(status)
    body = "<!doctype html><meta charset='utf-8'><title>ResolveOME comparison</title><h1>ResolveOME revised workflow</h1>"
    body += "<p>Genome-wide small variants; unchanged DeepVariant model; source somatic gates. SV clustering is exploratory, not ancestry.</p><ul>"
    for label, path in links:
        target = path if (root / path).is_file() else str(Path(path).parent / "status.json")
        body += f"<li><a href='{target}'>{html.escape(label)}</a></li>"
    for directory in ("basej_somatic", "basej_lineage", "cnv", "ginkgo", "coverage_qc"):
        body += f"<li><a href='{directory}/'>{directory} outputs</a></li>"
    completed = {key: value for key, value in state.items() if key != "final_dashboard"}
    body += "</ul><h2>Stages</h2><pre>" + html.escape(json.dumps(completed, indent=2)) + "</pre>"
    (root / "workflow_dashboard.html").write_text(body)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path, default=PROJECT / "results/dna_cells.tsv")
    parser.add_argument("--outdir", type=Path, default=PROJECT / "results/DNA_Illumina_parity_v2")
    parser.add_argument("--ref", default=str(TOOLS / "refs/genome/GRCh38.primary.fa"))
    parser.add_argument("--queue", default="standard")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--gender", choices=("male", "female"), default="male")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight", action="store_true", help="Verify tools/resources only (default)")
    mode.add_argument("--plan", action="store_true", help="Print stage commands without writing/submitting")
    mode.add_argument("--run", action="store_true", help="Execute stages and wait; use the LSF shell launcher")
    args = parser.parse_args(argv)
    args.cells = args.cells.absolute()
    args.outdir = args.outdir.absolute()
    args.ref = str(Path(args.ref).absolute())
    protected = (PROJECT / "results/DNA_Ilumina_output_1").resolve()
    resolved_outdir = args.outdir.resolve()
    if resolved_outdir == protected or protected in resolved_outdir.parents:
        parser.error("The previous validated results are protected; choose a new output root")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.queue):
        parser.error("Unsafe queue name")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.plan:
        for name, command in commands(args).items():
            print(name + ": " + shlex.join(map(str, command)))
        print("Remaining stages: cell QC, combined CNV, library dashboard, BaseJ somatic, "
              "integrated matrices, source lineage, both SV sensitivity analyses, final dashboard")
        return
    report = preflight(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "tool_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    print("[OK] Tool and reference preflight")
    if args.run:
        with (args.outdir / "workflow.lock").open("w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            perform(args)


if __name__ == "__main__":
    main()
