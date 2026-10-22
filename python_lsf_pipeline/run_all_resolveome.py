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
EXTRA = HERE.parent / "extra_methods"  # extra_methods/ lives at workflow_DNA root alongside resolveome_align/
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


def bjobs_status(job_ids):
    """Batched `bjobs` lookup (100 IDs/call, avoids an oversized argv) returning
    {job_id: LSF status string} for every id in job_ids. Raises if any id's
    status cannot be recovered at all (distinct from the id having an EXIT
    status, which is a normal, representable outcome here)."""
    ids = sorted({str(jid) for jid in job_ids})
    seen = {}
    for offset in range(0, len(ids), 100):
        batch = ids[offset:offset + 100]
        result = run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *batch],
                     capture_output=True)
        for line in result.stdout.splitlines():
            fields = line.split()
            if len(fields) != 2 or not fields[0].isdigit():
                raise RuntimeError(f"Unrecognized LSF status output: {line}")
            seen[fields[0]] = fields[1]
        missing = set(batch) - set(seen)
        if missing:
            raise RuntimeError(f"Cannot recover job status for {missing}; inspect LSF history before resuming")
    return seen


def wait_jobs(job_ids, poll_seconds=60):
    remaining = set(job_ids)
    if not remaining or any(not str(jid).isdigit() for jid in remaining):
        raise ValueError("No valid submitted LSF job IDs")
    while remaining:
        statuses = bjobs_status(remaining)
        pending = set()
        for jid, status in statuses.items():
            if status == "DONE":
                continue
            if status not in {"PEND", "RUN", "WAIT"}:
                raise RuntimeError(f"LSF job {jid}: {status}; dependent stages not started")
            pending.add(jid)
        remaining = pending
        if remaining:
            print(f"[WAIT] {len(remaining)} jobs outstanding", flush=True)
            time.sleep(poll_seconds)


def poll_multi_stage(pending, poll_seconds=60):
    """Tolerant, interleaved multi-stage job poller.

    pending: a SHARED, MUTABLE {stage_name: [job_ids]} dict for two or more
    LSF-job-backed stages submitted concurrently (e.g. the DeepVariant/VEP/
    SnpEff track and the CNVkit/Manta track). Every poll cycle checks ALL
    outstanding stages' jobs together in one batched bjobs sweep; the moment
    any single stage's *own* job set is fully terminal (every job DONE or
    EXIT), that stage is yielded immediately as (stage_name, failed_job_ids)
    and dropped from the active set -- the caller can act on it (e.g. fire
    that stage's own downstream dependents) without waiting for the other,
    still-running stage(s). Never raises on a job EXIT: unlike wait_jobs(),
    a failure is reported, not treated as fatal, so one stage's chunk
    failures can never block recognition of another concurrently-running
    stage's success.

    The caller MAY mutate `pending` (add brand-new stage_name -> job_ids
    entries) from inside its own per-yield handling code -- e.g. when one
    stage's resolution triggers submission of a newly-dependent LSF-job-
    backed stage (such as Beagle phasing firing only once DeepVariant's raw
    VCF exists). This generator re-scans `pending` for not-yet-seen entries
    at the top of every cycle, so a stage added between yields is picked up
    on its very next poll without the caller needing a second, separate
    polling loop.
    """
    seen = {}
    while True:
        for name, ids in list(pending.items()):
            if name in seen:
                continue
            ids = {str(jid) for jid in ids}
            if not ids or any(not jid.isdigit() for jid in ids):
                raise ValueError(f"No valid submitted LSF job IDs for stage {name!r}")
            seen[name] = ids
        if not pending:
            return
        statuses = bjobs_status(set().union(*seen.values()))
        newly_resolved = []
        for name in list(pending):
            ids = seen.get(name)
            if ids is None:
                continue
            these = {statuses[jid] for jid in ids}
            if these <= {"DONE", "EXIT"}:
                failed = sorted(jid for jid in ids if statuses[jid] == "EXIT")
                newly_resolved.append((name, failed))
        for name, failed in newly_resolved:
            del pending[name]
            del seen[name]
            yield name, failed
        if pending:
            outstanding = sum(len(ids) for ids in seen.values())
            print(f"[WAIT] {outstanding} jobs outstanding across {sorted(pending)}", flush=True)
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
    // basej-somatic's own nextflow.config gives SORT_INDEX_VEP_VCF only a
    // 1h time budget with no retry (unlike its sibling VEP_ANNOTATE, which
    // produces the same per-group VCF with a 4h budget) -- too tight for
    // some chromosomes' pooled multi-cell VCFs (confirmed: a real
    // "somatic_cells_X" run exceeded 1h and was killed by LSF, exit 140,
    // with no retry since this config's own `errorStrategy = 'terminate'`
    // wildcard above applies project-wide by design). Raised to match
    // VEP_ANNOTATE's own 4h budget for the same data.
    withName: 'SORT_INDEX_VEP_VCF' {{
        time = 4.h
    }}
    // VARIANT_PREFILTER_TABLE_PERCHR occasionally hangs indefinitely (0s of
    // CPU time accruing over 45+ minutes, confirmed via bjobs on a live run,
    // on 2 of 24 per-chromosome tasks -- not a data-size/timeout issue, since
    // the same task completes in well under a minute for every other
    // chromosome including larger ones). The upstream basej-somatic
    // nextflow.config already gives this process its own `errorStrategy =
    // 'retry'` with escalating memory on attempt 2 (15GB -> 60GB), but this
    // run's own project-wide `errorStrategy = 'terminate'` wildcard above
    // was overriding that, so one hung task (after being killed to unblock
    // it) took down the entire multi-hour run instead of just retrying.
    // Restoring a dedicated retry policy here so a future hang on this
    // process is handled automatically instead of requiring a manual
    // bkill + full relaunch.
    withName: 'VARIANT_PREFILTER_TABLE_PERCHR' {{
        errorStrategy = 'retry'
        maxRetries = 2
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


def extra_stage(name, *arguments):
    """Invoke a extra_methods/ script (the Beagle/Control-FREEC/MAPD/FDR/
    artifact-panel/phasing-&-LOH-summary/dashboard scripts) with the same
    interpreter this driver itself runs under."""
    return [sys.executable, EXTRA / name, *map(str, arguments)]


def commands(args):
    root = args.outdir
    passed = root / "manifests/passed_cells_only_manifest.tsv"
    gtf = TOOLS / "refs/annotation/gencode.v46.annotation.gtf.gz"
    model = ANNOTATION / "DeepVariant/Illumina/bioskryb-af-20241102"
    em = root / "extra_methods"  # matches extra_methods/organize_final_results.py's own hardcoded layout
    return {
        "alignment": python_stage("step1_submit_resolveome_alignment_IluminaOnly_lsf.py",
                                 "--cells", args.cells, "--outdir", root, "--ref", args.ref,
                                 "--queue", args.queue, "--threads", args.threads, "--env", "resolveome-align"),
        # The original single "variants" step3 invocation (DeepVariant+VEP+
        # SnpEff+CNVkit+Manta all in one) is split into two independent
        # invocations so a DeepVariant/VEP/SnpEff failure can never block or
        # crash the CNVkit/Manta track (or vice versa). Each gets its own
        # --job-manifest-name so the two concurrently-running step3
        # processes never truncate/append the same manifest file. See
        # perform() for how these two are submitted concurrently and polled
        # independently.
        "variants_small": python_stage("step3_submit_variants_cnv_sv_annot_lsf.py",
                                "--outdir", root, "--ref", args.ref, "--env", "resolveome-align",
                                "--queue", args.queue, "--threads", args.threads,
                                "--snv-tools", "deepvariant", "--cnv-tools", "none",
                                "--sv-tools", "none", "--annotators", "snpeff,vep", "--chroms", "all",
                                "--job-manifest-name", "step3_job_manifest.small_variants.tsv",
                                "--deepvariant-sif", TOOLS / "containers/deepvariant.sif",
                                "--deepvariant-model-dir", model / "model",
                                "--deepvariant-checkpoint-filename", "sc_af.model.ckpt",
                                "--deepvariant-population-vcfs", model / "population_vcfs_nochr",
                                "--dv-use-allele-frequency",
                                "--snpeff-db", "GRCh38.105", "--vep-sif", TOOLS / "containers/vep.r112.sif",
                                "--vep-cache-dir", ANNOTATION / "VEP", "--vep-fasta", args.ref,
                                "--vep-cache-version", "111"),
                                # --run-multiqc-if-missing deliberately NOT passed -- see the
                                # "multiqc" stage() call in perform() below for why.
        "variants_cnv_sv": python_stage("step3_submit_variants_cnv_sv_annot_lsf.py",
                                "--outdir", root, "--ref", args.ref, "--env", "resolveome-align",
                                "--queue", args.queue, "--threads", args.threads,
                                "--snv-tools", "none", "--cnv-tools", "cnvkit",
                                "--sv-tools", "manta", "--annotators", "none", "--chroms", "all",
                                "--job-manifest-name", "step3_job_manifest.cnv_sv.tsv",
                                "--manta-sif", TOOLS / "containers/manta.sif",
                                "--cnv-gene-bed", TOOLS / "refs/annotation/genes.bed",
                                "--cnv-reference-cnn", TOOLS / "step3_bundle/resources/cnvkit/reference.cnn"),
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

        # =====================================================================
        # extra_methods/ integration: these scripts already existed (written
        # for/shared with the WDL pipeline) but were never wired into this
        # production orchestrator. They are genuinely additive -- none
        # duplicates an existing production stage -- EXCEPT lineage_track.py/
        # verify_lineage.py, which independently rebuild a whole second
        # Sequoia/MPBoot lineage tree; those are deliberately NOT wired in
        # here since source_lineage (BaseJ lineage Nextflow) already performs
        # that exact analysis and running both would duplicate expensive
        # compute for no new information (see README's methods note).
        #
        # Output layout follows extra_methods/organize_final_results.py's own
        # existing, hardcoded expectations exactly (<outdir>/extra_methods/...)
        # so that script -- and build_extra_dashboard.py, which reads cohort
        # summary TSVs directly from that same directory -- work unmodified.
        # =====================================================================
        "loh_freec_fixed": extra_stage("loh_freec_lsf.py", "--bam-dir", root / "bam", "--manifest", passed,
                                "--ref-fasta", args.ref, "--chrom-fasta-dir", EXTRA / "refs/freec_chrom_fasta",
                                "--chr-len-file", EXTRA / "refs/freec_chrLen.fai",
                                "--snp-positions-vcf", EXTRA / "refs/freec_snp_positions/all_chrom.sites.vcf.gz",
                                "--freec-bin", "freec", "--sex", "XY" if args.gender == "male" else "XX",
                                "--ploidy", "2",
                                "--output-dir", em / "freec_loh", "--logs-dir", root / "logs_loh_freec_fixed",
                                "--queue", args.queue),
        "loh_freec_candidate": extra_stage("loh_freec_lsf.py", "--bam-dir", root / "bam", "--manifest", passed,
                                "--ref-fasta", args.ref, "--chrom-fasta-dir", EXTRA / "refs/freec_chrom_fasta",
                                "--chr-len-file", EXTRA / "refs/freec_chrLen.fai",
                                "--snp-positions-vcf", EXTRA / "refs/freec_snp_positions/all_chrom.sites.vcf.gz",
                                "--freec-bin", "freec", "--sex", "XY" if args.gender == "male" else "XX",
                                "--ploidy", "2,3,4,5,6,7,8",
                                "--output-dir", em / "freec_loh_ploidy", "--logs-dir", root / "logs_loh_freec_candidate",
                                "--queue", args.queue),
        "phase_variants": extra_stage("phase_variants_lsf.py", "--variants-dir", root / "variants",
                                     "--manifest", passed, "--beagle-jar", EXTRA / "tools/beagle/beagle.27Feb25.75f.jar",
                                     "--ref-panel-dir", EXTRA / "refs/phasing_panel_grch38",
                                     "--genetic-map-dir", EXTRA / "refs/genetic_maps",
                                     "--output-dir", em / "phasing", "--logs-dir", root / "logs_phasing",
                                     "--queue", args.queue),
        "mapd_qc": extra_stage("mapd_qc.py", "--cnv-dir", root / "cnv", "--manifest", passed,
                              "--output", em / "mapd_qc.tsv"),
        "cnv_overview": extra_stage("cnv_overview.py", "--cnv-dir", root / "cnv", "--manifest", passed,
                                   "--output", em / "cnv_overview.tsv"),
        "ploidy_determination": extra_stage("ploidy_determination.py",
                                           "--freec-logs-dir", root / "logs_loh_freec_candidate",
                                           "--manifest", passed,
                                           "--output-candidates", em / "ploidy_candidates.tsv",
                                           "--output-narrative", em / "ploidy_narrative.tsv"),
        "fdr_estimate": extra_stage("fdr_estimate.py", "--variants-dir", root / "variants",
                                   "--picard-dir", root / "qc/picard", "--manifest", passed,
                                   "--output", em / "fdr_estimate.tsv"),
        "artifact_panel": extra_stage("artifact_panel.py", "--variants-dir", root / "variants",
                                     "--manifest", passed,
                                     "--output-sites", em / "artifact_panel_sites.tsv",
                                     "--output-summary", em / "artifact_panel_summary.tsv"),
        "phasing_summary": extra_stage("phasing_summary.py", "--phasing-dir", em / "phasing",
                                      "--manifest", passed, "--output", em / "phasing_summary.tsv"),
        # loh_summary/ploidy_summary/loh_by_region all read the CANDIDATE-
        # ploidy FREEC run (freec_loh_ploidy), not the fixed-ploidy=2 run --
        # confirmed directly from each script's own source.
        "loh_summary": extra_stage("loh_summary.py", "--freec-dir", em / "freec_loh_ploidy",
                                  "--manifest", passed, "--output", em / "loh_summary.tsv"),
        "ploidy_summary": extra_stage("ploidy_summary.py", "--freec-dir", em / "freec_loh_ploidy",
                                     "--manifest", passed, "--output", em / "ploidy_summary.tsv"),
        "loh_by_region": extra_stage("loh_by_region.py", "--freec-dir", em / "freec_loh_ploidy",
                                    "--manifest", passed, "--chrom-len-file", EXTRA / "refs/freec_chrLen.fai",
                                    "--output-tsv", em / "loh_by_region.tsv",
                                    "--output-png", em / "loh_by_region_heatmap.png"),
        "qc_summary": extra_stage("qc_summary.py",
                                 "--multiqc-general-stats",
                                 root / "multiqc_step3/multiqc_data/multiqc_general_stats.txt",
                                 "--manifest", passed,
                                 "--output-tsv", em / "qc_summary.tsv",
                                 "--output-png", em / "qc_overview.png"),
        # Bare "Rscript" would resolve to the driver process's own active
        # env (resolveome-align, inherited from this whole orchestrator's
        # `conda run --prefix resolveome-align` launch), which has no
        # MutationalPatterns/BSgenome/VariantAnnotation installed -- those
        # packages only exist in the dedicated "mutsig" conda env (see
        # setup.sh's mutsig.yml). Must invoke via `conda run -n mutsig`
        # explicitly; confirmed this works correctly even nested inside an
        # already-conda-run'd resolveome-align process.
        "mutational_signatures": ["conda", "run", "-n", "mutsig", "Rscript",
                                 str(EXTRA / "mutational_signatures.R"),
                                 str(em / "artifact_panel_sites.tsv"), str(em / "mutational_signatures")],
        "build_extra_dashboard": extra_stage("build_extra_dashboard.py",
                                            "--extra-methods-dir", em,
                                            "--donor-id", root.name,
                                            "--multiqc-relpath", "../multiqc_step3/multiqc_report.html",
                                            "--output", em / "extra_dashboard.html"),
        "organize_final_results": extra_stage("organize_final_results.py",
                                             "--dna-results-dir", root, "--final-dir", root / "final"),
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

    def stage(name, action, tolerant=False):
        status = state.get(name, {}).get("status")
        resolved_statuses = ("DONE", "FAILED") if tolerant else ("DONE",)
        if status in resolved_statuses:
            print(f"[RESUME] {name} (previously {status})", flush=True)
            return
        if status == "STARTED":
            raise RuntimeError(f"Stage {name} was interrupted. Inspect its job manifest/logs; do not resubmit blindly")
        state[name] = {"status": "STARTED", "started_utc": datetime.now(timezone.utc).isoformat()}
        persist()
        if tolerant:
            try:
                action()
                state[name]["status"] = "DONE"
            except Exception as exc:  # noqa: BLE001 -- deliberately broad: this is the
                # resilience boundary the user asked for. A failure here (e.g. a
                # downstream analysis stage erroring because an upstream track it
                # consumes failed/was skipped) must be recorded and must NOT crash
                # the driver or block sibling stages; it is surfaced clearly in
                # workflow_state.json and on stdout instead.
                state[name]["status"] = "FAILED"
                state[name]["error"] = str(exc)
                print(f"[FAILED] {name}: {exc}", flush=True)
        else:
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

    def submit_only(name, command, manifest=None, manifest_format="tsv"):
        """Fire an LSF-submission-driving command and record its job IDs
        without waiting -- the counterpart of submitted(), split so several
        independent stages can be launched back-to-back before any of them
        is polled to completion. Pairs with poll_multi_stage()/
        finalize_stage() below.

        manifest_format="tsv" expects a TSV with a "job_id" column (step3's
        own manifests, ginkgo's job_manifest.tsv). manifest_format="lines"
        expects a plain newline-delimited job-ID file (phase_job_ids.txt/
        freec_job_ids.txt, as written by extra_methods/phase_variants_lsf.py
        and extra_methods/loh_freec_lsf.py)."""
        with (root / "logs_driver" / f"{name}.submission.log").open("w") as log:
            run(command, stdout=log, stderr=subprocess.STDOUT)
        if manifest and manifest_format == "lines":
            require_file(manifest)
            ids = [line.strip() for line in Path(manifest).read_text().splitlines() if line.strip()]
        elif manifest:
            ids = [row["job_id"] for row in read_rows(manifest)]
        else:
            text = (root / "logs_driver" / f"{name}.submission.log").read_text()
            ids = re.findall(r"Job <(\d+)>", text)
        state[name] = {"status": "STARTED", "started_utc": datetime.now(timezone.utc).isoformat(), "job_ids": ids}
        persist()
        return ids

    def finalize_stage(name, failed_ids):
        state[name]["status"] = "FAILED_PARTIAL" if failed_ids else "DONE"
        state[name]["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if failed_ids:
            state[name]["failed_job_ids"] = failed_ids
            print(f"[FAILED_PARTIAL] {name}: {len(failed_ids)} job(s) EXIT -- inspect "
                  f"logs_step3/ for these IDs: {failed_ids}", flush=True)
        else:
            print(f"[DONE] {name}", flush=True)
        persist()

    # Backward-compatibility migration: older runs (started before the single
    # "variants" stage -- DeepVariant+VEP+SnpEff+CNVkit+Manta submitted
    # together -- was split into independently-resolved "variants_small" and
    # "variants_cnv_sv" tracks) recorded completion only under the old name.
    # Without this shim, resuming such a run under the new code would see
    # neither new stage name present and resubmit the entire, already-
    # completed ~10,000-job variants population from scratch.
    if state.get("variants", {}).get("status") == "DONE":
        for migrated in ("variants_small", "variants_cnv_sv"):
            if migrated not in state:
                state[migrated] = {
                    "status": "DONE",
                    "migrated_from": "variants",
                    "note": ("Migrated from the pre-split single 'variants' stage, which "
                             "already completed successfully for both tracks combined "
                             "before the orchestrator was rewritten to submit/resolve the "
                             "small-variant and CNV/SV tracks independently."),
                    "job_ids": state["variants"].get("job_ids", []),
                }
        persist()

    cmd = commands(args)
    em = root / "extra_methods"
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

    # library_dashboard only scans qc/picard, qc/fastp, and qc/fastqc -- all
    # produced by alignment/cell_qc above -- so it has no real dependency on
    # any Tier-1 stage below. Run it now rather than waiting behind them.
    stage("library_dashboard", lambda: run(python_stage(
        "step2_library_qc_dashboard.py", "--outdir", root, "--ref", args.ref)), tolerant=True)

    # =========================================================================
    # Tier 1: four mutually-independent LSF-job-backed stages.
    #
    #   variants_small   -- DeepVariant + VEP + SnpEff (per-cell small variants)
    #   variants_cnv_sv  -- CNVkit + Manta (per-cell CNV/SV)
    #   ginkgo           -- Ginkgo CNV (BAM-only)
    #   coverage         -- coverage QC (BAM-only)
    #
    # None of these needs another's output, so all four are submitted
    # concurrently (no waiting between submissions) and polled together via
    # poll_multi_stage(): the moment any ONE stage's own job set is fully
    # terminal, it is resolved (DONE, or FAILED_PARTIAL if any chunk EXITed)
    # and its own Tier-2 dependent(s) fire immediately -- a DeepVariant/VEP/
    # SnpEff failure can never block or crash the CNVkit/Manta/Ginkgo/
    # coverage tracks (or vice versa).
    # =========================================================================
    tier1_specs = {
        "variants_small": (cmd["variants_small"], root / "logs_step3/step3_job_manifest.small_variants.tsv", "tsv"),
        "variants_cnv_sv": (cmd["variants_cnv_sv"], root / "logs_step3/step3_job_manifest.cnv_sv.tsv", "tsv"),
        "ginkgo": (cmd["ginkgo"], root / "ginkgo/job_manifest.tsv", "tsv"),
        "coverage": (cmd["coverage"], None, "tsv"),
        # Control-FREEC is BAM-only (same as ginkgo/coverage): both the
        # fixed-ploidy=2 and candidate-ploidy (2,3,4,5,6,7,8) runs are genuinely
        # Tier-1 eligible, independent of DeepVariant/CNVkit/Manta.
        "loh_freec_fixed": (cmd["loh_freec_fixed"], em / "freec_loh/freec_job_ids.txt", "lines"),
        "loh_freec_candidate": (cmd["loh_freec_candidate"], em / "freec_loh_ploidy/freec_job_ids.txt", "lines"),
    }
    tier1_resolved = {name for name in tier1_specs if state.get(name, {}).get("status") in ("DONE", "FAILED_PARTIAL")}
    tier1_outstanding = {}
    for name, (command, manifest, manifest_format) in tier1_specs.items():
        if name in tier1_resolved:
            print(f"[RESUME] {name} (previously {state[name]['status']})", flush=True)
            continue
        status = state.get(name, {}).get("status")
        if status == "STARTED":
            print(f"[RESUME] {name}: rejoining {len(state[name]['job_ids'])} already-submitted jobs", flush=True)
            tier1_outstanding[name] = state[name]["job_ids"]
            continue
        tier1_outstanding[name] = submit_only(name, command, manifest, manifest_format)

    tier2_triggered = set()

    def combined_cnv():
        for size, directory in ((1000000, "combined"), (200000, "combined_200kb")):
            prefix = root / "cnv" / directory / "combined_single_cell_cnv"
            run(python_stage("step4_combine_single_cell_cnv.py", "--cnv-dir", root / "cnv",
                             "--reference-fai", args.ref + ".fai", "--output-dir", prefix.parent,
                             "--bin-size", size))
            run(python_stage("step4_annotate_combined_cnv_genes.py", "--input-prefix", prefix,
                             "--gtf", TOOLS / "refs/annotation/gencode.v46.annotation.gtf.gz",
                             "--bin-size", size, "--recurrent-min-fraction", "0.20"))

    def cell_qc_gate_all():
        # BaseJumper-derived composite gate (preseq/chimeric/chrM/MAPD/skew),
        # run per cell against CNVkit's ALREADY-PRODUCED .cnr/.call.cns from
        # variants_cnv_sv -- no separate lightweight CNVkit pass is needed
        # (the WDL path runs one only because it gates BEFORE its own full
        # CNV stage exists yet; here variants_cnv_sv has already finished).
        # This is informational/additive: it reports PASS/Borderline/FAIL
        # per cell for the dashboard, but does NOT remove any cell from
        # downstream processing -- the existing flagstat-based cell_qc gate
        # remains the sole actual inclusion/exclusion gate, so adding this
        # richer, complementary report carries no risk of silently changing
        # which cells proceed.
        gate_dir = em / "cell_qc_gate"
        gate_dir.mkdir(parents=True, exist_ok=True)
        for row in read_rows(root / "manifests/passed_cells_only_manifest.tsv"):
            sample = row["sample_id"]
            picard_prefix = root / "qc/picard" / sample / sample
            cnv_dir = root / "cnv" / sample
            idxstats = root / "bam" / sample / f"{sample}.dedup.idxstats.txt"
            preseq_file = root / "qc/preseq" / sample / f"{sample}_preseq.txt"
            inputs = {
                "idxstats": idxstats,
                "alignment-summary-metrics": Path(f"{picard_prefix}.alignment_summary_metrics"),
                "preseq-file": preseq_file,
                "cnr": cnv_dir / f"{sample}.dedup.cnr",
                "call-cns": cnv_dir / f"{sample}.dedup.call.cns",
            }
            missing = [k for k, v in inputs.items() if not v.is_file()]
            if missing:
                print(f"[SKIP] cell_qc_gate for {sample}: missing {missing}", flush=True)
                continue
            run(extra_stage("cell_qc_gate.py", "--sample-id", sample,
                            "--idxstats", inputs["idxstats"],
                            "--alignment-summary-metrics", inputs["alignment-summary-metrics"],
                            "--preseq-file", inputs["preseq-file"],
                            "--cnr", inputs["cnr"], "--call-cns", inputs["call-cns"],
                            "--output-tsv", gate_dir / f"{sample}.qc_gate.tsv",
                            "--output-explanation", gate_dir / f"{sample}.qc_gate_explanation.txt",
                            "--output-pass-flag", gate_dir / f"{sample}.qc_gate_status.txt"))

    def multiqc():
        # Run MultiQC directly here (no bsub, no LSF dependency expression)
        # instead of via step3's own --run-multiqc-if-missing: by this point
        # every Tier-1 stage's own job set has reached a terminal state --
        # the same guarantee an LSF -w dependency would have provided, with
        # no risk of ever building an oversized dependency expression
        # (step3's own `--run-multiqc-if-missing` path built one job-ID-per-
        # prerequisite expression that, at this cohort's scale of ~10,447
        # jobs, was ~194 KB -- over Linux's 128 KB single-argument limit,
        # and crashed this orchestrator's own subprocess call with
        # `OSError: Argument list too long` on 2026-10-20). This conda env
        # (resolveome-align) is already active in this process, so
        # `multiqc` is directly callable.
        multiqc_dir = root / "multiqc_step3"
        if multiqc_dir.is_dir():
            print("[SKIP] multiqc_step3/ already exists", flush=True)
            return
        run(["multiqc", "-f", "-o", multiqc_dir, root])

    def on_tier1_resolved(name):
        tier1_resolved.add(name)
        if name == "variants_cnv_sv":
            if "combined_cnv" not in tier2_triggered:
                tier2_triggered.add("combined_cnv")
                stage("combined_cnv", combined_cnv, tolerant=True)
            if "cell_qc_gate" not in tier2_triggered:
                tier2_triggered.add("cell_qc_gate")
                stage("cell_qc_gate", cell_qc_gate_all, tolerant=True)
        if name == "variants_small" and "phase_variants" not in tier2_triggered:
            tier2_triggered.add("phase_variants")
            status = state.get("phase_variants", {}).get("status")
            if status in ("DONE", "FAILED_PARTIAL"):
                print(f"[RESUME] phase_variants (previously {status})", flush=True)
                on_tier1_resolved("phase_variants")
            elif status == "STARTED":
                tier1_outstanding["phase_variants"] = state["phase_variants"]["job_ids"]
            else:
                tier1_outstanding["phase_variants"] = submit_only(
                    "phase_variants", cmd["phase_variants"],
                    em / "phasing/phase_job_ids.txt", "lines")
        if name == "phase_variants" and "phasing_summary" not in tier2_triggered:
            tier2_triggered.add("phasing_summary")
            stage("phasing_summary", lambda: run(cmd["phasing_summary"]), tolerant=True)
        if name == "loh_freec_candidate" and "loh_summaries" not in tier2_triggered:
            tier2_triggered.add("loh_summaries")
            stage("loh_summary", lambda: run(cmd["loh_summary"]), tolerant=True)
            stage("ploidy_summary", lambda: run(cmd["ploidy_summary"]), tolerant=True)
            stage("loh_by_region", lambda: run(cmd["loh_by_region"]), tolerant=True)
            stage("ploidy_determination", lambda: run(cmd["ploidy_determination"]), tolerant=True)
        if tier1_resolved == set(tier1_specs) | {"phase_variants"} and "multiqc" not in tier2_triggered:
            tier2_triggered.add("multiqc")
            stage("multiqc", multiqc, tolerant=True)

    # Any Tier-1 stage already resolved from a prior driver invocation needs
    # its Tier-2 dependent(s) checked/fired now, before polling the rest.
    for name in sorted(tier1_resolved):
        on_tier1_resolved(name)
    for name, failed_ids in poll_multi_stage(tier1_outstanding):
        finalize_stage(name, failed_ids)
        on_tier1_resolved(name)

    # =========================================================================
    # Remaining extra_methods/ cohort-wide summary scripts: all are fast,
    # synchronous (pandas/R, seconds-to-minutes, no LSF submission), so no
    # real parallelism benefit from triggering them any earlier than "once
    # their single upstream prerequisite is terminal" -- which, by this
    # point, every Tier-1 member (including the dynamically-added
    # phase_variants) already is. Each is independently tolerant: a failure
    # in any one is recorded and does not block the others.
    # =========================================================================
    variants_small_ok = state.get("variants_small", {}).get("status") == "DONE"
    variants_cnv_sv_ok = state.get("variants_cnv_sv", {}).get("status") == "DONE"
    if variants_cnv_sv_ok:
        stage("mapd_qc", lambda: run(cmd["mapd_qc"]), tolerant=True)
        stage("cnv_overview", lambda: run(cmd["cnv_overview"]), tolerant=True)
    if variants_small_ok:
        stage("fdr_estimate", lambda: run(cmd["fdr_estimate"]), tolerant=True)
        stage("artifact_panel", lambda: run(cmd["artifact_panel"]), tolerant=True)
    if state.get("multiqc", {}).get("status") == "DONE":
        stage("qc_summary", lambda: run(cmd["qc_summary"]), tolerant=True)
    if state.get("artifact_panel", {}).get("status") == "DONE":
        stage("mutational_signatures", lambda: run(cmd["mutational_signatures"]), tolerant=True)
    stage("build_extra_dashboard", lambda: run(cmd["build_extra_dashboard"]), tolerant=True)

    # =========================================================================
    # Tier 2b / somatic_inputs onward: these stages read the RAW per-cell
    # DeepVariant VCF and Manta VCF for the FULL cohort (step5 hard-fails on
    # the first cell missing either file), so they only make sense once both
    # variant tracks have reached a terminal state. They are attempted
    # tolerantly: a failure here (e.g. because variants_small/variants_cnv_sv
    # was FAILED_PARTIAL for some cells) is recorded and does not crash the
    # driver, but the strictly-sequential stages built on top of its output
    # (somatic, integrated_matrices, source_lineage, sv_similarity*) are only
    # attempted if it actually produced its output.
    # =========================================================================
    stage("somatic_inputs", lambda: run(cmd["somatic_inputs"]), tolerant=True)
    somatic_inputs_ok = state.get("somatic_inputs", {}).get("status") == "DONE"
    output = root / "basej_somatic"
    index = output / "workflow_outputs/demo/TEST_WORKFLOW_ID/index"

    def not_computable(name, directory, reason):
        state[name] = {"status": "DONE", "note": reason}
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "status.json").write_text(json.dumps({
            "status": "NOT_COMPUTABLE", "reason": reason,
        }, indent=2) + "\n")
        print(f"[NOTE] {name}: {reason}", flush=True)
        persist()

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
    if somatic_inputs_ok:
        stage("somatic", somatic, tolerant=True)
    elif state.get("somatic", {}).get("status") not in ("DONE",):
        not_computable("somatic", output, "somatic_inputs did not complete; see its recorded error")
    somatic_ok = state.get("somatic", {}).get("status") == "DONE" and "note" not in state.get("somatic", {})

    if somatic_ok:
        stage("integrated_matrices", lambda: run(python_stage(
            "step5_integrate_sv_lineage.py", "--lineage-inputs", index / "lineage_inputs.csv",
            "--sv-matrix-dir", root / "somatic_inputs/structural_variants",
            "--output-dir", index / "sv_integrated", "--group", "somatic_cells")), tolerant=True)
    elif state.get("integrated_matrices", {}).get("status") not in ("DONE",):
        not_computable("integrated_matrices", index / "sv_integrated",
                       "somatic did not complete; cannot integrate SV/lineage matrices")
    integrated_ok = state.get("integrated_matrices", {}).get("status") == "DONE" and \
        "note" not in state.get("integrated_matrices", {})

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
    if integrated_ok:
        stage("source_lineage", lineage, tolerant=True)
    elif state.get("source_lineage", {}).get("status") not in ("DONE",):
        not_computable("source_lineage", root / "basej_lineage",
                       "integrated_matrices did not complete; cannot infer a source lineage")
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
        if somatic_inputs_ok:
            stage(name, sv_analysis, tolerant=True)
        elif state.get(name, {}).get("status") not in ("DONE",):
            not_computable(name, root / name, "somatic_inputs did not complete; cannot assess SV similarity")
    stage("final_dashboard", lambda: final_dashboard(root, state))
    # Symlink-based navigation tree over everything above; safe to attempt
    # even if some upstream stages are only partially complete (relink()
    # skips any target that doesn't exist yet rather than failing).
    stage("organize_final_results", lambda: run(cmd["organize_final_results"]), tolerant=True)


def final_dashboard(root, state):
    links = [
        ("Library QC and source thresholds", "qc/library_dashboard/source_threshold_dashboard.html", "library_dashboard"),
        ("FastQC / fastp / Picard MultiQC", "qc/library_dashboard/multiqc_library_qc.html", "library_dashboard"),
        ("SV similarity — all features", "sv_similarity/sv_similarity.report.html", "sv_similarity"),
        ("SV sensitivity — shared features only", "sv_similarity_shared/sv_similarity.report.html", "sv_similarity_shared"),
        ("Extra-methods dashboard (Beagle phasing, Control-FREEC LOH/ploidy, "
         "MAPD/FDR/artifact-panel QC, mutational signatures)",
         "extra_methods/extra_dashboard.html", "build_extra_dashboard"),
    ]
    for _, relative, owning_stage in links:
        if (root / relative).is_file():
            continue
        status = (root / relative).parent / "status.json"
        if status.is_file():
            continue
        # Unlike the original design (a hard require_file(status) crash),
        # write a fallback status.json here instead of raising: a stage that
        # feeds this link may have been attempted tolerantly (stage(...,
        # tolerant=True)) and genuinely failed or was skipped without ever
        # getting the chance to write its own status.json -- the final
        # dashboard itself must never crash just because one upstream
        # optional component did; it should report that clearly instead.
        status.parent.mkdir(parents=True, exist_ok=True)
        stage_state = state.get(owning_stage, {})
        status.write_text(json.dumps({
            "status": "NOT_COMPUTABLE",
            "reason": f"{relative} was not produced",
            "owning_stage_status": stage_state.get("status", "NOT_STARTED"),
            "owning_stage_error": stage_state.get("error"),
        }, indent=2) + "\n")
    body = "<!doctype html><meta charset='utf-8'><title>ResolveOME comparison</title><h1>ResolveOME revised workflow</h1>"
    body += "<p>Genome-wide small variants; unchanged DeepVariant model; source somatic gates. SV clustering is exploratory, not ancestry.</p><ul>"
    for label, path, _ in links:
        target = path if (root / path).is_file() else str(Path(path).parent / "status.json")
        body += f"<li><a href='{target}'>{html.escape(label)}</a></li>"
    for directory in ("basej_somatic", "basej_lineage", "cnv", "ginkgo", "coverage_qc", "extra_methods"):
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
