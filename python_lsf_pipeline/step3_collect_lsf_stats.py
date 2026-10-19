#!/usr/bin/env python3
"""
step3_collect_lsf_stats.py

Reads the Step3 manifest:
  <outdir>/logs_step3/step3_job_manifest.tsv

Queries LSF accounting for each job id (bacct -l preferred; bjobs -l fallback),
and writes:
  1) job-level stats:   <outdir>/logs_step3/step3_job_stats.tsv
  2) stage summary:     <outdir>/logs_step3/step3_stage_stats.tsv

Stage summary includes:
- per-sample stage totals
- SNV rolled up as complete stage timing (first SNV chunk start -> SNV merge end)
- cohort totals for SNV, CNV, SV

Usage:
python step3_collect_lsf_stats.py \
  --outdir /research_jude/rgs01_jude/groups/zhanggrp/projects/SJCloud/common/Zalman/Ultima_test/ResolveOME_illumina/results/DNA_Ilumina_output_1
"""

import argparse
import csv
import os
import re
import subprocess
from datetime import datetime
from typing import Dict, List, Optional, Tuple


def run_no_fail(cmd: str) -> Tuple[int, str, str]:
    p = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    return p.returncode, p.stdout, p.stderr


def parse_time(s: str) -> Optional[datetime]:
    s = s.strip()
    if not s:
        return None
    # Example format in LSF text:
    # Fri Sep  4 09:30:40 2026
    for fmt in ("%a %b %d %H:%M:%S %Y", "%a %b  %d %H:%M:%S %Y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


def fmt_dt(d: Optional[datetime]) -> str:
    return d.isoformat(sep=" ", timespec="seconds") if d else ""


def sec_to_hms(sec: Optional[int]) -> str:
    if sec is None:
        return ""
    h = sec // 3600
    m = (sec % 3600) // 60
    s = sec % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def parse_bacct_text(txt: str) -> Dict[str, str]:
    out = {
        "status": "",
        "exit_code": "",
        "start_time_raw": "",
        "end_time_raw": "",
        "cpu_time_seconds": "",
        "max_mem_mb": "",
    }

    # status / exit
    m = re.search(r"Completed <(\w+)>", txt)
    if m:
        out["status"] = m.group(1)
    m = re.search(r"Exited with exit code (\d+)", txt)
    if m:
        out["exit_code"] = m.group(1)
    elif out["status"] == "done":
        out["exit_code"] = "0"

    # start / end
    m = re.search(r"Started at\s+([A-Za-z]{3}\s+[A-Za-z]{3}\s+\d+\s+\d+:\d+:\d+\s+\d{4})", txt)
    if m:
        out["start_time_raw"] = m.group(1)
    m = re.search(r"Terminated at\s+([A-Za-z]{3}\s+[A-Za-z]{3}\s+\d+\s+\d+:\d+:\d+\s+\d{4})", txt)
    if m:
        out["end_time_raw"] = m.group(1)

    # cpu time
    m = re.search(r"CPU time\s*:\s*([0-9.]+)\s*sec", txt)
    if m:
        out["cpu_time_seconds"] = str(int(float(m.group(1))))

    # max memory (reported MB usually)
    m = re.search(r"Max Memory\s*:\s*([0-9.]+)\s*([A-Za-z]+)", txt)
    if m:
        val = float(m.group(1))
        unit = m.group(2).upper()
        if unit.startswith("GB"):
            val *= 1024
        elif unit.startswith("KB"):
            val /= 1024
        out["max_mem_mb"] = str(int(val))

    return out


def fetch_job_runtime(job_id: str) -> Dict[str, str]:
    """
    Try bacct first (best completed accounting), fallback to bjobs -l.
    Returns parsed dict.
    """
    rc, out, err = run_no_fail(f"bacct -l {job_id}")
    txt = (out or "") + "\n" + (err or "")
    if rc == 0 and ("Accounting information about this job" in txt or "CPU time" in txt):
        d = parse_bacct_text(txt)
        d["_source"] = "bacct"
        return d

    rc2, out2, err2 = run_no_fail(f"bjobs -l {job_id}")
    txt2 = (out2 or "") + "\n" + (err2 or "")
    d = parse_bacct_text(txt2)
    d["_source"] = "bjobs"
    return d


def read_manifest(path: str) -> List[Dict[str, str]]:
    rows = []
    with open(path, "r", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append(row)
    return rows


def write_tsv(path: str, rows: List[Dict[str, str]], columns: List[str]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in columns})


def main():
    ap = argparse.ArgumentParser(description="Collect LSF run stats for Step3 jobs")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    logs = os.path.join(args.outdir, "logs_step3")
    manifest = os.path.join(logs, "step3_job_manifest.tsv")
    if not os.path.isfile(manifest):
        raise SystemExit(f"[ERROR] manifest not found: {manifest}")

    rows = read_manifest(manifest)

    # -------- job-level stats --------
    job_stats = []
    for r in rows:
        jid = r["job_id"]
        rt = fetch_job_runtime(jid)

        submit_dt = parse_time(r["submission_ts"].replace("T", " "))
        start_dt = parse_time(rt.get("start_time_raw", ""))
        end_dt = parse_time(rt.get("end_time_raw", ""))

        wall_sec = ""
        if start_dt and end_dt:
            wall_sec = str(int((end_dt - start_dt).total_seconds()))

        job_stats.append(
            {
                **r,
                "lsf_source": rt.get("_source", ""),
                "status": rt.get("status", ""),
                "exit_code": rt.get("exit_code", ""),
                "submit_time": fmt_dt(submit_dt),
                "start_time": fmt_dt(start_dt),
                "end_time": fmt_dt(end_dt),
                "elapsed_wall_seconds": wall_sec,
                "elapsed_wall_hms": sec_to_hms(int(wall_sec)) if wall_sec else "",
                "cpu_time_seconds": rt.get("cpu_time_seconds", ""),
                "max_mem_mb": rt.get("max_mem_mb", ""),
            }
        )

    job_cols = [
        "sample", "stage", "job_name", "job_id", "status", "exit_code", "lsf_source",
        "queue", "requested_ncpus", "requested_mem_mb", "requested_wall",
        "submit_time", "start_time", "end_time",
        "elapsed_wall_seconds", "elapsed_wall_hms",
        "cpu_time_seconds", "max_mem_mb",
        "out_log", "err_log", "script_path", "depends_on"
    ]
    job_stats_path = os.path.join(logs, "step3_job_stats.tsv")
    write_tsv(job_stats_path, job_stats, job_cols)

    # -------- stage rollups --------
    # For SNV, user asked for complete part duration (not per chunk):
    # per sample SNV start = min(start of SNV chunk jobs)
    # per sample SNV end   = end of SNV merge job
    by_sample = {}
    for r in job_stats:
        by_sample.setdefault(r["sample"], []).append(r)

    stage_rows = []

    def add_simple_stage(sample: str, stage_name: str, filt_stage: str):
        subset = [x for x in by_sample.get(sample, []) if x["stage"] == filt_stage]
        if not subset:
            return
        starts = [parse_time(x["start_time"]) for x in subset if x["start_time"]]
        ends = [parse_time(x["end_time"]) for x in subset if x["end_time"]]
        start = min(starts) if starts else None
        end = max(ends) if ends else None
        wall = int((end - start).total_seconds()) if (start and end) else None
        req_cpus = sum(int(x["requested_ncpus"]) for x in subset if x["requested_ncpus"].isdigit())
        req_mem = sum(int(x["requested_mem_mb"]) for x in subset if x["requested_mem_mb"].isdigit())
        cpu_used = sum(int(x["cpu_time_seconds"]) for x in subset if x["cpu_time_seconds"].isdigit())
        max_mem = max([int(x["max_mem_mb"]) for x in subset if x["max_mem_mb"].isdigit()] or [0])
        statuses = ",".join(sorted(set(x["status"] for x in subset if x["status"])))

        stage_rows.append({
            "sample": sample,
            "stage": stage_name,
            "job_ids": ",".join(x["job_id"] for x in subset),
            "start_time": fmt_dt(start),
            "end_time": fmt_dt(end),
            "elapsed_wall_seconds": str(wall) if wall is not None else "",
            "elapsed_wall_hms": sec_to_hms(wall) if wall is not None else "",
            "requested_cpus_total": str(req_cpus),
            "requested_mem_mb_total": str(req_mem),
            "cpu_time_seconds_total": str(cpu_used),
            "max_mem_mb_peak": str(max_mem),
            "status_set": statuses,
        })

    for sample, items in by_sample.items():
        if sample == "ALL":
            continue

        # SNV complete stage rollup
        snv_chunks = [x for x in items if x["stage"] == "SNV_DEEPVARIANT_CHUNK"]
        snv_merge = [x for x in items if x["stage"] == "SNV_DEEPVARIANT_MERGE"]
        if snv_chunks and snv_merge:
            starts = [parse_time(x["start_time"]) for x in snv_chunks if x["start_time"]]
            ends = [parse_time(x["end_time"]) for x in snv_merge if x["end_time"]]
            start = min(starts) if starts else None
            end = max(ends) if ends else None
            wall = int((end - start).total_seconds()) if (start and end) else None

            snv_all = snv_chunks + snv_merge
            req_cpus = sum(int(x["requested_ncpus"]) for x in snv_all if x["requested_ncpus"].isdigit())
            req_mem = sum(int(x["requested_mem_mb"]) for x in snv_all if x["requested_mem_mb"].isdigit())
            cpu_used = sum(int(x["cpu_time_seconds"]) for x in snv_all if x["cpu_time_seconds"].isdigit())
            max_mem = max([int(x["max_mem_mb"]) for x in snv_all if x["max_mem_mb"].isdigit()] or [0])
            statuses = ",".join(sorted(set(x["status"] for x in snv_all if x["status"])))

            stage_rows.append({
                "sample": sample,
                "stage": "SNV_DEEPVARIANT_COMPLETE",
                "job_ids": ",".join(x["job_id"] for x in snv_all),
                "start_time": fmt_dt(start),
                "end_time": fmt_dt(end),
                "elapsed_wall_seconds": str(wall) if wall is not None else "",
                "elapsed_wall_hms": sec_to_hms(wall) if wall is not None else "",
                "requested_cpus_total": str(req_cpus),
                "requested_mem_mb_total": str(req_mem),
                "cpu_time_seconds_total": str(cpu_used),
                "max_mem_mb_peak": str(max_mem),
                "status_set": statuses,
            })

        add_simple_stage(sample, "CNV_CNVKIT_COMPLETE", "CNV_CNVKIT")
        add_simple_stage(sample, "SV_MANTA_COMPLETE", "SV_MANTA")
        add_simple_stage(sample, "ANN_SNPEFF_COMPLETE", "ANN_SNPEFF")
        add_simple_stage(sample, "ANN_VEP_COMPLETE", "ANN_VEP")

    # Cohort totals for SNV/CNV/SV
    def add_cohort_total(stage_prefix: str, out_name: str):
        subset = [x for x in stage_rows if x["stage"] == stage_prefix]
        if not subset:
            return
        starts = [parse_time(x["start_time"]) for x in subset if x["start_time"]]
        ends = [parse_time(x["end_time"]) for x in subset if x["end_time"]]
        start = min(starts) if starts else None
        end = max(ends) if ends else None
        wall = int((end - start).total_seconds()) if (start and end) else None

        stage_rows.append({
            "sample": "ALL_SAMPLES",
            "stage": out_name,
            "job_ids": ";".join(x["job_ids"] for x in subset),
            "start_time": fmt_dt(start),
            "end_time": fmt_dt(end),
            "elapsed_wall_seconds": str(wall) if wall is not None else "",
            "elapsed_wall_hms": sec_to_hms(wall) if wall is not None else "",
            "requested_cpus_total": str(sum(int(x["requested_cpus_total"]) for x in subset if x["requested_cpus_total"].isdigit())),
            "requested_mem_mb_total": str(sum(int(x["requested_mem_mb_total"]) for x in subset if x["requested_mem_mb_total"].isdigit())),
            "cpu_time_seconds_total": str(sum(int(x["cpu_time_seconds_total"]) for x in subset if x["cpu_time_seconds_total"].isdigit())),
            "max_mem_mb_peak": str(max([int(x["max_mem_mb_peak"]) for x in subset if x["max_mem_mb_peak"].isdigit()] or [0])),
            "status_set": ",".join(sorted(set(",".join(x["status_set"] for x in subset).split(",")))),
        })

    add_cohort_total("SNV_DEEPVARIANT_COMPLETE", "ALL_SAMPLES_SNV_DEEPVARIANT_TOTAL")
    add_cohort_total("CNV_CNVKIT_COMPLETE", "ALL_SAMPLES_CNV_CNVKIT_TOTAL")
    add_cohort_total("SV_MANTA_COMPLETE", "ALL_SAMPLES_SV_MANTA_TOTAL")

    stage_cols = [
        "sample", "stage", "job_ids",
        "start_time", "end_time",
        "elapsed_wall_seconds", "elapsed_wall_hms",
        "requested_cpus_total", "requested_mem_mb_total",
        "cpu_time_seconds_total", "max_mem_mb_peak", "status_set",
    ]
    stage_path = os.path.join(logs, "step3_stage_stats.tsv")
    write_tsv(stage_path, stage_rows, stage_cols)

    print(f"[OK] Wrote job-level stats:   {job_stats_path}")
    print(f"[OK] Wrote stage-level stats: {stage_path}")


if __name__ == "__main__":
    main()
