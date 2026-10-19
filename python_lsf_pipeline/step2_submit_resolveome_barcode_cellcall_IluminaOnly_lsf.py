#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path
from datetime import datetime

def parse_flagstat(flagstat_path: Path):
    total = mapped = proper = 0
    with flagstat_path.open() as f:
        for line in f:
            line = line.strip()
            if " in total " in line:
                m = re.match(r'^(\d+)\s+\+\s+\d+\s+in total', line)
                if m: total = int(m.group(1))
            elif " mapped (" in line and "primary mapped" not in line:
                m = re.match(r'^(\d+)\s+\+\s+\d+\s+mapped', line)
                if m: mapped = int(m.group(1))
            elif " properly paired (" in line:
                m = re.match(r'^(\d+)\s+\+\s+\d+\s+properly paired', line)
                if m: proper = int(m.group(1))
    return total, mapped, proper

def safe_rate(n, d): return (n / d) if d else 0.0

def main():
    ap = argparse.ArgumentParser(description="Step2: Barcode/cell-calling aggregation submitter (LSF).")
    ap.add_argument("--aligned-manifest", required=True, help="aligned_bams_*.tsv")
    ap.add_argument("--outdir", required=True, help="Output root for step2")
    ap.add_argument("--queue", default="standard")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--mem-mb", type=int, default=8000)
    ap.add_argument("--wall", default="145:00")
    ap.add_argument("--min-total-reads", type=int, default=1_000_000)
    ap.add_argument("--min-map-rate", type=float, default=0.70)
    ap.add_argument("--min-proper-pair-rate", type=float, default=0.60)
    ap.add_argument("--dependency", default=None, help='LSF dependency expression, e.g. "done(123)&&done(124)"')
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--local", action="store_true", help="Aggregate in the current process instead of submitting LSF")
    args = ap.parse_args()

    aligned_manifest = Path(args.aligned_manifest).resolve()
    outdir = Path(args.outdir).resolve()
    logs = outdir / "logs"
    scripts = outdir / "scripts"
    manifests = outdir / "manifests"
    for d in (outdir, logs, scripts, manifests):
        d.mkdir(parents=True, exist_ok=True)

    if not aligned_manifest.is_file():
        raise FileNotFoundError(f"Aligned manifest not found: {aligned_manifest}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    worker = scripts / f"run_barcode_cellcall_{run_id}.py"
    worker.write_text(f"""#!/usr/bin/env python3
import csv, re
from pathlib import Path

aligned_manifest = Path(r\"{aligned_manifest}\")
manifests = Path(r\"{manifests}\")
min_total_reads = {args.min_total_reads}
min_map_rate = {args.min_map_rate}
min_proper_pair_rate = {args.min_proper_pair_rate}

sample_re = re.compile(r'^(?P<cell_label>SY5Y_Cell\\d+)_DNA_(?P<c_global>SJNBL046426_C\\d+)_DNA_(?P<c_local>C\\d+)$')

def parse_flagstat(path: Path):
    total = mapped = proper = 0
    with path.open() as f:
        for line in f:
            line = line.strip()
            if " in total " in line:
                m = re.match(r'^(\\d+)\\s+\\+\\s+\\d+\\s+in total', line)
                if m: total = int(m.group(1))
            elif " mapped (" in line and "primary mapped" not in line:
                m = re.match(r'^(\\d+)\\s+\\+\\s+\\d+\\s+mapped', line)
                if m: mapped = int(m.group(1))
            elif " properly paired (" in line:
                m = re.match(r'^(\\d+)\\s+\\+\\s+\\d+\\s+properly paired', line)
                if m: proper = int(m.group(1))
    return total, mapped, proper

def rate(n,d): return (n/d) if d else 0.0

rows = list(csv.DictReader(aligned_manifest.open(), delimiter='\\t'))
if not rows or len({{row["sample_id"] for row in rows}}) != len(rows):
    raise ValueError("Aligned manifest is empty or contains duplicate sample IDs")
for row in rows:
    for key in ("bam", "bai", "flagstat", "idxstats"):
        if not Path(row[key]).is_file() or Path(row[key]).stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty {{key}} for {{row['sample_id']}}: {{row[key]}}")

whitelist = manifests / "whitelist_cells.tsv"
qc_table = manifests / "cell_qc_metrics.tsv"
cell_aware = manifests / "cell_aware_manifest.tsv"
passed = manifests / "passed_cells_only_manifest.tsv"
failed = manifests / "failed_cells_manifest.tsv"

with whitelist.open("w", newline="") as w:
    wr = csv.writer(w, delimiter='\\t')
    wr.writerow(["sample_id","cell_label","c_global","c_local","barcode_id"])
    for row in rows:
        sid = row["sample_id"]
        m = sample_re.match(sid)
        if m:
            cell_label, c_global, c_local = m.group("cell_label"), m.group("c_global"), m.group("c_local")
        else:
            cell_label, c_global, c_local = sid, "NA", "NA"
        wr.writerow([sid, cell_label, c_global, c_local, f"{{c_global}}_{{c_local}}"])

with qc_table.open("w", newline="") as qf, cell_aware.open("w", newline="") as cf, passed.open("w", newline="") as pf, failed.open("w", newline="") as ff:
    qwr = csv.writer(qf, delimiter='\\t')
    cwr = csv.writer(cf, delimiter='\\t')
    pwr = csv.writer(pf, delimiter='\\t')
    fwr = csv.writer(ff, delimiter='\\t')

    qwr.writerow(["sample_id","total_reads","mapped_reads","properly_paired","map_rate","proper_pair_rate","pass_qc"])
    hdr = ["sample_id","bam","bai","barcode_id","cell_label","c_global","c_local","pass_qc"]
    cwr.writerow(hdr); pwr.writerow(hdr); fwr.writerow(hdr)

    for row in rows:
        sid,bam,bai = row["sample_id"],row["bam"],row["bai"]
        m = sample_re.match(sid)
        if m:
            cell_label, c_global, c_local = m.group("cell_label"), m.group("c_global"), m.group("c_local")
        else:
            cell_label, c_global, c_local = sid, "NA", "NA"
        barcode_id = f"{{c_global}}_{{c_local}}"

        total,mapped,proper = parse_flagstat(Path(row["flagstat"]))
        map_rate = rate(mapped,total)
        proper_rate = rate(proper,total)
        pass_qc = total >= min_total_reads and map_rate >= min_map_rate and proper_rate >= min_proper_pair_rate

        qwr.writerow([sid,total,mapped,proper,f"{{map_rate:.6f}}",f"{{proper_rate:.6f}}",str(pass_qc).lower()])
        out = [sid,bam,bai,barcode_id,cell_label,c_global,c_local,str(pass_qc).lower()]
        cwr.writerow(out)
        (pwr if pass_qc else fwr).writerow(out)

print("DONE")
""")
    worker.chmod(0o755)
    if args.local:
        if args.dry_run:
            print(f"{sys.executable} {worker}")
        else:
            subprocess.run([sys.executable, str(worker)], check=True)
        return

    lsf = scripts / f"barcode_cellcall_{run_id}.lsf.sh"
    dep_line = f'#BSUB -w "{args.dependency}"\n' if args.dependency else ""
    lsf.write_text(f"""#!/bin/bash
#BSUB -J barcode_cellcall
#BSUB -q {args.queue}
#BSUB -n {args.threads}
#BSUB -M {args.mem_mb}
#BSUB -R "rusage[mem={args.mem_mb}]"
#BSUB -W {args.wall}
{dep_line}#BSUB -oo {logs}/barcode_cellcall.%J.out
#BSUB -eo {logs}/barcode_cellcall.%J.err
set -euo pipefail
"{sys.executable}" "{worker}"
""")
    lsf.chmod(0o755)

    cmd = f"bsub < {lsf}"
    if args.dry_run:
        print(cmd)
    else:
        proc = subprocess.run(["bash", "-lc", cmd], check=True, capture_output=True, text=True)
        print(proc.stdout.strip())
        print(f"[INFO] Step2 submitted. Outputs: {manifests}")

if __name__ == "__main__":
    main()
