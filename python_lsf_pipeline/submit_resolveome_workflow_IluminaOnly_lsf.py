#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
from pathlib import Path
from datetime import datetime


def parse_args():
    p = argparse.ArgumentParser(
        description="Submit ResolveOME Illumina workflow: step1 alignment + conditional step2 barcode/cellcall."
    )
    p.add_argument("--cells", required=True, help="TSV with header: sample_id\\tr1\\tr2")
    p.add_argument("--outdir", required=True, help="Output root directory")
    p.add_argument("--ref", required=True, help="Reference FASTA")
    p.add_argument("--queue", default="standard")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--mem-mb", type=int, default=32000)
    p.add_argument("--wall", default="24:00")
    p.add_argument("--env", default="resolveome-align")
    p.add_argument("--assay", default="DNA")
    p.add_argument("--platform", default="ILLUMINA", choices=["ILLUMINA"])
    p.add_argument("--sort-mem", default="3G", help="samtools sort memory per thread (forward to step1)")
    p.add_argument("--skip-fastqc", action="store_true", help="Skip FastQC in step1")

    # step2 knobs
    p.add_argument("--step2-threads", type=int, default=2)
    p.add_argument("--step2-mem-mb", type=int, default=8000)
    p.add_argument("--step2-wall", default="04:00")
    p.add_argument("--min-total-reads", type=int, default=1_000_000)
    p.add_argument("--min-map-rate", type=float, default=0.70)
    p.add_argument("--min-proper-pair-rate", type=float, default=0.60)

    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate_inputs(args):
    cells = Path(args.cells).resolve()
    ref = Path(args.ref).resolve()
    outdir = Path(args.outdir).resolve()

    if not cells.is_file():
        raise FileNotFoundError(f"--cells not found: {cells}")
    if not ref.is_file():
        raise FileNotFoundError(f"--ref not found: {ref}")
    if not Path(str(ref) + ".fai").is_file():
        raise FileNotFoundError(f"Missing FASTA index: {ref}.fai")
    outdir.mkdir(parents=True, exist_ok=True)
    return cells, ref, outdir


def expected_from_cells(cells_tsv: Path):
    expected = []
    with cells_tsv.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            sid = row["sample_id"]
            expected.append(sid)
    return expected


def sample_done(outdir: Path, sample_id: str):
    d = outdir / "bam" / sample_id
    return (
        (d / "alignment.complete").is_file() and
        any(d.glob("*.dedup.bam")) and
        any(d.glob("*.dedup.bam.bai")) and
        any(d.glob("*.dedup.flagstat.txt")) and
        any(d.glob("*.dedup.idxstats.txt"))
    )


def find_latest_aligned_manifest(outdir: Path):
    mdir = outdir / "manifests"
    pointer = mdir / "aligned_bams_latest.tsv"
    if pointer.is_file():
        target = Path(pointer.read_text().strip())
        if target.is_file():
            return target
    cands = sorted((p for p in mdir.glob("aligned_bams_*.tsv") if p != pointer),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0] if cands else None


def write_aligned_manifest(outdir: Path, sample_ids):
    mdir = outdir / "manifests"
    mdir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest = mdir / f"aligned_bams_{run_id}.tsv"

    with manifest.open("w", newline="") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(["sample_id", "bam", "bai", "flagstat", "idxstats"])
        for sid in sample_ids:
            d = outdir / "bam" / sid
            bam = sorted(d.glob("*.dedup.bam"))[0]
            bai = sorted(d.glob("*.dedup.bam.bai"))[0]
            flagstat = outdir / "qc" / "picard" / sid / f"{sid}.pre_dedup.flagstat.txt"
            if not flagstat.is_file():
                raise FileNotFoundError(flagstat)
            idxstats = sorted(d.glob("*.dedup.idxstats.txt"))[0]
            wr.writerow([sid, str(bam), str(bai), str(flagstat), str(idxstats)])

    return manifest


def run_cmd(cmd, dry_run=False):
    print("[CMD]", " ".join(cmd))
    if dry_run:
        return 0
    subprocess.run(cmd, check=True)
    return 0


def main():
    args = parse_args()
    cells, ref, outdir = validate_inputs(args)

    this_dir = Path(__file__).resolve().parent
    step1 = this_dir / "step1_submit_resolveome_alignment_IluminaOnly_lsf.py"
    step2 = this_dir / "step2_submit_resolveome_barcode_cellcall_IluminaOnly_lsf.py"

    if not step1.is_file():
        raise FileNotFoundError(f"Step1 script not found: {step1}")
    if not step2.is_file():
        raise FileNotFoundError(f"Step2 script not found: {step2}")

    expected = expected_from_cells(cells)
    done = [sid for sid in expected if sample_done(outdir, sid)]
    missing = [sid for sid in expected if sid not in done]

    print(f"[INFO] Step1 completion: {len(done)}/{len(expected)} samples complete.")
    if missing:
        print("[INFO] Missing samples:")
        for sid in missing:
            print("  -", sid)

    # If step1 incomplete, submit step1 and stop here.
    if missing:
        step1_cmd = [
            sys.executable, str(step1),
            "--cells", str(cells),
            "--outdir", str(outdir),
            "--ref", str(ref),
            "--queue", args.queue,
            "--threads", str(args.threads),
            "--mem-mb", str(args.mem_mb),
            "--wall", args.wall,
            "--env", args.env,
            "--assay", args.assay,
            "--platform", args.platform,
            "--sort-mem", args.sort_mem,
        ]
        if args.skip_fastqc:
            step1_cmd.append("--skip-fastqc")
        if args.dry_run:
            step1_cmd.append("--dry-run")

        print("[INFO] Submitting/running step1 because not all samples are complete.")
        run_cmd(step1_cmd, dry_run=False)
        print("[INFO] Step1 invoked. Re-run this wrapper after step1 jobs finish to trigger step2.")
        return

    # Step1 complete -> run step2
    aligned_manifest = find_latest_aligned_manifest(outdir)
    if aligned_manifest is None:
        print("[INFO] No aligned_bams_*.tsv found; creating one from completed BAM outputs.")
        aligned_manifest = write_aligned_manifest(outdir, expected)
        print(f"[INFO] Wrote manifest: {aligned_manifest}")

    step2_cmd = [
        sys.executable, str(step2),
        "--aligned-manifest", str(aligned_manifest),
        "--outdir", str(outdir),
        "--queue", args.queue,
        "--threads", str(args.step2_threads),
        "--mem-mb", str(args.step2_mem_mb),
        "--wall", args.step2_wall,
        "--min-total-reads", str(args.min_total_reads),
        "--min-map-rate", str(args.min_map_rate),
        "--min-proper-pair-rate", str(args.min_proper_pair_rate),
    ]
    if args.dry_run:
        step2_cmd.append("--dry-run")

    print("[INFO] Step1 complete for all samples; invoking step2.")
    run_cmd(step2_cmd, dry_run=args.dry_run)
    print("[INFO] Workflow wrapper done.")


if __name__ == "__main__":
    main()
