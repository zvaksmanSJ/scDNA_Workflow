#!/usr/bin/env python3
"""Calibrated per-cell false-discovery-rate (FDR) reporting for SNV/indel calls.

Applies a published synthetic-data-derived PTA amplification false-positive
rate (Luquette et al., used by the SCAN2 single-cell somatic caller and
reported in the companion lung/colon PTA paper this was modeled after:
0.013 SNV/Mb and 0.00073 indel/Mb of accessible genome) to each cell's own
callable genome size and call count, producing a per-cell calibrated FDR
rather than a bare PASS/FAIL gate.

FP_mut  = fp_rate_mut * callable_Mb
FDR_mut = FP_mut / N_mut   (capped at 1.0; NOT_COMPUTABLE if N_mut == 0)

"Callable Mb" is approximated as GENOME_TERRITORY * PCT_1X from this same
cell's own Picard CollectWgsMetrics output (already produced in Step 1),
i.e. the portion of the genome with at least 1x coverage in this cell --
no new upstream tool is required.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

# Luquette et al. synthetic PTA-artifact false-positive rate, per Mb of
# accessible genome (as reported/reused in the companion paper's Methods,
# "SNV and indel false discovery rate estimation").
SNV_FP_RATE_PER_MB = 0.013
INDEL_FP_RATE_PER_MB = 0.00073


def read_wgs_metrics(path: Path) -> dict:
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    return dict(zip(header, values))


def count_pass_calls(vcf_path: Path) -> tuple[int, int]:
    """Return (n_snv, n_indel) among PASS records (biallelic-aware count:
    each ALT allele counted separately by type, matching typical SNV/indel
    call-count conventions)."""
    result = subprocess.run(
        ["bcftools", "view", "-f", "PASS", "-H", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    n_snv = n_indel = 0
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        ref, alt_field = fields[3], fields[4]
        for alt in alt_field.split(","):
            if alt in (".", "<NON_REF>"):
                continue
            if len(ref) == 1 and len(alt) == 1:
                n_snv += 1
            else:
                n_indel += 1
    return n_snv, n_indel


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variants-dir", required=True, type=Path,
                    help="Step 3 variants/ directory")
    p.add_argument("--picard-dir", required=True, type=Path,
                    help="Step 1 qc/picard/ directory")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--final-vcf-suffix", default=".deepvariant.raw.vep.snpeff.vcf.gz")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    with args.manifest.open() as fh:
        samples = [row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")]

    rows = []
    for sample in samples:
        vcf_path = args.variants_dir / sample / f"{sample}{args.final_vcf_suffix}"
        metrics_path = args.picard_dir / sample / f"{sample}.wgs_metrics"
        if not vcf_path.is_file() or not metrics_path.is_file():
            rows.append({"sample_id": sample, "status": "MISSING_INPUT"})
            continue
        try:
            metrics = read_wgs_metrics(metrics_path)
            genome_territory = float(metrics["GENOME_TERRITORY"])
            pct_1x = float(metrics["PCT_1X"])
            callable_mb = genome_territory * pct_1x / 1e6

            n_snv, n_indel = count_pass_calls(vcf_path)
            fp_snv = SNV_FP_RATE_PER_MB * callable_mb
            fp_indel = INDEL_FP_RATE_PER_MB * callable_mb
            fdr_snv = min(1.0, fp_snv / n_snv) if n_snv else None
            fdr_indel = min(1.0, fp_indel / n_indel) if n_indel else None

            rows.append({
                "sample_id": sample, "callable_mb": f"{callable_mb:.2f}",
                "n_snv_pass": n_snv, "n_indel_pass": n_indel,
                "expected_fp_snv": f"{fp_snv:.3f}", "expected_fp_indel": f"{fp_indel:.4f}",
                "fdr_snv": f"{fdr_snv:.5f}" if fdr_snv is not None else "NOT_COMPUTABLE",
                "fdr_indel": f"{fdr_indel:.5f}" if fdr_indel is not None else "NOT_COMPUTABLE",
                "status": "OK",
            })
        except Exception as exc:  # noqa: BLE001 - report per-sample, keep going
            rows.append({"sample_id": sample, "status": f"ERROR: {exc}"})

    fieldnames = ["sample_id", "callable_mb", "n_snv_pass", "n_indel_pass",
                  "expected_fp_snv", "expected_fp_indel", "fdr_snv", "fdr_indel", "status"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    ok = [r for r in rows if r["status"] == "OK"]
    print(f"[OK] FDR computed for {len(ok)}/{len(rows)} samples -> {args.output}")


if __name__ == "__main__":
    main()
