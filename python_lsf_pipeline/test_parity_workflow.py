"""Focused regression tests for the revised launcher and alignment protocol."""
import csv
import gzip
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
import sys
from unittest.mock import patch

import run_all_resolveome as workflow
import step1_submit_resolveome_alignment_IluminaOnly_lsf as alignment


class ParityWorkflowTests(unittest.TestCase):
    def test_failed_captured_command_exposes_diagnostics(self):
        error = subprocess.CalledProcessError(
            1, ["probe"], output="probe report", stderr="underlying failure")
        with patch.object(workflow.subprocess, "run", side_effect=error), \
                patch.object(workflow.sys, "stderr", new_callable=io.StringIO) as stderr:
            with self.assertRaises(subprocess.CalledProcessError):
                workflow.run(["probe"], capture_output=True)
            self.assertIn("probe report", stderr.getvalue())
            self.assertIn("underlying failure", stderr.getvalue())

    def test_plan_is_non_submitting_and_preserves_deepvariant(self):
        with patch.object(workflow, "run") as execute:
            workflow.main(["--plan"])
            execute.assert_not_called()
        args = workflow.parse_args(["--plan"])
        command = workflow.commands(args)["variants"]
        self.assertIn("--dv-use-allele-frequency", command)
        self.assertIn("--dv-add-hp-channel", command)
        self.assertIn("sc_af.model.ckpt", command)
        self.assertNotIn("--snv-filter-pass-only", command)
        self.assertIn("reference", workflow.commands(args)["somatic_inputs"])

    def test_protects_validated_outputs(self):
        with self.assertRaises(SystemExit):
            workflow.parse_args(["--outdir", str(workflow.PROJECT / "results/DNA_Ilumina_output_1")])

    def test_worker_shell_and_required_commands(self):
        with tempfile.TemporaryDirectory() as temp:
            worker = Path(temp) / "worker.sh"
            alignment.write_worker_script(worker)
            subprocess.run(["bash", "-n", str(worker)], check=True)
            text = worker.read_text()
            self.assertIn("--detect_adapter_for_pe --trim_poly_g", text)
            self.assertIn("bwa-mem2 mem -M -Y -K 2500000000", text)
            self.assertIn("samtools markdup -r -s", text)
            self.assertIn("PROGRAM=CollectQualityYieldMetrics", text)
            self.assertIn("CollectWgsMetrics", text)
            self.assertIn("bam2mr -seg_len 100000", text)
            self.assertIn("preseq-bioskryb gc_extrap", text)
            self.assertLess(text.index("fastp --thread"), text.index("bwa-mem2 mem"))
            self.assertLess(text.index("picard -Xmx8g"), text.index('> "${completion}"'))
            self.assertNotIn("FastQC failed", text)
            self.assertIn('"${pre_flagstat}" "${idxstats}"', text)

    def test_duplicate_or_unsafe_sample_names_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            r1, r2 = root / "r1", root / "r2"
            r1.touch()
            r2.touch()
            table = root / "cells.tsv"
            for names in (["A", "A"], ["bad;id"]):
                with table.open("w") as stream:
                    writer = csv.writer(stream, delimiter="\t")
                    writer.writerow(["sample_id", "r1", "r2"])
                    writer.writerows([name, r1, r2] for name in names)
                with self.assertRaises(ValueError):
                    alignment.read_cells_tsv(table)

    def test_failed_job_stops_dependencies(self):
        result = subprocess.CompletedProcess([], 0, "123 EXIT\n", "")
        with patch.object(workflow, "run", return_value=result):
            with self.assertRaisesRegex(RuntimeError, "123: EXIT"):
                workflow.wait_jobs(["123"], 0)

    def test_successful_jobs_finish_without_polling(self):
        result = subprocess.CompletedProcess([], 0, "123 DONE\n124 DONE\n", "")
        with patch.object(workflow, "run", return_value=result), patch.object(workflow.time, "sleep") as sleep:
            workflow.wait_jobs(["123", "124"], 0)
            sleep.assert_not_called()

    def test_missing_job_status_is_an_error(self):
        result = subprocess.CompletedProcess([], 0, "", "")
        with patch.object(workflow, "run", return_value=result):
            with self.assertRaisesRegex(RuntimeError, "Cannot recover"):
                workflow.wait_jobs(["123"], 0)

    def test_reference_candidates_preserve_nonpass_and_noncoding_records(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sid = "Synthetic"
            bam = root / "synthetic.bam"
            bam.write_bytes(b"exists-only; this test does not analyze BAM contents")
            manifest = root / "passed.tsv"
            manifest.write_text(f"sample_id\tbam\tpass_qc\n{sid}\t{bam}\ttrue\n")
            raw = root / "variants" / sid / "support/merged" / f"{sid}.deepvariant.raw.vcf.gz"
            sv = root / "sv" / sid / "manta/results/variants/diploidSV.vcf.gz"
            header = (
                "##fileformat=VCFv4.2\n##contig=<ID=1,length=200000>\n"
                "##FILTER=<ID=LowQual,Description=\"Rejected caller candidate\">\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}\n"
            )
            for path, records in (
                (raw, "1\t100\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\n"
                      "1\t200\t.\tC\tA\t10\tLowQual\t.\tGT\t0/1\n"),
                (sv, ""),
            ):
                path.parent.mkdir(parents=True)
                with gzip.open(path, "wt") as stream:
                    stream.write(header + records)
            out = root / "inputs"
            result = subprocess.run([
                sys.executable, str(workflow.HERE / "step5_filter_variant_inputs.py"),
                "--manifest", str(manifest), "--manifest-base-dir", str(root),
                "--variants-dir", str(root / "variants"), "--sv-dir", str(root / "sv"),
                "--output-dir", str(out), "--output-csv", str(out / "inputs.csv"),
                "--chroms", "all", "--small-variant-method", "reference",
                "--bcftools", str(Path(sys.executable).parent / "bcftools"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with gzip.open(out / "small_variants/Synthetic.reference_candidates.vcf.gz", "rt") as stream:
                rows = [line for line in stream if not line.startswith("#")]
            self.assertEqual(len(rows), 2)
            self.assertIn("LowQual", rows[1])


if __name__ == "__main__":
    unittest.main()
