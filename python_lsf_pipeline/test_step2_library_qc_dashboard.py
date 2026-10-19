"""Focused synthetic tests; scratch fixtures stay inside the workflow directory."""

import json
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch
import uuid

import step2_library_qc_dashboard as qc


class LibraryQCTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parent / (".qc_test_" + uuid.uuid4().hex)
        self.root.mkdir()
        self.addCleanup(shutil.rmtree, self.root)
        self.sid = "cellA"
        self.prefix = self.root / "qc/picard/cellA/cellA"
        self.prefix.parent.mkdir(parents=True)
        (self.root / "qc/fastp/cellA").mkdir(parents=True)
        (self.root / "manifests").mkdir()
        self.ref = self.root / "reference.fa"
        self.ref.write_text(">chr1\nACGT\n")
        self.gate = dict(sample_id=self.sid, total_reads="100", mapped_reads="90",
                         properly_paired="80", map_rate="0.9", proper_pair_rate="0.8", pass_qc="true")
        qc.write_tsv(self.root / "manifests/cell_qc_metrics.tsv", [self.gate])
        self.config = json.loads(qc.DEFAULT_CONFIG.read_text())
        self.metric(".alignment_summary_metrics",
                    ["CATEGORY", "TOTAL_READS", "PF_READS", "PF_READS_ALIGNED",
                     "READS_ALIGNED_IN_PAIRS", "PCT_CHIMERAS", "SAMPLE", "LIBRARY", "READ_GROUP"],
                    [["FIRST_OF_PAIR", 50, 50, 45, 40, 0.1, "", "", ""],
                     ["SECOND_OF_PAIR", 50, 50, 45, 40, 0.1, "", "", ""],
                     ["PAIR", 100, 100, 90, 80, 0.1, "", "", ""]])
        self.metric(".insert_size_metrics",
                    ["MEDIAN_INSERT_SIZE", "MEAN_INSERT_SIZE", "STANDARD_DEVIATION", "READ_PAIRS", "PAIR_ORIENTATION"],
                    [[300, 310, 50, 40, "FR"]])
        self.metric(".gc_bias.summary_metrics",
                    ["ACCUMULATION_LEVEL", "TOTAL_CLUSTERS", "ALIGNED_READS", "AT_DROPOUT", "GC_DROPOUT"],
                    [["All Reads", 50, 90, 1.2, 2.3]])
        self.metric(".quality_yield_metrics", ["TOTAL_READS", "TOTAL_BASES", "Q20_BASES", "Q30_BASES"],
                    [[100, 10000, 9500, 9000]])
        Path(str(self.prefix) + ".markdup_stats.txt").write_text(
            "COMMAND: samtools markdup\nREAD: 110\nEXAMINED: 100\nDUPLICATE PAIR: 10\n"
            "DUPLICATE SINGLE: 2\nDUPLICATE TOTAL: 20\n")
        fastp = {
            "command": "fastp --in1 cellA_R1.fastq.gz --in2 cellA_R2.fastq.gz --out1 cellA_R1_trim.fastq.gz --out2 cellA_R2_trim.fastq.gz",
            "summary": {"before_filtering": dict(total_reads=110, total_bases=11000, q20_rate=.95,
                                                q30_rate=.9, gc_content=.4, read1_mean_length=100,
                                                read2_mean_length=100),
                        "after_filtering": dict(total_reads=100, q30_rate=.92)},
            "adapter_cutting": dict(adapter_trimmed_reads=3, adapter_trimmed_bases=20),
            "duplication": dict(rate=.1),
        }
        (self.root / "qc/fastp/cellA/cellA.fastp.json").write_text(json.dumps(fastp))

    def metric(self, suffix, headers, rows):
        path = Path(str(self.prefix) + suffix)
        classes = {
            ".alignment_summary_metrics": "AlignmentSummaryMetrics",
            ".insert_size_metrics": "InsertSizeMetrics",
            ".gc_bias.summary_metrics": "GcBiasSummaryMetrics",
            ".quality_yield_metrics": "QualityYieldMetrics",
        }
        path.write_text("## METRICS CLASS\tpicard.analysis." + classes[suffix] + "\n"
                        + "\t".join(headers) + "\n"
                        + "".join("\t".join(map(str, row)) + "\n" for row in rows)
                        + "\n## HISTOGRAM\tjava.lang.Integer\nvalue\tcount\n1\t100\n")
        return path

    def collect(self, extra=None):
        return qc.collect_cell(self.root, self.sid, self.gate, self.config, extra or {})

    def test_pair_row_trailing_columns_and_non_picard_dedup(self):
        row = self.collect()
        self.assertEqual(row["alignment_CATEGORY"], "PAIR")
        self.assertEqual(row["alignment_SAMPLE"], "")
        self.assertEqual(row["alignment_TOTAL_READS"], 100)
        self.assertAlmostEqual(row["markdup_primary_duplicate_fraction"], .12)
        self.assertEqual(row["picard_q30_rate"], .9)
        self.assertEqual(row["source_qc_informational_status"], "UNKNOWN")
        self.assertIsNone(row["source_composite_score"])
        self.assertEqual(row["source_components_known"], 1)
        self.assertEqual(row["flagstat_gate_pass_qc"], "true")

    def test_strict_source_thresholds_and_status(self):
        values = dict(preseq_count=3.5e9, PCT_CHIMERAS=.2, chrM=.2, MAPD_CNV_Log2=.25, SKEW_CNV=.25)
        _, score, status, _, _ = qc.source_scores(values, self.config)
        self.assertEqual((score, status), (0, "FAIL"))
        values.update(preseq_count=3.5e9 + 1, PCT_CHIMERAS=.1, chrM=.1, MAPD_CNV_Log2=.1)
        _, score, status, _, _ = qc.source_scores(values, self.config)
        self.assertEqual((score, status), (4, "PASS"))
        values["chrM"] = .3
        self.assertEqual(qc.source_scores(values, self.config)[1:3], (3, "Borderline"))

    def test_low_coverage_insert_header_only_is_na(self):
        self.metric(".insert_size_metrics",
                    ["MEDIAN_INSERT_SIZE", "MEAN_INSERT_SIZE", "STANDARD_DEVIATION", "READ_PAIRS"], [])
        row = self.collect()
        self.assertIsNone(row["insert_MEDIAN_INSERT_SIZE"])
        self.assertIn("header-only", row["insert_reason"])

    def test_missing_file_is_error_even_for_no_data(self):
        Path(str(self.prefix) + ".insert_size_metrics").unlink()
        with self.assertRaisesRegex(ValueError, "requires explicit"):
            self.collect()

    def test_missing_insert_requires_exact_matching_status(self):
        Path(str(self.prefix) + ".insert_size_metrics").unlink()
        status = Path(str(self.prefix) + ".insert_size.status.tsv")
        status.write_text("not_computable\tno_eligible_pairs\n")
        row = self.collect()
        self.assertIsNone(row["insert_MEDIAN_INSERT_SIZE"])
        self.assertIn("no_eligible_pairs", row["insert_reason"])
        status.write_text("not_computable\tlow_coverage\n")
        with self.assertRaisesRegex(ValueError, "expected not_computable"):
            self.collect()

    def test_status_conflicting_with_positive_insert_pairs_is_error(self):
        Path(str(self.prefix) + ".insert_size.status.tsv").write_text("not_computable\tno_eligible_pairs\n")
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.collect()

    def test_current_ginkgo_path_preferred_over_legacy(self):
        for directory, value in (("ginkgo", "0.1"), ("ginkgo_cnv", "0.3")):
            path = self.root / directory
            path.mkdir()
            (path / "source_metrics.tsv").write_text(
                "sample_name\tMAPD_CNV_Log2\tSKEW_CNV\ncellA\t" + value + "\t0.1\n")
        def fake_multiqc(command, **kwargs):
            output = Path(command[command.index("--outdir") + 1])
            (output / "multiqc_library_qc.html").write_text("<html>MultiQC</html>")
        with patch.object(qc.subprocess, "run", side_effect=fake_multiqc):
            report = qc.run(self.root, self.ref, multiqc="/existing/multiqc")
        provenance = json.loads((report / "provenance.json").read_text())
        self.assertEqual(provenance["source_metrics"], str(self.root / "ginkgo/source_metrics.tsv"))

    def test_missing_fields_and_bad_row_are_errors(self):
        path = Path(str(self.prefix) + ".quality_yield_metrics")
        path.write_text("TOTAL_READS\tTOTAL_BASES\n10\t100\n")
        with self.assertRaisesRegex(ValueError, "missing fields"):
            self.collect()
        path.write_text("TOTAL_READS\tTOTAL_BASES\tQ20_BASES\tQ30_BASES\n10\t100\t90\n")
        with self.assertRaisesRegex(ValueError, "field count"):
            self.collect()

    def test_unpaired_chimeras_are_unknown_not_false_pass(self):
        self.metric(".alignment_summary_metrics",
                    ["CATEGORY", "TOTAL_READS", "PF_READS", "PF_READS_ALIGNED", "READS_ALIGNED_IN_PAIRS", "PCT_CHIMERAS"],
                    [["UNPAIRED", 100, 100, 90, 0, 0]])
        row = self.collect()
        self.assertIsNone(row["PCT_CHIMERAS"])
        self.assertEqual(row["PCT_CHIMERAS_source_check"], "UNKNOWN")

    def test_explicit_upstream_na_and_full_score(self):
        extras = dict(preseq_count="4000000000", chrM="0.1", MAPD_CNV_Log2="0.1", SKEW_CNV="0.1")
        self.assertEqual(self.collect(extras)["source_composite_score"], 5)
        extras.update(SKEW_CNV="NA", SKEW_CNV_reason="Insufficient bins")
        row = self.collect(extras)
        self.assertEqual(row["SKEW_CNV_reason"], "Insufficient bins")
        self.assertIsNone(row["source_composite_score"])

    def test_na_needs_reason_and_picard_nan_has_reason(self):
        with self.assertRaisesRegex(ValueError, "needs SKEW_CNV_reason"):
            self.collect({"SKEW_CNV": "NA"})
        self.metric(".gc_bias.summary_metrics",
                    ["ACCUMULATION_LEVEL", "TOTAL_CLUSTERS", "ALIGNED_READS", "AT_DROPOUT", "GC_DROPOUT"],
                    [["All Reads", 50, 90, "NaN", 2.3]])
        row = self.collect()
        self.assertIsNone(row["gc_AT_DROPOUT"])
        self.assertIn("undefined", row["gc_AT_DROPOUT_reason"])

    def test_picard_320_gc_reads_used_enum_is_metadata(self):
        self.metric(".gc_bias.summary_metrics",
                    ["ACCUMULATION_LEVEL", "READS_USED", "TOTAL_CLUSTERS", "ALIGNED_READS",
                     "AT_DROPOUT", "GC_DROPOUT"],
                    [["All Reads", "ALL", 50, 90, 1.2, 2.3]])
        row = self.collect()
        self.assertEqual(row["gc_READS_USED"], "ALL")
        self.assertEqual(row["gc_AT_DROPOUT"], 1.2)
        self.assertEqual(row["gc_ALIGNED_READS"], 90)
        self.metric(".gc_bias.summary_metrics",
                    ["ACCUMULATION_LEVEL", "READS_USED", "TOTAL_CLUSTERS", "ALIGNED_READS",
                     "AT_DROPOUT", "GC_DROPOUT"],
                    [["All Reads", "ALL", 50, 90, "unexpected", 2.3]])
        with self.assertRaisesRegex(ValueError, "gc_AT_DROPOUT: expected a number"):
            self.collect()

    def test_pre_dedup_flagstat_preserved_separately(self):
        Path(str(self.prefix) + ".pre_dedup.flagstat.txt").write_text(
            "110 + 1 in total (QC-passed reads + QC-failed reads)\n"
            "99 + 1 mapped (90.00% : 100.00%)\n"
            "95 + 1 primary mapped (90.00% : 100.00%)\n"
            "88 + 0 properly paired (80.00% : N/A)\n")
        row = self.collect()
        self.assertEqual(row["pre_dedup_total_reads"], 110)
        self.assertEqual(row["pre_dedup_mapped_reads"], 99)
        self.assertEqual(row["pre_dedup_map_rate"], .9)
        self.assertEqual(row["total_reads"], 100)

    def test_fastp_missing_expected_metric_is_error(self):
        path = self.root / "qc/fastp/cellA/cellA.fastp.json"
        data = json.loads(path.read_text())
        del data["summary"]["before_filtering"]["q30_rate"]
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "missing fields"):
            self.collect()

    def test_optional_preseq_and_source_depth_fraction(self):
        directory = self.root / "qc/preseq/cellA"
        directory.mkdir(parents=True)
        (directory / "cellA_preseq.txt").write_text("4000000000\n")
        path = Path(str(self.prefix) + ".cov_sentieonmetrics.sample_interval_summary")
        path.write_text("Target\ttotal_coverage\nchr1\t900\nchrM\t100\n")
        row = self.collect()
        self.assertEqual(row["preseq_count"], 4e9)
        self.assertEqual(row["chrM"], .1)
        path.write_text("Target\ttotal_coverage\nchr1\t0\nchrM\t0\n")
        row = self.collect()
        self.assertIsNone(row["chrM"])
        self.assertIn("No aligned bases", row["chrM_reason"])

    def preseq_failure_fixture(self):
        directory = self.root / "qc/preseq/cellA"
        directory.mkdir(parents=True)
        (directory / "cellA.preseq.status.tsv").write_text(
            "status\treason\nNOT_COMPUTABLE\tpreseq model fit failed\n")
        (directory / "cellA.preseq.log").write_text("model fit failed: insufficient complexity\n")
        return directory

    def test_preseq_model_failure_is_unknown_with_reason(self):
        self.preseq_failure_fixture()
        row = self.collect()
        self.assertIsNone(row["preseq_count"])
        self.assertEqual(row["preseq_count_reason"], "preseq model fit failed")
        self.assertEqual(row["preseq_status"], "NOT_COMPUTABLE")
        self.assertEqual(row["preseq_count_source_check"], "UNKNOWN")

    def test_preseq_failure_marker_needs_raw_logs(self):
        directory = self.preseq_failure_fixture()
        (directory / "cellA.preseq.log").unlink()
        with self.assertRaisesRegex(ValueError, "requires nonempty raw"):
            self.collect()

    def test_preseq_failure_without_marker_is_error(self):
        directory = self.preseq_failure_fixture()
        (directory / "cellA.preseq.status.tsv").unlink()
        with self.assertRaisesRegex(ValueError, "requires NOT_COMPUTABLE"):
            self.collect()

    def test_preseq_noisy_status_is_error(self):
        directory = self.preseq_failure_fixture()
        (directory / "cellA.preseq.status.tsv").write_text(
            "status\treason\nNOT_COMPUTABLE\tpreseq model fit failed\nnoise\n")
        with self.assertRaisesRegex(ValueError, "expected one"):
            self.collect()

    def test_preseq_stale_marker_conflicts_with_success(self):
        directory = self.preseq_failure_fixture()
        (directory / "cellA_preseq.txt").write_text("4000000000\n")
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.collect()
        (directory / "cellA_preseq.txt").unlink()
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.collect({"preseq_count": "4000000000"})

    def test_preseq_invalid_estimate_needs_failure_marker(self):
        directory = self.preseq_failure_fixture()
        (directory / "cellA_preseq.txt").write_text("model-fit-noise\n")
        self.assertIsNone(self.collect()["preseq_count"])
        (directory / "cellA.preseq.status.tsv").unlink()
        with self.assertRaisesRegex(ValueError, "requires NOT_COMPUTABLE"):
            self.collect()

    def test_optional_wgs_is_informational_only(self):
        path = Path(str(self.prefix) + ".wgs_metrics")
        path.write_text("## METRICS CLASS\tpicard.analysis.WgsMetrics\n"
                        "GENOME_TERRITORY\tMEAN_COVERAGE\tPCT_1X\n1000\t0.5\t0.25\n\n")
        row = self.collect()
        self.assertEqual(row["wgs_MEAN_COVERAGE"], .5)
        self.assertEqual(row["source_qc_informational_status"], "UNKNOWN")

    def test_unknown_cells_rejected(self):
        extra = self.root / "extra.tsv"
        extra.write_text("sample_name\tchrM\notherCell\t0.1\n")
        with self.assertRaisesRegex(ValueError, "unknown cells"):
            qc.run(self.root, self.ref, source_metrics=extra)

    def test_report_is_self_contained_and_multiqc_called(self):
        (self.root / "qc/fastqc_raw/cellA").mkdir(parents=True)
        (self.root / "qc/fastqc_trimmed/cellA").mkdir(parents=True)
        def fake_multiqc(command, **kwargs):
            output = Path(command[command.index("--outdir") + 1])
            (output / "multiqc_library_qc.html").write_text("<html>MultiQC</html>")

        with patch.object(qc.subprocess, "run", side_effect=fake_multiqc) as execute:
            report = qc.run(self.root, self.ref, multiqc="/existing/multiqc")
        self.assertTrue(execute.call_args.kwargs["check"])
        dashboard = (report / "source_threshold_dashboard.html").read_text()
        self.assertIn("data:image/svg+xml;base64,", dashboard)
        self.assertIn("UNKNOWN", dashboard)
        self.assertNotIn('<script src=', dashboard)
        self.assertNotIn('<img src="http', dashboard)
        self.assertTrue((report / "cohort_checks.tsv").is_file())
        self.assertIn("NA", (report / "library_qc_metrics.tsv").read_text())
        provenance = json.loads((report / "provenance.json").read_text())
        self.assertEqual(len(provenance["fastqc"]), 2)
        self.assertIn(str(self.root / "qc/fastqc_raw"), execute.call_args.args[0])
        self.assertIn(str(self.root / "qc/fastqc_trimmed"), execute.call_args.args[0])

    def test_multiqc_failures_propagate(self):
        import subprocess
        with patch.object(qc, "render_dashboard"), patch.object(
            qc.subprocess, "run", side_effect=subprocess.CalledProcessError(1, ["multiqc"])
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                qc.run(self.root, self.ref, multiqc="/existing/multiqc")

    @unittest.skipUnless(qc.importlib.util.find_spec("multiqc"), "MultiQC is not installed in the selected interpreter")
    def test_real_multiqc_picard_html_smoke(self):
        report = qc.run(self.root, self.ref)
        self.assertGreater((report / "multiqc_library_qc.html").stat().st_size, 1000)
        log = (report / "multiqc.log").read_text()
        self.assertIn("AlignmentSummaryMetrics reports", log)
        provenance = json.loads((report / "provenance.json").read_text())
        self.assertEqual(provenance["multiqc_command"][:3], [qc.sys.executable, "-m", "multiqc"])


if __name__ == "__main__":
    unittest.main()
