import csv
import json
from pathlib import Path
import shutil
import subprocess
import unittest
from unittest.mock import patch
import uuid

import step7_ginkgo_cnv as stage


class GinkgoTests(unittest.TestCase):
    def setUp(self):
        self.work = stage.HERE / (".ginkgo_test_" + uuid.uuid4().hex)
        self.work.mkdir()
        self.addCleanup(shutil.rmtree, self.work)

    def manifest(self, text):
        path = self.work / "manifest.tsv"
        path.write_text(text)
        return path

    def test_standard_manifest_relative_to_results_root_and_pass_qc(self):
        root = self.work / "results"
        (root / "manifests").mkdir(parents=True)
        for name in ("B.bam", "B.bam.bai", "A.bam", "A.bam.bai"):
            (root / "run1" / name).parent.mkdir(exist_ok=True)
            (root / "run1" / name).write_text("x")
        manifest = root / "manifests/passed_cells_only_manifest.tsv"
        manifest.write_text("sample_id\tbam\tbai\tbarcode_id\tpass_qc\n"
                            "B\trun1/B.bam\trun1/B.bam.bai\tb\ttrue\n"
                            "C\trun1/missing.bam\trun1/missing.bai\tc\tFALSE\n"
                            "A\trun1/A.bam\trun1/A.bam.bai\ta\tTrue\n")
        cells, excluded = stage.load_manifest(manifest)
        self.assertEqual([c["sample_name"] for c in cells], ["B", "A"])
        self.assertEqual(cells[0]["bam_path"], str(root / "run1/B.bam"))
        self.assertEqual(cells[1]["bai_path"], str(root / "run1/A.bam.bai"))
        self.assertEqual(excluded, [{"sample_name": "C", "reason": "pass_qc=FALSE"}])
        cells, _ = stage.load_manifest(manifest, path_root=root)
        self.assertEqual(len(cells), 2)

    def test_manifest_errors(self):
        (self.work / "a.bam").write_text("x")
        for text, message in (
                ("foo\tbar\nx\ty\n", "requires sample_id"),
                ("sample_id\tbam\n", "no passing"),
                ("sample_id\tbam\nA\ta.bam\nA\ta.bam\n", "duplicate"),
                ("sample_id\tbam\n../A\ta.bam\n", "unsafe"),
                ("sample_id\tbam\nA\t\n", "empty bam"),
                ("sample_id\tbam\nA\tmissing.bam\n", "bam not found"),
                ("sample_id\tbam\tbai\nA\ta.bam\tmissing.bai\n", "bai not found"),
                ("sample_id\tbam\tpass_qc\nA\ta.bam\tmaybe\n", "invalid pass_qc"),
                ("sample_id\tbam\tpass_qc\nA\ta.bam\tfalse\n", "no passing"),
                ("sample_id\tbam\nA\ta.bam\textra\n", "column count")):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, message):
                stage.load_manifest(self.manifest(text))

    def test_real_50bp_resources_are_usable(self):
        paths, bins = stage.load_references(stage.REFERENCE, 1000000, 50)
        self.assertEqual(len(bins), 2467)
        self.assertEqual({c for c, _ in bins},
                         {f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]})
        self.assertTrue(all(Path(p).stat().st_size for p in paths.values()))

    def test_missing_reference_never_falls_back(self):
        with self.assertRaisesRegex(ValueError, "no read-length fallback"):
            stage.load_references(stage.REFERENCE, 1000000, 51)

    def test_chromosome_aliases_and_header_errors(self):
        fai = self.work / "ref.fai"
        fai.write_text("1\t100\t0\t0\t0\n2\t100\t0\t0\t0\nMT\t10\t0\t0\t0\n"
                       "chr1\t100\t0\t0\t0\n")
        bins = [("chr1", 50), ("chr2", 90)]
        header = "@SQ\tSN:1\tLN:100\n@SQ\tSN:2\tLN:100\n@SQ\tSN:MT\tLN:10\n"
        mapping, seq = stage.chromosome_map(header, bins, fai)
        self.assertEqual(mapping, {"1": "chr1", "2": "chr2"})
        self.assertIn("MT", seq)
        for bad, expected in ((header.replace("LN:100", "LN:99"), "disagrees"),
                              (header.replace("@SQ\tSN:2\tLN:100\n", ""), "missing modeled"),
                              (header + "@SQ\tSN:chr1\tLN:100\n", "Ambiguous"),
                              ("@HD\tVN:1.6\n", "no @SQ")):
            with self.subTest(expected=expected), self.assertRaisesRegex(ValueError, expected):
                stage.chromosome_map(bad, bins, fai)

    def test_real_binunsorted_smoke(self):
        binary = stage.SOURCE / "container/ginkgo/binUnsorted"
        paths, bins = stage.load_references(stage.REFERENCE, 1000000, 50)
        bed = self.work / "reads.bed"
        bed.write_text("chr2\t20\t30\tr3\t60\t+\nchr1\t3000000\t3000010\tr2\t60\t-\n"
                       "chr1\t20\t30\tr1\t60\t+\nchrY\t57227410\t57227415\tr4\t60\t+\n")
        output = self.work / "counts"
        subprocess.run([str(binary), paths["bins"], str(len(bins) + 1),
                        str(bed), "cell", str(output)],
                       check=True, capture_output=True, text=True)
        values = stage.read_counts(output, "cell", len(bins))
        self.assertEqual(sum(values), 4)
        self.assertEqual(values[:2], [1, 1])
        self.assertEqual(values[next(i for i, b in enumerate(bins) if b[0] == "chr2")], 1)
        self.assertEqual(values[-1], 1)

    def test_counts_header_and_number_errors(self):
        output = self.work / "counts"
        output.write_text("wrong\n1\n")
        with self.assertRaisesRegex(ValueError, "header"):
            stage.read_counts(output, "A", 1)
        output.write_text("A\n1\n")
        with self.assertRaisesRegex(ValueError, "expected 2"):
            stage.read_counts(output, "A", 2)

    def test_real_bam_bedtools_worker_with_numeric_chromosomes(self):
        tools = stage.tool_paths(Path(stage.sys.executable).parent)
        if not all(Path(tools[t]).exists() for t in ("samtools", "bedtools")):
            self.skipTest("Current env samtools/bedtools unavailable")
        resources, bins = stage.load_references(stage.REFERENCE, 1000000, 50)
        lengths = {}
        for chrom, end in bins:
            lengths[stage.canonical(chrom)] = end + 100
        lengths["MT"] = 16569
        ref = self.work / "synthetic.fa"
        ref.write_text(">synthetic\nA\n")
        Path(str(ref) + ".fai").write_text(
            "".join(f"{c}\t{n}\t0\t0\t0\n" for c, n in lengths.items()))
        sam = self.work / "synthetic.sam"
        sam.write_text("@HD\tVN:1.6\tSO:unsorted\n"
                       + "".join(f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in lengths.items())
                       + "".join(f"r{i}\t0\t{c}\t{pos}\t60\t50M\t*\t0\t0\t"
                                 + "A" * 50 + "\t" + "I" * 50 + "\n"
                                 for i, (c, pos) in enumerate(
                                     [("Y", 57227300), ("1", 30), ("2", 30), ("MT", 30)])))
        bam = self.work / "synthetic.bam"
        stage.command([tools["samtools"], "view", "-b", "-o", bam, sam], capture_output=True)
        (self.work / "cells/A").mkdir(parents=True)
        config = {"stage": str(self.work), "tools": tools, "source": str(stage.SOURCE),
                  "reference_dir": str(stage.REFERENCE), "resources": resources,
                  "bin_size": 1000000, "read_length": 50, "ref": str(ref),
                  "cells": [{"sample_name": "A", "bam_path": str(bam)}]}
        self.assertTrue(stage.bin_cell(config, "A"))
        status = json.loads((self.work / "cells/A/status.json").read_text())
        self.assertEqual(status["modeled_records"], 3)
        self.assertEqual(status["unmodeled_records"], 1)
        with (self.work / "cells/A/chromosome_audit.tsv").open() as stream:
            audit = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(audit[-1]["disposition"], "not_in_source_reference")
        self.assertFalse((self.work / "cells/A/reads.bed").exists())

    def test_failed_cohort_retains_every_cell(self):
        (self.work / "cells/A").mkdir(parents=True)
        (self.work / "cells/B").mkdir(parents=True)
        (self.work / "cells/A/status.json").write_text(
            json.dumps({"sample_name": "A", "status": "ERROR", "reason": "bad BAM header"}))
        config = {"stage": str(self.work), "cells": [{"sample_name": "A"}, {"sample_name": "B"}]}
        with self.assertRaisesRegex(ValueError, "no cells silently removed"):
            stage.merge_cohort(config)
        with (self.work / "cohort_status.tsv").open() as stream:
            statuses = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual([s["sample_name"] for s in statuses], ["A", "B"])
        self.assertIn("bad BAM header", statuses[0]["reason"])
        self.assertIn("did not produce", statuses[1]["reason"])
        self.assertFalse((self.work / "COMPLETE").exists())

    def test_lsf_dependency_and_manifest(self):
        for d in ("jobs", "logs"):
            (self.work / d).mkdir()
        config = {"stage": str(self.work), "tool_bin": "/env/bin",
                  "cells": [{"sample_name": "A"}, {"sample_name": "B"}]}
        responses = [subprocess.CompletedProcess([], 0, f"Job <{n}> submitted\n", "")
                     for n in (101, 102, 103)]
        with patch.object(stage, "command", side_effect=responses) as run:
            stage.submit(config, self.work / "config.json", "standard", 8000)
        call = run.call_args_list[-1].args[0]
        self.assertIn("ended(101) && ended(102)", call)
        with (self.work / "job_manifest.tsv").open() as stream:
            jobs = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual([j["job_id"] for j in jobs], ["101", "102", "103"])
        self.assertEqual(jobs[-1]["sample_name"], "COHORT")

    def test_undefined_metrics_are_explicit_na_with_reason(self):
        for value in ("NA", "NaN", "nan", "Inf", "-Inf", ""):
            result = stage.metric_values({"MAPD_CNV_Log2": value, "SKEW_CNV": "0.1"})
            self.assertEqual(result["MAPD_CNV_Log2"], "NA")
            self.assertTrue(result["MAPD_CNV_Log2_reason"])
            self.assertEqual(result["SKEW_CNV"], "0.1")
            self.assertEqual(result["SKEW_CNV_reason"], "")

if __name__ == "__main__":
    unittest.main()
