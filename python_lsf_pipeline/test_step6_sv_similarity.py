import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from step6_sv_similarity import (
    Matrix,
    align_matrix,
    clustering_newick,
    jaccard_counts,
    load_cohort,
    load_matrix,
)


class SvSimilarityTests(unittest.TestCase):
    def test_exact_jaccard_counts(self):
        binary = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.int64)
        distance, shared, union = jaccard_counts(binary)
        np.testing.assert_array_equal(shared, [[2, 1, 1], [1, 2, 0], [1, 0, 1]])
        np.testing.assert_array_equal(union, [[2, 3, 2], [3, 2, 3], [2, 3, 1]])
        np.testing.assert_allclose(distance, [[0, 2 / 3, 1 / 2], [2 / 3, 0, 1], [1 / 2, 1, 0]])

    def test_identical_cells_and_newick(self):
        distance, _, _ = jaccard_counts(np.array([[1, 1, 0], [0, 0, 1]], dtype=np.int64))
        tree = linkage(squareform(distance), method="average")
        text = clustering_newick(tree, ["Cell_A", "Cell_B", "Cell'C"])
        self.assertIn("'Cell_A':0", text)
        self.assertIn("'Cell_B':0", text)
        self.assertIn("'Cell''C':1", text)
        self.assertTrue(text.endswith(";\n"))

    def test_reject_uninformative_or_invalid_calls(self):
        for values in ([[1, 0]], [[0, 0], [1, 1]], [[1, 2]], [[1]]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                jaccard_counts(np.array(values, dtype=np.int64))

    def test_align_both_axes_without_filling_missing_cells(self):
        matrix = Matrix(["B", "A"], ["y", "x"], np.array([[4, 3], [2, 1]], dtype=np.int64))
        np.testing.assert_array_equal(align_matrix(matrix, ["A", "B"], ["x", "y"]), [[1, 2], [3, 4]])
        with self.assertRaises(ValueError):
            align_matrix(matrix, ["A", "C"], ["x", "y"])

    def test_matrix_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv"
            for text in (
                "VariantId\tA\tB\nx\t1\n",
                "VariantId\tA\tA\nx\t1\t0\n",
                "VariantId\tA\tB\nx\t1\t0\nx\t0\t1\n",
                "VariantId\tA\tB\nx\t-1\t0\n",
                "VariantId\tA\tB\n",
                "VariantId\tA\tB\nx\tNA\t0\n",
            ):
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    load_matrix(path)
            path.write_text("VariantId\tA\tB\nx\t1\t0\n")
            matrix = load_matrix(path)
            self.assertEqual(matrix.samples, ["A", "B"])
            self.assertEqual(matrix.features, ["x"])

    def test_cohort_requires_passed_unique_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.tsv"
            path.write_text("sample_id\tcell_label\tpass_qc\nA\tCell1\ttrue\nB\tCell2\ttrue\n")
            self.assertEqual(load_cohort(path), (["A", "B"], ["Cell1", "Cell2"]))
            path.write_text("sample_id\tcell_label\tpass_qc\nA\tCell1\ttrue\nB\tCell2\tfalse\n")
            with self.assertRaises(ValueError):
                load_cohort(path)

    def test_cli_writes_complete_reports_and_exact_distances(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sv = directory / "sv"
            sv.mkdir()
            for kind, values in (
                ("binary", "x\t1\t1\t0\ny\t1\t0\t1\nz\t0\t1\t0\n"),
                ("NR", "x\t10\t10\t0\ny\t10\t0\t10\nz\t0\t10\t0\n"),
                ("NV", "x\t5\t5\t0\ny\t5\t0\t5\nz\t0\t5\t0\n"),
            ):
                (sv / f"sv_{kind}.tsv").write_text("VariantId\tA\tB\tC\n" + values)
            (sv / "sv_features.tsv").write_text(
                "VariantId\tSVTYPE\tsupporting_cells\tbreakpoint_tolerance\n"
                "x\tBND\t2\t500\ny\tDEL\t2\t500\nz\tDUP\t1\t500\n"
            )
            manifest = directory / "manifest.tsv"
            manifest.write_text(
                "sample_id\tcell_label\tpass_qc\n"
                "A\tCell1\ttrue\nB\tCell2\ttrue\nC\tCell3\ttrue\n"
            )
            out = directory / "out"
            subprocess.run(
                [
                    sys.executable, str(Path(__file__).with_name("step6_sv_similarity.py")),
                    "--sv-matrix-dir", str(sv), "--manifest", str(manifest),
                    "--output-dir", str(out),
                ],
                check=True, capture_output=True, text=True,
            )
            expected = {
                "jaccard_distance.tsv", "pairwise_sv_similarity.tsv",
                "sample_sv_summary.tsv", "sv_feature_prevalence.tsv",
                "average_linkage.tsv", "sv_similarity.dendrogram.newick",
                "sv_similarity.dendrogram.png", "sv_similarity.distance_heatmap.png",
                "sv_similarity.report.pdf", "sv_similarity.report.html", "analysis_metadata.json",
            }
            self.assertEqual({path.name for path in out.iterdir()}, expected)
            self.assertTrue(all(path.stat().st_size > 0 for path in out.iterdir()))
            with (out / "jaccard_distance.tsv").open() as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(rows[0], ["sample_id", "A", "B", "C"])
            self.assertAlmostEqual(float(rows[1][2]), 2 / 3)
            self.assertAlmostEqual(float(rows[1][3]), 1 / 2)
            metadata = json.loads((out / "analysis_metadata.json").read_text())
            self.assertFalse(metadata["rooted_phylogeny"])
            self.assertEqual(metadata["zero_call_semantics"], "not_detected")
            self.assertEqual(metadata["feature_count"], 3)
            self.assertEqual(metadata["detected_cell_feature_pairs"], 5)
            self.assertIn("confirmed biological absence", (out / "sv_similarity.report.html").read_text())


if __name__ == "__main__":
    unittest.main()
