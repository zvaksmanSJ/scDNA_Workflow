"""Synthetic-only functional check of the installed source lineage methods."""
import csv
from pathlib import Path
import subprocess
import tempfile
import unittest

import run_all_resolveome as workflow


class SourceLineageSmoke(unittest.TestCase):
    def test_source_phylogeny_and_figures(self):
        image = workflow.TOOLS / "containers/basej-somatic/custom_snp_somatic_filter_sequoia_feb2026.sif"
        with tempfile.TemporaryDirectory(prefix="resolveome_lineage_smoke_") as temp:
            root = Path(temp)
            nr, nv, gt = root / "nr.tsv", root / "nv.tsv", root / "gt.tsv"
            names = ["Synthetic_A", "Synthetic_B", "Synthetic_C", "Synthetic_D"]
            patterns = [(1, 1, 0, 0), (0, 0, 1, 1), (1, 0, 0, 0),
                        (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            with nr.open("w") as r, nv.open("w") as v, gt.open("w") as g:
                rw, vw, gw = (csv.writer(s, delimiter="\t") for s in (r, v, g))
                rw.writerow(["VariantId", *names])
                vw.writerow(["VariantId", *names])
                for index in range(120):
                    variant = f"1_{1000+index}_C_T"
                    pattern = patterns[index % len(patterns)]
                    rw.writerow([variant, 40, 40, 40, 40])
                    vw.writerow([variant, *(20 * present for present in pattern)])
                    for name, present in zip(names, pattern):
                        gw.writerow([name, variant, "0/1" if present else "0/0"])
            output = root / "output"
            output.mkdir()
            def execute(argv):
                with (root / "lineage.log").open("a") as log:
                    result = subprocess.run(["apptainer", "exec", "-B", f"{root}:{root}", image,
                                             *map(str, argv)], stdout=log, stderr=subprocess.STDOUT, cwd=root)
                if result.returncode:
                    self.fail((root / "lineage.log").read_text()[-8000:])
            execute(["Rscript", "/usr/local/bin/scripts/rscript_sequoia_build_phylogeny_only.R",
                     "--donor_id", "Synthetic", "--input_nr", nr, "--input_nv", nv,
                     "--output_dir", str(output) + "/", "--only_snvs", "TRUE", "--gender", "male",
                     "--vaf_absent", "0.1", "--vaf_present", "0.3", "--create_multi_tree", "TRUE",
                     "--mpboot_path", "/usr/local/bin/scripts/", "--mpboot_bootstrap", "1000"])
            trees = list(output.glob("*tree_with_branch_length_selectedscheme.tree"))
            self.assertEqual(len(trees), 1, list(output.iterdir()))
            text = trees[0].read_text()
            for name in names:
                self.assertIn(name, text)
            placed = list(output.glob("*placed_variants_all.tsv"))
            self.assertEqual(len(placed), 1)
            execute(["Rscript", "/usr/local/bin/scripts/rscript_postprocess_vaf_tree_egg.R",
                     placed[0], nv, nr, trees[0], gt, "", "5", "false"])
            figures = list(root.rglob("*.png")) + list(root.rglob("*.pdf"))
            self.assertTrue(figures, (root / "lineage.log").read_text()[-6000:])


if __name__ == "__main__":
    unittest.main()
