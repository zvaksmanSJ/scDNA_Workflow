"""End-to-end CUSTOM Ginkgo smoke: synthetic BAMs -> source binUnsorted/R/CNVcaller/parsers."""
import csv
import os
import random
import shutil
import subprocess
import tarfile
import unittest
import uuid

import step7_ginkgo_cnv as stage

CELLS = ("Cell_A", "Cell_B_gain8", "Cell_C")


class RealGinkgoSmoke(unittest.TestCase):
    def setUp(self):
        self.tools = stage.tool_paths(stage.Path(stage.sys.executable).parent)
        _, errors = stage.environment_report(self.tools, stage.SOURCE)
        if errors:
            self.skipTest("Ginkgo environment unavailable: " + "; ".join(errors))
        self.work = stage.HERE / (".ginkgo_smoke_" + uuid.uuid4().hex)
        (self.work / "results/manifests").mkdir(parents=True)
        if not os.environ.get("GINKGO_SMOKE_KEEP"):
            self.addCleanup(shutil.rmtree, self.work)

    def make_bam(self, name, bins, lengths, header):
        rng = random.Random(name)
        bam = self.work / "results/bam" / f"{name}.bam"
        bam.parent.mkdir(exist_ok=True)
        sort = subprocess.Popen([self.tools["samtools"], "sort", "-o", str(bam), "-"],
                                stdin=subprocess.PIPE, text=True)
        sort.stdin.write(header)
        start, previous = 0, None
        for index, (chrom, end) in enumerate(bins):
            start = 0 if chrom != previous else start
            previous = chrom
            mean = 60 if name == "Cell_B_gain8" and chrom == "chr8" else 40
            reads = max(1, round(rng.gauss(mean, mean ** 0.5)))
            for read in range(reads):
                pos = rng.randint(start + 1, max(start + 1, min(end, lengths[chrom]) - 50))
                sort.stdin.write(f"{name}_{index}_{read}\t0\t{chrom}\t{pos}\t60\t50M\t*\t0\t0\t*\t*\n")
            start = end
        sort.stdin.write(f"{name}_mt\t0\tchrM\t100\t60\t50M\t*\t0\t0\t*\t*\n")
        sort.stdin.close()
        self.assertEqual(sort.wait(), 0)
        subprocess.run([self.tools["samtools"], "index", str(bam)], check=True)
        return bam

    def test_local_cohort_pipeline_outputs(self):
        _, bins = stage.load_references(stage.REFERENCE, 1000000, 50)
        lengths = {}
        for chrom, end in bins:
            lengths[chrom] = max(lengths.get(chrom, 0), end + 1000)
        lengths["chrM"] = 16569
        fasta = self.work / "synthetic.fa"
        fasta.write_text(">placeholder\nN\n")
        stage.Path(str(fasta) + ".fai").write_text(
            "".join(f"{c}\t{n}\t0\t0\t0\n" for c, n in lengths.items()))
        header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
            f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in lengths.items())
        rows = ["sample_id\tbam\tbai\tpass_qc"]
        for name in CELLS:
            self.make_bam(name, bins, lengths, header)
            rows.append(f"{name}\tbam/{name}.bam\tbam/{name}.bam.bai\ttrue")
        rows.append("Failed_cell\tbam/absent.bam\tbam/absent.bam.bai\tfalse")
        manifest = self.work / "results/manifests/passed_cells_only_manifest.tsv"
        manifest.write_text("\n".join(rows) + "\n")
        outdir = self.work / "results/DNA_Illumina_parity_smoke"
        self.assertEqual(stage.main(["--manifest", str(manifest), "--outdir", str(outdir),
                                     "--ref", str(fasta), "--mode", "local"]), 0)
        root = outdir / "ginkgo"
        cohort = root / "cohort"
        self.assertTrue((root / "COMPLETE").is_file())
        with (root / "excluded_cells.tsv").open() as stream:
            self.assertEqual([r["sample_name"] for r in csv.DictReader(stream, delimiter="\t")],
                             ["Failed_cell"])
        self.assertEqual((cohort / "ploidy.txt").read_text(),
                         "".join(f"{c}\t2\n" for c in CELLS))
        with (cohort / "SegCopy.binsize_1000000.tsv").open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["CHR", "START", "END", *CELLS])
            segcopy = list(reader)
        self.assertEqual(len(segcopy), len(bins))
        chr8 = [r for r in segcopy if r["CHR"] == "chr8"]
        median = lambda values: sorted(values)[len(values) // 2]
        self.assertEqual(median([int(r["Cell_B_gain8"]) for r in chr8]), 3)
        self.assertEqual(median([int(r["Cell_A"]) for r in chr8]), 2)
        with tarfile.open(cohort / "cnv_plots_binsize_1000000.tar.gz") as archive:
            self.assertEqual(sorted(archive.getnames()),
                             sorted(f"cnv_binsize_1000000_{c}_CN.jpeg" for c in CELLS))
        for name in CELLS:
            self.assertTrue((cohort / "per_sample" / f"{name}_clouds.tsv").stat().st_size)
            with (root / "cells" / name / "chromosome_audit.tsv").open() as stream:
                audit = {r["bam_chromosome"]: r for r in csv.DictReader(stream, delimiter="\t")}
            self.assertEqual(audit["chrM"]["disposition"], "not_in_source_reference")
            self.assertEqual(audit["chrM"]["bed_records"], "1")
        for name in ("CNV1_binsize_1000000.tsv", "CNV2_binsize_1000000.tsv",
                     "Legacy-GinkgoSegmentSummary.txt", "sampleploidy.tsv", "SegNorm"):
            self.assertTrue((cohort / name).stat().st_size, name)
        with (root / "source_metrics.tsv").open() as stream:
            metrics = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual([m["sample_name"] for m in metrics], list(CELLS))
        for row in metrics:
            for key in ("MAPD_CNV_Log2", "SKEW_CNV"):
                self.assertTrue(row[key] != "NA" or row[f"{key}_reason"])
            self.assertGreater(float(row["MAPD_CNV_Log2"]), 0)
            self.assertTrue(0 <= float(row["SKEW_CNV"]) <= 1)
        with (root / "ginkgo_metrics.tsv").open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            self.assertEqual(reader.fieldnames[0], "sample_name")
            self.assertIn("GenomePloidy", reader.fieldnames)
            self.assertEqual([r["sample_name"] for r in reader], list(CELLS))


if __name__ == "__main__":
    unittest.main()
