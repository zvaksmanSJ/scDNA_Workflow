"""Real synthetic-input smoke; never reads or submits the production cohort."""
import csv
import gzip
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import step1_submit_resolveome_alignment_IluminaOnly_lsf as alignment


class PreprocessingSmoke(unittest.TestCase):
    def test_real_fastp_bwa_picard_preseq_chain(self):
        env = Path(sys.executable).parent
        with tempfile.TemporaryDirectory(prefix="resolveome_preprocess_smoke_") as temp:
            root = Path(temp)
            rng = random.Random(317)
            genome = "".join(rng.choices("ACGT", k=200000))
            ref = root / "synthetic.fa"
            ref.write_text(">1\n" + genome + "\n")
            process_env = {**os.environ, "PATH": str(env) + os.pathsep + os.environ["PATH"]}
            def execute(command):
                with (root / "smoke.log").open("a") as log:
                    result = subprocess.run(list(map(str, command)), env=process_env,
                                            stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    self.fail((root / "smoke.log").read_text()[-7000:])
            execute([env / "samtools", "faidx", ref])
            execute([env / "bwa-mem2", "index", ref])
            execute([env / "picard", "CreateSequenceDictionary", f"R={ref}",
                     f"O={ref.with_suffix('.dict')}"])
            r1, r2 = root / "R1.fastq.gz", root / "R2.fastq.gz"
            complement = str.maketrans("ACGT", "TGCA")
            with gzip.open(r1, "wt") as left, gzip.open(r2, "wt") as right:
                for index in range(5000):
                    start = rng.randrange(100, len(genome) - 500)
                    insert = rng.randrange(240, 340)
                    seq1 = genome[start:start + 150]
                    seq2 = genome[start + insert - 150:start + insert].translate(complement)[::-1]
                    # Include known G tails; the synthetic reference is random.
                    if index % 3 == 0:
                        seq1, seq2 = seq1[:130] + "G" * 20, seq2[:130] + "G" * 20
                    left.write(f"@smoke:{index}/1\n{seq1}\n+\n{'I'*150}\n")
                    right.write(f"@smoke:{index}/2\n{seq2}\n+\n{'I'*150}\n")
            out = root / "out"
            (out / "manifests").mkdir(parents=True)
            manifest = out / "manifests/aligned.tsv"
            manifest.write_text("sample_id\tbam\tbai\tflagstat\tidxstats\n")
            bed = root / "intervals.bed"
            bed.write_text("1\t0\t200000\n")
            worker = root / "worker.sh"
            alignment.write_worker_script(worker)
            execute(["bash", worker, "Synthetic", r1, r2, out, ref, "2", "resolveome-align",
                     "ILLUMINA", "DNA", manifest, "100M", "0", bed])
            sample = out / "bam/Synthetic"
            self.assertTrue((sample / "alignment.complete").is_file())
            execute([env / "samtools", "quickcheck", sample / "Synthetic.dedup.bam"])
            prefix = out / "qc/picard/Synthetic/Synthetic"
            for suffix in ("alignment_summary_metrics", "gc_bias.summary_metrics",
                           "quality_yield_metrics", "wgs_metrics", "pre_dedup.flagstat.txt"):
                self.assertGreater(Path(str(prefix) + "." + suffix).stat().st_size, 0)
            with manifest.open() as stream:
                row = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(row["flagstat"], str(prefix) + ".pre_dedup.flagstat.txt")
            self.assertTrue(list((out / "qc/fastqc_raw").rglob("*_fastqc.zip")))
            self.assertTrue(list((out / "qc/fastqc_trimmed").rglob("*_fastqc.zip")))
            trimmed = out / "trimmed_fastq/Synthetic/Synthetic_R1_trim.fastq.gz"
            shortened_g_tails = 0
            with gzip.open(trimmed, "rt") as stream:
                for header in stream:
                    sequence = next(stream).strip()
                    next(stream)
                    next(stream)
                    index = int(header.split(":")[1].split("/")[0])
                    if index % 3 == 0 and len(sequence) <= 130:
                        shortened_g_tails += 1
            self.assertGreater(shortened_g_tails, 1500)
            self.assertTrue((out / "qc/preseq/Synthetic/Synthetic_preseq.txt").is_file()
                            or (out / "qc/preseq/Synthetic/Synthetic.preseq.status.tsv").is_file())
            here = Path(__file__).parent
            execute([sys.executable, here / "step2_submit_resolveome_barcode_cellcall_IluminaOnly_lsf.py",
                     "--aligned-manifest", manifest, "--outdir", out, "--local",
                     "--min-total-reads", "1"])
            execute([sys.executable, here / "step2_library_qc_dashboard.py",
                     "--outdir", out, "--ref", ref])
            self.assertGreater((out / "qc/library_dashboard/source_threshold_dashboard.html").stat().st_size, 0)
            self.assertGreater((out / "qc/library_dashboard/multiqc_library_qc.html").stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
