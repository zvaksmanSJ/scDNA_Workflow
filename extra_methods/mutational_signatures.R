#!/usr/bin/env Rscript
# Mutational signature refitting via MutationalPatterns (Bioconductor, GPL-3,
# free) -- the established R package for this task. Given this cohort has no
# matched-bulk sample, there is no curated "confirmed clonal somatic-only"
# call set; the best available proxy for cell-private (non-recurrent, likely
# PTA/clonal-acquisition) SNVs is the SINGLETON bucket from our own cross-cell
# artifact panel (sites observed as PASS in exactly one of the 15 cells).
# This is NOT equivalent to a true matched-bulk-validated private-somatic
# call set -- it will still contain some residual allele-dropout/amplification
# noise that a matched bulk would otherwise let us subtract -- and this
# caveat is carried through to the dashboard text.
#
# Because per-cell mutation counts here (several hundred to ~1,600) are too
# low to support de novo NMF signature discovery, this uses refitting
# (NNLS) against the published COSMIC reference signature catalog instead,
# via MutationalPatterns::fit_to_signatures_strict() (backward-elimination
# refit recommended by the package authors for exactly this low-count
# regime, to avoid overfitting many near-zero-contribution signatures).
#
# Usage: Rscript mutational_signatures.R <artifact_panel_sites.tsv> <output_dir>

suppressMessages({
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
sites_path <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
bases <- c("A", "C", "G", "T")

cat("[INFO] Reading", sites_path, "\n")
sites <- read.delim(sites_path, stringsAsFactors = FALSE)

snv <- sites[sites$n_cells == 1 &
             sites$classification == "SINGLETON" &
             nchar(sites$ref) == 1 &
             nchar(sites$alt) == 1 &
             sites$ref %in% bases &
             sites$alt %in% bases, ]
snv$sample_id <- sapply(strsplit(snv$cells, ","), `[`, 1)
cat("[INFO]", nrow(snv), "singleton SNV sites across", length(unique(snv$sample_id)), "cells\n")

samples <- sort(unique(snv$sample_id))
vcf_dir <- file.path(output_dir, "per_cell_vcfs")
dir.create(vcf_dir, showWarnings = FALSE, recursive = TRUE)

vcf_header <- c(
  "##fileformat=VCFv4.2",
  paste0("##contig=<ID=chr", c(1:22, "X", "Y"), ">"),
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

vcf_files <- character(length(samples))
for (i in seq_along(samples)) {
  s <- samples[i]
  d <- snv[snv$sample_id == s, ]
  d <- d[order(d$chrom, d$pos), ]
  records <- sprintf("chr%s\t%d\t.\t%s\t%s\t.\tPASS\t.", d$chrom, d$pos, d$ref, d$alt)
  vcf_path <- file.path(vcf_dir, paste0(s, ".vcf"))
  writeLines(c(vcf_header, records), vcf_path)
  vcf_files[i] <- vcf_path
}

cat("[INFO] Building 96-trinucleotide-context mutation matrix\n")
grl <- read_vcfs_as_granges(vcf_files, samples, genome = ref_genome, type = "snv")
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = get(ref_genome))
write.table(mut_mat, file.path(output_dir, "mutation_matrix_96.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("[INFO] Fetching COSMIC reference signatures (bundled with MutationalPatterns)\n")
signatures <- get_known_signatures(muttype = "snv", genome = "GRCh38")

cat("[INFO] Refitting via backward-elimination NNLS (fit_to_signatures_strict)\n")
fit_res <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
contrib <- fit_res$fit_res$contribution

contrib_df <- as.data.frame(t(contrib))
contrib_df$sample_id <- rownames(contrib_df)
contrib_df <- contrib_df[, c("sample_id", setdiff(colnames(contrib_df), "sample_id"))]
write.table(contrib_df, file.path(output_dir, "signature_contributions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rel_contrib <- sweep(contrib, 2, colSums(contrib), FUN = "/")
rel_contrib[is.na(rel_contrib)] <- 0
rel_df <- as.data.frame(t(rel_contrib))
rel_df$sample_id <- rownames(rel_df)
rel_df <- rel_df[, c("sample_id", setdiff(colnames(rel_df), "sample_id"))]
write.table(rel_df, file.path(output_dir, "signature_relative_contributions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] Cosine similarity between reconstructed and original profiles\n")
cos_sim <- diag(cos_sim_matrix(mut_mat, fit_res$fit_res$reconstructed))
cos_df <- data.frame(sample_id = colnames(mut_mat), cosine_similarity_reconstructed = cos_sim)
write.table(cos_df, file.path(output_dir, "reconstruction_cosine_similarity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] Writing per-cell 96-profile and signature-contribution plots\n")
png(file.path(output_dir, "signature_contribution_barplot.png"), width = 1400, height = 800, res = 130)
print(plot_contribution(contrib, mode = "relative", coord_flip = TRUE))
dev.off()

png(file.path(output_dir, "mutation_spectrum_96.png"), width = 1600, height = 1800, res = 130)
print(plot_96_profile(mut_mat, condensed = TRUE))
dev.off()

cat("[OK] Mutational signature analysis complete ->", output_dir, "\n")
