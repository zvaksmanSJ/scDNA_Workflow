#!/usr/bin/env Rscript
# Lineage verification using ape (Paradis & Schliep 2019, Bioinformatics),
# one of the most widely used phylogenetics packages in R -- already a
# dependency of the Sequoia lineage-building script itself, so no new
# software is introduced.
#
# Two checks, both established/standard rather than custom-invented:
#   1. Bootstrap support: ape::read.tree natively parses MPBoot's embedded
#      ultrafast-bootstrap internal-node labels into node.label. Any
#      internal branch below --threshold is flagged.
#   2. Topology stability: ape::dist.topo(method="PH85") computes the
#      standard Robinson-Foulds distance (Robinson & Foulds 1981) between
#      the primary tree and an independently re-subsampled tree built from
#      a different random seed. 0 = identical topology; higher = more
#      disagreement. This is the same class of metric (RF/quartet-based
#      tree distance) the reference paper itself used (TreeDist's
#      JaccardRobinsonFoulds) to compare bulk- vs single-cell-derived
#      phylogenies.

suppressWarnings(suppressMessages(library(ape)))
args <- commandArgs(trailingOnly = TRUE)
primary_path <- args[1]
threshold <- as.numeric(args[2])
bootstrap_out <- args[3]
summary_out <- args[4]
resample_path <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NA

primary <- read.tree(primary_path)

node_labels <- primary$node.label
support <- suppressWarnings(as.numeric(node_labels))
n_internal <- length(node_labels)
low_conf <- !is.na(support) & support < threshold

bootstrap_df <- data.frame(
  internal_node_index = seq_len(n_internal),
  bootstrap_support = ifelse(is.na(support), "NA", as.character(support)),
  low_confidence = low_conf,
  stringsAsFactors = FALSE
)
write.table(bootstrap_df, bootstrap_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[OK] Bootstrap support report -> %s\n", bootstrap_out))
cat(sprintf("[INFO] %d internal branches, %d below threshold %.0f\n",
            n_internal, sum(low_conf, na.rm = TRUE), threshold))

result <- list(
  n_internal_branches = n_internal,
  n_low_confidence_branches = as.integer(sum(low_conf, na.rm = TRUE)),
  bootstrap_threshold = threshold,
  min_bootstrap_support = if (all(is.na(support))) NA else min(support, na.rm = TRUE)
)

if (!is.na(resample_path)) {
  resample <- read.tree(resample_path)
  same_tips <- setequal(primary$tip.label, resample$tip.label)
  result$same_leaf_set <- same_tips
  if (same_tips) {
    resample <- keep.tip(resample, primary$tip.label)
    rf <- dist.topo(primary, resample, method = "PH85")
    bsd <- tryCatch(dist.topo(primary, resample, method = "score"), error = function(e) NA)
    result$robinson_foulds_distance <- as.numeric(rf)
    result$branch_score_distance <- if (is.na(bsd)) NA else as.numeric(bsd)
    cat(sprintf("[INFO] Robinson-Foulds distance (primary vs. independent re-subsample): %s (0 = identical topology)\n",
                as.character(as.numeric(rf))))
  } else {
    cat("[WARN] Leaf sets differ between primary and resample tree; RF distance not computed\n")
    result$robinson_foulds_distance <- NA
  }
}

json_lines <- c("{")
kv <- c()
for (name in names(result)) {
  val <- result[[name]]
  if (is.na(val)) {
    kv <- c(kv, sprintf('  "%s": null', name))
  } else if (is.logical(val)) {
    kv <- c(kv, sprintf('  "%s": %s', name, tolower(as.character(val))))
  } else if (is.numeric(val)) {
    kv <- c(kv, sprintf('  "%s": %s', name, as.character(val)))
  } else {
    kv <- c(kv, sprintf('  "%s": "%s"', name, val))
  }
}
json_lines <- c(json_lines, paste(kv, collapse = ",\n"), "}")
writeLines(json_lines, summary_out)
cat(sprintf("[OK] Verification summary -> %s\n", summary_out))
