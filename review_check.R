# Check sample-specific network perturbations and post-fit survival ranking against reference data.
# Run Rscript --vanilla review_check.R [--all]. Prints comparisons and writes no files.
library(qs2)
library(survival)
source("scripts/00_Functions_SSN.R")

edges <- qs_read("reference/lr_edges_ssn.qs2")
normal <- qs_read("reference/gtex_lr_expression_adjusted.qs2")
tumor <- qs_read("reference/tcga_tumor_lr_expression_adjusted.qs2")
reported <- qs_read("reference/ssn_tumor_spearman_perturbation.qs2")
stopifnot(ncol(normal) == 898L, identical(dim(reported), c(238L, 618L)),
          identical(colnames(tumor), colnames(reported)),
          identical(edges$ssn_edge_id, rownames(reported)))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--all"))
idx <- if ("--all" %in% args) seq_len(ncol(tumor)) else c(1L, 309L, 618L)
replay <- calculate_spearman_ssn_perturbation(normal, tumor[, idx, drop = FALSE], edges)
expected <- reported[, idx, drop = FALSE]
error <- max(abs(replay$perturbation - expected))
stopifnot(identical(dimnames(replay$perturbation), dimnames(expected)), error < 1e-12)
print(data.frame(edges = nrow(expected), tumors = ncol(expected), max_absolute_error = error))

# Each factor is split at its median in patients with observed OS, as in script 09.
scores <- read.delim("reference/mofa_factor_scores.tsv", check.names = FALSE)
clinical <- qs_read("data/tcga_overall_survival.qs2")
clinical <- clinical[complete.cases(clinical[, c("OS", "OS.time")]), ]
archived <- read.delim("reference/mofa_factor_os_association.tsv")
stopifnot(nrow(clinical) == 584L, sum(clinical$OS) == 121L, !anyDuplicated(clinical$ID),
          all(clinical$ID %in% scores$sample_id), length(archived$factor) == 10L)
scores <- scores[match(clinical$ID, scores$sample_id), ]
p <- vapply(archived$factor, function(f) {
  group <- scores[[f]] > median(scores[[f]])
  stat <- survdiff(Surv(OS.time, OS) ~ group, data = clinical)
  pchisq(stat$chisq, length(stat$n) - 1L, lower.tail = FALSE)
}, numeric(1))
adjusted <- p.adjust(p, "BH")
stopifnot(max(abs(p - archived$km_logrank_p)) < 1e-10,
          max(abs(adjusted - archived$km_logrank_fdr_bh)) < 1e-10)
print(data.frame(factor = names(which.min(adjusted)), logrank_BH = min(adjusted)))
