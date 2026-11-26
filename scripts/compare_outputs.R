# Compare outputs of a complete run with the study release. Run from the repository root.
# Reports numerical deviations and factor alignment. No outputs are overwritten.
library(qs2)
pairs <- qs_read("results/02_Define_LR_Edges/lr_pairs_directed.qs2")
ref_pairs <- read.delim("reference/lr_pairs_directed.tsv", check.names = FALSE)
stopifnot(identical(pairs, ref_pairs) || isTRUE(all.equal(pairs, ref_pairs, check.attributes = FALSE)))
edges <- qs_read("results/02_Define_LR_Edges/lr_edges_ssn.qs2")
stopifnot(identical(edges, qs_read("reference/lr_edges_ssn.qs2")))
ssn <- qs_read("results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2")
reference_ssn <- qs_read("reference/ssn_tumor_spearman_perturbation.qs2")
stopifnot(identical(dimnames(ssn), dimnames(reference_ssn)))
ssn_error <- max(abs(ssn - reference_ssn))
print(data.frame(edges = nrow(ssn), tumors = ncol(ssn), max_delta_rho_error = ssn_error))
stopifnot(ssn_error < 1e-12)

cohort <- qs_read("results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2")
paired <- qs_read("results/07_Validate_Paired_SSN/ssn_matched_pair_validation.qs2")
observed <- list(cohort = cohort, paired = paired)
expected <- list(cohort = read.delim("reference/ssn_cohort_significance.tsv"),
                 paired = read.delim("reference/ssn_matched_pair_validation.tsv"))
for (stage in names(observed)) {
  x <- observed[[stage]]
  y <- expected[[stage]]
  stopifnot(identical(x$ssn_edge_id, y$ssn_edge_id))
  numeric_fields <- names(x)[vapply(x, is.numeric, logical(1))]
  other_fields <- setdiff(names(x), numeric_fields)
  stopifnot(identical(names(x), names(y)),
            all(vapply(other_fields, function(f) identical(x[[f]], y[[f]]), logical(1))))
  difference <- max(abs(as.matrix(x[, numeric_fields]) - as.matrix(y[, numeric_fields])))
  print(data.frame(stage = stage, max_numeric_difference = difference))
  stopifnot(difference < 1e-10)
}
stability <- qs_read("results/06_Assess_SSN_Stability/ssn_reference_size_stability.qs2")
sizes <- aggregate(cbind(spearman_correlation_to_full_reference, direction_agreement_to_full_reference)
                   ~ reference_size, stability$resample_diagnostics, median)
ref_sizes <- read.delim("reference/reference_size_medians.tsv")
stopifnot(isTRUE(all.equal(sizes, ref_sizes, tolerance = 1e-10, check.attributes = FALSE)))

# Factor signs and numbering can vary by software environment. Report the full best match.
scores <- qs_read("results/09_Run_MOFA/mofa_factor_scores.qs2")
ref_scores <- read.delim("reference/mofa_factor_scores.tsv", check.names = FALSE)
weights <- qs_read("results/09_Run_MOFA/mofa_edge_loadings.qs2")
ref_weights <- read.delim("reference/mofa_edge_loadings.tsv", check.names = FALSE)
stopifnot(identical(scores$sample_id, ref_scores$sample_id),
          identical(weights$ssn_edge_id, ref_weights$ssn_edge_id))
factors <- setdiff(names(scores), "sample_id")
rho <- cor(ref_scores[, -1], scores[, factors])
match_index <- apply(abs(rho), 1, which.max)
index <- cbind(seq_len(nrow(rho)), match_index)
weight_rho <- cor(ref_weights[, -1], weights[, factors])
stopifnot(!anyDuplicated(match_index))
os <- qs_read("results/09_Run_MOFA/mofa_factor_os_association.qs2")
ref_os <- read.delim("reference/mofa_factor_os_association.tsv")
os <- os[match(factors[match_index], os$factor), ]
ref_os <- ref_os[match(rownames(rho), ref_os$factor), ]
comparison <- data.frame(study = rownames(rho), refit = factors[match_index],
                         score_correlation = rho[index], loading_correlation = weight_rho[index],
                         refit_BH = os$km_logrank_fdr_bh, study_BH = ref_os$km_logrank_fdr_bh)
print(comparison, row.names = FALSE)
study_top <- comparison$study[which.min(comparison$study_BH)]
refit_top <- comparison$study[which.min(comparison$refit_BH)]
cat("Study factor prioritised in this refit:", refit_top, "\n")
variance <- qs_read("results/09_Run_MOFA/mofa_factor_variance.qs2")
ref_variance <- read.delim("reference/mofa_factor_variance.tsv")
refit_largest <- as.character(variance$factor[which.max(variance$explained_variance)])
study_largest <- ref_variance$factor[which.max(ref_variance$explained_variance)]
# Release-agreement checks on the fitted outputs. All ten factors enter the analysis.
stopifnot(all(abs(comparison$score_correlation) > 0.999),
          all(abs(comparison$loading_correlation) > 0.999),
          all(sign(comparison$score_correlation) == sign(comparison$loading_correlation)),
          max(abs(comparison$refit_BH - comparison$study_BH)) < 1e-8,
          identical(refit_top, study_top),
          identical(comparison$study[comparison$refit == refit_largest], study_largest))
