# Purpose: independently validate discovery-cohort Spearman SSN perturbations in 50 matched TCGA
# tumour-normal pairs using patient-level paired contrasts.
# Inputs: data/tcga_matched_pairs.qs2,
#         results/02_Define_LR_Edges/lr_edges_ssn.qs2,
#         results/04_Adjust_LR_Expression/*_lr_expression_adjusted.qs2,
#         results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2, and
#         results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2.
# Outputs: results/07_Validate_Paired_SSN/ssn_matched_pair_perturbation.qs2,
#          results/07_Validate_Paired_SSN/ssn_matched_pair_bootstrap_mean.qs2,
#          results/07_Validate_Paired_SSN/ssn_matched_pair_validation.qs2, and
#          results/07_Validate_Paired_SSN/ssn_matched_pair_provenance.qs2.
#
# Run from the repository root with Rscript scripts/07_Validate_Paired_SSN.R
# Each patient-level contrast is Delta rho(tumour) - Delta rho(matched normal), with both samples
# evaluated against the same composition-adjusted GTEx normal reference.

library(qs2)

source("scripts/00_Functions_SSN.R")

dir.create("results/07_Validate_Paired_SSN", recursive = TRUE, showWarnings = FALSE)

set.seed(20260731)
n_pair_bootstrap <- 1000L
n_paired_sign_permutation <- 10000L

matched_pairs <- qs2::qs_read("data/tcga_matched_pairs.qs2")
lr_edges_ssn <- qs2::qs_read("results/02_Define_LR_Edges/lr_edges_ssn.qs2")
gtex_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2")
tcga_normal_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/tcga_normal_lr_expression_adjusted.qs2")
ssn_tumor_perturbation <- qs2::qs_read("results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2")
ssn_cohort_significance <- qs2::qs_read("results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2")

matched_pairs <- matched_pairs[
  matched_pairs$tumor_id %in% colnames(ssn_tumor_perturbation) &
    matched_pairs$normal_id %in% colnames(tcga_normal_expression_adjusted),
  , drop = FALSE
]
matched_pairs <- matched_pairs[order(matched_pairs$patient_id), , drop = FALSE]
lr_edges_ssn <- lr_edges_ssn[match(rownames(ssn_tumor_perturbation), lr_edges_ssn$ssn_edge_id), , drop = FALSE]
ssn_cohort_significance <- ssn_cohort_significance[
  match(rownames(ssn_tumor_perturbation), ssn_cohort_significance$ssn_edge_id), , drop = FALSE
]
stopifnot(
  nrow(matched_pairs) == 50L,
  !anyDuplicated(matched_pairs$patient_id),
  !anyDuplicated(matched_pairs$tumor_id),
  !anyDuplicated(matched_pairs$normal_id),
  identical(lr_edges_ssn$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  identical(ssn_cohort_significance$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  identical(rownames(gtex_expression_adjusted), rownames(tcga_normal_expression_adjusted)),
  all(lr_edges_ssn$Gene1 %in% rownames(gtex_expression_adjusted)),
  all(lr_edges_ssn$Gene2 %in% rownames(gtex_expression_adjusted))
)

paired_tumor_perturbation <- ssn_tumor_perturbation[, matched_pairs$tumor_id, drop = FALSE]
colnames(paired_tumor_perturbation) <- matched_pairs$patient_id
paired_normal_expression <- tcga_normal_expression_adjusted[, matched_pairs$normal_id, drop = FALSE]
colnames(paired_normal_expression) <- matched_pairs$patient_id

message("Computing matched-normal perturbations for ", ncol(paired_normal_expression), " TCGA pairs.")
paired_normal_result <- calculate_spearman_ssn_perturbation(
  gtex_expression_adjusted, paired_normal_expression, lr_edges_ssn
)
stopifnot(
  identical(paired_normal_result$reference_spearman$ssn_edge_id, rownames(paired_tumor_perturbation)),
  nrow(paired_normal_result$dropped_edges) == 0L,
  all(is.finite(paired_normal_result$perturbation))
)
paired_normal_perturbation <- paired_normal_result$perturbation
paired_contrast <- paired_tumor_perturbation - paired_normal_perturbation
stopifnot(all(is.finite(paired_contrast)))

observed_mean_paired_contrast <- rowMeans(paired_contrast)
paired_bootstrap_mean <- vapply(seq_len(n_pair_bootstrap), function(i) {
  rowMeans(paired_contrast[, sample.int(ncol(paired_contrast), replace = TRUE), drop = FALSE])
}, numeric(nrow(paired_contrast)))
rownames(paired_bootstrap_mean) <- rownames(paired_contrast)
paired_bootstrap_ci <- t(apply(paired_bootstrap_mean, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
paired_bootstrap_direction_consistency <- rowMeans(
  sweep(paired_bootstrap_mean, 1, sign(observed_mean_paired_contrast), "*") > 0
)

paired_permutation_signs <- matrix(
  sample(c(-1, 1), ncol(paired_contrast) * n_paired_sign_permutation, replace = TRUE),
  nrow = ncol(paired_contrast)
)
paired_sign_permutation_mean <- t(paired_permutation_signs) %*% t(paired_contrast) / ncol(paired_contrast)
paired_empirical_p_value <- (colSums(abs(paired_sign_permutation_mean) >= abs(observed_mean_paired_contrast)) + 1) /
  (n_paired_sign_permutation + 1)

primary_direction <- sign(ssn_cohort_significance$mean_delta_rho)
paired_direction <- sign(observed_mean_paired_contrast)
direction_matches_primary <- primary_direction == paired_direction
direction_matches_primary[primary_direction == 0 | paired_direction == 0] <- NA

ssn_matched_pair_validation <- cbind(
  ssn_cohort_significance[, c(
    "Gene1", "Gene2", "ssn_edge_id", "direction_support", "n_directional_entries",
    "directional_lr_entries", "mean_delta_rho", "empirical_sign_permutation_fdr_bh"
  )],
  n_matched_pairs = ncol(paired_contrast),
  mean_tumor_delta_rho = rowMeans(paired_tumor_perturbation),
  mean_matched_normal_delta_rho = rowMeans(paired_normal_perturbation),
  mean_paired_contrast = observed_mean_paired_contrast,
  median_paired_contrast = apply(paired_contrast, 1, stats::median),
  paired_contrast_ci_lower = paired_bootstrap_ci[, 1],
  paired_contrast_ci_upper = paired_bootstrap_ci[, 2],
  paired_bootstrap_direction_consistency = paired_bootstrap_direction_consistency,
  paired_empirical_sign_permutation_p = paired_empirical_p_value,
  paired_empirical_sign_permutation_fdr_bh = stats::p.adjust(paired_empirical_p_value, method = "BH"),
  paired_direction = paired_direction,
  direction_matches_primary = direction_matches_primary,
  paired_direction_consistency = rowMeans(sweep(paired_contrast, 1, paired_direction, "*") > 0)
)
rownames(ssn_matched_pair_validation) <- ssn_cohort_significance$ssn_edge_id

ssn_matched_pair_perturbation <- list(
  matched_pairs = matched_pairs,
  tumor_delta_rho = paired_tumor_perturbation,
  matched_normal_delta_rho = paired_normal_perturbation,
  paired_contrast = paired_contrast
)
paired_provenance <- data.frame(
  item = c(
    "random_seed", "correlation_method", "reference_cohort", "reference_sample_count",
    "matched_pair_count", "paired_contrast", "tumor_perturbation_source", "normal_perturbation_method",
    "pair_bootstrap_iterations", "pair_bootstrap_unit", "paired_sign_permutation_iterations",
    "paired_sign_permutation_null", "multiple_testing", "fdr_filtering"
  ),
  value = c(
    "20260731", "Spearman rank correlation", "GTEx colon normal", as.character(ncol(gtex_expression_adjusted)),
    as.character(ncol(paired_contrast)), "Delta rho(tumour) - Delta rho(matched normal)",
    "subset of results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2",
    "same GTEx reference; calculate_spearman_ssn_perturbation() on composition-adjusted matched-normal expression",
    as.character(n_pair_bootstrap), "patient-level matched tumour-normal pair",
    as.character(n_paired_sign_permutation), "random sign reversal of each patient-level paired contrast",
    "Benjamini-Hochberg across retained SSN edges", "none; validation supports interpretation only"
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(ssn_matched_pair_perturbation, "results/07_Validate_Paired_SSN/ssn_matched_pair_perturbation.qs2")
qs2::qs_save(paired_bootstrap_mean, "results/07_Validate_Paired_SSN/ssn_matched_pair_bootstrap_mean.qs2")
qs2::qs_save(ssn_matched_pair_validation, "results/07_Validate_Paired_SSN/ssn_matched_pair_validation.qs2")
qs2::qs_save(paired_provenance, "results/07_Validate_Paired_SSN/ssn_matched_pair_provenance.qs2")

message("Saved paired SSN validation for ", nrow(ssn_matched_pair_validation), " edges across ",
        ncol(paired_contrast), " matched tumour-normal pairs.")
