# Purpose: prepare the complete tumour SSN perturbation matrix, aligned tumour metadata, and
# evidence annotations for unsupervised MOFA factorisation.
# Inputs: data/tcga_selected_sample_map.qs2,
#         results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2,
#         results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2,
#         results/07_Validate_Paired_SSN/ssn_matched_pair_validation.qs2, and
#         data/tcga_overall_survival.qs2.
# Outputs: results/08_Prepare_MOFA_Input/mofa_ssn_input_scaled.qs2,
#          results/08_Prepare_MOFA_Input/mofa_tumor_metadata.qs2,
#          results/08_Prepare_MOFA_Input/mofa_edge_annotation.qs2,
#          results/08_Prepare_MOFA_Input/mofa_input_qc.qs2, and
#          results/08_Prepare_MOFA_Input/mofa_input_provenance.qs2.
#
# Run from the repository root with Rscript scripts/08_Prepare_MOFA_Input.R
# All retained SSN edges are used. Cohort and paired statistics are attached solely for later
# interpretation of MOFA loadings and never determine whether an edge enters the factorisation.

library(qs2)

dir.create("results/08_Prepare_MOFA_Input", recursive = TRUE, showWarnings = FALSE)

ssn_tumor_perturbation <- qs2::qs_read("results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2")
tcga_sample_map <- qs2::qs_read("data/tcga_selected_sample_map.qs2")
ssn_cohort_significance <- qs2::qs_read("results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2")
ssn_matched_pair_validation <- qs2::qs_read("results/07_Validate_Paired_SSN/ssn_matched_pair_validation.qs2")
tcga_clinical <- qs_read("data/tcga_overall_survival.qs2")

tcga_tumor_map <- tcga_sample_map[tcga_sample_map$sample_type == "01", , drop = FALSE]
tcga_tumor_map <- tcga_tumor_map[match(colnames(ssn_tumor_perturbation), tcga_tumor_map$sample_id), , drop = FALSE]
ssn_cohort_significance <- ssn_cohort_significance[
  match(rownames(ssn_tumor_perturbation), ssn_cohort_significance$ssn_edge_id), , drop = FALSE
]
ssn_matched_pair_validation <- ssn_matched_pair_validation[
  match(rownames(ssn_tumor_perturbation), ssn_matched_pair_validation$ssn_edge_id), , drop = FALSE
]
stopifnot(
  identical(tcga_tumor_map$sample_id, colnames(ssn_tumor_perturbation)),
  !anyDuplicated(tcga_tumor_map$sample_id),
  !anyDuplicated(tcga_tumor_map$patient_id),
  identical(ssn_cohort_significance$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  identical(ssn_matched_pair_validation$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  all(is.finite(ssn_tumor_perturbation))
)

sample_sd_before_scaling <- apply(ssn_tumor_perturbation, 2, stats::sd)
edge_sd_before_scaling <- apply(ssn_tumor_perturbation, 1, stats::sd)
stopifnot(
  all(is.finite(sample_sd_before_scaling)), all(sample_sd_before_scaling > 0),
  all(is.finite(edge_sd_before_scaling)), all(edge_sd_before_scaling > 0)
)

mofa_ssn_input_scaled <- scale(ssn_tumor_perturbation)
dimnames(mofa_ssn_input_scaled) <- dimnames(ssn_tumor_perturbation)
stopifnot(all(is.finite(mofa_ssn_input_scaled)))

tcga_clinical <- tcga_clinical[match(tcga_tumor_map$sample_id, tcga_clinical$ID), , drop = FALSE]
mofa_tumor_metadata <- cbind(tcga_tumor_map, tcga_clinical[, setdiff(names(tcga_clinical), "ID"), drop = FALSE])
mofa_tumor_metadata$has_os_annotation <- complete.cases(mofa_tumor_metadata[, c("OS", "OS.time")])
stopifnot(identical(mofa_tumor_metadata$sample_id, colnames(mofa_ssn_input_scaled)))

mofa_edge_annotation <- cbind(
  ssn_cohort_significance,
  ssn_matched_pair_validation[, c(
    "mean_paired_contrast", "paired_contrast_ci_lower", "paired_contrast_ci_upper",
    "paired_bootstrap_direction_consistency", "paired_empirical_sign_permutation_p",
    "paired_empirical_sign_permutation_fdr_bh", "direction_matches_primary"
  ), drop = FALSE]
)
rownames(mofa_edge_annotation) <- mofa_edge_annotation$ssn_edge_id
stopifnot(identical(rownames(mofa_edge_annotation), rownames(mofa_ssn_input_scaled)))

mofa_input_qc <- list(
  feature_summary = data.frame(
    n_ssn_edges = nrow(mofa_ssn_input_scaled),
    n_zero_variance_edges_before_scaling = sum(edge_sd_before_scaling == 0),
    n_nonfinite_input_values = sum(!is.finite(ssn_tumor_perturbation))
  ),
  sample_summary = data.frame(
    n_tcga_primary_tumors = ncol(mofa_ssn_input_scaled),
    n_zero_variance_samples_before_scaling = sum(sample_sd_before_scaling == 0),
    n_samples_with_os_annotation = sum(mofa_tumor_metadata$has_os_annotation),
    n_samples_without_os_annotation = sum(!mofa_tumor_metadata$has_os_annotation)
  ),
  sample_scaling = data.frame(
    sample_id = colnames(mofa_ssn_input_scaled),
    mean_after_scaling = colMeans(mofa_ssn_input_scaled),
    sd_after_scaling = apply(mofa_ssn_input_scaled, 2, stats::sd)
  )
)

mofa_input_provenance <- data.frame(
  item = c(
    "mofa_input", "n_ssn_edges", "n_tcga_primary_tumors", "matrix_orientation",
    "scaling_method", "edge_filtering", "cohort_statistics_role", "paired_statistics_role",
    "clinical_metadata_source", "clinical_annotation_handling"
  ),
  value = c(
    "composition-adjusted Spearman SSN tumour perturbation matrix", as.character(nrow(mofa_ssn_input_scaled)),
    as.character(ncol(mofa_ssn_input_scaled)), "rows = undirected LR edges; columns = TCGA primary tumours",
    "scale(matrix): per-tumour centering and unit-SD scaling across all retained LR edges, matching the original MOFA recipe",
    "none; all retained SSN edges enter MOFA regardless of confidence interval, P value, FDR, or paired result",
    "loading annotation only; not used for feature selection or scaling",
    "loading annotation only; not used for feature selection or scaling",
    "data/tcga_overall_survival.qs2", "all 618 tumours retained for MOFA; 584 samples with complete OS/OS.time retained for later prognosis analysis"
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(mofa_ssn_input_scaled, "results/08_Prepare_MOFA_Input/mofa_ssn_input_scaled.qs2")
qs2::qs_save(mofa_tumor_metadata, "results/08_Prepare_MOFA_Input/mofa_tumor_metadata.qs2")
qs2::qs_save(mofa_edge_annotation, "results/08_Prepare_MOFA_Input/mofa_edge_annotation.qs2")
qs2::qs_save(mofa_input_qc, "results/08_Prepare_MOFA_Input/mofa_input_qc.qs2")
qs2::qs_save(mofa_input_provenance, "results/08_Prepare_MOFA_Input/mofa_input_provenance.qs2")

message("Prepared MOFA input with ", nrow(mofa_ssn_input_scaled), " SSN edges and ",
        ncol(mofa_ssn_input_scaled), " TCGA primary tumours.")
