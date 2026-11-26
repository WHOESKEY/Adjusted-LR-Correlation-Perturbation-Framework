# Purpose: calculate Spearman correlation perturbation for each undirected LR edge and TCGA tumour.
# Inputs: results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2,
#         results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2, and
#         results/02_Define_LR_Edges/lr_edges_ssn.qs2.
# Outputs: results/05_Calculate_SSN_Perturbation/ssn_reference_spearman.qs2,
#          results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2,
#          results/05_Calculate_SSN_Perturbation/ssn_dropped_edges.qs2, and
#          results/05_Calculate_SSN_Perturbation/ssn_perturbation_provenance.qs2.
#
# Run from the repository root with Rscript scripts/05_Calculate_SSN_Perturbation.R
# A tumour-specific perturbation is Delta rho = rho(GTEx reference + tumour) - rho(GTEx reference),
# where rho is Spearman rank correlation. The Pearson-SSN analytic Z statistic is deliberately not
# calculated because it is not valid for Spearman correlation.

library(qs2)

dir.create("results/05_Calculate_SSN_Perturbation", recursive = TRUE, showWarnings = FALSE)

gtex_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2")
tcga_tumor_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2")
lr_edges_ssn <- qs2::qs_read("results/02_Define_LR_Edges/lr_edges_ssn.qs2")

stopifnot(
  identical(rownames(gtex_expression_adjusted), rownames(tcga_tumor_expression_adjusted)),
  all(lr_edges_ssn$Gene1 %in% rownames(gtex_expression_adjusted)),
  all(lr_edges_ssn$Gene2 %in% rownames(gtex_expression_adjusted)),
  all(is.finite(gtex_expression_adjusted)),
  all(is.finite(tcga_tumor_expression_adjusted))
)

reference_rho <- vapply(seq_len(nrow(lr_edges_ssn)), function(i) {
  suppressWarnings(stats::cor(
    gtex_expression_adjusted[lr_edges_ssn$Gene1[i], ],
    gtex_expression_adjusted[lr_edges_ssn$Gene2[i], ],
    method = "spearman"
  ))
}, numeric(1))

valid_edges <- is.finite(reference_rho)
ssn_dropped_edges <- lr_edges_ssn[!valid_edges, , drop = FALSE]
ssn_dropped_edges$drop_reason <- rep("nonfinite_reference_spearman_correlation", nrow(ssn_dropped_edges))
lr_edges_ssn <- lr_edges_ssn[valid_edges, , drop = FALSE]
reference_rho <- reference_rho[valid_edges]
stopifnot(length(reference_rho) > 0)

reference_spearman <- cbind(lr_edges_ssn, reference_spearman_rho = reference_rho)
tumor_ids <- colnames(tcga_tumor_expression_adjusted)
message("Computing Spearman perturbation for ", length(tumor_ids), " tumours and ", nrow(lr_edges_ssn), " LR edges.")

ssn_tumor_spearman_perturbation <- vapply(tumor_ids, function(tumor_id) {
  vapply(seq_len(nrow(lr_edges_ssn)), function(i) {
    case_rho <- suppressWarnings(stats::cor(
      c(gtex_expression_adjusted[lr_edges_ssn$Gene1[i], ], tcga_tumor_expression_adjusted[lr_edges_ssn$Gene1[i], tumor_id]),
      c(gtex_expression_adjusted[lr_edges_ssn$Gene2[i], ], tcga_tumor_expression_adjusted[lr_edges_ssn$Gene2[i], tumor_id]),
      method = "spearman"
    ))
    case_rho - reference_rho[i]
  }, numeric(1))
}, numeric(nrow(lr_edges_ssn)))
rownames(ssn_tumor_spearman_perturbation) <- lr_edges_ssn$ssn_edge_id
colnames(ssn_tumor_spearman_perturbation) <- tumor_ids
stopifnot(all(is.finite(ssn_tumor_spearman_perturbation)))

perturbation_provenance <- data.frame(
  item = c(
    "correlation_method", "perturbation_definition", "reference_cohort", "reference_sample_count",
    "case_cohort", "case_sample_count", "n_input_ssn_edges", "n_retained_ssn_edges",
    "n_dropped_nonfinite_reference_edges", "analytic_ssn_z_statistic"
  ),
  value = c(
    "Spearman rank correlation", "rho(reference + one tumour) - rho(reference)", "GTEx colon normal",
    as.character(ncol(gtex_expression_adjusted)), "TCGA primary tumour", as.character(length(tumor_ids)),
    as.character(length(valid_edges)), as.character(nrow(lr_edges_ssn)), as.character(nrow(ssn_dropped_edges)),
    "not calculated; Pearson-specific formula is not valid for Spearman"
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(reference_spearman, "results/05_Calculate_SSN_Perturbation/ssn_reference_spearman.qs2")
qs2::qs_save(ssn_tumor_spearman_perturbation, "results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2")
qs2::qs_save(ssn_dropped_edges, "results/05_Calculate_SSN_Perturbation/ssn_dropped_edges.qs2")
qs2::qs_save(perturbation_provenance, "results/05_Calculate_SSN_Perturbation/ssn_perturbation_provenance.qs2")

message("Saved a ", nrow(ssn_tumor_spearman_perturbation), " x ", ncol(ssn_tumor_spearman_perturbation), " Spearman perturbation matrix.")
