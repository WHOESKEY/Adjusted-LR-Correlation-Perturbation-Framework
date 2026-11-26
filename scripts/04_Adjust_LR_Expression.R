# Purpose: generate composition-adjusted expression for the undirected LR genes used by SSN.
# Inputs: results/01_Bulk_Processed/*_expression_bc_log2tpm.qs2,
#         results/02_Define_LR_Edges/lr_edges_ssn.qs2, and
#         results/03_Estimate_Composition_Purity/composition_covariates.qs2.
# Outputs: results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2,
#          results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2,
#          results/04_Adjust_LR_Expression/tcga_normal_lr_expression_adjusted.qs2,
#          results/04_Adjust_LR_Expression/lr_expression_adjustment_summary.qs2, and
#          results/04_Adjust_LR_Expression/lr_expression_adjustment_provenance.qs2.
#
# Run from the repository root with Rscript scripts/04_Adjust_LR_Expression.R
# Each gene is adjusted by QR residualization. The original within-group gene mean is then added
# back to the residual, preserving a log2(TPM + 1)-like scale for subsequent correlation analysis.

library(qs2)

dir.create("results/04_Adjust_LR_Expression", recursive = TRUE, showWarnings = FALSE)

gtex_bc <- qs2::qs_read("results/01_Bulk_Processed/gtex_colon_expression_bc_log2tpm.qs2")
tcga_tumor_bc <- qs2::qs_read("results/01_Bulk_Processed/tcga_tumor_expression_bc_log2tpm.qs2")
tcga_normal_bc <- qs2::qs_read("results/01_Bulk_Processed/tcga_normal_expression_bc_log2tpm.qs2")
lr_edges_ssn <- qs2::qs_read("results/02_Define_LR_Edges/lr_edges_ssn.qs2")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--archived-covariates"))
covariate_file <- if ("--archived-covariates" %in% args) {
  "reference/composition_covariates.qs2"
} else {
  "results/03_Estimate_Composition_Purity/composition_covariates.qs2"
}
composition_covariates <- qs2::qs_read(covariate_file)

lr_genes <- sort(unique(c(lr_edges_ssn$Gene1, lr_edges_ssn$Gene2)))
stopifnot(
  all(lr_genes %in% rownames(gtex_bc)),
  all(lr_genes %in% rownames(tcga_tumor_bc)),
  all(lr_genes %in% rownames(tcga_normal_bc))
)

gtex_lr_expression <- gtex_bc[lr_genes, , drop = FALSE]
tcga_tumor_lr_expression <- tcga_tumor_bc[lr_genes, , drop = FALSE]
tcga_normal_lr_expression <- tcga_normal_bc[lr_genes, , drop = FALSE]

gtex_covariates <- composition_covariates[
  match(colnames(gtex_lr_expression), composition_covariates$sample_id),
  c("sample_id", "composition_pc1", "composition_pc2", "composition_pc3", "epic_converged"),
  drop = FALSE
]
tcga_tumor_covariates <- composition_covariates[
  match(colnames(tcga_tumor_lr_expression), composition_covariates$sample_id),
  c("sample_id", "composition_pc1", "composition_pc2", "composition_pc3", "tumor_purity", "epic_converged"),
  drop = FALSE
]
tcga_normal_covariates <- composition_covariates[
  match(colnames(tcga_normal_lr_expression), composition_covariates$sample_id),
  c("sample_id", "composition_pc1", "composition_pc2", "composition_pc3", "epic_converged"),
  drop = FALSE
]
stopifnot(
  identical(gtex_covariates$sample_id, colnames(gtex_lr_expression)),
  identical(tcga_tumor_covariates$sample_id, colnames(tcga_tumor_lr_expression)),
  identical(tcga_normal_covariates$sample_id, colnames(tcga_normal_lr_expression)),
  all(complete.cases(gtex_covariates[, c("composition_pc1", "composition_pc2", "composition_pc3")])),
  all(complete.cases(tcga_tumor_covariates[, c("composition_pc1", "composition_pc2", "composition_pc3", "tumor_purity")])),
  all(complete.cases(tcga_normal_covariates[, c("composition_pc1", "composition_pc2", "composition_pc3")]))
)

gtex_design <- model.matrix(~ composition_pc1 + composition_pc2 + composition_pc3, data = gtex_covariates)
gtex_residual <- t(qr.resid(qr(gtex_design), t(gtex_lr_expression)))
gtex_expression_adjusted <- sweep(gtex_residual, 1, rowMeans(gtex_lr_expression), "+")
dimnames(gtex_expression_adjusted) <- dimnames(gtex_lr_expression)

tcga_tumor_design <- model.matrix(
  ~ composition_pc1 + composition_pc2 + composition_pc3 + tumor_purity,
  data = tcga_tumor_covariates
)
tcga_tumor_residual <- t(qr.resid(qr(tcga_tumor_design), t(tcga_tumor_lr_expression)))
tcga_tumor_expression_adjusted <- sweep(tcga_tumor_residual, 1, rowMeans(tcga_tumor_lr_expression), "+")
dimnames(tcga_tumor_expression_adjusted) <- dimnames(tcga_tumor_lr_expression)

tcga_normal_design <- model.matrix(~ composition_pc1 + composition_pc2 + composition_pc3, data = tcga_normal_covariates)
tcga_normal_residual <- t(qr.resid(qr(tcga_normal_design), t(tcga_normal_lr_expression)))
tcga_normal_expression_adjusted <- sweep(tcga_normal_residual, 1, rowMeans(tcga_normal_lr_expression), "+")
dimnames(tcga_normal_expression_adjusted) <- dimnames(tcga_normal_lr_expression)
stopifnot(
  all(is.finite(gtex_expression_adjusted)),
  all(is.finite(tcga_tumor_expression_adjusted)),
  all(is.finite(tcga_normal_expression_adjusted))
)

adjustment_summary <- rbind(
  data.frame(
    matrix = "GTEx_normal", n_sample = ncol(gtex_lr_expression), n_lr_gene = nrow(gtex_lr_expression),
    covariates = "composition_pc1 + composition_pc2 + composition_pc3",
    n_epic_nonconverged = sum(!gtex_covariates$epic_converged),
    mean_gene_variance_before = mean(apply(gtex_lr_expression, 1, var)),
    mean_gene_variance_after = mean(apply(gtex_expression_adjusted, 1, var))
  ),
  data.frame(
    matrix = "TCGA_tumor", n_sample = ncol(tcga_tumor_lr_expression), n_lr_gene = nrow(tcga_tumor_lr_expression),
    covariates = "composition_pc1 + composition_pc2 + composition_pc3 + tumor_purity",
    n_epic_nonconverged = sum(!tcga_tumor_covariates$epic_converged),
    mean_gene_variance_before = mean(apply(tcga_tumor_lr_expression, 1, var)),
    mean_gene_variance_after = mean(apply(tcga_tumor_expression_adjusted, 1, var))
  ),
  data.frame(
    matrix = "TCGA_matched_normal", n_sample = ncol(tcga_normal_lr_expression), n_lr_gene = nrow(tcga_normal_lr_expression),
    covariates = "composition_pc1 + composition_pc2 + composition_pc3",
    n_epic_nonconverged = sum(!tcga_normal_covariates$epic_converged),
    mean_gene_variance_before = mean(apply(tcga_normal_lr_expression, 1, var)),
    mean_gene_variance_after = mean(apply(tcga_normal_expression_adjusted, 1, var))
  )
)

adjustment_provenance <- data.frame(
  item = c(
    "ssn_edge_input", "n_ssn_edges", "n_lr_genes", "expression_input", "adjustment_method",
    "normal_covariates", "tumor_covariates", "sample_exclusion", "epic_nonconvergence_handling"
  ),
  value = c(
    "results/02_Define_LR_Edges/lr_edges_ssn.qs2", as.character(nrow(lr_edges_ssn)),
    as.character(length(lr_genes)), "ComBat-corrected log2(TPM + 1)",
    "gene-wise QR residualization with original within-group gene mean added back",
    "composition_pc1 + composition_pc2 + composition_pc3",
    "composition_pc1 + composition_pc2 + composition_pc3 + tumor_purity",
    "none", "all samples retained; epic_converged remains available in prior covariate object"
  ),
  stringsAsFactors = FALSE
)

adjustment_provenance <- rbind(adjustment_provenance,
                               data.frame(item = "composition_covariates", value = covariate_file))
qs2::qs_save(gtex_expression_adjusted, "results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2")
qs2::qs_save(tcga_tumor_expression_adjusted, "results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2")
qs2::qs_save(tcga_normal_expression_adjusted, "results/04_Adjust_LR_Expression/tcga_normal_lr_expression_adjusted.qs2")
qs2::qs_save(adjustment_summary, "results/04_Adjust_LR_Expression/lr_expression_adjustment_summary.qs2")
qs2::qs_save(adjustment_provenance, "results/04_Adjust_LR_Expression/lr_expression_adjustment_provenance.qs2")

message("Adjusted ", length(lr_genes), " LR genes for ", ncol(gtex_expression_adjusted), " GTEx normals, ",
        ncol(tcga_tumor_expression_adjusted), " TCGA tumours, and ",
        ncol(tcga_normal_expression_adjusted), " TCGA matched normals.")
