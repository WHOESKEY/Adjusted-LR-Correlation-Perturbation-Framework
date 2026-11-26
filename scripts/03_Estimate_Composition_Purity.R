# Purpose: estimate bulk cell composition and TCGA tumour purity for SSN covariate control.
# Inputs: results/00_Input/gtex_colon_expression_tpm.qs2,
#         results/00_Input/tcga_selected_expression_tpm.qs2,
#         data/tcga_selected_sample_map.qs2, and data/tcga_matched_pairs.qs2.
# Outputs: results/03_Estimate_Composition_Purity/composition_covariates.qs2,
#          results/03_Estimate_Composition_Purity/composition_pca.qs2,
#          results/03_Estimate_Composition_Purity/tcga_tumor_purity.qs2,
#          results/03_Estimate_Composition_Purity/estimate_scores_tcga_tumors.qs2,
#          results/03_Estimate_Composition_Purity/epic_fit_quality.qs2, and
#          results/03_Estimate_Composition_Purity/composition_purity_provenance.qs2;
#          figures/03_Estimate_Composition_Purity/composition_pca.pdf.
#
# Run from the repository root with Rscript scripts/03_Estimate_Composition_Purity.R
# EPIC is run on raw TPM expression with identical settings for GTEx normal, TCGA primary tumour,
# and TCGA matched-normal samples. ESTIMATE-derived tumour purity is calculated for TCGA primary
# tumours only; normal samples retain missing purity rather than an artificial value of zero.

library(EPIC)
library(ggplot2)
library(IOBR)
library(qs2)

dir.create("results/03_Estimate_Composition_Purity", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/03_Estimate_Composition_Purity", recursive = TRUE, showWarnings = FALSE)

gtex_tpm <- qs2::qs_read("results/00_Input/gtex_colon_expression_tpm.qs2")
tcga_tpm <- qs2::qs_read("results/00_Input/tcga_selected_expression_tpm.qs2")
tcga_map <- qs2::qs_read("data/tcga_selected_sample_map.qs2")
tcga_pairs <- qs2::qs_read("data/tcga_matched_pairs.qs2")

tumor_ids <- tcga_map$sample_id[tcga_map$sample_type == "01"]
normal_ids <- tcga_pairs$normal_id
stopifnot(length(tumor_ids) == 618, length(normal_ids) == 50)
stopifnot(all(tumor_ids %in% colnames(tcga_tpm)), all(normal_ids %in% colnames(tcga_tpm)))
stopifnot(all(is.finite(gtex_tpm)), all(is.finite(tcga_tpm)))

gtex_epic <- EPIC::EPIC(bulk = gtex_tpm, withOtherCells = TRUE)
tcga_tumor_epic <- EPIC::EPIC(bulk = tcga_tpm[, tumor_ids, drop = FALSE], withOtherCells = TRUE)
tcga_normal_epic <- EPIC::EPIC(bulk = tcga_tpm[, normal_ids, drop = FALSE], withOtherCells = TRUE)

gtex_fractions <- as.data.frame(gtex_epic$cellFractions)
tcga_tumor_fractions <- as.data.frame(tcga_tumor_epic$cellFractions)
tcga_normal_fractions <- as.data.frame(tcga_normal_epic$cellFractions)
stopifnot(
  identical(colnames(gtex_fractions), colnames(tcga_tumor_fractions)),
  identical(colnames(gtex_fractions), colnames(tcga_normal_fractions))
)

epic_columns <- paste0("epic_", make.names(colnames(gtex_fractions)))
colnames(gtex_fractions) <- epic_columns
colnames(tcga_tumor_fractions) <- epic_columns
colnames(tcga_normal_fractions) <- epic_columns

gtex_epic_converged <- setNames(gtex_epic$fit.gof$convergeCode == 0, rownames(gtex_epic$fit.gof))
tcga_tumor_epic_converged <- setNames(tcga_tumor_epic$fit.gof$convergeCode == 0, rownames(tcga_tumor_epic$fit.gof))
tcga_normal_epic_converged <- setNames(tcga_normal_epic$fit.gof$convergeCode == 0, rownames(tcga_normal_epic$fit.gof))
stopifnot(
  identical(names(gtex_epic_converged), rownames(gtex_fractions)),
  identical(names(tcga_tumor_epic_converged), rownames(tcga_tumor_fractions)),
  identical(names(tcga_normal_epic_converged), rownames(tcga_normal_fractions))
)

gtex_covariates <- data.frame(
  sample_id = rownames(gtex_fractions), cohort = "GTEx", tissue = "Normal", sample_type = "GTEx_normal",
  tumor_purity = NA_real_, epic_converged = gtex_epic_converged, gtex_fractions, check.names = FALSE
)
tcga_tumor_covariates <- data.frame(
  sample_id = rownames(tcga_tumor_fractions), cohort = "TCGA", tissue = "Tumor", sample_type = "01",
  tumor_purity = NA_real_, epic_converged = tcga_tumor_epic_converged, tcga_tumor_fractions, check.names = FALSE
)
tcga_normal_covariates <- data.frame(
  sample_id = rownames(tcga_normal_fractions), cohort = "TCGA", tissue = "Normal", sample_type = "11",
  tumor_purity = NA_real_, epic_converged = tcga_normal_epic_converged, tcga_normal_fractions, check.names = FALSE
)
composition_covariates <- rbind(gtex_covariates, tcga_tumor_covariates, tcga_normal_covariates)
stopifnot(nrow(composition_covariates) == 1566, !anyDuplicated(composition_covariates$sample_id))

# deconvo_estimate() requires the IOBR namespace to be attached and uses temporary GCT files.
# The scratch directory prevents those files from being written into this analysis workspace.
estimate_scratch <- tempfile("estimate_")
dir.create(estimate_scratch)
original_working_directory <- getwd()
setwd(estimate_scratch)
estimate_scores <- tryCatch(IOBR::deconvo_estimate(
  eset = as.data.frame(log2(tcga_tpm[, tumor_ids, drop = FALSE] + 1)),
  project = NULL,
  platform = "affymetrix"
), finally = setwd(original_working_directory))
unlink(estimate_scratch, recursive = TRUE, force = TRUE)

tumor_purity <- data.frame(sample_id = estimate_scores$ID, tumor_purity = estimate_scores$TumorPurity_estimate)
stopifnot(setequal(tumor_purity$sample_id, tumor_ids))
tumor_purity <- tumor_purity[match(tumor_ids, tumor_purity$sample_id), , drop = FALSE]
if (any(!is.finite(tumor_purity$tumor_purity) | tumor_purity$tumor_purity <= 0 | tumor_purity$tumor_purity >= 1)) {
  stop("ESTIMATE-derived tumour purity must be finite and strictly between 0 and 1.")
}
composition_covariates$tumor_purity[match(tumor_purity$sample_id, composition_covariates$sample_id)] <- tumor_purity$tumor_purity
stopifnot(sum(is.finite(composition_covariates$tumor_purity)) == length(tumor_ids))

epic_matrix <- as.matrix(composition_covariates[, epic_columns, drop = FALSE])
stopifnot(all(is.finite(epic_matrix)))
keep_epic <- apply(epic_matrix, 2, function(x) is.finite(sd(x)) && sd(x) > 0)
stopifnot(sum(keep_epic) >= 3)
composition_pca_fit <- prcomp(epic_matrix[, keep_epic, drop = FALSE], center = TRUE, scale. = TRUE)
selected_pcs <- 3
composition_scores <- composition_pca_fit$x[, seq_len(selected_pcs), drop = FALSE]
colnames(composition_scores) <- paste0("composition_pc", seq_len(selected_pcs))
composition_covariates <- cbind(composition_covariates, composition_scores)

composition_pca <- list(
  rotation = composition_pca_fit$rotation,
  sdev = composition_pca_fit$sdev,
  variance_explained = composition_pca_fit$sdev^2 / sum(composition_pca_fit$sdev^2),
  selected_components = colnames(composition_scores),
  scores = data.frame(sample_id = composition_covariates$sample_id, composition_scores, check.names = FALSE)
)

gtex_fit_quality <- as.data.frame(gtex_epic$fit.gof)
gtex_fit_quality$sample_id <- rownames(gtex_epic$fit.gof)
gtex_fit_quality$cohort <- "GTEx"
gtex_fit_quality$tissue <- "Normal"
tcga_tumor_fit_quality <- as.data.frame(tcga_tumor_epic$fit.gof)
tcga_tumor_fit_quality$sample_id <- rownames(tcga_tumor_epic$fit.gof)
tcga_tumor_fit_quality$cohort <- "TCGA"
tcga_tumor_fit_quality$tissue <- "Tumor"
tcga_normal_fit_quality <- as.data.frame(tcga_normal_epic$fit.gof)
tcga_normal_fit_quality$sample_id <- rownames(tcga_normal_epic$fit.gof)
tcga_normal_fit_quality$cohort <- "TCGA"
tcga_normal_fit_quality$tissue <- "Normal"
epic_fit_quality <- rbind(gtex_fit_quality, tcga_tumor_fit_quality, tcga_normal_fit_quality)

composition_pca_plot <- ggplot(
  composition_covariates,
  aes(composition_pc1, composition_pc2, color = cohort, shape = tissue)
) +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_color_manual(values = c(GTEx = "#66C2A5", TCGA = "#B2182B")) +
  theme_classic(base_size = 12) +
  labs(
    title = "EPIC cell-composition principal components",
    x = paste0("PC1 (", sprintf("%.1f", 100 * composition_pca$variance_explained[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", 100 * composition_pca$variance_explained[2]), "%)")
  )
ggsave("figures/03_Estimate_Composition_Purity/composition_pca.pdf", composition_pca_plot, width = 7, height = 5)

composition_purity_provenance <- data.frame(
  item = c(
    "epic_package_version", "iobr_package_version", "epic_input_scale", "epic_withOtherCells",
    "tumor_purity_method", "tumor_purity_input_scale", "estimate_platform", "n_gtex_normal",
    "n_tcga_tumor", "n_tcga_matched_normal", "n_epic_components", "composition_pc_count",
    "n_epic_converged", "n_epic_nonconverged"
  ),
  value = c(
    as.character(getNamespaceVersion("EPIC")), as.character(getNamespaceVersion("IOBR")), "TPM",
    "TRUE", "IOBR::deconvo_estimate (ESTIMATE score transform)", "log2(TPM + 1)", "affymetrix",
    as.character(ncol(gtex_tpm)), as.character(length(tumor_ids)), as.character(length(normal_ids)),
    as.character(sum(keep_epic)), as.character(selected_pcs), as.character(sum(composition_covariates$epic_converged)),
    as.character(sum(!composition_covariates$epic_converged))
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(composition_covariates, "results/03_Estimate_Composition_Purity/composition_covariates.qs2")
qs2::qs_save(composition_pca, "results/03_Estimate_Composition_Purity/composition_pca.qs2")
qs2::qs_save(tumor_purity, "results/03_Estimate_Composition_Purity/tcga_tumor_purity.qs2")
qs2::qs_save(estimate_scores, "results/03_Estimate_Composition_Purity/estimate_scores_tcga_tumors.qs2")
qs2::qs_save(epic_fit_quality, "results/03_Estimate_Composition_Purity/epic_fit_quality.qs2")
qs2::qs_save(composition_purity_provenance, "results/03_Estimate_Composition_Purity/composition_purity_provenance.qs2")

message("EPIC composition covariates saved for ", nrow(composition_covariates), " samples.")
message("EPIC convergence: ", sum(composition_covariates$epic_converged), " converged; ",
        sum(!composition_covariates$epic_converged), " nonconverged (flagged, not removed).")
message("ESTIMATE tumour purity saved for ", nrow(tumor_purity), " TCGA primary tumours.")
