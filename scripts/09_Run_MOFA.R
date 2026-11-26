# Purpose: fit a single unsupervised MOFA model to the complete composition-adjusted SSN perturbation
# matrix, then annotate every active factor with OS associations and edge stability data.
# Inputs: results/08_Prepare_MOFA_Input/mofa_ssn_input_scaled.qs2,
#         results/08_Prepare_MOFA_Input/mofa_tumor_metadata.qs2, and
#         results/08_Prepare_MOFA_Input/mofa_edge_annotation.qs2.
# Outputs: results/09_Run_MOFA/ssn_mofa_model.hdf5,
#          results/09_Run_MOFA/mofa_factor_scores.qs2,
#          results/09_Run_MOFA/mofa_factor_variance.qs2,
#          results/09_Run_MOFA/mofa_edge_loadings.qs2,
#          results/09_Run_MOFA/mofa_factor_os_association.qs2,
#          results/09_Run_MOFA/mofa_primary_factor.qs2, and
#          results/09_Run_MOFA/mofa_run_provenance.qs2.
#
# Run from the repository root with Rscript scripts/09_Run_MOFA.R
# The fit uses no survival data. Variance rank and subsequent OS prioritisation are separate.

library(qs2)
library(MOFA2)
library(reticulate)
library(survival)

dir.create("results/09_Run_MOFA", recursive = TRUE, showWarnings = FALSE)

set.seed(20260801)
mofa_seed <- 20260801L
mofa_max_factors <- 10L

mofa_ssn_input_scaled <- qs2::qs_read("results/08_Prepare_MOFA_Input/mofa_ssn_input_scaled.qs2")
mofa_tumor_metadata <- qs2::qs_read("results/08_Prepare_MOFA_Input/mofa_tumor_metadata.qs2")
mofa_edge_annotation <- qs2::qs_read("results/08_Prepare_MOFA_Input/mofa_edge_annotation.qs2")
stopifnot(
  identical(colnames(mofa_ssn_input_scaled), mofa_tumor_metadata$sample_id),
  identical(rownames(mofa_ssn_input_scaled), rownames(mofa_edge_annotation)),
  all(is.finite(mofa_ssn_input_scaled))
)

# Set RETICULATE_PYTHON to an interpreter containing mofapy2 before starting R.
if (!reticulate::py_module_available("mofapy2")) {
  stop("Install mofapy2 and set RETICULATE_PYTHON to that environment's Python executable.")
}
mofa_python_version <- as.character(reticulate::import("importlib.metadata")$version("mofapy2"))

mofa_views <- list(SSN_LR = mofa_ssn_input_scaled)
mofa_object <- MOFA2::create_mofa(mofa_views)
mofa_model_options <- MOFA2::get_default_model_options(mofa_object)
mofa_model_options$num_factors <- mofa_max_factors
mofa_training_options <- MOFA2::get_default_training_options(mofa_object)
mofa_training_options$seed <- mofa_seed
mofa_training_options$maxiter <- 1000L
mofa_training_options$verbose <- FALSE
mofa_object <- MOFA2::prepare_mofa(
  mofa_object,
  data_options = MOFA2::get_default_data_options(mofa_object),
  model_options = mofa_model_options,
  training_options = mofa_training_options
)

message("Fitting MOFA with ", nrow(mofa_ssn_input_scaled), " SSN edges and ",
        ncol(mofa_ssn_input_scaled), " tumours.")
mofa_model <- MOFA2::run_mofa(
  mofa_object,
  outfile = "results/09_Run_MOFA/ssn_mofa_model.hdf5",
  save_data = TRUE,
  use_basilisk = FALSE
)

mofa_factor_scores <- MOFA2::get_factors(mofa_model, as.data.frame = FALSE)
if (is.list(mofa_factor_scores)) mofa_factor_scores <- mofa_factor_scores[[1]]
mofa_factor_scores <- as.matrix(mofa_factor_scores)
mofa_edge_weights <- MOFA2::get_weights(mofa_model, views = "SSN_LR", as.data.frame = FALSE)
if (is.list(mofa_edge_weights)) mofa_edge_weights <- mofa_edge_weights[[1]]
mofa_edge_weights <- as.matrix(mofa_edge_weights)
stopifnot(
  identical(rownames(mofa_factor_scores), colnames(mofa_ssn_input_scaled)),
  identical(rownames(mofa_edge_weights), rownames(mofa_ssn_input_scaled))
)

for (factor_index in seq_len(ncol(mofa_factor_scores))) {
  orientation_edge <- which.max(abs(mofa_edge_weights[, factor_index]))
  if (mofa_edge_weights[orientation_edge, factor_index] < 0) {
    mofa_factor_scores[, factor_index] <- -mofa_factor_scores[, factor_index]
    mofa_edge_weights[, factor_index] <- -mofa_edge_weights[, factor_index]
  }
}

mofa_variance_raw <- MOFA2::get_variance_explained(mofa_model, views = "SSN_LR", as.data.frame = TRUE)$r2_per_factor
mofa_factor_variance <- stats::aggregate(value ~ factor, data = mofa_variance_raw, FUN = sum)
colnames(mofa_factor_variance) <- c("factor", "explained_variance")
mofa_factor_variance <- mofa_factor_variance[order(-mofa_factor_variance$explained_variance), , drop = FALSE]
mofa_factor_variance$variance_rank <- seq_len(nrow(mofa_factor_variance))
stopifnot(all(mofa_factor_variance$factor %in% colnames(mofa_factor_scores)))

mofa_factor_scores <- data.frame(sample_id = rownames(mofa_factor_scores), mofa_factor_scores, check.names = FALSE)
mofa_edge_loadings <- cbind(mofa_edge_annotation, as.data.frame(mofa_edge_weights, check.names = FALSE))
rownames(mofa_edge_loadings) <- mofa_edge_loadings$ssn_edge_id

os_metadata <- mofa_tumor_metadata[mofa_tumor_metadata$has_os_annotation, c("sample_id", "OS", "OS.time"), drop = FALSE]
os_factor_scores <- mofa_factor_scores[match(os_metadata$sample_id, mofa_factor_scores$sample_id), , drop = FALSE]
stopifnot(identical(os_factor_scores$sample_id, os_metadata$sample_id), nrow(os_metadata) == 584L)

mofa_factor_os_association <- data.frame()
for (factor_name in colnames(mofa_factor_scores)[-1]) {
  os_data <- data.frame(OS = os_metadata$OS, OS.time = os_metadata$OS.time, score = os_factor_scores[[factor_name]])
  cox_result <- summary(survival::coxph(survival::Surv(OS.time, OS) ~ score, data = os_data))
  median_cutoff <- stats::median(os_data$score)
  os_data$group <- ifelse(os_data$score > median_cutoff, "High", "Low")
  logrank_result <- survival::survdiff(survival::Surv(OS.time, OS) ~ group, data = os_data)
  mofa_factor_os_association <- rbind(
    mofa_factor_os_association,
    data.frame(
      factor = factor_name,
      n_os_annotated = nrow(os_data),
      explained_variance = mofa_factor_variance$explained_variance[match(factor_name, mofa_factor_variance$factor)],
      cox_hazard_ratio = cox_result$conf.int[1, "exp(coef)"],
      cox_ci_lower = cox_result$conf.int[1, "lower .95"],
      cox_ci_upper = cox_result$conf.int[1, "upper .95"],
      cox_p = cox_result$coefficients[1, "Pr(>|z|)"],
      km_cutoff = median_cutoff,
      km_logrank_p = stats::pchisq(logrank_result$chisq, length(logrank_result$n) - 1, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  )
}
mofa_factor_os_association$cox_fdr_bh <- stats::p.adjust(mofa_factor_os_association$cox_p, method = "BH")
mofa_factor_os_association$km_logrank_fdr_bh <- stats::p.adjust(mofa_factor_os_association$km_logrank_p, method = "BH")
mofa_factor_os_association <- mofa_factor_os_association[
  match(mofa_factor_variance$factor, mofa_factor_os_association$factor), , drop = FALSE
]

primary_factor_name <- as.character(mofa_factor_variance$factor[1])
os_order <- order(mofa_factor_os_association$km_logrank_fdr_bh,
                  mofa_factor_os_association$km_logrank_p)
mofa_primary_factor <- data.frame(
  criterion = c("Largest explained variance", "Smallest BH-adjusted median-split log-rank P"),
  factor = c(primary_factor_name, mofa_factor_os_association$factor[os_order[1]])
)

mofa_run_provenance <- data.frame(
  item = c(
    "random_seed", "mofapy2_version", "mofapy2_available", "view_name", "n_input_edges",
    "n_input_tumors", "maximum_requested_factors", "active_factors", "input_scaling",
    "edge_filtering", "factor_orientation", "primary_factor_selection", "os_association",
    "km_grouping", "multiple_testing"
  ),
  value = c(
    as.character(mofa_seed), mofa_python_version, "TRUE", "SSN_LR", as.character(nrow(mofa_ssn_input_scaled)),
    as.character(ncol(mofa_ssn_input_scaled)), as.character(mofa_max_factors),
    as.character(ncol(mofa_factor_scores) - 1L), "precomputed per-tumour scaling from script 08",
    "none; all 238 SSN edges entered MOFA", "largest absolute loading positive for each factor",
    "maximum SSN-view explained variance", "univariate Cox model across all active factors after unsupervised fit",
    "median factor score", "BH correction separately across active-factor Cox P values and KM log-rank P values"
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(mofa_factor_scores, "results/09_Run_MOFA/mofa_factor_scores.qs2")
qs2::qs_save(mofa_factor_variance, "results/09_Run_MOFA/mofa_factor_variance.qs2")
qs2::qs_save(mofa_edge_loadings, "results/09_Run_MOFA/mofa_edge_loadings.qs2")
qs2::qs_save(mofa_factor_os_association, "results/09_Run_MOFA/mofa_factor_os_association.qs2")
qs2::qs_save(mofa_primary_factor, "results/09_Run_MOFA/mofa_primary_factor.qs2")
qs2::qs_save(mofa_model@status, "results/09_Run_MOFA/mofa_model_status.qs2")
qs2::qs_save(mofa_model@training_stats, "results/09_Run_MOFA/mofa_training_stats.qs2")
qs2::qs_save(mofa_run_provenance, "results/09_Run_MOFA/mofa_run_provenance.qs2")

message("Saved MOFA results with ", ncol(mofa_factor_scores) - 1L, " active factors; largest-variance factor is ",
        primary_factor_name, ".")
