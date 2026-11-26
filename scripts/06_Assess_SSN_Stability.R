# Purpose: quantify cohort-level uncertainty, empirical significance, and normal-reference stability
# for the composition-adjusted Spearman SSN perturbation scores.
# Inputs: results/02_Define_LR_Edges/lr_edges_ssn.qs2,
#         results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2,
#         results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2, and
#         results/05_Calculate_SSN_Perturbation/ssn_*_spearman*.qs2.
# Outputs: results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2,
#          results/06_Assess_SSN_Stability/ssn_tumor_bootstrap_mean.qs2,
#          results/06_Assess_SSN_Stability/ssn_reference_bootstrap.qs2, and
#          results/06_Assess_SSN_Stability/ssn_reference_size_stability.qs2, and
#          results/06_Assess_SSN_Stability/ssn_stability_provenance.qs2;
#          figures/06_Assess_SSN_Stability/reference_size_stability.{png,pdf}.
#
# Run from the repository root with Rscript scripts/06_Assess_SSN_Stability.R
# This script deliberately uses nonparametric resampling because the analytic Pearson-SSN Z statistic
# is not valid after replacing Pearson correlation with Spearman rank correlation.

library(qs2)
library(Rcpp)
library(ggplot2)

source("scripts/00_Functions_SSN.R")

dir.create("results/06_Assess_SSN_Stability", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/06_Assess_SSN_Stability", recursive = TRUE, showWarnings = FALSE)

set.seed(20260730)
n_tumor_bootstrap <- 1000L
n_sign_permutation <- 1000L
n_reference_bootstrap <- 1000L
reference_sizes <- c(50L, 100L, 250L, 500L, 898L)
n_reference_size_resamples <- 1000L

lr_edges_ssn <- qs2::qs_read("results/02_Define_LR_Edges/lr_edges_ssn.qs2")
gtex_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/gtex_lr_expression_adjusted.qs2")
tcga_tumor_expression_adjusted <- qs2::qs_read("results/04_Adjust_LR_Expression/tcga_tumor_lr_expression_adjusted.qs2")
reference_spearman <- qs2::qs_read("results/05_Calculate_SSN_Perturbation/ssn_reference_spearman.qs2")
ssn_tumor_perturbation <- qs2::qs_read("results/05_Calculate_SSN_Perturbation/ssn_tumor_spearman_perturbation.qs2")

lr_edges_ssn <- lr_edges_ssn[match(rownames(ssn_tumor_perturbation), lr_edges_ssn$ssn_edge_id), , drop = FALSE]
reference_spearman <- reference_spearman[match(rownames(ssn_tumor_perturbation), reference_spearman$ssn_edge_id), , drop = FALSE]
stopifnot(
  identical(lr_edges_ssn$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  identical(reference_spearman$ssn_edge_id, rownames(ssn_tumor_perturbation)),
  identical(rownames(gtex_expression_adjusted), rownames(tcga_tumor_expression_adjusted)),
  all(lr_edges_ssn$Gene1 %in% rownames(gtex_expression_adjusted)),
  all(lr_edges_ssn$Gene2 %in% rownames(gtex_expression_adjusted)),
  all(is.finite(ssn_tumor_perturbation))
)

observed_mean_delta_rho <- rowMeans(ssn_tumor_perturbation)
tumor_bootstrap_mean <- vapply(seq_len(n_tumor_bootstrap), function(i) {
  rowMeans(ssn_tumor_perturbation[, sample.int(ncol(ssn_tumor_perturbation), replace = TRUE), drop = FALSE])
}, numeric(nrow(ssn_tumor_perturbation)))
rownames(tumor_bootstrap_mean) <- rownames(ssn_tumor_perturbation)

tumor_bootstrap_ci <- t(apply(tumor_bootstrap_mean, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
tumor_bootstrap_direction_consistency <- rowMeans(
  sweep(tumor_bootstrap_mean, 1, sign(observed_mean_delta_rho), "*") > 0
)

permutation_signs <- matrix(
  sample(c(-1, 1), ncol(ssn_tumor_perturbation) * n_sign_permutation, replace = TRUE),
  nrow = ncol(ssn_tumor_perturbation)
)
sign_permutation_mean <- t(permutation_signs) %*% t(ssn_tumor_perturbation) / ncol(ssn_tumor_perturbation)
empirical_p_value <- (colSums(abs(sign_permutation_mean) >= abs(observed_mean_delta_rho)) + 1) /
  (n_sign_permutation + 1)

# The reference-size analysis needs only the mean perturbation across tumours for each resample.
# This compiled implementation returns the exact Spearman result while avoiding storage of a large
# edge-by-tumour matrix for every resampled reference cohort.
Rcpp::cppFunction(code = '
Rcpp::List calculate_reference_summary_cpp(
    Rcpp::NumericMatrix reference_expression,
    Rcpp::NumericMatrix case_expression,
    Rcpp::IntegerVector gene1_index,
    Rcpp::IntegerVector gene2_index) {
  const int n_edge = gene1_index.size();
  const int n_reference = reference_expression.ncol();
  const int n_case = case_expression.ncol();
  Rcpp::NumericVector reference_spearman_rho(n_edge);
  Rcpp::NumericVector mean_delta_rho(n_edge);

  for (int edge_index = 0; edge_index < n_edge; ++edge_index) {
    const int gene1 = gene1_index[edge_index] - 1;
    const int gene2 = gene2_index[edge_index] - 1;
    std::vector<double> reference_gene1(n_reference), reference_gene2(n_reference);
    std::vector<int> order_gene1(n_reference), order_gene2(n_reference);
    std::iota(order_gene1.begin(), order_gene1.end(), 0);
    std::iota(order_gene2.begin(), order_gene2.end(), 0);

    for (int reference_index = 0; reference_index < n_reference; ++reference_index) {
      reference_gene1[reference_index] = reference_expression(gene1, reference_index);
      reference_gene2[reference_index] = reference_expression(gene2, reference_index);
    }
    std::sort(order_gene1.begin(), order_gene1.end(), [&](int left, int right) {
      return reference_gene1[left] < reference_gene1[right];
    });
    std::sort(order_gene2.begin(), order_gene2.end(), [&](int left, int right) {
      return reference_gene2[left] < reference_gene2[right];
    });

    std::vector<double> rank_gene1(n_reference), rank_gene2(n_reference);
    for (int start = 0; start < n_reference;) {
      int end = start + 1;
      while (end < n_reference &&
             reference_gene1[order_gene1[end]] == reference_gene1[order_gene1[start]]) ++end;
      const double rank = (static_cast<double>(start + 1) + static_cast<double>(end)) / 2.0;
      for (int rank_index = start; rank_index < end; ++rank_index) rank_gene1[order_gene1[rank_index]] = rank;
      start = end;
    }
    for (int start = 0; start < n_reference;) {
      int end = start + 1;
      while (end < n_reference &&
             reference_gene2[order_gene2[end]] == reference_gene2[order_gene2[start]]) ++end;
      const double rank = (static_cast<double>(start + 1) + static_cast<double>(end)) / 2.0;
      for (int rank_index = start; rank_index < end; ++rank_index) rank_gene2[order_gene2[rank_index]] = rank;
      start = end;
    }

    double sum_rank_gene1_sq = 0.0, sum_rank_gene2_sq = 0.0, sum_rank_product = 0.0;
    double sum_rank_gene1 = 0.0, sum_rank_gene2 = 0.0;
    for (int reference_index = 0; reference_index < n_reference; ++reference_index) {
      sum_rank_gene1_sq += rank_gene1[reference_index] * rank_gene1[reference_index];
      sum_rank_gene2_sq += rank_gene2[reference_index] * rank_gene2[reference_index];
      sum_rank_product += rank_gene1[reference_index] * rank_gene2[reference_index];
      sum_rank_gene1 += rank_gene1[reference_index];
      sum_rank_gene2 += rank_gene2[reference_index];
    }
    const double reference_mean_gene1 = sum_rank_gene1 / n_reference;
    const double reference_mean_gene2 = sum_rank_gene2 / n_reference;
    const double reference_covariance = sum_rank_product -
      n_reference * reference_mean_gene1 * reference_mean_gene2;
    const double reference_variance_gene1 = sum_rank_gene1_sq -
      n_reference * reference_mean_gene1 * reference_mean_gene1;
    const double reference_variance_gene2 = sum_rank_gene2_sq -
      n_reference * reference_mean_gene2 * reference_mean_gene2;
    const double reference_rho = reference_covariance /
      std::sqrt(reference_variance_gene1 * reference_variance_gene2);
    reference_spearman_rho[edge_index] = reference_rho;

    std::vector<double> gene1_sorted(n_reference), gene2_sorted(n_reference);
    std::vector<double> gene1_cumulative_rank_gene1(n_reference + 1, 0.0);
    std::vector<double> gene1_cumulative_rank_gene2(n_reference + 1, 0.0);
    std::vector<double> gene2_cumulative_rank_gene1(n_reference + 1, 0.0);
    std::vector<double> gene2_cumulative_rank_gene2(n_reference + 1, 0.0);
    for (int sorted_index = 0; sorted_index < n_reference; ++sorted_index) {
      const int gene1_reference_index = order_gene1[sorted_index];
      const int gene2_reference_index = order_gene2[sorted_index];
      gene1_sorted[sorted_index] = reference_gene1[gene1_reference_index];
      gene2_sorted[sorted_index] = reference_gene2[gene2_reference_index];
      gene1_cumulative_rank_gene1[sorted_index + 1] =
        gene1_cumulative_rank_gene1[sorted_index] + rank_gene1[gene1_reference_index];
      gene1_cumulative_rank_gene2[sorted_index + 1] =
        gene1_cumulative_rank_gene2[sorted_index] + rank_gene2[gene1_reference_index];
      gene2_cumulative_rank_gene1[sorted_index + 1] =
        gene2_cumulative_rank_gene1[sorted_index] + rank_gene1[gene2_reference_index];
      gene2_cumulative_rank_gene2[sorted_index + 1] =
        gene2_cumulative_rank_gene2[sorted_index] + rank_gene2[gene2_reference_index];
    }

    std::vector<int> case_order(n_case);
    std::iota(case_order.begin(), case_order.end(), 0);
    std::sort(case_order.begin(), case_order.end(), [&](int left, int right) {
      return case_expression(gene1, left) > case_expression(gene1, right);
    });
    std::vector<int> reference_order_descending = order_gene1;
    std::reverse(reference_order_descending.begin(), reference_order_descending.end());
    std::vector<double> gene2_unique = gene2_sorted;
    gene2_unique.erase(std::unique(gene2_unique.begin(), gene2_unique.end()), gene2_unique.end());
    std::vector<int> fenwick_tree(gene2_unique.size() + 1, 0);
    std::vector<int> joint_greater_count(n_case, 0);
    int inserted_reference_count = 0;
    int reference_pointer = 0;
    for (int case_rank = 0; case_rank < n_case; ++case_rank) {
      const int case_index = case_order[case_rank];
      const double case_gene1 = case_expression(gene1, case_index);
      const double case_gene2 = case_expression(gene2, case_index);
      while (reference_pointer < n_reference &&
             reference_gene1[reference_order_descending[reference_pointer]] > case_gene1) {
        const double reference_gene2_value = reference_gene2[reference_order_descending[reference_pointer]];
        const int fenwick_index = std::lower_bound(gene2_unique.begin(), gene2_unique.end(), reference_gene2_value) -
          gene2_unique.begin() + 1;
        for (int update_index = fenwick_index; update_index < static_cast<int>(fenwick_tree.size());
             update_index += update_index & -update_index) fenwick_tree[update_index] += 1;
        ++inserted_reference_count;
        ++reference_pointer;
      }
      const int less_or_equal_gene2 = std::upper_bound(gene2_unique.begin(), gene2_unique.end(), case_gene2) -
        gene2_unique.begin();
      int fenwick_sum = 0;
      for (int query_index = less_or_equal_gene2; query_index > 0; query_index -= query_index & -query_index) {
        fenwick_sum += fenwick_tree[query_index];
      }
      joint_greater_count[case_index] = inserted_reference_count - fenwick_sum;
    }

    double delta_sum = 0.0;
    for (int case_index = 0; case_index < n_case; ++case_index) {
      const double case_gene1 = case_expression(gene1, case_index);
      const double case_gene2 = case_expression(gene2, case_index);
      const int gene1_less_or_equal = std::upper_bound(gene1_sorted.begin(), gene1_sorted.end(), case_gene1) -
        gene1_sorted.begin();
      const int gene2_less_or_equal = std::upper_bound(gene2_sorted.begin(), gene2_sorted.end(), case_gene2) -
        gene2_sorted.begin();
      const int gene1_greater_count = n_reference - gene1_less_or_equal;
      const int gene2_greater_count = n_reference - gene2_less_or_equal;
      const double gene1_rank_sum_greater = sum_rank_gene1 - gene1_cumulative_rank_gene1[gene1_less_or_equal];
      const double gene2_rank_sum_greater = sum_rank_gene2 - gene2_cumulative_rank_gene2[gene2_less_or_equal];
      const double gene2_rank_sum_by_gene1_greater = sum_rank_gene2 - gene1_cumulative_rank_gene2[gene1_less_or_equal];
      const double gene1_rank_sum_by_gene2_greater = sum_rank_gene1 - gene2_cumulative_rank_gene1[gene2_less_or_equal];
      const double case_rank_gene1 = n_reference - gene1_greater_count + 1.0;
      const double case_rank_gene2 = n_reference - gene2_greater_count + 1.0;
      const double case_sum_rank_product = sum_rank_product + gene1_rank_sum_by_gene2_greater +
        gene2_rank_sum_by_gene1_greater + joint_greater_count[case_index] +
        case_rank_gene1 * case_rank_gene2;
      const double case_sum_rank_gene1_sq = sum_rank_gene1_sq + 2.0 * gene1_rank_sum_greater +
        gene1_greater_count + case_rank_gene1 * case_rank_gene1;
      const double case_sum_rank_gene2_sq = sum_rank_gene2_sq + 2.0 * gene2_rank_sum_greater +
        gene2_greater_count + case_rank_gene2 * case_rank_gene2;
      const double case_mean_rank = (n_reference + 2.0) / 2.0;
      const double case_covariance = case_sum_rank_product - (n_reference + 1.0) * case_mean_rank * case_mean_rank;
      const double case_variance_gene1 = case_sum_rank_gene1_sq -
        (n_reference + 1.0) * case_mean_rank * case_mean_rank;
      const double case_variance_gene2 = case_sum_rank_gene2_sq -
        (n_reference + 1.0) * case_mean_rank * case_mean_rank;
      const double case_rho = case_covariance /
        std::sqrt(case_variance_gene1 * case_variance_gene2);
      delta_sum += case_rho - reference_rho;
    }
    mean_delta_rho[edge_index] = delta_sum / n_case;
  }
  return Rcpp::List::create(
    Rcpp::Named("reference_spearman_rho") = reference_spearman_rho,
    Rcpp::Named("mean_delta_rho") = mean_delta_rho
  );
}
')

gene1_index <- match(lr_edges_ssn$Gene1, rownames(gtex_expression_adjusted))
gene2_index <- match(lr_edges_ssn$Gene2, rownames(gtex_expression_adjusted))
stopifnot(
  all(!is.na(gene1_index)), all(!is.na(gene2_index)),
  !any(vapply(rownames(gtex_expression_adjusted), function(gene) {
    any(gtex_expression_adjusted[gene, ] %in% tcga_tumor_expression_adjusted[gene, ])
  }, logical(1)))
)

validation_reference <- gtex_expression_adjusted[, sample.int(ncol(gtex_expression_adjusted), 50L, replace = TRUE), drop = FALSE]
validation_cases <- tcga_tumor_expression_adjusted[, seq_len(5L), drop = FALSE]
validation_reference_result <- calculate_spearman_ssn_perturbation(
  validation_reference, validation_cases, lr_edges_ssn
)
validation_cpp_result <- calculate_reference_summary_cpp(
  validation_reference, validation_cases, gene1_index, gene2_index
)
stopifnot(
  max(abs(validation_cpp_result$reference_spearman_rho -
          validation_reference_result$reference_spearman$reference_spearman_rho)) < 1e-10,
  max(abs(validation_cpp_result$mean_delta_rho - rowMeans(validation_reference_result$perturbation))) < 1e-10
)

reference_bootstrap_mean <- matrix(
  NA_real_, nrow = nrow(ssn_tumor_perturbation), ncol = n_reference_bootstrap,
  dimnames = list(rownames(ssn_tumor_perturbation), paste0("bootstrap_", seq_len(n_reference_bootstrap)))
)
reference_bootstrap_rho <- reference_bootstrap_mean
message("Recalculating perturbations for ", n_reference_bootstrap, " bootstrap GTEx reference cohorts.")
for (bootstrap_index in seq_len(n_reference_bootstrap)) {
  bootstrap_reference <- gtex_expression_adjusted[, sample.int(ncol(gtex_expression_adjusted), replace = TRUE), drop = FALSE]
  bootstrap_result <- calculate_reference_summary_cpp(
    bootstrap_reference, tcga_tumor_expression_adjusted, gene1_index, gene2_index
  )
  reference_bootstrap_mean[, bootstrap_index] <- bootstrap_result$mean_delta_rho
  reference_bootstrap_rho[, bootstrap_index] <- bootstrap_result$reference_spearman_rho
  if (bootstrap_index %% 100L == 0L) message("Completed ", bootstrap_index, " / ", n_reference_bootstrap, " GTEx bootstrap references.")
}

reference_bootstrap_delta_ci <- t(apply(reference_bootstrap_mean, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
reference_bootstrap_rho_ci <- t(apply(reference_bootstrap_rho, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
reference_bootstrap_direction_consistency <- rowMeans(
  sweep(reference_bootstrap_mean, 1, sign(observed_mean_delta_rho), "*") > 0
)

reference_size_bootstrap_mean <- list(`898` = reference_bootstrap_mean)
reference_size_bootstrap_rho <- list(`898` = reference_bootstrap_rho)
for (reference_size in reference_sizes[reference_sizes < ncol(gtex_expression_adjusted)]) {
  size_bootstrap_mean <- matrix(
    NA_real_, nrow = nrow(ssn_tumor_perturbation), ncol = n_reference_size_resamples,
    dimnames = list(rownames(ssn_tumor_perturbation), paste0("resample_", seq_len(n_reference_size_resamples)))
  )
  size_bootstrap_rho <- size_bootstrap_mean
  message("Assessing stability with normal-reference size n = ", reference_size, ".")
  for (resample_index in seq_len(n_reference_size_resamples)) {
    size_reference <- gtex_expression_adjusted[, sample.int(ncol(gtex_expression_adjusted), reference_size, replace = TRUE), drop = FALSE]
    size_result <- calculate_reference_summary_cpp(
      size_reference, tcga_tumor_expression_adjusted, gene1_index, gene2_index
    )
    size_bootstrap_mean[, resample_index] <- size_result$mean_delta_rho
    size_bootstrap_rho[, resample_index] <- size_result$reference_spearman_rho
    if (resample_index %% 100L == 0L) {
      message("Completed ", resample_index, " / ", n_reference_size_resamples,
              " resamples at n = ", reference_size, ".")
    }
  }
  reference_size_bootstrap_mean[[as.character(reference_size)]] <- size_bootstrap_mean
  reference_size_bootstrap_rho[[as.character(reference_size)]] <- size_bootstrap_rho
}

reference_size_stability_summary <- do.call(rbind, lapply(names(reference_size_bootstrap_mean), function(reference_size) {
  size_bootstrap_mean <- reference_size_bootstrap_mean[[reference_size]]
  size_bootstrap_rho <- reference_size_bootstrap_rho[[reference_size]]
  size_delta_ci <- t(apply(size_bootstrap_mean, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
  size_rho_ci <- t(apply(size_bootstrap_rho, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE))
  data.frame(
    ssn_edge_id = rownames(size_bootstrap_mean),
    reference_size = as.integer(reference_size),
    n_resamples = ncol(size_bootstrap_mean),
    mean_delta_rho = observed_mean_delta_rho,
    mean_delta_rho_ci_lower = size_delta_ci[, 1],
    mean_delta_rho_ci_upper = size_delta_ci[, 2],
    delta_direction_consistency = rowMeans(sweep(size_bootstrap_mean, 1, sign(observed_mean_delta_rho), "*") > 0),
    reference_spearman_rho_ci_lower = size_rho_ci[, 1],
    reference_spearman_rho_ci_upper = size_rho_ci[, 2],
    stringsAsFactors = FALSE
  )
}))

reference_size_resample_diagnostics <- do.call(rbind, lapply(names(reference_size_bootstrap_mean), function(reference_size) {
  size_bootstrap_mean <- reference_size_bootstrap_mean[[reference_size]]
  data.frame(
    reference_size = as.integer(reference_size),
    resample = colnames(size_bootstrap_mean),
    spearman_correlation_to_full_reference = apply(size_bootstrap_mean, 2, function(x) {
      stats::cor(x, observed_mean_delta_rho, method = "spearman")
    }),
    direction_agreement_to_full_reference = apply(size_bootstrap_mean, 2, function(x) {
      mean(sign(x) == sign(observed_mean_delta_rho))
    }),
    mean_absolute_deviation_from_full_reference = colMeans(abs(size_bootstrap_mean - observed_mean_delta_rho)),
    stringsAsFactors = FALSE
  )
}))

reference_size_stability <- list(
  mean_delta_rho = reference_size_bootstrap_mean,
  reference_spearman_rho = reference_size_bootstrap_rho,
  edge_summary = reference_size_stability_summary,
  resample_diagnostics = reference_size_resample_diagnostics
)

source("scripts/06_Plot_SSN_Stability.R")

ssn_cohort_significance <- cbind(
  reference_spearman,
  mean_delta_rho = observed_mean_delta_rho,
  tumor_bootstrap_ci_lower = tumor_bootstrap_ci[, 1],
  tumor_bootstrap_ci_upper = tumor_bootstrap_ci[, 2],
  tumor_bootstrap_direction_consistency = tumor_bootstrap_direction_consistency,
  empirical_sign_permutation_p = empirical_p_value,
  empirical_sign_permutation_fdr_bh = stats::p.adjust(empirical_p_value, method = "BH"),
  reference_bootstrap_mean_delta_ci_lower = reference_bootstrap_delta_ci[, 1],
  reference_bootstrap_mean_delta_ci_upper = reference_bootstrap_delta_ci[, 2],
  reference_bootstrap_delta_direction_consistency = reference_bootstrap_direction_consistency,
  reference_bootstrap_rho_ci_lower = reference_bootstrap_rho_ci[, 1],
  reference_bootstrap_rho_ci_upper = reference_bootstrap_rho_ci[, 2]
)

reference_bootstrap <- list(
  mean_delta_rho = reference_bootstrap_mean,
  reference_spearman_rho = reference_bootstrap_rho
)
stability_provenance <- data.frame(
  item = c(
    "random_seed", "correlation_method", "tumor_bootstrap_iterations", "tumor_bootstrap_unit",
    "sign_permutation_iterations", "sign_permutation_null", "reference_bootstrap_iterations",
    "reference_bootstrap_unit", "reference_bootstrap_score", "reference_size_sensitivity",
    "reference_size_resamples", "multiple_testing", "fdr_filtering"
  ),
  value = c(
    "20260730", "Spearman rank correlation", as.character(n_tumor_bootstrap), "TCGA primary tumour",
    as.character(n_sign_permutation), "independent random sign assigned to each tumour-level Delta rho",
    as.character(n_reference_bootstrap), "GTEx normal sample; resampled with replacement at original n = 898",
    "recomputed mean Delta rho across all 618 tumours and reference Spearman rho",
    paste(reference_sizes, collapse = "; "), as.character(n_reference_size_resamples),
    "Benjamini-Hochberg across retained SSN edges", "none; full perturbation matrix remains input to MOFA"
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(ssn_cohort_significance, "results/06_Assess_SSN_Stability/ssn_cohort_significance.qs2")
qs2::qs_save(tumor_bootstrap_mean, "results/06_Assess_SSN_Stability/ssn_tumor_bootstrap_mean.qs2")
qs2::qs_save(reference_bootstrap, "results/06_Assess_SSN_Stability/ssn_reference_bootstrap.qs2")
qs2::qs_save(reference_size_stability, "results/06_Assess_SSN_Stability/ssn_reference_size_stability.qs2")
qs2::qs_save(stability_provenance, "results/06_Assess_SSN_Stability/ssn_stability_provenance.qs2")

message("Saved stability statistics for ", nrow(ssn_cohort_significance), " SSN edges.")
