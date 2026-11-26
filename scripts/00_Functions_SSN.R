# Purpose: shared calculation used when recalculating Spearman SSN perturbations with an alternative
# normal reference. This helper does not read or write project objects.

calculate_spearman_ssn_perturbation <- function(reference_expression, case_expression, lr_edges_ssn) {
  stopifnot(
    identical(rownames(reference_expression), rownames(case_expression)),
    all(lr_edges_ssn$Gene1 %in% rownames(reference_expression)),
    all(lr_edges_ssn$Gene2 %in% rownames(reference_expression)),
    all(is.finite(reference_expression)),
    all(is.finite(case_expression))
  )

  reference_rho <- vapply(seq_len(nrow(lr_edges_ssn)), function(i) {
    suppressWarnings(stats::cor(
      reference_expression[lr_edges_ssn$Gene1[i], ],
      reference_expression[lr_edges_ssn$Gene2[i], ],
      method = "spearman"
    ))
  }, numeric(1))
  valid_edges <- is.finite(reference_rho)
  retained_edges <- lr_edges_ssn[valid_edges, , drop = FALSE]
  reference_rho <- reference_rho[valid_edges]

  perturbation <- vapply(colnames(case_expression), function(case_id) {
    vapply(seq_len(nrow(retained_edges)), function(i) {
      case_rho <- suppressWarnings(stats::cor(
        c(reference_expression[retained_edges$Gene1[i], ], case_expression[retained_edges$Gene1[i], case_id]),
        c(reference_expression[retained_edges$Gene2[i], ], case_expression[retained_edges$Gene2[i], case_id]),
        method = "spearman"
      ))
      case_rho - reference_rho[i]
    }, numeric(1))
  }, numeric(nrow(retained_edges)))
  rownames(perturbation) <- retained_edges$ssn_edge_id
  colnames(perturbation) <- colnames(case_expression)

  list(
    reference_spearman = cbind(retained_edges, reference_spearman_rho = reference_rho),
    perturbation = perturbation,
    dropped_edges = lr_edges_ssn[!valid_edges, , drop = FALSE]
  )
}
