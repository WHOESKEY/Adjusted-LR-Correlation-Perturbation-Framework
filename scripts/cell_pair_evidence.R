# Assess spatial co-abundance and inferred interactions for fibroblasts and macrophages.
# Run at the repository root. Replaces two tables under results/cell_pair_evidence/.
library(qs2)
spatial <- qs_read("data/spatial_rctd_weights.qs2")
counts <- read.delim("data/cellchat_interaction_counts.tsv", row.names = 1, check.names = FALSE)
stopifnot(identical(names(spatial$total_umi), rownames(spatial$weights)))
bounds <- quantile(spatial$total_umi, c(0.05, 0.98))
keep <- spatial$total_umi > bounds[1] & spatial$total_umi < bounds[2]
# Correlate the stored RCTD weights in spots strictly between the two UMI quantiles.
weights <- spatial$weights[keep, , drop = FALSE]
rho <- cor(weights, method = "spearman")
source_cells <- c("Fibroblast", "Macrophage")
stopifnot(all(colnames(weights) %in% colnames(counts)), all(source_cells %in% rownames(counts)))
pairs <- expand.grid(source = source_cells, target = colnames(weights), stringsAsFactors = FALSE)
pairs$spatial_spearman <- rho[cbind(match(pairs$source, rownames(rho)),
                                   match(pairs$target, colnames(rho)))]
pairs$interaction_count <- as.matrix(counts)[cbind(match(pairs$source, rownames(counts)),
                                                    match(pairs$target, colnames(counts)))]
dir.create("results/cell_pair_evidence", recursive = TRUE, showWarnings = FALSE)
write.table(pairs, "results/cell_pair_evidence/source_target.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
qc <- data.frame(spots = length(keep), retained = sum(keep), lower = bounds[[1]], upper = bounds[[2]])
write.table(qc, "results/cell_pair_evidence/spatial_qc.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(pairs[pairs$source != pairs$target & pairs$target %in% source_cells, ])

# Sensitivity check for the two nominated cell types when spot weights sum to one.
weight_total <- rowSums(weights)
stopifnot(all(weight_total > 0))
row_fraction <- weights[, source_cells, drop = FALSE] / weight_total
print(data.frame(stored_weight_rho = rho["Fibroblast", "Macrophage"],
                 row_fraction_rho = cor(row_fraction[, 1], row_fraction[, 2], method = "spearman")))
