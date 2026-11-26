# Optional provenance step for holders of the original processed study objects.
# Arguments: cohort marker RDS, spatial RCTD Seurat RDS, CellChat RData (object "cellchat").
# Writes to results/rebuilt_cell_inputs and verifies the three bundled inputs.
library(qs2)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Expected three paths: marker RDS, spatial RCTD RDS, CellChat RData.")
}
dir.create("results/rebuilt_cell_inputs", recursive = TRUE, showWarnings = FALSE)

markers <- readRDS(args[1])
for (cohort in names(markers)) {
  for (lineage in c("Fibroblast", "Macrophage")) {
    table <- as.data.frame(markers[[cohort]][[lineage]])
    table <- table[, c("feature", "pct.in", "pct.out", "padj")]
    rownames(table) <- NULL
    markers[[cohort]][[lineage]] <- table
  }
}
stopifnot(identical(markers, qs_read("data/lineage_marker_statistics.qs2")))
qs_save(markers, "results/rebuilt_cell_inputs/lineage_marker_statistics.qs2")

spatial <- readRDS(args[2])
total <- Matrix::colSums(spatial@assays$Spatial@counts)
weights <- as.matrix(spatial@tools$RCTD@results$weights)
stopifnot(identical(rownames(weights), names(total)))
rctd <- list(total_umi = total, weights = weights)
stopifnot(identical(rctd, qs_read("data/spatial_rctd_weights.qs2")))
qs_save(rctd, "results/rebuilt_cell_inputs/spatial_rctd_weights.qs2")

objects <- new.env()
load(args[3], envir = objects)
counts <- objects$cellchat@net$count
bundled <- as.matrix(read.delim("data/cellchat_interaction_counts.tsv", row.names = 1,
                                check.names = FALSE))
stopifnot(isTRUE(all.equal(counts, bundled, check.attributes = FALSE)),
          identical(rownames(counts), rownames(bundled)), identical(colnames(counts), colnames(bundled)))
write.table(counts, "results/rebuilt_cell_inputs/cellchat_interaction_counts.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
