# Purpose: define the undirected CAF--macrophage ligand-receptor edge universe for SSN analysis.
# Inputs: results/*_expression_bc_log2tpm.qs2; data/lineage_marker_statistics.qs2;
#         results/00_Input/omni_resources.rds.
# Outputs: results/02_Define_LR_Edges/lr_pairs_directed.qs2,
#          results/02_Define_LR_Edges/lr_edges_ssn.qs2,
#          results/02_Define_LR_Edges/lr_marker_selection.qs2,
#          results/02_Define_LR_Edges/lr_dropped_edges.qs2, and
#          results/02_Define_LR_Edges/lr_edge_provenance.qs2.
#
# Run from the repository root with Rscript scripts/02_Define_LR_Edges.R
# The study marker rule used pct.in > 20, pct.out < 20,
# adjusted P < 0.05, and presence in at least 8 of all 16 single-cell cohorts.
# Both CAF-ligand -> macrophage-receptor and macrophage-ligand -> CAF-receptor entries are
# retained as directional annotations, then collapsed to unordered gene pairs for SSN, whose
# Spearman-correlation perturbation statistic is intrinsically undirected.

library(qs2)

dir.create("results/02_Define_LR_Edges", recursive = TRUE, showWarnings = FALSE)

gtex_bc <- qs2::qs_read("results/01_Bulk_Processed/gtex_colon_expression_bc_log2tpm.qs2")
tcga_tumor_bc <- qs2::qs_read("results/01_Bulk_Processed/tcga_tumor_expression_bc_log2tpm.qs2")
tcga_normal_bc <- qs2::qs_read("results/01_Bulk_Processed/tcga_normal_expression_bc_log2tpm.qs2")

shared_genes <- Reduce(intersect, list(rownames(gtex_bc), rownames(tcga_tumor_bc), rownames(tcga_normal_bc)))
stopifnot(length(shared_genes) > 10000)

fmdea_list <- qs_read("data/lineage_marker_statistics.qs2")
lr_resource <- readRDS("results/00_Input/omni_resources.rds")
resource_names <- c("Consensus", "Baccin2019", "CellCall", "CellChatDB", "Cellinker",
                    "CellPhoneDB", "CellTalkDB", "connectomeDB2020", "EMBRACE", "Guide2Pharma",
                    "HPMR", "ICELLNET", "iTALK", "Kirouac2010", "LRdb", "Ramilowski2015", "OmniPath")
stopifnot(all(resource_names %in% names(lr_resource)))
lr_resource <- lr_resource[resource_names]

cohort_names <- c(
  "CNP0004138", "EMTAB8107", "GSE132465", "GSE144735", "GSE161277", "GSE164522",
  "GSE166555", "GSE178318", "GSE178341", "GSE188711", "GSE200997", "GSE205506",
  "GSE236698", "GSE245552", "JIRGSU", "Inhouse"
)
stopifnot(identical(sort(cohort_names), sort(names(fmdea_list))))

pct_cut <- 20
n_cohort <- 8
padj_cut <- 0.05

fib_marker_by_cohort <- lapply(fmdea_list[cohort_names], function(cohort) {
  fib <- cohort$Fibroblast
  fib$feature[fib$pct.in > pct_cut & fib$pct.out < pct_cut & fib$padj < padj_cut]
})
mac_marker_by_cohort <- lapply(fmdea_list[cohort_names], function(cohort) {
  mac <- cohort$Macrophage
  mac$feature[mac$pct.in > pct_cut & mac$pct.out < pct_cut & mac$padj < padj_cut]
})

fib_marker_count <- sort(table(unlist(fib_marker_by_cohort)), decreasing = TRUE)
mac_marker_count <- sort(table(unlist(mac_marker_by_cohort)), decreasing = TRUE)
fib_markers <- names(fib_marker_count)[fib_marker_count >= n_cohort]
mac_markers <- names(mac_marker_count)[mac_marker_count >= n_cohort]
stopifnot(length(fib_markers) > 0, length(mac_markers) > 0)

marker_selection <- rbind(
  data.frame(gene = names(fib_marker_count), cell_type = "Fibroblast", n_cohort = as.integer(fib_marker_count),
             selected = names(fib_marker_count) %in% fib_markers, stringsAsFactors = FALSE),
  data.frame(gene = names(mac_marker_count), cell_type = "Macrophage", n_cohort = as.integer(mac_marker_count),
             selected = names(mac_marker_count) %in% mac_markers, stringsAsFactors = FALSE)
)

lr_blocks <- lapply(lr_resource, function(resource) {
  resource <- as.data.frame(resource[, c("source_genesymbol", "target_genesymbol")], stringsAsFactors = FALSE)
  colnames(resource) <- c("Ligand", "Receptor")
  resource
})
lr_pairs_raw <- do.call(rbind, lr_blocks)
lr_pairs_raw$Ligand <- trimws(as.character(lr_pairs_raw$Ligand))
lr_pairs_raw$Receptor <- trimws(as.character(lr_pairs_raw$Receptor))
lr_pairs_raw$LR <- paste(lr_pairs_raw$Ligand, lr_pairs_raw$Receptor, sep = " - ")

invalid_lr <- is.na(lr_pairs_raw$Ligand) | is.na(lr_pairs_raw$Receptor) |
  lr_pairs_raw$Ligand == "" | lr_pairs_raw$Receptor == "" | lr_pairs_raw$Ligand == lr_pairs_raw$Receptor
lr_invalid <- lr_pairs_raw[invalid_lr, c("Ligand", "Receptor", "LR"), drop = FALSE]
lr_invalid$drop_reason <- rep("missing_or_self_pair", nrow(lr_invalid))
lr_pairs_raw <- lr_pairs_raw[!invalid_lr, , drop = FALSE]
lr_pairs_raw <- lr_pairs_raw[!duplicated(lr_pairs_raw$LR), , drop = FALSE]

fib_to_mac <- lr_pairs_raw[
  lr_pairs_raw$Ligand %in% fib_markers & lr_pairs_raw$Receptor %in% mac_markers,
  c("Ligand", "Receptor", "LR"), drop = FALSE
]
fib_to_mac$L_source <- "Fibroblast"
fib_to_mac$R_source <- "Macrophage"
fib_to_mac$direction <- "CAF_ligand_to_macrophage_receptor"

mac_to_fib <- lr_pairs_raw[
  lr_pairs_raw$Ligand %in% mac_markers & lr_pairs_raw$Receptor %in% fib_markers,
  c("Ligand", "Receptor", "LR"), drop = FALSE
]
mac_to_fib$L_source <- "Macrophage"
mac_to_fib$R_source <- "Fibroblast"
mac_to_fib$direction <- "Macrophage_ligand_to_CAF_receptor"

lr_pairs_directed <- rbind(fib_to_mac, mac_to_fib)
lr_pairs_directed <- lr_pairs_directed[!duplicated(paste(lr_pairs_directed$LR, lr_pairs_directed$direction)), , drop = FALSE]
stopifnot(nrow(lr_pairs_directed) > 0)

edge_in_bulk <- lr_pairs_directed$Ligand %in% shared_genes & lr_pairs_directed$Receptor %in% shared_genes
lr_not_in_bulk <- lr_pairs_directed[!edge_in_bulk, , drop = FALSE]
lr_not_in_bulk$drop_reason <- rep("gene_absent_from_one_or_more_corrected_bulk_matrices", nrow(lr_not_in_bulk))
lr_pairs_directed <- lr_pairs_directed[edge_in_bulk, , drop = FALSE]
stopifnot(nrow(lr_pairs_directed) >= 15)

lr_pairs_directed$ssn_edge_id <- vapply(
  seq_len(nrow(lr_pairs_directed)),
  function(i) {
    genes <- c(lr_pairs_directed$Ligand[i], lr_pairs_directed$Receptor[i])
    paste(sort(genes, method = "radix"), collapse = "_")
  },
  character(1)
)

edge_groups <- split(lr_pairs_directed, lr_pairs_directed$ssn_edge_id)
lr_edges_ssn <- do.call(rbind, lapply(edge_groups, function(edge_group) {
  edge_genes <- sort(c(edge_group$Ligand[1], edge_group$Receptor[1]), method = "radix")
  direction_support <- sort(unique(edge_group$direction))
  data.frame(
    Gene1 = edge_genes[1],
    Gene2 = edge_genes[2],
    ssn_edge_id = edge_group$ssn_edge_id[1],
    direction_support = paste(direction_support, collapse = "; "),
    n_directional_entries = nrow(edge_group),
    directional_lr_entries = paste(edge_group$LR[order(edge_group$direction)], collapse = "; "),
    stringsAsFactors = FALSE
  )
}))
# Sort the two gene symbols separately so a symbol such as CD4 precedes CD40 on every platform.
edge_order <- order(lr_edges_ssn$Gene1, lr_edges_ssn$Gene2, method = "radix")
lr_edges_ssn <- lr_edges_ssn[edge_order, , drop = FALSE]
rownames(lr_edges_ssn) <- NULL

for (field in c("L_source", "R_source", "direction")) {
  lr_invalid[[field]] <- rep(NA_character_, nrow(lr_invalid))
}
lr_dropped_edges <- rbind(lr_invalid[, names(lr_not_in_bulk)], lr_not_in_bulk)
provenance <- data.frame(
  item = c(
    "marker_source", "lr_source", "directional_candidate_rule", "ssn_edge_rule", "marker_rule", "pct_cut",
    "padj_cut", "n_cohort", "n_sc_cohorts", "n_shared_bulk_genes", "n_raw_lr_pairs",
    "n_unique_valid_lr_pairs", "n_fibroblast_markers", "n_macrophage_markers",
    "n_caf_to_macrophage_candidates", "n_macrophage_to_caf_candidates", "n_directed_edges",
    "n_ssn_undirected_edges", "n_directed_entries_collapsed", "n_invalid_lr_pairs", "n_edges_missing_from_bulk"
  ),
  value = c(
    "data/lineage_marker_statistics.qs2", "results/00_Input/omni_resources.rds",
    "union of CAF-ligand_to_macrophage-receptor and macrophage-ligand_to_CAF-receptor",
    "unordered gene pair; retain all directional annotations", "pct.in > pct_cut; pct.out < pct_cut; padj < padj_cut",
    as.character(pct_cut), as.character(padj_cut), as.character(n_cohort), as.character(length(cohort_names)),
    as.character(length(shared_genes)), as.character(sum(vapply(lr_blocks, nrow, integer(1)))),
    as.character(nrow(lr_pairs_raw)), as.character(length(fib_markers)), as.character(length(mac_markers)),
    as.character(nrow(fib_to_mac)), as.character(nrow(mac_to_fib)), as.character(nrow(lr_pairs_directed)),
    as.character(nrow(lr_edges_ssn)), as.character(nrow(lr_pairs_directed) - nrow(lr_edges_ssn)),
    as.character(nrow(lr_invalid)), as.character(nrow(lr_not_in_bulk))
  ),
  stringsAsFactors = FALSE
)

qs2::qs_save(lr_pairs_directed, "results/02_Define_LR_Edges/lr_pairs_directed.qs2")
qs2::qs_save(lr_edges_ssn, "results/02_Define_LR_Edges/lr_edges_ssn.qs2")
qs2::qs_save(marker_selection, "results/02_Define_LR_Edges/lr_marker_selection.qs2")
qs2::qs_save(lr_dropped_edges, "results/02_Define_LR_Edges/lr_dropped_edges.qs2")
qs2::qs_save(provenance, "results/02_Define_LR_Edges/lr_edge_provenance.qs2")

message("CAF markers retained: ", length(fib_markers))
message("Macrophage markers retained: ", length(mac_markers))
message("Directional CAF--macrophage LR entries retained: ", nrow(lr_pairs_directed))
message("Undirected SSN edges retained: ", nrow(lr_edges_ssn))
