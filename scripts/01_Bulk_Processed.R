# Purpose: batch-correct GTEx colon and TCGA-CRC expression while preserving tissue state.
# Inputs: the two TPM matrices from prepare_inputs.R and data/tcga_matched_pairs.qs2.
# Outputs: corrected matrices in results/01_Bulk_Processed/ and a PCA diagnostic in
#          figures/01_Bulk_Processed/.

library(ggplot2)
library(sva)

dir.create("figures/01_Bulk_Processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/01_Bulk_Processed", recursive = TRUE, showWarnings = FALSE)

tcga <- qs2::qs_read("results/00_Input/tcga_selected_expression_tpm.qs2")
gtex <- qs2::qs_read("results/00_Input/gtex_colon_expression_tpm.qs2")
tcga_pairs <- qs2::qs_read("data/tcga_matched_pairs.qs2")

all_tumor_ids <- colnames(tcga)[substr(colnames(tcga), 14, 16) == "01A"]
stopifnot(all(substr(all_tumor_ids, 14, 16) == "01A"))
message("Retained ", length(all_tumor_ids), " TCGA 01A primary-tumour samples.")
gtex_ids <- colnames(gtex)
shared_genes <- intersect(rownames(tcga), rownames(gtex))

# Include matched TCGA normal samples and protect Tissue in ComBat so the GTEx/TCGA cohort effect
# can be estimated without treating tumour-versus-normal biology as batch variation.
normal_ids_bc <- tcga_pairs$normal_id
combat_ids <- c(gtex_ids, all_tumor_ids, normal_ids_bc)
combat_tissue <- c(rep("Normal", length(gtex_ids)), rep("Tumor", length(all_tumor_ids)), rep("Normal", length(normal_ids_bc)))
combat_cohort <- c(rep("GTEx", length(gtex_ids)), rep("TCGA", length(all_tumor_ids) + length(normal_ids_bc)))
combat_clin <- data.frame(ID = combat_ids, Tissue = combat_tissue, Cohort = combat_cohort, row.names = combat_ids)
table(combat_clin$Cohort, combat_clin$Tissue)

combat_expr_raw <- cbind(log2(gtex[shared_genes, gtex_ids] + 1), log2(tcga[shared_genes, c(all_tumor_ids, normal_ids_bc)] + 1))
stopifnot(identical(colnames(combat_expr_raw), combat_ids))

mod_combat <- model.matrix(~Tissue, data = combat_clin)
combat_expr <- ComBat(dat = combat_expr_raw, batch = combat_clin$Cohort, mod = mod_combat)

# --- diagnostic: does correction converge GTEx-normal onto TCGA-normal while keeping TCGA-tumour
# separated, as expected if batch (not biology) was removed? -----------------------------------
pca_before <- prcomp(t(combat_expr_raw[apply(combat_expr_raw, 1, sd) > 0, ]), scale. = TRUE, rank. = 2)
pca_after <- prcomp(t(combat_expr[apply(combat_expr, 1, sd) > 0, ]), scale. = TRUE, rank. = 2)
group <- interaction(combat_clin$Cohort, combat_clin$Tissue, drop = TRUE)
centroid_dist <- function(pca) {
  cc <- aggregate(pca$x[, 1:2], list(group = group), mean)
  list(
    gtex_normal_to_tcga_normal = dist(rbind(cc[cc$group == "GTEx.Normal", 2:3], cc[cc$group == "TCGA.Normal", 2:3])),
    gtex_normal_to_tcga_tumor = dist(rbind(cc[cc$group == "GTEx.Normal", 2:3], cc[cc$group == "TCGA.Tumor", 2:3])),
    tcga_normal_to_tcga_tumor = dist(rbind(cc[cc$group == "TCGA.Normal", 2:3], cc[cc$group == "TCGA.Tumor", 2:3]))
  )
}
message("PC1-2 centroid distances BEFORE correction:"); print(centroid_dist(pca_before))
message("PC1-2 centroid distances AFTER correction:"); print(centroid_dist(pca_after))

plot_df <- rbind(
  data.frame(pca_before$x[, 1:2], Cohort = combat_clin$Cohort, Tissue = combat_clin$Tissue, stage = "Before ComBat"),
  data.frame(pca_after$x[, 1:2], Cohort = combat_clin$Cohort, Tissue = combat_clin$Tissue, stage = "After ComBat")
)
plot_df$stage <- factor(plot_df$stage, levels = c("Before ComBat", "After ComBat"))
p_combat_pca <- ggplot(plot_df, aes(PC1, PC2, color = Cohort, shape = Tissue)) +
  geom_point(alpha = 0.55, size = 1.5) +
  facet_wrap(~stage, scales = "free") +
  scale_color_manual(values = c(GTEx = "#66C2A5", TCGA = "#B2182B")) +
  labs(
    title = "PCA before and after ComBat correction",
    # subtitle = "GTEx normal, TCGA tumour, and TCGA matched normal; model preserves tissue state",
    x = "PC1",
    y = "PC2",
    color = "Cohort",
    shape = "Tissue"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30", margin = margin(b = 10)),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    strip.text = element_text(size = 15, face = "bold", margin = margin(t = 5, b = 5)),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = grid::unit(0.65, "cm"),
    legend.key.width = grid::unit(0.65, "cm"),
    panel.spacing.x = grid::unit(0.55, "cm"),
    plot.margin = margin(10, 12, 8, 10)
  )
ggsave("figures/01_Bulk_Processed/combat_pca_diagnostic.pdf", p_combat_pca, width = 10, height = 5.5)
ggsave("figures/01_Bulk_Processed/combat_pca_diagnostic.png", p_combat_pca, width = 10, height = 5.5, dpi = 300)

gtex_bc <- combat_expr[, gtex_ids]
tcga_tumor_bc <- combat_expr[, all_tumor_ids]
tcga_normal_bc <- combat_expr[, normal_ids_bc]

qs2::qs_save(gtex_bc, "results/01_Bulk_Processed/gtex_colon_expression_bc_log2tpm.qs2")
qs2::qs_save(tcga_tumor_bc, "results/01_Bulk_Processed/tcga_tumor_expression_bc_log2tpm.qs2")
qs2::qs_save(tcga_normal_bc, "results/01_Bulk_Processed/tcga_normal_expression_bc_log2tpm.qs2")
