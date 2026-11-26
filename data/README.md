# Inputs and attribution

| Files | Contents and provenance |
|---|---|
| `tpm/GTEx/` | The study's GTEx v10 colon TPM matrix, 57,853 symbols × 898 samples. Columns are consecutive blocks in the original order, with full numerical precision and original public sample IDs. The project used the v10 RNASeQC 2.4.2 gene TPM GCT and v10 sample attributes. The historical preparation selected colon specimens, grouped duplicate gene descriptions using the per-sample maximum and retained gene symbols. [GTEx expression downloads](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression) |
| `tpm/TCGA/`, `tcga_selected_sample_map.qs2`, `tcga_matched_pairs.qs2` | The study's TCGA COAD/READ TPM input, 38,624 symbols × 668 samples. The selected matrix contains 618 primary tumour aliquots (01A) and 50 patient-matched adjacent normal aliquots (11A). The pair map identifies the exact tumour and normal used per patient. The release supplies the actual analysed values; the original GDC download manifest and release were not recorded with this matrix. [TCGA-COAD](https://portal.gdc.cancer.gov/projects/TCGA-COAD) and [TCGA-READ](https://portal.gdc.cancer.gov/projects/TCGA-READ) |
| `tcga_overall_survival.qs2` | Public TCGA tumour sample ID, OS status (1 = event) and OS time in years for 584 patients. The other 34 tumours remain in the factor fit and have missing OS after the sample join. |
| `lineage_marker_statistics.qs2` | Fibroblast and macrophage marker statistics for 16 single-cell cohorts. Each table retains all archived rows with gene symbol, within-lineage and outside-lineage expression percentages, and adjusted P. Step 02 lists the cohort identifiers and applies the recurrence and expression thresholds. |
| `spatial_rctd_weights.qs2` | Total UMIs per spot and archived RCTD weights for 12 cell types in all 2,054 spots of GSM7058756. The cell-pair evidence script applies the strict 5th/98th percentile filter to spot totals and calculates Spearman correlations using the stored weights. [GSE225857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225857) |
| `cellchat_interaction_counts.tsv` | The full sender × receiver count matrix from the archived CellChat object for the in-house colorectal cancer cohort. Rows are senders and columns are receivers. Entries count inferred LR interactions. The evidence script reports both directions and preserves all targets. [Processed sequencing data for the in-house cohort](https://doi.org/10.5281/zenodo.18795278) |

## LR resource snapshot

Step 02 uses the gene-symbol fields from 17 human collections in LIANA. The repository resolves their original resource tables from a pinned upstream snapshot. During release preparation the ligand and receptor symbols and their row order were verified across all 17 collections against the study archive. Field names matter because the OmniPath table places its gene symbols in different positions in the package snapshot. The snapshot MD5 is `b807279afdb1b2f34fa83a3dd51bbcb5`. For an offline full run, place the matching upstream file at `results/00_Input/omni_resources.rds` before running `prepare_inputs.R`. Upstream data remain subject to their source terms. [LIANA source revision](https://github.com/saezlab/liana/tree/6cab46c54234f861ea176c3de77c4b8aa45ecb3d)

## Single-cell starting point

The original marker script read the annotated Seurat cohort list and called `RunDEA(sc, group.by = "Major.type3", return.all = TRUE)` when that field existed, using `Major.type2` otherwise. It retained rows with group `Fibroblast` or `Macrophage`. The historical wrapper version was not recorded in that script. The supplied marker tables are the starting point for reapplying every threshold and reconstructing the edge universe. Raw read processing, original single-cell annotation, original RCTD fitting and original CellChat inference precede the release entry points.

Holders of the original marker RDS, spatial Seurat RDS and CellChat RData can also run `Rscript --vanilla scripts/extract_cell_inputs.R MARKER_RDS SPATIAL_RDS CELLCHAT_RDATA` from the repository root. This optional script checks the published summaries against the processed objects and writes a copy under `results/rebuilt_cell_inputs/`.

## Archived comparisons

The sibling `reference/` folder contains study outputs used to check a rerun. The adjusted matrices contain only the LR genes needed for the SSN calculation. The edge, cohort and paired summaries retain all tested edges. The TSV files provide the oriented factor scores and loadings, explained variance and reference-size diagnostics for direct inspection. Recomputed files go under `results/`; the reference files are kept as comparison targets.
