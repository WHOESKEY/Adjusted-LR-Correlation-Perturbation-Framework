# Adjusted LR-correlation perturbation framework

The framework measures how much an individual tumour changes the expression correlation of a ligand–receptor (LR) gene pair relative to a normal tissue reference. It defines candidate pairs from cell-type markers and LR resources, adjusts bulk expression for cohort and cellular composition, and accounts for purity in tumour samples. For each pair, it calculates the change in Spearman correlation when one tumour is added to the reference. The resulting sample-specific network (SSN) matrix is the input to Multi-Omics Factor Analysis (MOFA), which estimates how edge perturbations covary across tumours.

The bundled colorectal cancer analysis examines cancer-associated fibroblasts (CAFs) and macrophages. It starts with processed single-cell marker statistics, spatial cell-composition estimates, cell–cell interaction counts and bulk expression in TPM. It yields 238 LR perturbations per tumour for 618 primary tumours, using 898 normal colon specimens as the reference. The scripts include the cohort labels and sample checks for this application.

## Verify the example results

From the repository root, run the following in a shell with R, `qs2` and `survival` installed. This check uses the bundled study data and runs offline.

```sh
Rscript --vanilla review_check.R
Rscript --vanilla review_check.R --all
```

The first command recalculates every edge in three fixed tumours. The second recalculates all 147,084 perturbations. Both test agreement with the archived matrix (absolute tolerance 1e−12) and reproduce the survival-based factor ranking from the oriented study scores. The expected maximum perturbation error in the release check is zero. The smallest BH-adjusted median-split log-rank P is 0.0462873 for Factor 5. To inspect the spatial and communication evidence that motivated the CAF–macrophage pair, run `Rscript --vanilla scripts/cell_pair_evidence.R`.

## Reproduce the workflow

File digests for all bundled numerical inputs and reference tables are in [SHA256SUMS](SHA256SUMS). From a shell providing `sha256sum`, run `sha256sum -c SHA256SUMS` to verify the downloaded files.

Run each command from the repository root in the listed order. Scripts write to `results/` and `figures/` and replace their own previous outputs. The `reference/` directory contains the numerical results used to verify a run. The full run requires R with `qs2`, `sva`, `ggplot2`, `EPIC`, `IOBR`, `Rcpp`, `MOFA2`, `reticulate` and `survival`, a working R C++ toolchain for the reference bootstrap, and a Python environment with `mofapy2`. Set `RETICULATE_PYTHON` to that environment's Python executable before starting step 09. The original composition analysis recorded EPIC 1.1.7 and IOBR 2.2.1. See [release validation](VALIDATION.md) for the versions and comparisons tested here.

```sh
Rscript --vanilla scripts/cell_pair_evidence.R
Rscript --vanilla scripts/prepare_inputs.R
Rscript --vanilla scripts/01_Bulk_Processed.R
Rscript --vanilla scripts/02_Define_LR_Edges.R
Rscript --vanilla scripts/03_Estimate_Composition_Purity.R
Rscript --vanilla scripts/04_Adjust_LR_Expression.R
Rscript --vanilla scripts/05_Calculate_SSN_Perturbation.R
Rscript --vanilla scripts/06_Assess_SSN_Stability.R
Rscript --vanilla scripts/07_Validate_Paired_SSN.R
Rscript --vanilla scripts/08_Prepare_MOFA_Input.R
Rscript --vanilla scripts/09_Run_MOFA.R
Rscript --vanilla scripts/compare_outputs.R
```

Step 06 calls the plotting script automatically. Steps 06 and 07 load `00_Functions_SSN.R`. For a replay using the original EPIC/ESTIMATE covariates, run step 04 with `--archived-covariates` in place of steps 03 and 04. Its output records which covariate file was used. Step 01 and the resampling in step 06 are the most computationally demanding stages. Input preparation uses the local LIANA snapshot if its digest matches, otherwise it downloads the pinned upstream file. [Data sources and inputs](data/README.md)

## Methods and results in colorectal cancer

| Stage | Calculation and study result |
|---|---|
| Cell-pair evidence | In Visium section GSM7058756, 1,909 of 2,054 spots passed the library-size filter. The fibroblast and macrophage RCTD weights had Spearman ρ = 0.805 across these spots. Dividing each spot's weights by their sum gave ρ = 0.795. CellChat inferred 106 interactions from fibroblasts to macrophages and 50 in the reverse direction. These quantities motivated the choice of the two broad cell types. |
| Gene pairing | In each of 16 cohorts, lineage markers required expression in >20% of cells inside the lineage and <20% outside, with adjusted P < 0.05. Markers present in ≥8 cohorts were matched against 17 human LR resources in both directions. After matching bulk gene symbols, 280 directional entries collapsed to 238 unordered pairs. |
| Bulk adjustment | ComBat was fitted to log2(TPM + 1) with cohort as batch and tumour/normal tissue protected. All groups were adjusted for three joint EPIC composition PCs. Primary tumours were additionally adjusted for ESTIMATE-derived purity. Residualisation was performed separately within each group and the original within-group gene mean was restored. All 1,566 samples were retained, including 65 with flagged EPIC convergence (33 GTEx, 29 primary tumours, 3 adjacent normals). |
| Perturbation | For each edge and each tumour, Δρ is the Spearman correlation after adding that tumour to the 898-sample GTEx normal reference minus the correlation in the reference alone. This estimates a change in bulk expression association. Physical binding and signalling activity require independent evidence. |
| Stability and paired analysis | The original analysis used 1,000 tumour bootstraps, 1,000 tumour sign flips and 1,000 normal-reference draws with replacement per reference size. Reference sizes were 50, 100, 250, 500 and 898. Normal samples were resampled as specimens (898 specimens from 623 donors). The 50 paired TCGA cases were resampled by patient, with 1,000 bootstraps and 10,000 patient-level sign flips. Their tumours also contribute to the discovery cohort. Reference-size statistics describe the **mean edge profile** across tumours. |
| Factorisation | Each tumour column was centred and scaled across the 238 edges. All edges entered a single-view MOFA fit, including edges with nonsignificant individual tests. The loading with largest absolute value orients each factor to a positive sign. Factor 1 explained most variance. After fitting, median-split log-rank tests in 584 patients with 121 events prioritised Factor 5 across ten factors. This last choice uses survival outcomes and is exploratory. |

Exploratory analyses considered alternative marker thresholds, cohort recurrence requirements and factor counts. The values above specify the configuration used for this dataset. Factor labels can change between fits. To compare results across environments, align edges and tumour IDs, then match loading and score profiles allowing for differences in factor order and sign.
