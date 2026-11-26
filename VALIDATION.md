# Validation of the adjusted LR-correlation perturbation framework

The release was tested from the bundled GTEx and TCGA TPM blocks and the processed single-cell and spatial inputs. The full replay estimated the cell composition and tumour purity anew, then ran steps 04–09. The study's archived covariates were used for comparison. The separate offline check starts with adjusted expression and recalculates the sample-specific perturbations directly.

## Results of the full replay

| Comparison | Observed result |
|---|---|
| Cell-pair evidence | The original processed objects yielded the bundled lineage marker statistics, RCTD weights and CellChat counts. Of 2,054 spatial spots, 1,909 passed the specified filter. The Spearman correlation between fibroblast and macrophage weights was 0.8049. Interaction counts were 106 from fibroblasts to macrophages and 50 in the reverse direction. |
| LR features | 1,581 fibroblast markers and 1,003 macrophage markers passed the recurrence rule. Two directed entries lacked a bulk gene. The 280 remaining directional entries and all 238 unordered edges matched the archive, including the final feature order. |
| Composition and purity | All 1,566 sample IDs and convergence flags matched. Maximum absolute differences were 2.23 × 10⁻¹⁶ or less for the eight EPIC fractions, 4.27 × 10⁻¹⁴ or less for the three composition PCs and 8.01 × 10⁻¹⁴ or less for tumour purity. Among the 65 flagged fits, 33 were GTEx normals, 29 were TCGA tumours and 3 were TCGA adjacent normals. All contributed to the original and replayed analysis. |
| Adjusted expression | The adjusted LR gene matrices had identical gene and sample labels. Maximum absolute differences from the archive were 8.93 × 10⁻¹⁴ for GTEx and 2.55 × 10⁻¹³ for the primary tumours. |
| SSN perturbation | All 147,084 values in the 238-edge × 618-tumour ΔSpearman matrix matched the study matrix exactly. Maximum absolute difference was zero. |
| Resampling and paired contrasts | Across the numeric columns of the cohort and 50-pair result tables, maximum absolute differences from the archived TSVs were below 5 × 10⁻¹⁶. The median correlation of the resampled mean edge profile with the full reference rose from 0.94584 at 50 normal specimens to 0.99675 at 898. In both the cohort and paired analyses, zero individual edges met BH-adjusted permutation P < 0.05. All 238 edges entered MOFA. |
| MOFA fit | All ten factors matched one to one, in the same order and orientation in this replay. The smallest corresponding score and loading correlations both exceeded 0.99999999998. The largest absolute differences were 1.53 × 10⁻⁵ for scores and 2.08 × 10⁻⁶ for loadings. Factor 1 explained 14.50% of variance and Factor 5 explained 3.51%. |
| Survival prioritisation | Among 584 patients with observed survival and 121 deaths, Factor 5 had the smallest median-split log-rank P across the ten factors. Raw P = 0.00462873 and BH-adjusted P = 0.04628730. The fit used all 618 tumours before survival testing. |

To inspect the tables, use [the cohort permutation results](reference/ssn_cohort_significance.tsv), [paired validation](reference/ssn_matched_pair_validation.tsv), [reference-size medians](reference/reference_size_medians.tsv), [factor variance](reference/mofa_factor_variance.tsv) and [all ten factor–survival tests](reference/mofa_factor_os_association.tsv).

## Tested environment

The replay used Windows, R 4.5.1, Rtools45 for step 06, qs2 0.1.6, sva 3.56.0, ggplot2 3.5.2, Rcpp 1.1.0, MOFA2 1.18.0, reticulate 1.43.0 and survival 3.8-3. The MOFA interpreter used Python 3.9.25, mofapy2 0.7.4, numpy 2.0.2, scipy 1.13.1 and h5py 3.14.0. A single BLAS thread was used for the replay.

The original analysis recorded EPIC 1.1.7 and IOBR 2.2.1. The independently tested installation used EPIC 1.1.5, IOBR 0.99.99 and estimate 1.0.13. The composition comparisons above quantify agreement for these installations. The tested source revisions are linked here. [EPIC](https://github.com/GfellerLab/EPIC/tree/dcdcbc5515ed72552f06ff6c1dce8a83332c949b) and [IOBR](https://github.com/IOBR/IOBR/tree/4ce31d009f2937f99f72aabd323e40a67fcb57c1)

The default input preparation retrieves a pinned LR database snapshot over the network if a matching local copy is unavailable. The quick and full offline checks in `review_check.R` require only R, qs2, survival and the bundled files. To repeat the full workflow and run the final output comparison, follow the commands in [the repository README](README.md). See [input provenance](data/README.md) for the entry points and the upstream raw-data processing outside this release.
