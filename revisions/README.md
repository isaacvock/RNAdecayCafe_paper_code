# RNAdecayCafe Paper Code (Revisions)

This directory contains curated and polished R scripts to reproduce the main figures (**Figures 1–2**) and supplemental figures (**Supplemental Figures 1–10**) of the *RNAdecayCafe* manuscript.

## Scripts Overview

| Script | Reproduces | Output Panels |
| :--- | :--- | :--- |
| `Figure_1.R` | **Figure 1A–F** | Sample-to-sample $\log(t_{1/2})$ correlation heatmaps for Agarwal et al. 2022 non-NR-seq data (1A), pulse-chase NR-seq (1B), pulse-label NR-seq (1C), RNAdecayCafe (1D), RNAdecayCafe vs. high-concordance non-NR-seq subset (1E), and violin/boxplot of correlation distributions (1F). |
| `Figure_2.R` | **Figure 2A–D** | X-chromosome vs. autosome $\log(t_{1/2})$ density and cell-line boxplots (2A), within-cell $\log(t_{1/2})$ vs. $\log(\text{RPKM})$ K562 scatter plot and cell-line $R^2$ bar plot (2B), `NR6A1` and `KIF1B` across-cell-line regressions in 3'-end and full-length data (2C), and histograms of $R^2$ and $\min(m, |1/m|)$ across DE genes (2D). |
| `Supplemental_Figure_1.R` | **Supplemental Figure 1** | Pulse-label NR-seq correlation heatmap and density scatter plots for well-correlated (`SRR15654895` vs. `SRR19981794`) and poorly correlated (`SRR14475312` vs. `SRR11458920`) sample pairs. |
| `Supplemental_Figure_2.R` | **Supplemental Figure 2A–B** | Combined correlation heatmap (2A) and violin/boxplot (2B) directly comparing pulse-label and pulse-chase NR-seq datasets. |
| `Supplemental_Figure_3.R` | **Supplemental Figure 3** | Violin/boxplot comparing sample-to-sample $\log(t_{1/2})$ Pearson correlations across the 4 datasets from Figure 1A–D (`Agarwal`, `PC`, `PL`, `RNAdecayCafe`). |
| `Supplemental_Figure_4.R` | **Supplemental Figure 4A–B** | RNAdecayCafe correlation heatmap annotated by inclusion of $-\text{s}^4\text{U}$ controls (4A) and dodged boxplot of cross-dataset correlations stratified by label time and $-\text{s}^4\text{U}$ control inclusion (4B). |
| `Supplemental_Figure_5.R` | **Supplemental Figure 5** | Heatmap of LASSO regression coefficients relating 9 transcript sequence features to z-scored $\log(k_{\text{deg}})$ across RNAdecayCafe cell lines and Agarwal et al. PC1 estimates. |
| `Supplemental_Figure_6.R` | **Supplemental Figure 6** | Bar plots of the number of genes passing read count filtering (left) and $\log_{10}(\text{max sequencing depth})$ (right) per cell line in RNAdecayCafe. |
| `Supplemental_Figure_7.R` | **Supplemental Figure 7** | Sample-to-sample Pearson correlation heatmap of normalized read counts ($\log_{10}\text{RPKM}$) across RNAdecayCafe, annotated by cell line and sequencing strategy (full-length vs. 3'-end). |
| `Supplemental_Figure_8.R` | **Supplemental Figure 8** | Additional examples of cell-type-specific expression patterns across MV411, MOLM-13, and K562 in 3'-end and full-length data ($t_{1/2}$-driven: `GMNN`, `IFRD1`, `GPAA1`/`SBDS`; not $t_{1/2}$-driven: `USP45`, `BLVRB`, `PLP2`). |
| `Supplemental_Figure_9.R` | **Supplemental Figure 9A–B** | GO enrichment dot plot (`enrichGO`) (9A) and HOMER 3'-UTR FASTA generation + motif enrichment bar plot (9B) for the 99 post-transcriptionally regulated genes from Figure 2D. |
| `Supplemental_Figure_10.R` | **Supplemental Figure 10A–C** | *In vivo* mouse SLAM-ITseq (Matsushima et al. 2018) analysis in Adipose and Brain tissues: X-chromosome vs. autosome $\log(t_{1/2})$ boxplot (10A), within-tissue $t_{1/2}$ vs. TPM $R^2$ bar plot (10B), and across-tissue $\text{LFD}(\text{TPM})$ vs. $\text{LFD}(t_{1/2})$ density scatter plot (10C). |
| `helper_functions.R` | **Helper functions** | Utility functions for cell-line differential expression and stability correlation analysis (`cell_line_analysis`), dropout correction, and kinetic estimation. |

## Note on Input Files (Avoiding Null-Byte Corrupted CSVs)

Several large CSV exports in `kdeg_curation/Data/` (`Data/Analyses/All_NRseq_kdegs_refined.csv`, `Data/RNAdegDB_OneTable_filtered_refined.csv`, and `Data/Tables/AvgKdegs_genes_v3.csv`) contain embedded null bytes (`0x00`) that corrupt a small number of rows when parsed with `readr::read_csv()`. All scripts in this directory are configured to read from the uncorrupted equivalents:
- **All NR-seq kinetic estimates**: `kdeg_curation/Data/Analyses/Revisions/All_NRseq_kdegs_refined.csv` (0 null bytes)
- **RNAdecayCafe complete table**: Reconstructed directly in R from `kdeg_curation/Data/RNAdegDB_filtered_refined.rds` (`kdegs %>% inner_join(feature_metadata, by = "feature_ID") %>% inner_join(sample_metadata, by = "sample")`, matching `Build_database_refined.R`)
- **RNAdecayCafe cell-line averages**: `kdeg_curation/Data/RNAdegDB_averages_filtered_refined.csv` (0 null bytes)
