### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 7:
## - Sample-to-sample Pearson correlation heatmap of normalized read counts
##   (log RPKM) across samples in RNAdecayCafe, annotated by cell line and
##   sequencing strategy (Full-length vs. 3'-end; ThreePseq_vs_RNAseq_heatmap.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# File paths -------------------------------------------------------------------

path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"
figure_dir               <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels_v2/"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Color palettes & helper functions --------------------------------------------

cell_cols <- c(
  "Bcell" = "#F2B5A0",
  "H1ESC" = "#EAC3FC",
  "HEK293" = "#242455",
  "HeLa" = "#9E65AA",
  "HepG2" = "#A0C2F2",
  "K562" = "#F3EA15",
  "LCL" = "#78F4AA",
  "MCF7" = "#976F2A",
  "RPE" = "#50F703",
  "Calu3" = "#962255",
  "CH22" = "#838A71",
  "MDAMB" = "#551E0D",
  "MOLM" = "#7C4292",
  "MV411" = "#7BA56F",
  "Nalm6" = "#663730",
  "U2OS" = "#D22027",
  "THP1" = "#0080FF",
  "HEL" = "black",
  "BC1" = "gray",
  "BE2C" = "blue",
  "MCC" = "magenta",
  "OCI/AML-3" = "orange"
)

col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),
  rev(brewer.pal(11, "RdYlBu"))
)

relu <- function(x) {
  ifelse(x < 0, 0, x)
}


# Compute normalized read count correlations & plot heatmap --------------------

database_list <- readRDS(path_to_rnadecaycafe_rds)
rnadegrad_db <- database_list$kdegs %>%
  dplyr::inner_join(
    database_list$feature_metadata,
    by = "feature_ID"
  ) %>%
  dplyr::inner_join(
    database_list$sample_metadata,
    by = "sample"
  )

cnt_mat <- rnadegrad_db %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    RPKM = ifelse(
      threePseq,
      (donorm_reads + 1) / (sum(donorm_reads) / 1000000),
      (donorm_reads + 1) / (exon_length / 1000) / (sum(donorm_reads) / 1000000)
    )
  ) %>%
  dplyr::mutate(
    log_RPKM = log(RPKM)
  ) %>%
  dplyr::select(
    sample, XF, log_RPKM
  ) %>%
  tidyr::pivot_wider(
    values_from = log_RPKM,
    names_from = sample
  ) %>%
  dplyr::select(-XF) %>%
  na.omit() %>%
  as.matrix()

corr_mat <- cor(cnt_mat)

annot_tbl <- rnadegrad_db %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::select(
    sample, cell_line, threePseq
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell_line = case_when(
      grepl("^MOLM", cell_line) ~ "MOLM",
      grepl("^HEK293", cell_line) ~ "HEK293",
      .default = cell_line
    )
  )

annot_tbl <- annot_tbl[match(rownames(corr_mat), annot_tbl$sample), ]

row_ha <- rowAnnotation(
  cell = annot_tbl$cell_line,
  seq = annot_tbl$threePseq,
  col = list(
    seq = c("TRUE" = "black", "FALSE" = "white"),
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_cnts <- ComplexHeatmap::Heatmap(
  relu(corr_mat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  right_annotation = row_ha,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  border = FALSE,
  heatmap_legend_param = list(
    title = "|r|",
    at = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  show_heatmap_legend = FALSE
)

pdf(
  file.path(figure_dir, "ThreePseq_vs_RNAseq_heatmap.pdf"),
  width = 2.4,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_cnts)
dev.off()

# Save version with legend
row_ha_leg <- rowAnnotation(
  cell = annot_tbl$cell_line,
  seq = annot_tbl$threePseq,
  col = list(
    seq = c("TRUE" = "black", "FALSE" = "white"),
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_cnts_leg <- ComplexHeatmap::Heatmap(
  relu(corr_mat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  right_annotation = row_ha_leg,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  border = FALSE,
  heatmap_legend_param = list(
    title = "|r|",
    at = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  show_heatmap_legend = TRUE
)

pdf(
  file.path(figure_dir, "ThreePseq_vs_RNAseq_heatmap_legend.pdf"),
  width = 6.4,
  height = 6.1,
  useDingbats = FALSE
)
draw(ht_cnts_leg)
dev.off()
