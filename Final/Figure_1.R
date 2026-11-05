### PURPOSE OF THIS SCRIPT
## Reproduce Figure 1 panels (A-F):
## - Figure 1A: Sample-to-sample log(half-life) Pearson correlation heatmap for
##              non-NR-seq datasets curated in Agarwal et al., 2022.
## - Figure 1B: Sample-to-sample log(half-life) Pearson correlation heatmap for
##              pulse-chase NR-seq datasets.
## - Figure 1C: Sample-to-sample log(half-life) Pearson correlation heatmap for
##              pulse-label NR-seq datasets.
## - Figure 1D: Sample-to-sample log(half-life) Pearson correlation heatmap for
##              RNAdecayCafe (filtered pulse-label NR-seq datasets).
## - Figure 1E: Sample-to-sample log(half-life) Pearson correlation heatmap
##              comparing RNAdecayCafe and high-concordance non-NR-seq datasets
##              from Agarwal et al., 2022.
## - Figure 1F: Violin/boxplot comparing Pearson correlations within and between
##              RNAdecayCafe (NR) and non-NR-seq (Non) datasets from Figure 1E.


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# File paths -------------------------------------------------------------------

# Input paths
# Note: Table S2 of Agarwal et al., 2022 (human sheet).
# Using processed_half_lives_TableS2.xlsx (skip = 2 to skip the 2 header note rows,
# and dropping the summary `half-life (PC1)` column to retain the 54 individual samples).
path_to_agarwal_data     <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Saluki_data/processed_half_lives_TableS2.xlsx"
path_to_all_nrseq        <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/All_NRseq_kdegs_refined.csv"
# Note: Using RNAdegDB_filtered_refined.rds (the uncorrupted R object from Build_database_refined.R)
# rather than RNAdegDB_OneTable_filtered_refined.csv, which has 32 rows corrupted by embedded null bytes.
path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"

# Output directory for Figure 1 panels
figure_dir <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels_v2"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Color palettes ---------------------------------------------------------------

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


# Figure 1A: Agarwal et al., 2022 (Saluki) database ----------------------------

saluki_hls <- readxl::read_excel(
  path_to_agarwal_data,
  sheet = "human",
  skip = 2
)

M_sal <- saluki_hls %>%
  dplyr::select(-any_of(c("Ensembl Gene Id", "Gene name", "half-life (PC1)"))) %>%
  as.matrix()

class(M_sal) <- "numeric"

new_corrmat <- cor(
  M_sal,
  method = "pearson",
  use = "pairwise.complete.obs"
)

ht_sal <- ComplexHeatmap::Heatmap(
  abs(new_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
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
  file.path(figure_dir, "Saluki_database_heatmap.pdf"),
  width = 2.2,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_sal)
dev.off()

# Save version with heatmap legend
ht_sal2 <- ComplexHeatmap::Heatmap(
  abs(new_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
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
  file.path(figure_dir, "Heatmap_with_legend.pdf"),
  width = 3,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_sal2)
dev.off()

# Identify high-concordance subset (ok_samples) and filter out NR-seq samples
ht_sal <- draw(ht_sal)
ord_idx_sal <- row_order(ht_sal)
sal_col_order <- colnames(new_corrmat)[ord_idx_sal]

bad_start <- min(which(grepl("^Gejman_", sal_col_order)))
ok_samples <- sal_col_order[1:(bad_start - 1)]

ord_idx_sal_nonrseq <- ord_idx_sal[which(!grepl("^Simon_", sal_col_order) & !grepl("^Shendure", sal_col_order))]

good_tbl <- tibble(
  sample = rownames(new_corrmat)[ord_idx_sal_nonrseq]
) %>%
  dplyr::mutate(
    good = factor(ifelse(
      sample %in% ok_samples,
      "Yes",
      "No"
    ),
    levels = c("Yes", "No")),
    cell = sapply(
      str_split(sample, "_"), function(x) x[3]
    ),
    method = sapply(
      str_split(sample, "_"), function(x) x[2]
    )
  ) %>%
  dplyr::mutate(
    cell = ifelse(
      grepl("^GM", cell),
      "LCL",
      cell
    ),
    method = case_when(
      method %in% c("ActD", "Aman") ~ "Transcription inhibition",
      .default = "Metabolic labeling (enrichment)"
    )
  )

row_ha <- rowAnnotation(
  cell = good_tbl$cell,
  method = good_tbl$method,
  col = list(
    cell = cell_cols,
    method = c(
      "Transcription inhibition" = "#D40D52",
      "Metabolic labeling (enrichment)" = "#6315B5"
    )
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_sal_nonr <- ComplexHeatmap::Heatmap(
  abs(new_corrmat[ord_idx_sal_nonrseq, ord_idx_sal_nonrseq]),
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
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
  file.path(figure_dir, "Saluki_with_annotation.pdf"),
  width = 2.5,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_sal_nonr)
dev.off()

# Save version with annotation legends
row_ha_legend <- rowAnnotation(
  cell = good_tbl$cell,
  method = good_tbl$method,
  col = list(
    cell = cell_cols,
    method = c(
      "Transcription inhibition" = "#D40D52",
      "Metabolic labeling (enrichment)" = "#6315B5"
    )
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_sal_nonr_legend <- ComplexHeatmap::Heatmap(
  abs(new_corrmat[ord_idx_sal_nonrseq, ord_idx_sal_nonrseq]),
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  right_annotation = row_ha_legend,
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
  file.path(figure_dir, "Saluki_with_annotation_legends2.pdf"),
  width = 5.5,
  height = 4.1,
  useDingbats = FALSE
)
draw(ht_sal_nonr_legend)
dev.off()


# Figure 1B: Pulse-chase NR-seq data -------------------------------------------

all_kinetics <- read_csv(path_to_all_nrseq)

pulse_chase_ests <- all_kinetics %>%
  dplyr::filter(
    pulse_chase & pnew >= 0.01 &
      tchase >= 1
  )

pc_wide <- pulse_chase_ests %>%
  dplyr::select(sample, XF, log_kdeg) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "log_kdeg"
  )

M_pc <- pc_wide %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_pc) <- "numeric"

pc_corrmat <- cor(
  M_pc,
  method = "pearson",
  use = "pairwise.complete.obs"
)

ht_pc <- ComplexHeatmap::Heatmap(
  abs(pc_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
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
  file.path(figure_dir, "NRseq_pulsechase_heatmap.pdf"),
  width = 2.2,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_pc)
dev.off()

# Annotate by cell line
annot_tbl <- pulse_chase_ests %>%
  dplyr::select(sample, dataset) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell = sapply(
      str_split(dataset, "_"), function(x) x[4]
    )
  ) %>%
  dplyr::mutate(
    cell = ifelse(
      cell == "HEK293T",
      "HEK293",
      cell
    )
  )

annot_tbl <- annot_tbl[match(rownames(pc_corrmat), annot_tbl$sample), ]

row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_pc_annot <- ComplexHeatmap::Heatmap(
  abs(pc_corrmat),
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
  file.path(figure_dir, "NRseq_pulsechase_heatmap_annotate.pdf"),
  width = 2.3,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_pc_annot)
dev.off()


# Figure 1C: Pulse-label NR-seq data -------------------------------------------

pulse_label_ests <- all_kinetics %>%
  dplyr::filter(
    !pulse_chase
  )

get_muhar_cell <- function(samples, cells) {

  cell_tbl <- tibble(
    sample = c(
      "SRR5806774", "SRR5806775", "SRR5806776",
      "SRR5806780", "SRR5806781", "SRR5806782",
      "SRR5806786", "SRR5806787", "SRR5806788",
      "SRR5806809", "SRR5806810"
    ),
    tl = 1,
    cell = c(
      "K562", "K562", "K562",
      "MV4-11", "MV4-11", "MV4-11",
      "MOLM-13", "MOLM-13", "MOLM-13",
      "OCI/AML-3", "OCI/AML-3"
    )
  )

  this_cell <- rep("", times = length(samples))
  for (s in seq_along(samples)) {
    if (cells[s] == "many") {
      this_samp <- samples[s]
      this_cell[s] <- cell_tbl$cell[cell_tbl$sample == this_samp]
    } else {
      this_cell[s] <- cells[s]
    }
  }

  return(this_cell)
}

annot_tbl <- pulse_label_ests %>%
  dplyr::select(sample, dataset) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell = sapply(
      str_split(dataset, "_"), function(x) x[4]
    )
  ) %>%
  dplyr::mutate(
    cell = case_when(
      cell == "HEK293T" ~ "HEK293",
      cell == "many" ~ get_muhar_cell(sample, cell),
      .default = cell
    )
  ) %>%
  dplyr::mutate(
    cell = case_when(
      grepl("^MOLM", cell) ~ "MOLM",
      cell == "MV4-11" ~ "MV411",
      .default = cell
    )
  )

pl_wide <- pulse_label_ests %>%
  dplyr::select(sample, XF, log_kdeg) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "log_kdeg"
  )

M_pl <- pl_wide %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_pl) <- "numeric"

pl_corrmat <- cor(
  M_pl,
  method = "pearson",
  use = "pairwise.complete.obs"
)

annot_tbl <- annot_tbl[match(rownames(pl_corrmat), annot_tbl$sample), ]

row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_pl <- ComplexHeatmap::Heatmap(
  abs(pl_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  border = FALSE,
  right_annotation = row_ha,
  heatmap_legend_param = list(
    title = "|r|",
    at = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  show_heatmap_legend = FALSE
)

pdf(
  file.path(figure_dir, "NRseq_pulselabel_heatmap_annotated.pdf"),
  width = 2.3,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_pl)
dev.off()

# Save version with legend
row_ha_legend <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_pl_legend <- ComplexHeatmap::Heatmap(
  abs(pl_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  border = FALSE,
  right_annotation = row_ha_legend,
  heatmap_legend_param = list(
    title = "|r|",
    at = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  show_heatmap_legend = TRUE
)

pdf(
  file.path(figure_dir, "NRseq_pulselabel_heatmap_annotated_legend.pdf"),
  width = 5.3,
  height = 5.1,
  useDingbats = FALSE
)
draw(ht_pl_legend)
dev.off()


# Figure 1D: RNAdecayCafe (filtered pulse-label NR-seq) ------------------------

database_list <- readRDS(path_to_rnadecaycafe_rds)
RNAdegDB_complete <- database_list$kdegs %>%
  dplyr::inner_join(
    database_list$feature_metadata,
    by = "feature_ID"
  ) %>%
  dplyr::inner_join(
    database_list$sample_metadata,
    by = "sample"
  )

pl_filt_wide <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::select(
    sample, XF, log_kdeg
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "log_kdeg"
  )

M_plf <- pl_filt_wide %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_plf) <- "numeric"

plf_corrmat <- cor(
  M_plf,
  method = "pearson",
  use = "pairwise.complete.obs"
)

ht_rdb <- ComplexHeatmap::Heatmap(
  abs(plf_corrmat),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
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
  file.path(figure_dir, "RDB_heatmap.pdf"),
  width = 2.2,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_rdb)
dev.off()

# Annotate by cell line and label time
label_times <- RNAdegDB_complete %>%
  dplyr::select(sample, label_time, cell_line, total_reads) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    label_time = factor(case_when(
      label_time < 2 ~ "< 2hr",
      label_time <= 4 ~ "[2, 4]hr",
      .default = "> 4hr"
    ),
    levels = c("< 2hr", "[2, 4]hr", "> 4hr")),
    cell_line = case_when(
      grepl("^MOLM", cell_line) ~ "MOLM",
      grepl("^HEK293", cell_line) ~ "HEK293",
      .default = cell_line
    ),
    total_reads = log10(total_reads)
  )

ann_cols <- c(
  "< 2hr" = "deepskyblue",
  "[2, 4]hr" = "darkgray",
  "> 4hr" = "darkred"
)

reads_cols <- colorRamp2(
  seq(from = min(label_times$total_reads),
      to = max(label_times$total_reads),
      length.out = 9),
  brewer.pal(9, "Greens")
)

label_times <- label_times[match(rownames(plf_corrmat), label_times$sample), ]

row_ha <- rowAnnotation(
  cell = label_times$cell_line,
  tl = label_times$label_time,
  col = list(
    tl = ann_cols,
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_rdb_col <- ComplexHeatmap::Heatmap(
  abs(plf_corrmat),
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
  file.path(figure_dir, "RDB_heatmap_annot.pdf"),
  width = 2.4,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_rdb_col)
dev.off()

# Save version with legend
row_ha_leg <- rowAnnotation(
  cell = label_times$cell_line,
  tl = label_times$label_time,
  col = list(
    tl = ann_cols,
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_rdb_col_leg <- ComplexHeatmap::Heatmap(
  abs(plf_corrmat),
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
  file.path(figure_dir, "RDB_heatmap_annot_legend.pdf"),
  width = 6.4,
  height = 6.1,
  useDingbats = FALSE
)
draw(ht_rdb_col_leg)
dev.off()


# Figure 1E: RNAdecayCafe vs. high-concordance Agarwal et al. subset -----------

nrseq_ests <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::select(sample, XF, log_kdeg) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = log_kdeg
  )

saluki_Compare <- saluki_hls %>%
  dplyr::select(-`Ensembl Gene Id`, -any_of("half-life (PC1)")) %>%
  dplyr::rename(XF = `Gene name`) %>%
  dplyr::select(XF, !!ok_samples) %>%
  dplyr::select(-dplyr::starts_with("Simon"),
                -dplyr::starts_with("Shendure")) %>%
  dplyr::mutate(across(-XF, ~ -as.numeric(.)))

non_nr_samps <- colnames(saluki_Compare)[colnames(saluki_Compare) != "XF"]
nr_samps <- colnames(nrseq_ests)[colnames(nrseq_ests) != "XF"]

combined_ests <- nrseq_ests %>%
  dplyr::inner_join(
    saluki_Compare,
    by = "XF"
  )

M_comb <- combined_ests %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_comb) <- "numeric"

combined_corrmat <- cor(
  M_comb,
  method = "pearson",
  use = "pairwise.complete.obs"
)

ht_rdb <- draw(ht_rdb)
ord_idx <- row_order(ht_rdb)

order_to_use <- c(
  rownames(plf_corrmat)[ord_idx],
  ok_samples[!grepl("^Simon_", ok_samples) & !grepl("^Shendure", ok_samples)]
)

new_id <- match(order_to_use, rownames(combined_corrmat))

final_mat <- abs(combined_corrmat[new_id, new_id])

annot_rdb <- RNAdegDB_complete %>%
  dplyr::select(
    sample, cell_line
  ) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    cell = cell_line
  ) %>%
  dplyr::mutate(
    cell = case_when(
      cell == "MOLM-13" ~ "MOLM",
      cell == "MOLM14" ~ "MOLM",
      cell == "HEK293T" ~ "HEK293",
      .default = cell
    )
  )

saluki_tbl <- tibble(
  sample = rownames(new_corrmat)[ord_idx_sal_nonrseq]
) %>%
  dplyr::mutate(
    cell = sapply(
      str_split(sample, "_"), function(x) x[3]
    )
  ) %>%
  dplyr::mutate(
    cell = ifelse(
      grepl("^GM", cell),
      "LCL",
      cell
    )
  )

cell_tbl <- bind_rows(
  annot_rdb,
  saluki_tbl
)
cell_tbl <- cell_tbl[match(rownames(final_mat), cell_tbl$sample), ]

row_ha <- rowAnnotation(
  cell = cell_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- ComplexHeatmap::Heatmap(
  final_mat,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
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
  file.path(figure_dir, "RDB_vs_Saluki_heatmap_annot.pdf"),
  width = 2.3,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht)
dev.off()

# Save version without annotation
ht_noannot <- ComplexHeatmap::Heatmap(
  final_mat,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
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
  file.path(figure_dir, "RDB_vs_Saluki_heatmap.pdf"),
  width = 2.2,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht_noannot)
dev.off()


# Figure 1F: Violin/boxplot of correlation distributions -----------------------

nr_vs_nr <- rownames(final_mat) %in% nr_samps
nr_vs_nr_corrs <- unique(as.vector(final_mat[nr_vs_nr, nr_vs_nr]))
nr_vs_nr_corrs <- nr_vs_nr_corrs[nr_vs_nr_corrs < 1]

non_vs_non <- rownames(final_mat) %in% non_nr_samps
non_vs_non_corrs <- unique(as.vector(final_mat[non_vs_non, non_vs_non]))
non_vs_non_corrs <- non_vs_non_corrs[non_vs_non_corrs < 1]

non_vs_nr_corrs <- unique(as.vector(final_mat[non_vs_non, nr_vs_nr]))
non_vs_nr_corrs <- non_vs_nr_corrs[non_vs_nr_corrs < 1]

corr_df <- tibble(
  corrs = c(nr_vs_nr_corrs,
            non_vs_non_corrs,
            non_vs_nr_corrs),
  comparison = factor(c(
    rep("NR-NR", times = length(nr_vs_nr_corrs)),
    rep("Non-Non", times = length(non_vs_non_corrs)),
    rep("Non-NR", times = length(non_vs_nr_corrs))
  ),
  levels = c("Non-Non", "Non-NR", "NR-NR"))
)

gcors <- corr_df %>%
  ggplot(
    aes(x = comparison,
        y = corrs)
  ) +
  geom_violin(
    fill = "darkgray"
  ) +
  geom_boxplot(
    width = 0.2,
    outliers = FALSE,
    fill = "white",
    notch = TRUE
  ) +
  theme_classic() +
  xlab("Comparison") +
  ylab("Pearson r") +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  plot = gcors,
  filename = file.path(figure_dir, "Corr_violin_box.pdf"),
  width = 2.25,
  height = 2.1
)
