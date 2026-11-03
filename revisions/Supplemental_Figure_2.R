### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 2 panels (A-B):
## - Supplemental Figure 2A: Sample-to-sample log(half-life) Pearson correlation
##                           heatmap directly comparing pulse-label and pulse-chase
##                           NR-seq datasets (Pulse_chase_vs_Pulse_label.pdf).
## - Supplemental Figure 2B: Violin/boxplot of Pearson correlations within and
##                           between pulse-chase (pc) and pulse-label (pl) datasets
##                           (PC_vs_PL_Corr_violin_box.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# File paths -------------------------------------------------------------------

path_to_all_nrseq <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/All_NRseq_kdegs_refined.csv"
figure_dir        <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels_v2"

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


# Cluster pulse-chase and pulse-label datasets individually --------------------

all_kinetics <- read_csv(path_to_all_nrseq)

### 1) Pulse-chase ordering
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
  relu(pc_corrmat),
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

ht_pc <- draw(ht_pc)
ord_idx_pc <- row_order(ht_pc)


### 2) Pulse-label ordering
pulse_label_ests <- all_kinetics %>%
  dplyr::filter(
    !pulse_chase
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

ht_pl <- ComplexHeatmap::Heatmap(
  relu(pl_corrmat),
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

ht_pl <- draw(ht_pl)
ord_idx_pl <- row_order(ht_pl)


# Supplemental Figure 2A: Combined pulse-label vs. pulse-chase heatmap ---------

nrseq_ests <- all_kinetics %>%
  dplyr::filter(
    !pulse_chase |
      (pnew >= 0.01 & tchase >= 1)
  ) %>%
  dplyr::select(sample, XF, log_kdeg) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = log_kdeg
  )

pl_samps <- unique(all_kinetics$sample[!all_kinetics$pulse_chase])
pc_samps <- unique(all_kinetics$sample[all_kinetics$pulse_chase & all_kinetics$pnew >= 0.01 & all_kinetics$tchase >= 1])

M_comb <- nrseq_ests %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_comb) <- "numeric"

combined_corrmat <- cor(
  M_comb,
  method = "pearson",
  use = "pairwise.complete.obs"
)

order_to_use <- c(
  rownames(pl_corrmat)[ord_idx_pl],
  rownames(pc_corrmat)[ord_idx_pc]
)

new_id <- match(order_to_use, rownames(combined_corrmat))

final_mat <- relu(combined_corrmat[new_id, new_id])

annot_tbl <- all_kinetics %>%
  dplyr::filter(
    !pulse_chase |
      (pnew >= 0.01 & tchase >= 1)
  ) %>%
  dplyr::select(sample, dataset, pulse_chase) %>%
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

annot_tbl <- annot_tbl[match(rownames(final_mat), annot_tbl$sample), ]

row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  method = annot_tbl$pulse_chase,
  col = list(
    cell = cell_cols,
    method = c("TRUE" = "black", "FALSE" = "white")
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
  file.path(figure_dir, "Pulse_chase_vs_Pulse_label.pdf"),
  width = 2.5,
  height = 2.1,
  useDingbats = FALSE
)
draw(ht)
dev.off()

# Save version with legend
row_ha_leg <- rowAnnotation(
  cell = annot_tbl$cell,
  method = annot_tbl$pulse_chase,
  col = list(
    cell = cell_cols,
    method = c("TRUE" = "black", "FALSE" = "white")
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_leg <- ComplexHeatmap::Heatmap(
  final_mat,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
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
  file.path(figure_dir, "Pulse_chase_vs_Pulse_label_with_legend.pdf"),
  width = 6.5,
  height = 6.1,
  useDingbats = FALSE
)
draw(ht_leg)
dev.off()


# Supplemental Figure 2B: Violin/boxplot of correlation distributions ----------

nr_vs_nr <- rownames(final_mat) %in% pl_samps
nr_vs_nr_corrs <- unique(as.vector(final_mat[nr_vs_nr, nr_vs_nr]))
nr_vs_nr_corrs <- nr_vs_nr_corrs[nr_vs_nr_corrs < 1]

non_vs_non <- rownames(final_mat) %in% pc_samps
non_vs_non_corrs <- unique(as.vector(final_mat[non_vs_non, non_vs_non]))
non_vs_non_corrs <- non_vs_non_corrs[non_vs_non_corrs < 1]

non_vs_nr_corrs <- unique(as.vector(final_mat[non_vs_non, nr_vs_nr]))
non_vs_nr_corrs <- non_vs_nr_corrs[non_vs_nr_corrs < 1]

corr_df <- tibble(
  corrs = c(nr_vs_nr_corrs,
            non_vs_non_corrs,
            non_vs_nr_corrs),
  comparison = factor(c(
    rep("pl-pl", times = length(nr_vs_nr_corrs)),
    rep("pc-pc", times = length(non_vs_non_corrs)),
    rep("pc-pl", times = length(non_vs_nr_corrs))
  ),
  levels = c("pl-pl", "pc-pc", "pc-pl"))
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
    fill = "white"
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
  filename = file.path(figure_dir, "PC_vs_PL_Corr_violin_box.pdf"),
  width = 2.25,
  height = 2.1
)
