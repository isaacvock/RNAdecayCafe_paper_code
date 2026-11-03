### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 1 panels:
## - Schematic depicting correlation analysis:
##   1) Pulse-label NR-seq sample-to-sample log(half-life) Pearson correlation
##      heatmap (NRseq_pulselabel_heatmap_annotated.pdf).
##   2) Density-colored scatter plot of two well-correlated samples:
##      SRR15654895 vs. SRR19981794 (r = 0.89; Good_corr_example.pdf).
##   3) Density-colored scatter plot of two poorly correlated samples:
##      SRR11458920 vs. SRR14475312 (r = 0.09; Bad_corr_example.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(EZbakR)


# File paths -------------------------------------------------------------------

path_to_all_nrseq <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/All_NRseq_kdegs_refined.csv"
figure_dir        <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels"

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


# Pulse-label NR-seq correlation heatmap ---------------------------------------

all_kinetics <- read_csv(path_to_all_nrseq)

pulse_label_ests <- all_kinetics %>%
  dplyr::filter(
    !pulse_chase
  )

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


# Example 1: Well-correlated pair (SRR15654895 vs. SRR19981794) ----------------

good_corr_ex <- pl_wide %>%
  dplyr::select(
    SRR15654895, SRR19981794
  ) %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = SRR15654895,
      y = SRR19981794,
      n = 100
    )
  ) %>%
  ggplot(
    aes(
      x = SRR15654895,
      y = SRR19981794,
      color = density
    )
  ) +
  geom_point(
    size = 0.1
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  coord_cartesian(
    xlim = c(-5, 2.25),
    ylim = c(-5, 2.25)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  filename = file.path(figure_dir, "Good_corr_example.pdf"),
  plot = good_corr_ex,
  width = 2.4,
  height = 2.1
)


# Example 2: Poorly correlated pair (SRR11458920 vs. SRR14475312) --------------

bad_corr_ex <- pl_wide %>%
  dplyr::select(
    SRR11458920, SRR14475312
  ) %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = SRR14475312,
      y = SRR11458920,
      n = 100
    )
  ) %>%
  ggplot(
    aes(
      x = SRR14475312,
      y = SRR11458920,
      color = density
    )
  ) +
  geom_point(
    size = 0.1
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  coord_cartesian(
    xlim = c(-9, 2),
    ylim = c(-9, 2)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  filename = file.path(figure_dir, "Bad_corr_example.pdf"),
  plot = bad_corr_ex,
  width = 2.4,
  height = 2.1
)
