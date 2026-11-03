### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 6 panels:
## - Left: Bar plot of the number of genes passing read count filtering in each
##         cell line in RNAdecayCafe (Ngenes_hist.pdf).
## - Right: Bar plot of log10(maximum sequencing depth) for samples from each
##          cell line in RNAdecayCafe (Nreads_hist.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)


# File paths -------------------------------------------------------------------

path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"
figure_dir               <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Supplemental_Panels_v2"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Load RNAdecayCafe database ---------------------------------------------------

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


# Left panel: Number of genes represented in each cell line --------------------

gg_ngenes <- RNAdegDB_complete %>%
  dplyr::group_by(
    cell_line
  ) %>%
  dplyr::summarise(
    ngene = length(unique(XF))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cell_line = factor(cell_line,
                       levels = cell_line[order(ngene)])
  ) %>%
  ggplot(
    aes(
      x = cell_line,
      y = ngene
    )
  ) +
  geom_bar(
    stat = "identity",
    fill = "darkgray"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab("Number of genes") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )


# Right panel: Maximum sequencing depth per cell line --------------------------

gg_nreads <- RNAdegDB_complete %>%
  dplyr::group_by(
    cell_line,
    sample
  ) %>%
  dplyr::summarise(
    reads = sum(reads)
  ) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(
    max_reads = max(reads)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cell_line = factor(cell_line,
                       levels = cell_line[order(max_reads)])
  ) %>%
  ggplot(
    aes(
      x = cell_line,
      y = log10(max_reads)
    )
  ) +
  geom_bar(
    stat = "identity",
    fill = "darkgray"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab("log10(max reads)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  coord_cartesian(ylim = c(6, 8))

ggsave(
  filename = file.path(figure_dir, "Ngenes_hist.pdf"),
  plot = gg_ngenes,
  width = 4,
  height = 2.8
)

ggsave(
  filename = file.path(figure_dir, "Nreads_hist.pdf"),
  plot = gg_nreads,
  width = 4,
  height = 2.8
)
