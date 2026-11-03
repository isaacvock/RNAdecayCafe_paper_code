### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 8 panels:
## - Additional examples of genes with cell-type-specific expression patterns
##   across MV411, MOLM-13/14, and K562 in 3'-end and Full-length data:
##   1) t1/2 driven: GMNN, IFRD1, and GPAA1 (saved as SBDS_post_scatter_*.pdf).
##   2) Not t1/2 driven: USP45, BLVRB, and PLP2.


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)


# File paths -------------------------------------------------------------------

path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"
figure_dir               <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Supplemental_Panels"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Normalize 3'-end and full-length data across MV411, MOLM, and K562 -----------

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

### Muhar data (3'-end)
db_3pend <- RNAdegDB_complete %>%
  dplyr::filter(
    dataset == "Muhar_etal_2018_many"
  )

### Full-length RNA data
datasets_to_compare <- c(
  "Ietswaart_etal_2024_K562",
  "Swartzel_etal_2022_MOLM14",
  "Harada_etal_2022_MV411"
)

db_full_len_subset <- RNAdegDB_complete %>%
  dplyr::filter(
    dataset %in% datasets_to_compare &
      label_time <= 2
  )

full_len_data <- db_full_len_subset %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    RPKM = ifelse(
      threePseq,
      donorm_reads / (sum(donorm_reads) / 1000000),
      donorm_reads / (exon_length / 1000) / (sum(donorm_reads) / 1000000)
    )
  ) %>%
  ### Median of ratios normalization
  dplyr::group_by(XF) %>%
  dplyr::mutate(
    KiR = prod(RPKM)^(1 / dplyr::n())
  ) %>%
  dplyr::mutate(
    ratio = RPKM / KiR
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    scale = median(ratio)
  ) %>%
  dplyr::mutate(
    RPKM = RPKM / scale
  ) %>%
  dplyr::mutate(
    log_RPKM = log(RPKM)
  )

threep_data <- db_3pend %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    RPKM = ifelse(
      threePseq,
      donorm_reads / (sum(donorm_reads) / 1000000),
      donorm_reads / (exon_length / 1000) / (sum(donorm_reads) / 1000000)
    )
  ) %>%
  ### Median of ratios normalization
  dplyr::group_by(XF) %>%
  dplyr::mutate(
    KiR = prod(RPKM)^(1 / dplyr::n())
  ) %>%
  dplyr::mutate(
    ratio = RPKM / KiR
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    scale = median(ratio)
  ) %>%
  dplyr::mutate(
    RPKM = RPKM / scale
  ) %>%
  dplyr::mutate(
    log_RPKM = log(RPKM)
  )


# Post-transcriptionally driven examples (full-length) -------------------------

full_len_scatter_post_S1 <- full_len_data %>%
  dplyr::filter(
    XF == "GMNN"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(3.5, 6),
    xlim = c(0, 2.5)
  )

full_len_scatter_post_S2 <- full_len_data %>%
  dplyr::filter(
    XF == "IFRD1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(0.5, 3.5),
    xlim = c(-1, 2)
  )

full_len_scatter_post_S3 <- full_len_data %>%
  dplyr::filter(
    XF == "GPAA1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2.5, 5),
    xlim = c(1.5, 4)
  )


# Post-transcriptionally driven examples (3'-end) ------------------------------

threep_scatter_post_S1 <- threep_data %>%
  dplyr::filter(
    XF == "GMNN"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(4.5, 6),
    xlim = c(0, 1.5)
  )

threep_scatter_post_S2 <- threep_data %>%
  dplyr::filter(
    XF == "IFRD1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2, 5),
    xlim = c(-1, 2)
  )

threep_scatter_post_S3 <- threep_data %>%
  dplyr::filter(
    XF == "GPAA1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2.5, 7),
    xlim = c(1.5, 6)
  )


# Not post-transcriptionally driven examples (full-length) ---------------------

full_len_scatter_nopost_S1 <- full_len_data %>%
  dplyr::filter(
    XF == "USP45"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(-0.5, 2.5),
    xlim = c(0, 3)
  )

full_len_scatter_nopost_S2 <- full_len_data %>%
  dplyr::filter(
    XF == "BLVRB"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2.5, 6),
    xlim = c(1.5, 5)
  )

full_len_scatter_nopost_S3 <- full_len_data %>%
  dplyr::filter(
    XF == "PLP2"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2.5, 6),
    xlim = c(1, 4.5)
  )


# Not post-transcriptionally driven examples (3'-end) --------------------------

threep_scatter_nopost_S1 <- threep_data %>%
  dplyr::filter(
    XF == "USP45"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(1.5, 4.5),
    xlim = c(-0.5, 2.5)
  )

threep_scatter_nopost_S2 <- threep_data %>%
  dplyr::filter(
    XF == "BLVRB"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(2, 6.5),
    xlim = c(0.5, 5)
  )

threep_scatter_nopost_S3 <- threep_data %>%
  dplyr::filter(
    XF == "PLP2"
  ) %>%
  ggplot(
    aes(
      x = log(log(2) / exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_point(
    size = 0.75,
    aes(
      color = cell_line
    )
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("deepskyblue",
               "forestgreen",
               "darkred")
  ) +
  xlab("log(kdeg)") +
  ylab("log(RPKM)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(
    ylim = c(3.5, 6),
    xlim = c(1.5, 4)
  )


# Save all figures -------------------------------------------------------------

# Post-transcriptionally driven
ggsave(
  filename = file.path(figure_dir, "GMNN_post_scatter_full_len.pdf"),
  plot = full_len_scatter_post_S1,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "GMNN_post_scatter_threep.pdf"),
  plot = threep_scatter_post_S1,
  width = 2.4,
  height = 1.75
)

ggsave(
  filename = file.path(figure_dir, "IFRD1_post_scatter_full_len.pdf"),
  plot = full_len_scatter_post_S2,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "IFRD1_post_scatter_threep.pdf"),
  plot = threep_scatter_post_S2,
  width = 2.4,
  height = 1.75
)

ggsave(
  filename = file.path(figure_dir, "SBDS_post_scatter_full_len.pdf"),
  plot = full_len_scatter_post_S3,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "SBDS_post_scatter_threep.pdf"),
  plot = threep_scatter_post_S3,
  width = 2.4,
  height = 1.75
)

# Not post-transcriptionally driven
ggsave(
  filename = file.path(figure_dir, "USP45_nopost_scatter_full_len.pdf"),
  plot = full_len_scatter_nopost_S1,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "USP45_nopost_scatter_threep.pdf"),
  plot = threep_scatter_nopost_S1,
  width = 2.4,
  height = 1.75
)

ggsave(
  filename = file.path(figure_dir, "BLVRB_nopost_scatter_full_len.pdf"),
  plot = full_len_scatter_nopost_S2,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "BLVRB_nopost_scatter_threep.pdf"),
  plot = threep_scatter_nopost_S2,
  width = 2.4,
  height = 1.75
)

ggsave(
  filename = file.path(figure_dir, "PLP2_nopost_scatter_full_len.pdf"),
  plot = full_len_scatter_nopost_S3,
  width = 2.4,
  height = 1.75
)
ggsave(
  filename = file.path(figure_dir, "PLP2_nopost_scatter_threep.pdf"),
  plot = threep_scatter_nopost_S3,
  width = 2.4,
  height = 1.75
)
