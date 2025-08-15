### PURPOSE OF THIS SCRIPT
## Focused and polished Figure 2 panels


# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(devtools)
load_all("C:/Users/isaac/Documents/Simon_Lab/EZbakR/")
library(readr)
library(arrow)
library(corrplot)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(purrr)
library(broom)
library(intervals)
library(msigdbr)
library(msigdbdf)
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ggbreak)

# Helper functions
source("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Scripts/helper_functions.R")


# Database
RNAdegDB_complete <- read_csv(
  "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/RNAdegDB_OneTable_filtered_refined.csv"
)

# Within-cell gene expression vs. kdeg -----------------------------------------

##### K562 correlation #####

k562_corr <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      cell_line == "K562" &
      !threePseq &
      seqnames %in% paste0("chr", 1:22)
  ) %>%
  dplyr::mutate(
    log_RPK = case_when(
      threePseq ~ log(reads),
      .default = log(reads / (exon_length / 1000))
    )
  ) %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(log(2)/exp(donorm_log_kdeg)),
      y = log_RPK,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(log(2)/exp(donorm_log_kdeg)),
        y = log_RPK)
  ) +
  geom_point(
    aes(color = density),
    size = 0.2
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(t1/2)") +
  ylab("log(expression)") +
  stat_smooth(
    method = "lm",
    linewidth = 0.4,
    color = "darkred"
  )

k562_corr

ggsave(
  plot = k562_corr +
    theme(legend.position = "none"),
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/K562_thalf_vs_expression.pdf",
  width = 2.5,
  height = 2
)
ggsave(
  plot = k562_corr,
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/K562_thalf_vs_expression_legend.pdf",
  width = 2.5,
  height = 2.0875
)


##### Cell type R^2's #####


### How well does kdeg explain gene expression
kdeg_vs_RNA_df <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      seqnames %in% paste0("chr", 1:22)
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    log_RPKM = case_when(
      threePseq ~ log(reads / (sum(reads) / 1000000) ),
      .default = log((reads / (sum(reads) / 1000000)) / (exon_length / 1000))
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    log_thalf = log(log(2)/exp(donorm_log_kdeg))
  ) %>%
  dplyr::group_by(
    sample, dataset, cell_line, threePseq, pnew, pold, total_reads, median_halflife
  ) %>%
  dplyr::summarise(
    R2_slope1 = 1 / (1 + var(log_RPKM + donorm_log_kdeg) / var(log_RPKM)),
    kdeg_vs_RNA_corr = abs(cor(log(log(2)/exp(donorm_log_kdeg)), log_RPKM)),
    ksyn_vs_RNA_corr = abs(cor(log_RPKM + donorm_log_kdeg, log_RPKM)),
    R2 = cor(log(log(2)/exp(donorm_log_kdeg)), log_RPKM) ^ 2,
    slope = lm(log_RPKM ~ log_thalf, .)$coefficients[2],
    CCC = (2*mean((log_thalf - mean(log_thalf))*(log_RPKM - mean(log_RPKM))) ) / (var(log_thalf) + var(log_RPKM)) 
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    kdeg_RNA_R2_adj = kdeg_vs_RNA_corr^2 / (kdeg_vs_RNA_corr^2 + ksyn_vs_RNA_corr^2)
  )



### Figure 2A right
within_cell_R2_withprev <- kdeg_vs_RNA_df %>%
  dplyr::filter(
    cell_line != "BC1"
  ) %>%
  dplyr::mutate(
    cell_line = ifelse(
      grepl("^MOLM", cell_line),
      "MOLM",
      cell_line
    )
  ) %>%
  dplyr::group_by(
    cell_line
  ) %>%
  dplyr::summarise(
    max_corr = max(R2),
    avg_corr = mean(kdeg_vs_RNA_corr),
    max_R2_adj = max(kdeg_RNA_R2_adj),
    avg_R2_adj = mean(kdeg_RNA_R2_adj),
    threePseq = threePseq[kdeg_vs_RNA_corr == max(kdeg_vs_RNA_corr)]
  ) %>%
  dplyr::mutate(
    cell_line = factor(cell_line,
                       levels = cell_line[order(max_R2_adj)])
  ) %>%
  ggplot(
    aes(
      x = cell_line,
      y = max_R2_adj
    )
  ) +
  geom_bar(
    stat = "identity"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("R^2") +
  coord_cartesian(
    ylim = c(0, 1)
  ) + 
  geom_hline(
    yintercept = 0.06,
    color = "orange",
    linetype = "dotted",
    linewidth = 0.4
  ) + 
  xlab("Cell line") +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10))


within_cell_R2_withprev

ggsave(
  plot = within_cell_R2_withprev +
    theme(
      legend.position = "none"
    ),
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/thalf_vs_expression_R2s_previouscol.pdf",
  width = 4,
  height = 2.4
)


# X-chromosome RNA stability ---------------------------------------------------

##### All cell lines together (half-life) #####

chrX_thalfs <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::mutate(
    chrX = factor(seqnames == "chrX",
                  levels = c(TRUE, FALSE))
  ) %>%
  ggplot(
    aes(
      x = log(log(2)/exp(donorm_log_kdeg)),
      color = chrX
    )
  ) +
  geom_density(
    linewidth = 0.4
  ) +
  theme_classic() + 
  scale_color_manual(
    values = c("deepskyblue", "darkgray")
  ) +
  xlab("log(t1/2)") +
  ylab("density") + 
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")

ggsave(
  plot = chrX_thalfs,
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/chrX_thalfs.pdf",
  width = 2.5,
  height = 2
)


##### Each cell line separately #####

chrX_thalfs_cell <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" &
      cell_line != "BC1"
  ) %>%
  dplyr::mutate(
    chrX = factor(seqnames == "chrX",
                  levels = c(FALSE, TRUE)),
    cell_line = ifelse(
      cell_line == "MOLM-13",
      "MOLM14",
      cell_line
    )
  ) %>%
  ggplot(
    aes(
      x = cell_line,
      y = log(log(2)/exp(donorm_log_kdeg)),
      fill = chrX
    )
  ) +
  #geom_violin() +
  geom_boxplot(
    color = "black",
    outliers = FALSE,
    notch = TRUE,
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c("darkgray", "deepskyblue")
  ) +
  theme_classic() + 
  xlab("Cell line") +
  ylab("log(t1/2)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")

chrX_thalfs_cell

ggsave(
  plot = chrX_thalfs_cell,
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/chrX_thalfs_bycell.pdf",
  width = 4,
  height = 2.4
)


# cell-line expression variance ------------------------------------------------


##### Correlation analysis #####


### Muhar data
db_3pend <- RNAdegDB_complete %>%
  dplyr::filter(
    dataset == "Muhar_etal_2018_many"
  )

kdeg_corr_3pend <- cell_line_analysis(db_3pend,
                                      FDR_cutoff = 0.01)


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

kdeg_corr_full_len_subset <- cell_line_analysis(db_full_len_subset,
                                                FDR_cutoff = 0.01)


##### Data for plots #####

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
    KiR = prod(RPKM)^(1/dplyr::n())
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
    KiR = prod(RPKM)^(1/dplyr::n())
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


##### Compile correlation analyses #####


combined_tbl <- kdeg_corr_3pend %>%
  dplyr::mutate(
    threePseq = TRUE
  ) %>%
  dplyr::bind_rows(
    kdeg_corr_full_len_subset %>%
      dplyr::mutate(
        threePseq = FALSE
      )
  ) %>%
  na.omit() %>%
  ### Only consider genes differentially expressed by both
  dplyr::group_by(
    XF
  ) %>%
  dplyr::filter(
    all(c(TRUE, FALSE) %in% threePseq)
  ) %>%
  dplyr::mutate(
    slope = slope_samp
  )

combined_tbl <- combined_tbl %>%
  dplyr::mutate(
    R2_samp = ifelse(
      slope_samp <= 0,
      0,
      R2_samp
    ),
    slope = pmin(
      slope_samp, abs(1/slope_samp)
    )
  )


R2_cutoff <- 0.25
slope_cutoff <- 0.25

combined_tbl <- combined_tbl %>%
  dplyr::group_by(
    XF
  ) %>%
  dplyr::mutate(
    both = all(R2_samp > R2_cutoff) & all(slope > slope_cutoff) #& all(slope_samp < 2)
  )


##### Specific gene scatters #####


### NR6A1 (Full-length)
full_len_scatter_post <- full_len_data %>%
  dplyr::filter(
    XF == "NR6A1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2)/exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  )  +
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
  ylab("log(RPKM)")  +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none") +
  coord_cartesian(
    ylim = c(
      -0.5,
      3.5
    ),
    xlim = c(
      -1.5,
      2.5
    )
  )


### NR6A1 (3'-end)
threep_scatter_post <- threep_data %>%
  dplyr::filter(
    XF == "NR6A1"
  ) %>%
  ggplot(
    aes(
      x = log(log(2)/exp(donorm_log_kdeg)),
      y = log_RPKM
    )
  ) +
  stat_smooth(
    method = "lm",
    color = "black",
    linetype = "dotted",
    linewidth = 0.4
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
  ylab("log(RPKM)")  +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none") +
  coord_cartesian(
    xlim = c(-2, 1.5),
    ylim = c(1.5, 5)
  )


### NR6A1 (Full-length)
full_len_scatter_none <- full_len_data %>%
  dplyr::filter(
    XF == "KIF1B"
  ) %>%
  ggplot(
    aes(
      x = log(log(2)/exp(donorm_log_kdeg)),
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
  ylab("log(RPKM)")  +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none") +
  coord_cartesian(
    ylim = c(
      -1.5,
      2
    ),
    xlim = c(
      0,
      3.5
    )
  )


### NR6A1 (3'-end)
threep_scatter_none <- threep_data %>%
  dplyr::filter(
    XF == "KIF1B"
  ) %>%
  ggplot(
    aes(
      x = log(log(2)/exp(donorm_log_kdeg)),
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
  ylab("log(RPKM)")  +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none") +
  coord_cartesian(
    ylim = c(
      0.5,
      4.5
    ),
    xlim = c(
      0,
      4
    )
  )


threep_scatter_post
full_len_scatter_post

threep_scatter_none
full_len_scatter_none


ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/FullRNA_cell_NR6A1_plot_regression.pdf",
  plot = full_len_scatter_post,
  width = 2.4,
  height = 1.25
)
ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/ThreepSeq_cell_NR6A1_plot_regression.pdf",
  plot = threep_scatter_post,
  width = 2.4,
  height = 1.25
)


ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/FullRNA_cell_KIF1B_plot_regression.pdf",
  plot = full_len_scatter_none,
  width = 2.4,
  height = 1.25
)
ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/ThreepSeq_cell_KIF1B_plot_regression.pdf",
  plot = threep_scatter_none,
  width = 2.4,
  height = 1.25
)



##### Figure 2C right: R2 and slope distributions #####

slope_full_len <- combined_tbl %>%
  dplyr::mutate(
    slope = ifelse(
      slope < 0,
      0,
      slope
    )
  ) %>%
  dplyr::filter(
    !threePseq
  ) %>%
  ggplot(
    aes(x = slope,
        y = after_stat(count) / sum(after_stat(count)),
        fill = both)) +
  geom_histogram(binwidth = 0.05,
                 position = "stack",
                 colour   = "black",
                 linewidth = .2) +
  scale_fill_manual(values = col_vals) +
  theme_classic() +
  xlab("min(m, |1/m|)") +
  scale_y_break(
    breaks = c(0.1, 0.3),   # gap from 100-1500
    scales = 0.5,            # upper part gets 30 % of the height
    space  = 0.1            # width of the zig-zag
  ) +
  scale_y_continuous(
    limits = c(0, 0.5)
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  ) +
  ylab("Fraction of DE genes") +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")



R2_full_len <- combined_tbl %>%
  dplyr::filter(
    !threePseq
  ) %>%
  ggplot(
    aes(x = R2_samp,
        y = after_stat(count) / sum(after_stat(count)),
        fill = both)) +
  geom_histogram(binwidth = 0.05,
                 position = "stack",
                 colour   = "black",
                 linewidth = .2) +
  scale_fill_manual(values = col_vals) +
  theme_classic() +
  xlab("R2") +
  scale_y_break(
    breaks = c(0.1, 0.3),   # gap from 100-1500
    scales = 0.5,            # upper part gets 30 % of the height
    space  = 0.1            # width of the zig-zag
  ) +
  scale_y_continuous(
    limits = c(0, 0.5)
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  ) +
  ylab("Fraction of DE genes") +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")





slope_threep <- combined_tbl %>%
  dplyr::mutate(
    slope = ifelse(
      slope < 0,
      0,
      slope
    )
  ) %>%
  dplyr::filter(
    threePseq
  ) %>%
  ggplot(
    aes(x = slope,
        y = after_stat(count) / sum(after_stat(count)),
        fill = both)) +
  geom_histogram(binwidth = 0.05,
                 position = "stack",
                 colour   = "black",
                 linewidth = .2) +
  scale_fill_manual(values = col_vals) +
  theme_classic() +
  xlab("min(m, |1/m|)") +
  scale_y_break(
    breaks = c(0.1, 0.3),   # gap from 100-1500
    scales = 0.5,            # upper part gets 30 % of the height
    space  = 0.1            # width of the zig-zag
  ) +
  scale_y_continuous(
    limits = c(0, 0.5)
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  ) +
  ylab("Fraction of DE genes") +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")


R2_threep <- combined_tbl %>%
  dplyr::filter(
    threePseq
  ) %>%
  ggplot(
    aes(x = R2_samp,
        y = after_stat(count) / sum(after_stat(count)),
        fill = both)) +
  geom_histogram(binwidth = 0.05,
                 position = "stack",
                 colour   = "black",
                 linewidth = .2) +
  scale_fill_manual(values = col_vals) +
  theme_classic() +
  xlab("R2") +
  scale_y_break(
    breaks = c(0.1, 0.3),   # gap from 100-1500
    scales = 0.5,            # upper part gets 30 % of the height
    space  = 0.1            # width of the zig-zag
  ) +
  scale_y_continuous(
    limits = c(0, 0.5)
  ) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
  ) +
  ylab("Fraction of DE genes") +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=8), #change font size of legend text
        legend.title=element_text(size=10),
        legend.position = "none")



ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/slopes_full_len.pdf",
  plot = slope_full_len,
  width = 2.2,
  height = 2.5
)
ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/R2s_full_len.pdf",
  plot = R2_full_len,
  width = 2.2,
  height = 2.5
)

ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/slopes_threep.pdf",
  plot = slope_threep,
  width = 2.2,
  height = 2.5
)
ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure2_panels/R2s_threep.pdf",
  plot = R2_threep,
  width = 2.2,
  height = 2.5
)

