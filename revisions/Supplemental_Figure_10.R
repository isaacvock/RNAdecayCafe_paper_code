### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 10 panels (A-C):
## - In vivo mouse SLAM-seq (SLAM-ITseq; Matsushima et al., 2018) data from
##   Adipose and Brain tissues processed with fastq2EZbakR and EZbakR:
##   1) Supplemental Figure 10A: Notched boxplot of log(t1/2) for X-chromosomal
##                               vs. autosomal genes in Adipose and Brain
##                               (SLAMIT_chrX_kdegs.pdf).
##   2) Supplemental Figure 10B: Bar plot of R^2 when correlating log(t1/2)
##                               and log(TPM) in Adipose and Brain
##                               (SLAMIT_kdeg_vs_TPM.pdf).
##   3) Supplemental Figure 10C: Density-colored scatter plot of log-fold
##                               differences in TPM vs. log-fold differences
##                               in t1/2 between Adipose and Brain
##                               (SLAMIT_celltype_diffs.pdf).


# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(arrow)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(EZbakR)


# File paths -------------------------------------------------------------------

dataset_str       <- "Matsushima_etal_2018_adipocyte_and_brain"
path_to_mouse_dir <- "C:/Users/isaac/Yale University/Simon Lab (General) - RNAdegDB_curation/fastq2EZbakR_output/Mouse_data/"
path_to_mouse_gtf <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Genomes/Mmus/Mus_musculus.GRCm39.112_filtered.gtf"
figure_dir        <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Load SLAM-ITseq dataset & estimate kinetics with EZbakR ----------------------

metadf <- tibble(
  sample = c(
    "Adipocyte_1",
    "Adipocyte_2",
    "Adipocyte_3",
    "Brain_1",
    "Brain_2",
    "Brain_3"
  ),
  tl = 1,
  celltype = rep(c("Adipocyte", "Brain"),
                 each = 3)
)

### Run settings
use_arrow <- TRUE
dropout_correct <- FALSE
dropout_normalize <- TRUE
normalize_across_tl <- FALSE

### Analysis
if (use_arrow) {
  cB <- arrow::open_dataset(paste0(path_to_mouse_dir, dataset_str, "/arrow_dataset"))
  ezbdo <- EZbakRArrowData(cB, metadf)
} else {
  cB <- fread(paste0(path_to_mouse_dir, dataset_str, "/cB/cB.csv.gz"))
  ezbdo <- EZbakRData(cB, metadf)
}

ezbdo <- EstimateFractions(
  ezbdo,
  features = "XF",
  pold_from_nolabel = FALSE
)

if (dropout_correct) {
  ezbdo <- CorrectDropout(ezbdo)
}

if (dropout_normalize) {
  ezbdo <- NormalizeForDropout(
    ezbdo,
    normalize_across_tls = normalize_across_tl
  )
}

ezbdo <- EstimateKinetics(ezbdo)


# Annotate chromosome & exonic length from mouse GTF ---------------------------

gtf <- rtracklayer::import(path_to_mouse_gtf)

chrdf <- as_tibble(gtf) %>%
  dplyr::filter(
    type %in% c("gene", "transcript")
  ) %>%
  dplyr::select(
    seqnames, gene_id
  ) %>%
  dplyr::distinct() %>%
  na.omit() %>%
  dplyr::rename(
    XF = gene_id
  )

txdb <- txdbmaker::makeTxDbFromGFF(path_to_mouse_gtf)

ex_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
ex_by_gene_reduced <- GenomicRanges::reduce(ex_by_gene, ignore.strand = FALSE)

len_by_gene <- sum(GenomicRanges::width(ex_by_gene_reduced))
df <- tibble::tibble(
  XF = names(len_by_gene),
  length = as.integer(unlist(len_by_gene))
)


# Supplemental Figure 10B: Correlation of degradation with expression ----------

kdeg_vs_exp <- EZget(
  ezbdo,
  type = "kinetics"
) %>%
  dplyr::inner_join(
    df,
    by = "XF"
  ) %>%
  dplyr::mutate(
    RPK = n / (length / 1000)
  ) %>%
  dplyr::group_by(
    sample
  ) %>%
  dplyr::mutate(
    TPM = RPK / (sum(RPK) / 1000000)
  ) %>%
  dplyr::filter(
    TPM > 1 & n > 50
  ) %>%
  dplyr::inner_join(
    metadf,
    by = "sample"
  ) %>%
  dplyr::group_by(
    celltype
  ) %>%
  dplyr::summarise(
    R2 = cor(log_kdeg, log10(TPM))^2
  ) %>%
  ggplot(
    aes(
      x = celltype,
      y = R2
    )
  ) +
  geom_bar(
    fill = "darkgray",
    color = "black",
    stat = "identity",
    linewidth = 0.4
  ) +
  theme_classic() +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )


# Supplemental Figure 10A: ChrX stability --------------------------------------

chrX_kdegs <- EZget(
  ezbdo,
  type = "kinetics"
) %>%
  dplyr::inner_join(
    chrdf,
    by = "XF"
  ) %>%
  dplyr::inner_join(
    metadf,
    by = "sample"
  ) %>%
  dplyr::filter(
    n > 50
  ) %>%
  dplyr::mutate(
    chrX = factor(seqnames == "X",
                  levels = c(FALSE, TRUE)),
    log_thalf = log(log(2) / exp(log_kdeg))
  ) %>%
  ggplot(
    aes(
      x = celltype,
      y = log_thalf,
      fill = chrX
    )
  ) +
  geom_boxplot(
    color = "black",
    outliers = FALSE,
    notch = TRUE,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c("darkgray", "deepskyblue")
  ) +
  theme_classic() +
  xlab("Cell line") +
  ylab("log(t1/2)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )


# Supplemental Figure 10C: Stability vs. expression differences ----------------

kinetics <- EZget(
  ezbdo,
  type = "kinetics"
) %>%
  dplyr::inner_join(
    metadf,
    by = "sample"
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    TPM = n / (sum(n) / 1000000),
    log_thalf = log(log(2) / kdeg)
  ) %>%
  dplyr::filter(
    n > 50
  ) %>%
  dplyr::group_by(
    celltype, XF
  ) %>%
  dplyr::summarise(
    avg_lthalf = mean(log_thalf),
    avg_lTPM = mean(log10(TPM))
  ) %>%
  dplyr::group_by(
    XF
  ) %>%
  dplyr::filter(
    dplyr::n() == 2
  ) %>%
  dplyr::summarise(
    LFC_thalf = avg_lthalf[celltype == "Adipocyte"] - avg_lthalf[celltype == "Brain"],
    LFC_TPM = avg_lTPM[celltype == "Adipocyte"] - avg_lTPM[celltype == "Brain"]
  )

celltype_diffs <- kinetics %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = LFC_thalf,
      y = LFC_TPM,
      n = 100
    )
  ) %>%
  ggplot(
    aes(
      x = LFC_thalf,
      y = LFC_TPM
    )
  ) +
  geom_point(
    size = 0.2,
    aes(
      color = density
    )
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("LFC(t1/2)") +
  ylab("LFC(TPM)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  coord_cartesian(
    xlim = c(-8, 6.5),
    ylim = c(-8, 6.5)
  ) +
  stat_smooth(
    method = "lm",
    color = "darkred",
    linewidth = 0.4
  )


# Save plots -------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "SLAMIT_kdeg_vs_TPM.pdf"),
  plot = kdeg_vs_exp,
  width = 2.5,
  height = 2.2
)
ggsave(
  filename = file.path(figure_dir, "SLAMIT_chrX_kdegs.pdf"),
  plot = chrX_kdegs,
  width = 2.5,
  height = 2.2
)
ggsave(
  filename = file.path(figure_dir, "SLAMIT_celltype_diffs.pdf"),
  plot = celltype_diffs,
  width = 2.5,
  height = 2.2
)
