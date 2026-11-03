### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 9 panels (A-B):
## - Supplemental Figure 9A: Gene Ontology (GO) enrichment dot plot for the 99
##                           differentially expressed genes meeting both the R^2
##                           (> 0.25) and slope (min(m, |1/m|) > 0.25) criteria
##                           in 3'-end and full-length datasets from Figure 2D
##                           (Figure2CD_GO_analysis.pdf).
## - Supplemental Figure 9B: Generation of target and background 3'-UTR FASTA
##                           files for HOMER motif discovery, and bar plot of
##                           the percentage of target vs. background sequences
##                           containing the top 4 enriched 3'-UTR motifs
##                           (HOMER_motif_barplot.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(Biostrings)
library(clusterProfiler)
library(org.Hs.eg.db)


# File paths -------------------------------------------------------------------

path_to_fig2_celltype_data <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/RDB_celltype_specificity.csv"
path_to_seq_and_kd         <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/MANE_sequence_and_kdegs.csv"
path_to_hits_fasta         <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/MANE_threepUTRs_posttxn_hits.fasta"
path_to_bg_fasta           <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/MANE_threepUTRs_background_DE.fasta"
figure_dir                 <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels_v2/"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Supplemental Figure 9A: GO enrichment dot plot -------------------------------

celltype_kdeg <- read_csv(path_to_fig2_celltype_data)

hit_genes <- celltype_kdeg %>%
  dplyr::filter(
    both
  ) %>%
  dplyr::select(
    XF
  ) %>%
  unlist() %>%
  unname() %>%
  unique()

universe_of_genes <- celltype_kdeg$XF %>%
  unique()

ego <- enrichGO(
  gene = hit_genes,
  universe = universe_of_genes,
  keyType = "SYMBOL",
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

# Pick top terms by FDR
top_n <- 10
df <- as_tibble(ego) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  mutate(Description = str_wrap(Description, width = 45)) %>%
  dplyr::mutate(
    Description = factor(Description,
                         levels = Description[order(-p.adjust)])
  )

p_dot <- ggplot(df, aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(
    low = "darkgray",
    high = "darkred"
  ) +
  labs(
    x = "Fold Enrichment",
    y = NULL,
    color = "-log10(FDR)",
    size = "Gene count"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    plot.margin = margin(5.5, 25, 5.5, 5.5)
  ) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  plot = p_dot,
  filename = file.path(figure_dir, "Figure2CD_GO_analysis.pdf"),
  width = 8,
  height = 4
)


# Supplemental Figure 9B: HOMER 3'-UTR FASTA export & motif bar plot -----------

### 1) Generate target and background 3'-UTR FASTAs for HOMER
sequence_info <- read_csv(path_to_seq_and_kd)

regulated_genes <- celltype_kdeg %>%
  dplyr::filter(
    both
  ) %>%
  dplyr::select(
    XF
  ) %>%
  unlist() %>%
  unname() %>%
  unique()

gene_universe <- sequence_info %>%
  dplyr::filter(
    cell_line %in% c("K562", "MV411", "MOLM")
  ) %>%
  dplyr::count(
    gene_id
  ) %>%
  dplyr::filter(
    n == 3
  ) %>%
  dplyr::select(
    gene_id
  ) %>%
  unlist() %>%
  unname()

hits <- sequence_info %>%
  dplyr::filter(
    gene_id %in% regulated_genes &
      !is.na(threepseq)
  ) %>%
  dplyr::select(
    gene_id,
    threepseq
  ) %>%
  dplyr::distinct()

background <- sequence_info %>%
  dplyr::filter(
    gene_id %in% gene_universe &
      !(gene_id %in% regulated_genes) &
      !is.na(threepseq)
  ) %>%
  dplyr::select(
    gene_id,
    threepseq
  ) %>%
  dplyr::distinct()

background_DE <- sequence_info %>%
  dplyr::filter(
    gene_id %in% unique(celltype_kdeg$XF) &
      !(gene_id %in% regulated_genes) &
      !is.na(threepseq)
  ) %>%
  dplyr::select(
    gene_id,
    threepseq
  ) %>%
  dplyr::distinct()

dna_hits <- DNAStringSet(hits$threepseq)
names(dna_hits) <- hits$gene_id

dna_bg_DE <- DNAStringSet(background_DE$threepseq)
names(dna_bg_DE) <- background_DE$gene_id

writeXStringSet(
  dna_hits,
  filepath = path_to_hits_fasta
)

writeXStringSet(
  dna_bg_DE,
  filepath = path_to_bg_fasta
)

### 2) Plot HOMER 3'-UTR motif enrichment results
motif_df <- tibble(
  motif = rep(c(
    "m1",
    "m2",
    "m3",
    "m4"
  ),
  each = 2),
  perc = c(
    5.05,
    0.00,
    13.13,
    1.73,
    4.04,
    0.04,
    6.06,
    0.29
  ),
  source = rep(
    c(
      "target",
      "background"
    ),
    times = 4
  )
)

motif_plot <- motif_df %>%
  ggplot(
    aes(
      x = motif,
      y = perc,
      fill = source
    )
  ) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      background = "darkgray",
      target = "darkred"
    )
  ) +
  theme_classic() +
  xlab("Motif") +
  ylab("% of sequences with motif") +
  labs(fill = "sequence type") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(
  plot = motif_plot,
  filename = file.path(figure_dir, "HOMER_motif_barplot.pdf"),
  width = 7,
  height = 5
)
