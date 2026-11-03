### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 5:
## - Heatmap of LASSO regression coefficients fit to z-score normalized
##   stability estimates in each cell type in RNAdecayCafe or PC1 stability
##   estimates from Agarwal et al., 2022 across 9 transcript sequence features
##   (LASSO_correlations.pdf).


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)
library(glmnet)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)


# File paths -------------------------------------------------------------------

path_to_avg_kdegs   <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_averages_filtered_refined.csv"
path_to_saluki_xlsx <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/Saluki_kdegs.xlsx"
mane_gtf_path       <- "C:/Users/isaac/Main_Docs/Kundaje_Lab/Data/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz"
path_to_seq_and_kd  <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/MANE_sequence_and_kdegs.csv"
figure_dir          <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels_v2"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Steps 1-5: Extract MANE sequence features & combine with kdeg estimates ------
# Note: If MANE.GRCh38.v1.4.refseq_genomic.gtf.gz is available, regenerates
# normalized_df directly from the uncorrupted averages table and GTF; otherwise
# loads MANE_sequence_and_kdegs.csv (filtering to valid cell_line rows).

valid_cell_lines <- c(
  "Agarwal", "Calu3", "CH22", "HEK293", "HeLa", "K562",
  "MCF7", "MDAMB", "MOLM", "MV411", "Nalm6", "U2OS"
)

if (file.exists(mane_gtf_path)) {

  Avg_kdegs <- read_csv(path_to_avg_kdegs) %>%
    dplyr::select(
      gene_id = feature_ID,
      cell_line,
      avg_donorm_log_kdeg
    ) %>%
    dplyr::group_by(cell_line) %>%
    dplyr::mutate(
      log_kdeg_z = (avg_donorm_log_kdeg - mean(avg_donorm_log_kdeg)) / sd(avg_donorm_log_kdeg)
    )

  Saluki_kdegs <- read_xlsx(
    path_to_saluki_xlsx,
    skip = 2
  ) %>%
    dplyr::select(
      `Gene name`,
      `half-life (PC1)`
    ) %>%
    dplyr::rename(
      gene_id = `Gene name`,
      thalf = `half-life (PC1)`
    ) %>%
    dplyr::mutate(
      log_kdeg_z = -thalf / sd(thalf)
    )

  MANE_gtf <- rtracklayer::import(mane_gtf_path)

  filter <- MANE_gtf$tag != "MANE Plus Clinical"
  filter <- ifelse(
    is.na(filter),
    TRUE,
    filter
  )

  MANE_gtf <- MANE_gtf[filter]
  txdb <- txdbmaker::makeTxDbFromGRanges(MANE_gtf)

  ### Number of exons
  num_exon_df <- MANE_gtf %>%
    dplyr::as_tibble() %>%
    dplyr::filter(
      type == "exon"
    ) %>%
    dplyr::count(gene_id) %>%
    dplyr::rename(
      exon_count = n
    )

  ### 0) gene ID-to-transcript ID mapping
  gene_to_tx <- MANE_gtf %>%
    dplyr::as_tibble() %>%
    dplyr::filter(
      type == "exon"
    ) %>%
    dplyr::select(
      gene_id, transcript_id
    ) %>%
    dplyr::distinct()

  ### 1) 5' UTR
  fiveprimeutr <- fiveUTRsByTranscript(
    txdb,
    use.names = TRUE
  )

  fiveprimeutr_seq <- getSeq(
    Hsapiens,
    fiveprimeutr
  ) %>%
    unlist()

  fiveprimeutr_str <- as.character(fiveprimeutr_seq)

  fiveprimeutr_df <- tibble(
    seq = fiveprimeutr_str,
    transcript_id = names(fiveprimeutr_str)
  ) %>%
    group_by(transcript_id) %>%
    summarise(seq = paste(seq, collapse = ""))

  ### 2) CDS
  cds <- cdsBy(
    txdb,
    by = "tx",
    use.names = TRUE
  )

  cds_seq <- getSeq(
    Hsapiens,
    cds
  ) %>%
    unlist()

  cds_str <- as.character(cds_seq)

  cds_df <- tibble(
    seq = cds_str,
    transcript_id = names(cds_str)
  ) %>%
    group_by(transcript_id) %>%
    summarise(seq = paste(seq, collapse = ""))

  ### 3) 3' UTR
  threeprimeutr <- threeUTRsByTranscript(
    txdb,
    use.names = TRUE
  )

  threeprimeutr_seq <- getSeq(
    Hsapiens,
    unlist(threeprimeutr)
  )

  threeprimeutr_str <- as.character(threeprimeutr_seq)

  threeprimeutr_df <- tibble(
    seq = threeprimeutr_str,
    transcript_id = names(threeprimeutr_str)
  ) %>%
    group_by(transcript_id) %>%
    summarise(seq = paste(seq, collapse = ""))

  ### 4) Combine
  combined_df <- threeprimeutr_df %>%
    dplyr::rename(
      threepseq = seq
    ) %>%
    dplyr::inner_join(
      gene_to_tx,
      by = "transcript_id"
    ) %>%
    dplyr::inner_join(
      fiveprimeutr_df %>%
        dplyr::rename(
          fivepseq = seq
        ),
      by = "transcript_id"
    ) %>%
    dplyr::inner_join(
      cds_df %>%
        dplyr::rename(
          CDSseq = seq
        ),
      by = "transcript_id"
    ) %>%
    dplyr::select(
      gene_id, everything()
    )

  ### 5) Compute on these features
  add_3utr_features <- function(df) {
    dna <- DNAStringSet(toupper(df$threepseq))
    utr_len <- width(dna)
    drach_hits <- vmatchPattern("DRACH", dna, fixed = FALSE) |> elementNROWS()
    drach_per_kb <- ifelse(utr_len > 0, drach_hits / (utr_len / 1000), NA_real_)

    df %>%
      mutate(
        threep_drach_norm = drach_per_kb,
        threep_len = utr_len
      )
  }

  combined_annot <- add_3utr_features(
    combined_df
  ) %>%
    dplyr::mutate(
      fivep_len = str_length(fivepseq),
      cds_len = str_length(CDSseq),
      threep_gc_content = ((str_count(threepseq, "G") + str_count(threepseq, "C")) + 0.5) / (threep_len + 1),
      fivep_gc_content = ((str_count(fivepseq, "G") + str_count(fivepseq, "C")) + 0.5) / (fivep_len + 1),
      cds_gc_content = ((str_count(CDSseq, "G") + str_count(CDSseq, "C") + 0.5) / (cds_len + 1))
    )

  logit <- function(x) log(x / (1 - x))

  normalized_df <- combined_annot %>%
    dplyr::mutate(
      threep_gc_z = (logit(threep_gc_content) - mean(logit(threep_gc_content))) / sd(logit(threep_gc_content)),
      fivep_gc_z = (logit(fivep_gc_content) - mean(logit(fivep_gc_content))) / sd(logit(fivep_gc_content)),
      cds_gc_z = (logit(cds_gc_content) - mean(logit(cds_gc_content))) / sd(logit(cds_gc_content)),
      threep_drach_z = (log(threep_drach_norm + 1) - mean(log(threep_drach_norm + 1))) / sd(log(threep_drach_norm + 1)),
      threep_len_z = (log(threep_len) - mean(log(threep_len))) / sd(log(threep_len)),
      fivep_len_z = (log(fivep_len) - mean(log(fivep_len))) / sd(log(fivep_len)),
      cds_len_z = (log(cds_len) - mean(log(cds_len))) / sd(log(cds_len))
    ) %>%
    dplyr::inner_join(
      Avg_kdegs %>%
        dplyr::select(
          gene_id, cell_line, log_kdeg_z
        ) %>%
        dplyr::bind_rows(
          Saluki_kdegs %>%
            dplyr::select(
              gene_id, log_kdeg_z
            ) %>%
            dplyr::mutate(
              cell_line = "Agarwal"
            )
        ),
      by = "gene_id"
    ) %>%
    dplyr::inner_join(
      num_exon_df,
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      exon_cnt_z = (log(exon_count) - mean(log(exon_count))) / sd(log(exon_count))
    ) %>%
    dplyr::select(
      cell_line,
      gene_id,
      log_kdeg_z,
      everything()
    )

} else {
  normalized_df <- read_csv(path_to_seq_and_kd) %>%
    dplyr::filter(cell_line %in% valid_cell_lines) %>%
    dplyr::filter(!is.na(CDSseq) & !is.na(threepseq) & !is.na(fivepseq))
}


# Step 6: LASSO regression for all cell lines & plot heatmap -------------------

normalized_df <- normalized_df %>%
  dplyr::mutate(
    cds_drach_cnt = vmatchPattern(
      "DRACH", DNAStringSet(toupper(CDSseq)), fixed = FALSE
    ) %>% elementNROWS(),
    length = cds_len + fivep_len + threep_len,
    exon_density = exon_count / cds_len
  ) %>%
  dplyr::mutate(
    cds_drach_norm = cds_drach_cnt / (cds_len / 1000),
    cds_drach_z = (log(cds_drach_norm + 1) - mean(log(cds_drach_norm + 1), na.rm = TRUE)) / sd(log(cds_drach_norm + 1), na.rm = TRUE),
    length_z = (log(length) - mean(log(length), na.rm = TRUE)) / sd(log(length), na.rm = TRUE),
    exon_dens_z = (log(exon_density) - mean(log(exon_density), na.rm = TRUE)) / sd(log(exon_density), na.rm = TRUE)
  )

cell_lines <- unique(normalized_df$cell_line)
all_coefs <- tibble()

for (c in seq_along(cell_lines)) {

  sub_df <- normalized_df %>%
    dplyr::filter(
      cell_line == cell_lines[c]
    ) %>%
    dplyr::select(
      log_kdeg_z,
      threep_gc_z,
      fivep_gc_z,
      cds_gc_z,
      threep_drach_z,
      threep_len_z,
      fivep_len_z,
      cds_len_z,
      exon_dens_z,
      cds_drach_z
    ) %>%
    na.omit()

  X <- sub_df %>%
    dplyr::select(-log_kdeg_z) %>%
    as.matrix()

  y <- sub_df %>%
    dplyr::select(log_kdeg_z) %>%
    unlist() %>%
    unname()

  # Cross validation to tune hyperparameter
  cv_lasso <- cv.glmnet(
    X,
    y,
    alpha = 1
  )

  best_lambda <- cv_lasso$lambda.min

  lasso_fit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = best_lambda
  )

  coefs <- tibble(
    name = coef(lasso_fit) %>% rownames(),
    coef = coef(lasso_fit)[, 1] %>% unname(),
    cell_line = cell_lines[c]
  ) %>%
    dplyr::filter(
      name != "(Intercept)"
    )

  all_coefs <- bind_rows(
    all_coefs,
    coefs
  )
}

LASSO_coefs <- all_coefs %>%
  ggplot(aes(x = name, y = cell_line, fill = coef)) +
  geom_tile() +
  scale_fill_gradient2(
    name = "coef",
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(figure_dir, "LASSO_correlations.pdf"),
  plot = LASSO_coefs,
  width = 5.5,
  height = 4.1
)
