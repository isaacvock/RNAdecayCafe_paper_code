### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 3:
## - Violin/boxplot comparing distributions of sample-to-sample half-life
##   estimate Pearson correlations across the 4 datasets shown in Figure 1A-D:
##   1) Agarwal et al., 2022 non-NR-seq datasets ("Agarwal")
##   2) Pulse-chase NR-seq datasets ("PC")
##   3) Pulse-label NR-seq datasets ("PL")
##   4) RNAdecayCafe filtered pulse-label NR-seq datasets ("RNAdecayCafe")


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)


# File paths -------------------------------------------------------------------

path_to_agarwal_data     <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Saluki_data/processed_half_lives_TableS2.xlsx"
path_to_all_nrseq        <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/All_NRseq_kdegs_refined.csv"
path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"
figure_dir               <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Supplemental_Panels_v2"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# 1) Agarwal et al., 2022 correlations -----------------------------------------

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


# 2) Pulse-chase NR-seq correlations -------------------------------------------

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


# 3) Pulse-label NR-seq correlations -------------------------------------------

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


# 4) RNAdecayCafe correlations -------------------------------------------------

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


# Combine correlation distributions & plot -------------------------------------

cor_tbl <- tibble(
  corrs = c(abs(unique(as.vector(new_corrmat))),
            unique(as.vector(pc_corrmat)),
            unique(as.vector(pl_corrmat)),
            unique(as.vector(plf_corrmat))),
  source = factor(c(
    rep("Agarwal", times = length(unique(as.vector(new_corrmat)))),
    rep("PC", times = length(unique(as.vector(pc_corrmat)))),
    rep("PL", times = length(unique(as.vector(pl_corrmat)))),
    rep("RDC", times = length(unique(as.vector(plf_corrmat))))
  ),
  levels = c("Agarwal", "PC", "PL", "RDC"))
) %>%
  dplyr::mutate(
    corrs = ifelse(
      corrs < 0,
      0,
      corrs
    )
  ) %>%
  dplyr::filter(corrs < 1)

cor_comp <- cor_tbl %>%
  ggplot(
    aes(x = source,
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
  xlab("Data source") +
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
  plot = cor_comp,
  filename = file.path(figure_dir, "Correlations.pdf"),
  width = 2.25,
  height = 2.1
)
