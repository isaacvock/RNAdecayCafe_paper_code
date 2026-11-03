### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 4 panels (A-B):
## - Supplemental Figure 4A: Sample-to-sample log(half-life) Pearson correlation
##                           heatmap for RNAdecayCafe with right annotations for
##                           cell line, label time, and inclusion of a no-s4U
##                           control (RDC_correlations_with_ctl_color.pdf).
## - Supplemental Figure 4B: Dodged notched boxplot comparing sample-to-sample
##                           Pearson correlations for short (< 2hr) and long
##                           (> 4hr) label times with and without no-s4U controls
##                           (Ctl_vs_noctl_pearsons.pdf).


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

path_to_rnadecaycafe_rds <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/RNAdegDB_filtered_refined.rds"
path_to_ctl_lookup       <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Data/Analyses/Revisions/nos4U_ctl_lookup.csv"
figure_dir               <- "C:/Users/isaac/Main_Docs/Simon_Lab/kdeg_curation/Manuscript/Figures/Revision_Supplemental_Panels"

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

logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) exp(x) / (1 + exp(x))


# Load or create no-s4U control lookup table -----------------------------------

if (file.exists(path_to_ctl_lookup)) {
  ctl_lookup <- read_csv(path_to_ctl_lookup)
} else {
  ctl_lookup <- tibble(
    sample = c(
      "SRR10887642", "SRR10887643", "SRR10887644", "SRR10887645",
      "SRR10887648", "SRR10887649", "SRR10887650", "SRR10887651",
      "SRR11067765", "SRR11067766", "SRR11067767",
      "SRR11142916", "SRR11142918", "SRR11142920", "SRR11142922",
      "SRR11294627", "SRR11294628", "SRR11294629", "SRR11294630",
      "SRR11294635", "SRR11294636", "SRR11294637", "SRR11294638",
      "SRR11458917", "SRR11458918", "SRR11458919", "SRR11458920",
      "SRR11458921", "SRR11458922",
      "SRR12168632", "SRR12168633", "SRR12168634", "SRR12168635",
      "SRR12168636", "SRR12168637",
      "SRR12539883", "SRR12539884", "SRR12539885", "SRR12539886",
      "SRR12539887", "SRR12539888", "SRR12539890", "SRR12539891",
      "SRR12539892", "SRR12539895", "SRR12539896",
      "SRR13159432", "SRR13159433", "SRR13159434", "SRR13159438",
      "SRR13159439", "SRR13159440", "SRR13159444", "SRR13159445",
      "SRR13159446",
      "SRR13573505", "SRR13573506", "SRR13573507", "SRR13573508",
      "SRR13573509", "SRR13573510", "SRR13573511", "SRR13573512",
      "SRR13573513", "SRR13573514", "SRR13573515", "SRR13573516",
      "SRR13573529", "SRR13573530", "SRR13573531", "SRR13573532",
      "SRR13573533", "SRR13573534", "SRR13573535", "SRR13573536",
      "SRR13573537", "SRR13573538", "SRR13573539", "SRR13573540",
      "SRR14475305", "SRR14475306", "SRR14475307", "SRR14475308",
      "SRR14475311", "SRR14475312", "SRR14475313", "SRR14475314",
      "SRR15654895", "SRR15654896", "SRR15654897", "SRR15655142",
      "SRR17803609", "SRR17803610", "SRR17803611", "SRR17803612",
      "SRR18583413", "SRR18583414", "SRR18583415", "SRR18583420",
      "SRR18583421", "SRR18583422", "SRR18583423", "SRR18583424",
      "SRR18583425", "SRR18583431", "SRR18583432", "SRR18583433",
      "SRR18583434", "SRR18583435", "SRR18583436", "SRR18583445",
      "SRR18583446", "SRR18583447", "SRR18583448", "SRR18583449",
      "SRR18583450", "SRR18583455", "SRR18583456", "SRR18583457",
      "SRR18583458", "SRR18583459", "SRR18583460",
      "SRR19165056", "SRR19165057", "SRR19165058",
      "SRR19981794", "SRR19981795",
      "SRR20080260", "SRR20080261", "SRR20080283", "SRR20080284",
      "SRR29164313", "SRR29164314", "SRR29164315", "SRR29164316",
      "SRR29164317", "SRR29164318", "SRR29164319", "SRR29164320",
      "SRR29164321", "SRR29164322", "SRR29164323", "SRR29164324",
      "SRR30007013", "SRR30007014", "SRR30007015", "SRR30007019",
      "SRR30007020", "SRR30007021",
      "SRR30252738", "SRR30252739", "SRR30252740", "SRR30252741",
      "SRR30252742", "SRR30252743", "SRR30252744", "SRR30252745",
      "SRR30252746",
      "SRR31699594", "SRR31699597", "SRR31699605",
      "SRR5806774", "SRR5806775", "SRR5806776", "SRR5806780",
      "SRR5806781", "SRR5806782", "SRR5806786", "SRR5806787",
      "SRR5806788", "SRR5806809", "SRR5806810",
      "SRR6262128", "SRR6262129",
      "SRR7993823", "SRR7993824", "SRR7993825", "SRR7993826",
      "SRR7993827", "SRR7993828",
      "SRR8244477", "SRR8244478", "SRR8244479", "SRR8244480",
      "SRR8244481", "SRR8244482",
      "SRR8571241", "SRR8571242", "SRR8571243", "SRR8571244",
      "SRR8571245", "SRR8571246", "SRR8571247", "SRR8571248",
      "SRR8571249", "SRR8571250", "SRR8571251", "SRR8571252",
      "SRR8571253", "SRR8571254", "SRR8571255", "SRR8571256",
      "SRR8571257", "SRR8571258",
      "SRR10316942", "SRR10316943", "SRR10316944", "SRR10316945",
      "SRR10316946", "SRR10316947", "SRR10316948",
      "DMSO_1", "DMSO_2", "DMSO_3"
    ),
    includes_nos4U = c(
      rep(TRUE, 8),
      rep(TRUE, 3),
      rep(TRUE, 4),
      rep(FALSE, 8),
      rep(FALSE, 6),
      rep(TRUE, 6),
      rep(TRUE, 11),
      rep(TRUE, 9),
      rep(TRUE, 24),
      rep(FALSE, 8),
      rep(TRUE, 4),
      rep(FALSE, 4),
      rep(FALSE, 26),
      rep(FALSE, 3),
      rep(TRUE, 2),
      rep(TRUE, 4),
      rep(FALSE, 12),
      rep(FALSE, 6),
      rep(FALSE, 9),
      rep(TRUE, 3),
      rep(TRUE, 11),
      rep(TRUE, 2),
      rep(TRUE, 6),
      rep(FALSE, 6),
      rep(TRUE, 18),
      rep(TRUE, 7),
      rep(TRUE, 3)
    )
  )
}


# Supplemental Figure 4A: RNAdecayCafe heatmap with no-s4U annotation ----------

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
  ) %>%
  dplyr::inner_join(
    ctl_lookup,
    by = "sample"
  )

ann_cols <- c(
  "< 2hr" = "deepskyblue",
  "[2, 4]hr" = "darkgray",
  "> 4hr" = "darkred"
)

ctl_cols <- c(
  "TRUE" = "seagreen3",
  "FALSE" = "sienna4"
)

label_times <- label_times[match(rownames(plf_corrmat), label_times$sample), ]

row_ha <- rowAnnotation(
  cell = label_times$cell_line,
  tl = label_times$label_time,
  ctl = label_times$includes_nos4U,
  col = list(
    tl = ann_cols,
    cell = cell_cols,
    ctl = ctl_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_rdb_col <- ComplexHeatmap::Heatmap(
  relu(plf_corrmat),
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
  file.path(figure_dir, "RDC_correlations_with_ctl_color.pdf"),
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
  ctl = label_times$includes_nos4U,
  col = list(
    tl = ann_cols,
    cell = cell_cols,
    ctl = ctl_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_rdb_col_leg <- ComplexHeatmap::Heatmap(
  relu(plf_corrmat),
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
  file.path(figure_dir, "RDC_correlations_with_ctl_color_legend.pdf"),
  width = 6.4,
  height = 6.1,
  useDingbats = FALSE
)
draw(ht_rdb_col_leg)
dev.off()


# Supplemental Figure 4B: Boxplot by label time and no-s4U control -------------

corr_df <- tibble(
  row_name = rep(rownames(plf_corrmat), each = ncol(plf_corrmat)),
  col_name = rep(colnames(plf_corrmat), times = nrow(plf_corrmat)),
  correlation = as.vector(t(plf_corrmat))
) %>%
  dplyr::inner_join(
    label_times %>%
      dplyr::select(sample, label_time, includes_nos4U, cell_line),
    by = c("row_name" = "sample")
  ) %>%
  dplyr::filter(
    row_name != col_name
  ) %>%
  dplyr::filter(
    label_time %in% c("< 2hr", "> 4hr")
  )

ctl_pearsons <- corr_df %>%
  ggplot(
    aes(
      x = label_time,
      y = inv_logit(logit(correlation)),
      fill = includes_nos4U
    )
  ) +
  geom_boxplot(
    outliers = FALSE,
    notch = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "seagreen3",
      "FALSE" = "sienna4"
    )
  ) +
  theme_classic() +
  xlab("Label time") +
  ylab("Pearson r") +
  coord_cartesian(
    ylim = c(0.4, 1)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  filename = file.path(figure_dir, "Ctl_vs_noctl_pearsons.pdf"),
  plot = ctl_pearsons,
  width = 2.5,
  height = 2.1
)
