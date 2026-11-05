### PURPOSE OF THIS SCRIPT
## Code to make final panels for Figure 1


# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


##### Helper functions #####

source(
  "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Scripts/helper_functions.R"
)




##### Color palettes #####

# Cell annotation
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



# Remake  plots with ComplexHeatmap --------------------------------------------


##### Saluki database #####

saluki_hls <- readxl::read_excel(
  'C:/Users/isaac/Documents/Simon_Lab/Sandbox/Kdeg_estimation/13059_2022_2811_MOESM2_ESM.xlsx',
  sheet = 1)

M_sal <- saluki_hls %>%
  dplyr::select(-`Ensembl Gene Id`, -`Gene name`) %>%
  as.matrix()

class(M_sal) <- "numeric"

new_corrmat <- cor(M_sal,
                   method = "pearson",
                   use = "pairwise.complete.obs")


col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),         
  rev(brewer.pal(11, "RdYlBu"))         
)

ht_sal <- ComplexHeatmap::Heatmap(abs(new_corrmat),
                                  col               = col_fun,
                                  cluster_rows      = TRUE,
                                  cluster_columns   = TRUE,
                                  show_row_names    = FALSE,
                                  show_column_names = FALSE,
                                  show_row_dend = FALSE,
                                  show_column_dend = FALSE,
                                  row_names_gp      = gpar(fontsize = 6),
                                  column_names_gp   = gpar(fontsize = 6),
                                  border            = FALSE,            # border around the whole heatmap
                                  heatmap_legend_param = list(
                                    title  = "|r|",
                                    at     = seq(0, 1, by = 0.2),
                                    labels = seq(0, 1, by = 0.2)
                                  ),
                                  show_heatmap_legend = FALSE
)


pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/Saluki_database_heatmap.pdf",
    width  = 2.2,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_sal)                  
dev.off()


### Save one with legend
ht_sal2 <- ComplexHeatmap::Heatmap(abs(new_corrmat),
                                   col               = col_fun,
                                   cluster_rows      = TRUE,
                                   cluster_columns   = TRUE,
                                   show_row_names    = FALSE,
                                   show_column_names = FALSE,
                                   show_row_dend = FALSE,
                                   show_column_dend = FALSE,
                                   row_names_gp      = gpar(fontsize = 6),
                                   column_names_gp   = gpar(fontsize = 6),
                                   border            = FALSE,            # border around the whole heatmap
                                   heatmap_legend_param = list(
                                     title  = "|r|",
                                     at     = seq(0, 1, by = 0.2),
                                     labels = seq(0, 1, by = 0.2)
                                   ),
                                   show_heatmap_legend = TRUE
)

pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/Heatmap_with_legend.pdf",
    width  = 3,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_sal2)                  
dev.off()

### Get decent samples

ht_sal <- draw(ht_sal)
ord_idx_sal <- row_order(ht_sal)
sal_col_order <- colnames(new_corrmat)[ord_idx_sal]

bad_start <- min(which(grepl("^Gejman_", sal_col_order)))

ok_samples <- sal_col_order[1:(bad_start - 1)]


ord_idx_sal_nonrseq <- ord_idx_sal[which(!grepl("^Simon_", sal_col_order) & !grepl("^Shendure", sal_col_order))]


### Redraw with NR-seq samples filtered out and good concordance subset used later noted

good_tbl <- tibble(
  sample = rownames(new_corrmat)[ord_idx_sal_nonrseq]
) %>%
  dplyr::mutate(
    good = factor(as.character(sample %in% ok_samples)),
    method = sapply(str_split(
      sample, "_"
    ),
    function(x) x[2]
    ),
    cell = sapply(str_split(
      sample, "_"
    ),
    function(x) x[3]
    ),
    dataset = sapply(str_split(
      sample, "_"
    ),
    function(x) x[1]
    )
  ) %>%
  dplyr::mutate(
    method = case_when(
      method %in% c("ActD", "Aman") ~ "txni",
      method %in% c("4sU", "BrU", "BrU4sU") ~ "enrich"
    ),
    cell = case_when(
      grepl("^GM", cell) ~ "LCL",
      .default = cell
    )
  )

good_cols <- c("TRUE" = "#d95f02", 
               "FALSE" = "darkgray")



row_ha <- rowAnnotation(
  good = good_tbl$good,
  col = list(good = good_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

set.seed(42)
row_ha_right <- rowAnnotation(
  cell = good_tbl$cell,
  method = good_tbl$method,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)


ht_sal_annot <- ComplexHeatmap::Heatmap(abs(new_corrmat)[ord_idx_sal_nonrseq, ord_idx_sal_nonrseq],
                                        col               = col_fun,
                                        cluster_rows      = FALSE,
                                        cluster_columns   = FALSE,
                                        show_row_names    = FALSE,
                                        show_column_names = FALSE,
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE,
                                        right_annotation = row_ha_right,
                                        row_names_gp      = gpar(fontsize = 6),
                                        column_names_gp   = gpar(fontsize = 6),
                                        border            = FALSE,            # border around the whole heatmap
                                        heatmap_legend_param = list(
                                          title  = "|r|",
                                          at     = seq(0, 1, by = 0.2),
                                          labels = seq(0, 1, by = 0.2)
                                        ),
                                        show_heatmap_legend = FALSE
)


draw(ht_sal_annot)

pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/Saluki_with_annotation.pdf",
    width  = 2.5,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_sal_annot)                  
dev.off()



### With legends

set.seed(42)
row_ha_right <- rowAnnotation(
  cell = good_tbl$cell,
  method = good_tbl$method,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)


ht_sal_annot <- ComplexHeatmap::Heatmap(abs(new_corrmat)[ord_idx_sal_nonrseq, ord_idx_sal_nonrseq],
                                        col               = col_fun,
                                        cluster_rows      = FALSE,
                                        cluster_columns   = FALSE,
                                        show_row_names    = FALSE,
                                        show_column_names = FALSE,
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE,
                                        left_annotation = row_ha,
                                        right_annotation = row_ha_right,
                                        row_names_gp      = gpar(fontsize = 6),
                                        column_names_gp   = gpar(fontsize = 6),
                                        border            = FALSE,            # border around the whole heatmap
                                        heatmap_legend_param = list(
                                          title  = "|r|",
                                          at     = seq(0, 1, by = 0.2),
                                          labels = seq(0, 1, by = 0.2)
                                        ),
                                        show_heatmap_legend = FALSE
)


draw(ht_sal_annot)

pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/Saluki_with_annotation_legends2.pdf",
    width  = 3,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_sal_annot)                  
dev.off()



##### Pulse-chase #####
all_kinetics <- read_csv("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/Analyses/All_NRseq_kdegs_refined.csv")


pulse_chase_ests <- all_kinetics %>%
  dplyr::filter(pulse_chase)

pc_wide <- pulse_chase_ests %>%
  dplyr::select(
    sample, XF, log_kdeg
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "log_kdeg"
  )


M_pc <- pc_wide %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_pc) <- "numeric"

pc_corrmat <- cor(M_pc,
                  method = "pearson",
                  use = "pairwise.complete.obs")


col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),         
  rev(brewer.pal(11, "RdYlBu"))         
)



annot_tbl <- pulse_chase_ests %>%
  dplyr::select(
    sample,
    dataset
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell = sapply(
      str_split(dataset,
                "_"),
      function(x) x[4]
    )
  ) %>%
  dplyr::mutate(
    cell = ifelse(
      cell == "HEK293T",
      "HEK293",
      cell
    )
  )


row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_pc <- ComplexHeatmap::Heatmap(abs(pc_corrmat),
                                 col               = col_fun,
                                 cluster_rows      = TRUE,
                                 cluster_columns   = TRUE,
                                 show_row_names    = FALSE,
                                 show_column_names = FALSE,
                                 show_row_dend = FALSE,
                                 show_column_dend = FALSE,
                                 right_annotation = row_ha,
                                 row_names_gp      = gpar(fontsize = 6),
                                 column_names_gp   = gpar(fontsize = 6),
                                 border            = FALSE,            # border around the whole heatmap
                                 heatmap_legend_param = list(
                                   title  = "|r|",
                                   at     = seq(0, 1, by = 0.2),
                                   labels = seq(0, 1, by = 0.2)
                                 ),
                                 show_heatmap_legend = FALSE
)

pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/NRseq_pulsechase_heatmap_annotate.pdf",
    width  = 2.3,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_pc)                  
dev.off()

##### Pulse-label data #####


all_kinetics <- read_csv("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/Analyses/All_NRseq_kdegs_refined.csv")

pulse_label_ests <- all_kinetics %>%
  dplyr::filter(!pulse_chase & 
                  tchase >= 1 & # No crazy short label times
                  pnew > 0.01 # s4U needs to have gotten incorporated
  )


pl_wide <- pulse_label_ests %>%
  dplyr::select(
    sample, XF, log_kdeg
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "log_kdeg"
  ) 


M_pl <- pl_wide %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_pl) <- "numeric"

pl_corrmat <- cor(M_pl,
                  method = "pearson",
                  use = "pairwise.complete.obs")


col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),         
  rev(brewer.pal(11, "RdYlBu"))         
)

get_muhar_cell <- function(samples, cells){
  
  
  cell_tbl <- tibble(
    sample = c("SRR5806774", "SRR5806775", "SRR5806776",
               "SRR5806780", "SRR5806781", "SRR5806782",
               "SRR5806786", "SRR5806787", "SRR5806788",
               "SRR5806809", "SRR5806810"),
    tl = 1,
    cell = c("K562", "K562", "K562",
             "MV4-11", "MV4-11", "MV4-11",
             "MOLM-13", "MOLM-13", "MOLM-13",
             "OCI/AML-3", "OCI/AML-3")
  )
  
  this_cell <- rep("", times = length(samples))
  for(s in seq_along(samples)){
    
    if(cells[s] == "many"){
      this_samp <- samples[s]
      
      this_cell[s] <- cell_tbl$cell[cell_tbl$sample == this_samp]
      
    }else{
      this_cell[s] <- cells[s]
    }
    
  }
  
  return(this_cell)
  
  
}

annot_tbl <- pulse_label_ests %>%
  dplyr::select(
    sample,
    dataset
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell = sapply(
      str_split(dataset,
                "_"),
      function(x) x[4]
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


row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht_pl <- ComplexHeatmap::Heatmap(abs(pl_corrmat),
                                 col               = col_fun,
                                 cluster_rows      = TRUE,
                                 cluster_columns   = TRUE,
                                 show_row_names    = FALSE,
                                 show_column_names = FALSE,
                                 show_row_dend = FALSE,
                                 show_column_dend = FALSE,
                                 right_annotation = row_ha,
                                 row_names_gp      = gpar(fontsize = 6),
                                 column_names_gp   = gpar(fontsize = 6),
                                 border            = FALSE,            # border around the whole heatmap
                                 heatmap_legend_param = list(
                                   title  = "|r|",
                                   at     = seq(0, 1, by = 0.2),
                                   labels = seq(0, 1, by = 0.2)
                                 ),
                                 show_heatmap_legend = FALSE
)


pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/NRseq_pulselabel_heatmap_annot.pdf",
    width  = 2.3,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_pl)                  
dev.off()


### Get one with legend

row_ha <- rowAnnotation(
  cell = annot_tbl$cell,
  col = list(
    cell = cell_cols
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

ht_pl <- ComplexHeatmap::Heatmap(abs(pl_corrmat),
                                 col               = col_fun,
                                 cluster_rows      = TRUE,
                                 cluster_columns   = TRUE,
                                 show_row_names    = FALSE,
                                 show_column_names = FALSE,
                                 show_row_dend = FALSE,
                                 show_column_dend = FALSE,
                                 right_annotation = row_ha,
                                 row_names_gp      = gpar(fontsize = 6),
                                 column_names_gp   = gpar(fontsize = 6),
                                 border            = FALSE,            # border around the whole heatmap
                                 heatmap_legend_param = list(
                                   title  = "|r|",
                                   at     = seq(0, 1, by = 0.2),
                                   labels = seq(0, 1, by = 0.2)
                                 ),
                                 show_heatmap_legend = FALSE
)


pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/NRseq_pulselabel_heatmap_annot_legend.pdf",
    width  = 4,          
    height = 4,
    useDingbats = FALSE)          
draw(ht_pl)                  
dev.off()

##### RNAdecayCafe #####

RNAdegDB_complete <- read_csv("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/RNAdegDB_OneTable_filtered_refined.csv")

pulse_label_filt_ests <- RNAdegDB_complete %>%
  dplyr::filter(
    type == "gene" # No isoform-level estimates
    &
      cell_line != "BC1" # bad data to filter later
  )

pl_filt_wide <- pulse_label_filt_ests %>%
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

plf_corrmat <- cor(M_plf,
                   method = "pearson",
                   use = "pairwise.complete.obs")


col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),         
  rev(brewer.pal(11, "RdYlBu"))         
)

label_times <- RNAdegDB_complete %>%
  dplyr::filter(
    cell_line != "BC1"
  ) %>%
  dplyr::select(
    sample,
    label_time,
    cell_line,
    total_reads
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    label_factor = factor(
      case_when(
        label_time < 2 ~ "< 2hr",
        label_time >= 2 & label_time <= 4 ~ "[2, 4]hr",
        label_time > 4 ~ "> 4hr"
      ),
      levels = c(
        "< 2hr",
        "[2, 4]hr",
        "> 4hr"
      )
    ),
    cell_line = factor(
      case_when(
        grepl("^MOLM", cell_line) ~ "MOLM",
        grepl("^HEK293", cell_line) ~ "HEK293",
        .default = cell_line
      )
    ),
    seqdepth_factor = factor( # Nothing useful came from this one
      case_when(
        total_reads < 1000000 ~ "< 1 mil",
        total_reads < 10000000 ~ "[1, 10] mil",
        total_reads < 50000000 ~ "[10, 50] mil",
        total_reads > 50000000 ~ "> 50 mil"
      )
    )
  )


# Label time colors
ann_cols <- c("< 2hr" = "deepskyblue", 
              "[2, 4]hr" = "darkgray",
              "> 4hr" = "darkred")



row_ha <- rowAnnotation(
  label_time = label_times$label_factor,
  cell_line = label_times$cell_line,
  col = list(label_time = ann_cols,
             cell_line = cell_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

row_ha_legend <- rowAnnotation(
  label_time = label_times$label_factor,
  cell_line = label_times$cell_line,
  col = list(label_time = ann_cols,
             cell_line = cell_cols),
  show_annotation_name = FALSE,
  show_legend = TRUE
)



# With annotations
ht_rdb_col <- ComplexHeatmap::Heatmap(abs(plf_corrmat),
                                      col               = col_fun,
                                      cluster_rows      = TRUE,
                                      cluster_columns   = TRUE,
                                      show_row_names    = FALSE,
                                      show_column_names = FALSE,
                                      show_row_dend = FALSE,
                                      show_column_dend = FALSE,
                                      right_annotation = row_ha,
                                      row_names_gp      = gpar(fontsize = 6),
                                      column_names_gp   = gpar(fontsize = 6),
                                      border            = FALSE,            # border around the whole heatmap
                                      heatmap_legend_param = list(
                                        title  = "|r|",
                                        at     = seq(0, 1, by = 0.2),
                                        labels = seq(0, 1, by = 0.2)
                                      ),
                                      show_heatmap_legend = FALSE
)


ht_rdb_col


# With annotations and legends
ht_rdb_col_leg <- ComplexHeatmap::Heatmap(abs(plf_corrmat),
                                          col               = col_fun,
                                          cluster_rows      = TRUE,
                                          cluster_columns   = TRUE,
                                          show_row_names    = FALSE,
                                          show_column_names = FALSE,
                                          show_row_dend = FALSE,
                                          show_column_dend = FALSE,
                                          right_annotation = row_ha_legend,
                                          row_names_gp      = gpar(fontsize = 6),
                                          column_names_gp   = gpar(fontsize = 6),
                                          border            = FALSE,            # border around the whole heatmap
                                          heatmap_legend_param = list(
                                            title  = "|r|",
                                            at     = seq(0, 1, by = 0.2),
                                            labels = seq(0, 1, by = 0.2)
                                          ),
                                          show_heatmap_legend = TRUE
)



pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/RDB_heatmap_annot.pdf",
    width  = 2.5,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht_rdb_col)                  
dev.off()

pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/RDB_heatmap_annot_legend.pdf",
    width  = 5,          
    height = 7,
    useDingbats = FALSE)          
draw(ht_rdb_col_leg)                  
dev.off()



##### Saluki database (high quality) vs. RNAdegDB ######
### Need to run everything above this before running this section


### Make combined matrix
ht_rdb <- draw(ht_rdb)
ord_idx <- row_order(ht_rdb)

coln_saluki <- colnames(saluki_hls)
ok_samples <- ok_samples[!grepl("^Simon", ok_samples) &
                           !grepl("^Shendure", ok_samples)]
ok_cols <- which(coln_saluki %in% ok_samples)

saluki_ok <- saluki_hls[,c(2, ok_cols)]

saluki_ok %>%
  dplyr::rename(
    XF = `Gene name`
  ) %>%
  dplyr::inner_join(
    pl_filt_wide,
    by = "XF"
  )


M_comp <- saluki_ok %>%
  dplyr::rename(
    XF = `Gene name`
  ) %>%
  dplyr::inner_join(
    pl_filt_wide,
    by = "XF"
  ) %>%
  dplyr::select(-XF) %>%
  as.matrix()

class(M_comp) <- "numeric"

comp_corrmat <- cor(M_comp,
                    method = "pearson",
                    use = "pairwise.complete.obs")


### Figure out manual ordering
new_order <- c(rownames(plf_corrmat)[ord_idx], ok_samples)
new_id <- match(new_order, rownames(comp_corrmat))



### Matrix to plot
mat <- abs(comp_corrmat[new_id, new_id])


## Color scale
col_fun <- colorRamp2(
  seq(0, 1, length.out = 11),
  rev(brewer.pal(11, "RdYlBu"))
)


##### Construct cell tibbles #####

annot_rdb <- RNAdegDB_complete %>%
  dplyr::filter(
    cell_line != "BC1"
  ) %>%
  dplyr::select(
    sample,
    cell_line
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    cell_line = factor(
      case_when(
        grepl("^MOLM", cell_line) ~ "MOLM",
        grepl("^HEK293", cell_line) ~ "HEK293",
        .default = cell_line
      )
    )
  ) %>%
  dplyr::rename(
    cell = cell_line
  )


# Sketchy b/c have to run top of script first
saluki_tbl <- tibble(
  sample = rownames(new_corrmat)[ord_idx_sal_nonrseq]
) %>%
  dplyr::mutate(
    good = factor(as.character(sample %in% ok_samples)),
    method = sapply(str_split(
      sample, "_"
    ),
    function(x) x[2]
    ),
    cell = sapply(str_split(
      sample, "_"
    ),
    function(x) x[3]
    ),
    dataset = sapply(str_split(
      sample, "_"
    ),
    function(x) x[1]
    )
  ) %>%
  dplyr::mutate(
    method = case_when(
      method %in% c("ActD", "Aman") ~ "txni",
      method %in% c("4sU", "BrU", "BrU4sU") ~ "enrich"
    ),
    cell = case_when(
      grepl("^GM", cell) ~ "LCL",
      .default = cell
    )
  ) %>%
  dplyr::filter(
    good == "TRUE"
  ) %>%
  dplyr::select(
    sample, cell
  )

cell_tbl <- bind_rows(
  annot_rdb,
  saluki_tbl
)

row_ha <- rowAnnotation(
  cell = cell_tbl$cell,
  col    = list(cell = cell_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  mat,
  col               = col_fun,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 6),
  column_names_gp   = gpar(fontsize = 6),
  border            = FALSE,            # border around the whole heatmap
  right_annotation   = row_ha,
  #top_annotation    = top_ha,
  heatmap_legend_param = list(
    title  = "|r|",
    at     = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  show_heatmap_legend = FALSE
)



pdf("C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/RDB_vs_Saluki_heatmap_annot.pdf",
    width  = 2.3,          
    height = 2.1,
    useDingbats = FALSE)          
draw(ht)                  
dev.off()


##### Make box plots of correlations (1E) #####

nr_nrow <- length(ord_idx)
nonnr_nrow <- length(ok_samples)

nr_vs_nr_idx <- expand.grid(1:nr_nrow, 1:nr_nrow) %>%
  as_tibble() %>%
  dplyr::filter(
    Var2 < Var1
  )

nr_vs_nr_corrs <- mat[cbind(nr_vs_nr_idx$Var1, nr_vs_nr_idx$Var2)]


non_vs_non_idx <- expand.grid((nr_nrow+1):nrow(mat), (nr_nrow+1):nrow(mat)) %>%
  as_tibble() %>%
  dplyr::filter(
    Var2 < Var1
  )

non_vs_non_corrs <- mat[cbind(non_vs_non_idx$Var1, non_vs_non_idx$Var2)]


non_vs_nr_idx <- expand.grid((nr_nrow+1):nrow(mat), 1:nr_nrow) %>%
  as_tibble() %>%
  dplyr::filter(
    Var2 < Var1
  )

non_vs_nr_corrs <- mat[cbind(non_vs_nr_idx$Var1, non_vs_nr_idx$Var2)]


corr_tbl <- tibble(
  corr = c(nr_vs_nr_corrs,
           non_vs_non_corrs,
           non_vs_nr_corrs),
  type = factor(rep(
    c("NR-NR",
      "Non-Non",
      "Non-NR"),
    times = c(
      length(nr_vs_nr_corrs),
      length(non_vs_non_corrs),
      length(non_vs_nr_corrs)
    )
  ),
  levels = c("Non-Non",
             "Non-NR",
             "NR-NR")
  )
)

gcors <- corr_tbl %>%
  ggplot(
    aes(
      x = type,
      y = corr
    )
  ) +
  geom_violin(
    fill = "darkgray",
    color = "black",
    linewidth = 0.4
  ) +
  geom_boxplot(
    outliers = FALSE,
    notch = FALSE,
    fill = "white",
    color = "black",
    linewidth = 0.2,
    width = 0.3
  ) +
  theme_classic() +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  xlab(
    "Comparison"
  ) +
  ylab(
    "Pearson r"
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10))  #change font size of legend title

ggsave(
  filename = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Manuscript/Figures/Figure1_panels/Corr_violin_box.pdf",
  plot = gcors,
  width = 2,
  height = 2.1
)
