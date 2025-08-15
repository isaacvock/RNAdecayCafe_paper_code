### PURPOSE OF THIS SCRIPT
## Compile some functions used across scripts and analyses

# Assessing kdeg vs. cell line correlation -------------------------------------

cell_line_analysis <- function(db, FDR_cutoff = 0.01){
  
  # Conservative, home-spun ANOVA
  # Replaces mean within-group variance with
  # max within-group variance.
  expression_anova <- conservative_anova(
    db,
    type = "expression"
  )
  
  expression_hits <- expression_anova %>%
    dplyr::filter(
      FDR < FDR_cutoff
    ) %>%
    dplyr::select(
      XF
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      !grepl("-", XF) # Remove fusion genes that can be a problem
    )
  
  # Step 2: Correlate half-life with cell-type gene expression ###
  
  kdeg_corrs_list <- calc_kdeg_corr(db)
  
  kdeg_corrs <- kdeg_corrs_list$corr
  
  # Step 3: Only care about expression hits
  
  kdeg_corrs <- kdeg_corrs %>%
    dplyr::mutate(
      R2 = ifelse(
        kdeg_corr > 0,
        0,
        kdeg_corr^2
      )
    ) %>%
    right_join(
      expression_hits,
      by = "XF"
    ) %>%
    tidyr::replace_na(
      list(
        R2 = 0
      )
    )
  
  return(kdeg_corrs)
}

conservative_anova <- function(db, type = c("expression", "kdeg")){
  
  type <- match.arg(type)
  
  if(type == "expression"){
    
    anova_res <- db %>%
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
      ) %>%
      ### End of normalization
      dplyr::group_by(
        XF
      ) %>%
      dplyr::mutate(
        log_RPKM_z = (log_RPKM - mean(log_RPKM)) / sd(log_RPKM)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(XF, cell_line) %>%
      # Only keep cell lines with multiple replicates of filtered data
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::group_by(XF) %>%
      # Only keep genes with multiple cell lines of filtered data
      dplyr::filter(length(unique(cell_line)) > 1) %>%
      dplyr::group_by(XF, cell_line) %>%
      dplyr::mutate(
        within_group_var = var(log_RPKM_z)
      ) %>%
      dplyr::group_by(XF) %>%
      dplyr::mutate(
        within_group_var = max(within_group_var),
        between_group_var = var(log_RPKM_z)
      ) %>%
      dplyr::summarise(
        Fstat = mean(between_group_var) / mean(within_group_var),
        ntot = dplyr::n(),
        ncell = length(unique(cell_line))
      ) %>%
      dplyr::mutate(
        pval = pf(Fstat, ncell - 1, ntot - ncell, lower.tail = FALSE)
      ) %>%
      dplyr::mutate(
        FDR = p.adjust(pval, method = "BH")
      ) %>%
      dplyr::select(
        XF, Fstat, pval, FDR
      )
    
  }else{
    
    
    anova_res <- db %>%
      dplyr::group_by(XF) %>%
      dplyr::mutate(
        # Maybe don't actually want to z-score normalize, and certainly not by sample,
        # could introduce kdeg variance that doesn't exist. Using dropout normalized
        # values should be enough to put these all on similar scales
        log_kdeg_z = donorm_log_kdeg 
        #log_kdeg_z = (donorm_log_kdeg - mean(donorm_log_kdeg)) / sd(donorm_log_kdeg)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(XF, cell_line) %>%
      # Only keep cell lines with multiple replicates of filtered data
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::group_by(XF) %>%
      # Only keep genes with multiple cell lines of filtered data
      dplyr::filter(length(unique(cell_line)) > 1) %>%
      dplyr::group_by(XF, cell_line) %>%
      dplyr::mutate(
        within_group_var = var(log_kdeg_z)
      ) %>%
      dplyr::group_by(XF) %>%
      dplyr::mutate(
        within_group_var = max(within_group_var),
        between_group_var = var(log_kdeg_z)
      ) %>%
      dplyr::summarise(
        Fstat = mean(between_group_var) / mean(within_group_var),
        ntot = dplyr::n(),
        ncell = length(unique(cell_line))
      ) %>%
      dplyr::mutate(
        pval = pf(Fstat, ncell - 1, ntot - ncell, lower.tail = FALSE)
      ) %>%
      dplyr::mutate(
        FDR = p.adjust(pval, method = "BH")
      ) %>%
      dplyr::select(
        XF, Fstat, pval, FDR
      )
  }
  

  return(anova_res)
  
  
}


get_nsamp <- function(db){
  
  nsamps <- db %>%
    dplyr::select(sample) %>%
    dplyr::distinct() %>%
    unlist() %>% unname() %>% length()
  
  return(nsamps)
}

get_ncell <- function(db){
  
  ncells <- db %>%
    dplyr::select(cell_line) %>%
    dplyr::distinct() %>%
    unlist() %>% unname() %>% length()
  
  
  return(ncells)
  
}


calc_kdeg_corr <- function(db){
  
  nsamps <- get_nsamp(db)
  ncells <- get_ncell(db)
    
  
  
  # Correlation sample-to-sample
  kdeg_corrs_samps <- db %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      RPKM = ifelse(
        threePseq,
        reads / (sum(reads) /  1000000),
        (reads / (exon_length / 1000)) / (sum(reads) /  1000000)
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
    ) %>%
    ### End of normalization
    dplyr::ungroup() %>%
    dplyr::mutate(
      log_thalf = log(log(2)/exp(donorm_log_kdeg))
    ) %>%
    dplyr::group_by(XF) %>%
    dplyr::filter(dplyr::n() == nsamps) %>%
    dplyr::summarise(
      kdeg_corr_pearson = cor(donorm_log_kdeg, log_RPKM),
      range_comp =  diff(range(donorm_log_kdeg)) / diff(range(log_RPKM)),
      R2_identity = 1 - var(log_RPKM + donorm_log_kdeg)/var(log_RPKM),
      slope = cov(log_thalf, log_RPKM) / var(log_thalf),
      R2 = cor(donorm_log_kdeg, log_RPKM)^2
    ) %>%
    dplyr::mutate(
      range_comp = ifelse(
        range_comp > 1 & kdeg_corr_pearson < 0,
        1,
        ifelse(kdeg_corr_pearson > 0,
               0,
               range_comp)
      ),
      R2_identity = ifelse(
        kdeg_corr_pearson > 0,
        0,
        R2_identity
      )
    )
  
  # Correlation cell-type to cell-type
  
  kdeg_corrs_cells <- db %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      RPKM = ifelse(
        threePseq,
        reads / (sum(reads) /  1000000),
        (reads / (exon_length / 1000)) / (sum(reads) /  1000000)
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
    ) %>%
    ### End of normalization
    dplyr::ungroup() %>%
    dplyr::group_by(XF, cell_line) %>%
    dplyr::summarise(
      donorm_log_kdeg = mean(donorm_log_kdeg),
      log_RPKM = mean(log_RPKM)
    ) %>%
    dplyr::mutate(
      log_thalf = log(log(2)/exp(donorm_log_kdeg))
    ) %>%
    dplyr::filter(dplyr::n() == ncells) %>%
    dplyr::summarise(
      kdeg_corr_pearson = cor(donorm_log_kdeg, log_RPKM),
      range_comp =  diff(range(donorm_log_kdeg)) / diff(range(log_RPKM)),
      R2_identity = 1 - var(log_RPKM + donorm_log_kdeg)/var(log_RPKM),
      slope = cov(log_RPKM, log_thalf)/var(log_thalf),
      R2 = cor(donorm_log_kdeg, log_RPKM)^2
    ) %>%
    dplyr::mutate(
      range_comp = ifelse(
        range_comp > 1 & kdeg_corr_pearson < 0,
        1,
        ifelse(kdeg_corr_pearson > 0,
               0,
               range_comp)
      ),
      R2_identity = ifelse(
        kdeg_corr_pearson > 0,
        0,
        R2_identity
      )
    )
  
  ##### NEED TO COMPLETE RANGE CALCULATION
  
  # Correlation of interest is minimum in magnitude
  kdeg_corrs <- kdeg_corrs_samps %>%
    dplyr::inner_join(
      kdeg_corrs_cells,
      by = "XF",
      suffix = c("_samp", "_cell")
    ) %>%
    dplyr::mutate(
      kdeg_corr = pmin(abs(kdeg_corr_pearson_samp), abs(kdeg_corr_pearson_cell)),
      kdeg_range_fxn = pmin(abs(range_comp_samp), abs(range_comp_cell)),
      R2_identity = pmin(R2_identity_samp, R2_identity_cell)
    ) %>%
    dplyr::mutate(
      kdeg_corr = ifelse(
        kdeg_corr == abs(kdeg_corr_pearson_samp),
        kdeg_corr_pearson_samp,
        kdeg_corr_pearson_cell
      )
    )
  
  return(
    list(
      samp_corr = kdeg_corrs_samps,
      cell_corr = kdeg_corrs_cells,
      corr = kdeg_corrs
    )
  )
  
}

RNA_aov <- function(db){
  
  RNA_aov <- db  %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      RPKM = ifelse(
        threePseq,
        reads / (sum(reads) /  1000000),
        (reads / (exon_length / 1000)) / (sum(reads) /  1000000)
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
    ) %>%
    ### End of normalization
    dplyr::ungroup() %>%
    dplyr::group_by(XF) %>%
    dplyr::filter(dplyr::n() > 1 & length(unique(cell_line)) > 1) %>%
    dplyr::group_modify(
      ~ broom::tidy(
        aov(log_RPKM ~ cell_line, data = .x)
      )
    ) %>%
    ungroup() %>%
    filter(term == "cell_line")
  
  
  RNA_aov_table <- RNA_aov %>%
    transmute(
      XF,
      Fstat = statistic,
      pval = p.value,
      FDR = p.adjust(p.value, method = "BH")
    ) %>%
    na.omit()
  
  
  return(RNA_aov_table)
  
}

id_RNA_var <- function(RNA_aov_table,
                       log10_FDR_cutoff = 4){
  
  RNA_that_varies <- RNA_aov_table %>%
    dplyr::filter((-log10(FDR)) > log10_FDR_cutoff) %>%
    dplyr::select(XF)
  
  return(RNA_that_varies)
  
}


make_heatmap_df <- function(db, kdeg_corrs, variable_RNA){
  
  
  rownorm_heatmap_cell <- db %>%
    dplyr::ungroup() %>%
    ### Be honest about dynamic range
    # For this analysis, the conservative dynamic range is a half-life ranging
    # from much lower than the longest label time to much higher
    # than the shortest label time
    dplyr::mutate(
      donorm_log_kdeg = case_when(
        exp(donorm_log_kdeg) > -log(1 - (reads - 1)/reads)/label_time ~ log(-log(1 - (reads - 1)/reads)/min(label_time)),
        exp(donorm_log_kdeg) < -log(1 - (1)/reads)/label_time ~ log(-log(1 - (5)/reads)/max(label_time)),
        .default = donorm_log_kdeg
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      RPKM = ifelse(
        threePseq,
        reads / (sum(reads) /  1000000),
        (reads / (exon_length / 1000)) / (sum(reads) /  1000000)
      )
    )  %>%
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
    ) %>%
    ### End of normalization
    dplyr::inner_join(
      RNA_that_varies,
      by = c("XF")
    ) %>%
    dplyr::group_by(XF, cell_line) %>%
    dplyr::summarise(
      cell_line_avg = mean( log(log(2)/ exp(donorm_log_kdeg)) ),
      log_RPKM = mean(log_RPKM)
    ) %>%
    dplyr::inner_join(
      kdeg_corrs,
      by = "XF"
    ) %>%
    dplyr::group_by(XF) %>%
    dplyr::mutate(
      lkdeg_sd = sd(cell_line_avg),
      lkdeg_range = max(cell_line_avg) - min(cell_line_avg)
    ) %>%
    dplyr::arrange(
      kdeg_corr
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      avg_sd = mean(lkdeg_sd),
      med_sd = median(lkdeg_sd),
      max_sd = quantile(lkdeg_sd)[4],
      max_range = max(lkdeg_range)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(XF) %>%
    dplyr::mutate(
      rownorm_log_kdeg = (cell_line_avg - mean(cell_line_avg)) / max_range
      #rownorm_log_kdeg = cell_line_avg - mean(cell_line_avg)
    ) %>%
    dplyr::mutate(
      RNA_order = factor(paste0("col", rank(log_RPKM)),
                         levels = paste0("col", 1:dplyr::n()))
    ) %>%
    dplyr::mutate(
      min_expression_cell = cell_line[RNA_order == "col1"],
      max_expression_cell = cell_line[RNA_order == paste0("col", dplyr::n()) ]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      XF, kdeg_corr, RNA_order, rownorm_log_kdeg, 
      min_expression_cell, max_expression_cell
    ) %>%
    tidyr::pivot_wider(
      values_from = rownorm_log_kdeg,
      names_from = RNA_order
    )
  
  # kdeg's are actually half-lifes for the paper figure
  rownorm_heatmap_samp <- db %>%
    dplyr::ungroup() %>%
    ### Be honest about dynamic range
    # For this analysis, the conservative dynamic range is a half-life ranging
    # from much lower than the longest label time to much higher
    # than the shortest label time
    dplyr::mutate(
      donorm_log_kdeg = case_when(
        exp(donorm_log_kdeg) > -log(1 - (reads - 1)/reads)/label_time ~ log(-log(1 - (reads - 1)/reads)/min(label_time)),
        exp(donorm_log_kdeg) < -log(1 - (1)/reads)/label_time ~ log(-log(1 - (5)/reads)/max(label_time)),
        .default = donorm_log_kdeg
      )
    ) %>%
    # Convert to half-life
    dplyr::mutate(
      donorm_log_kdeg = log(log(2)/exp(donorm_log_kdeg))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      RPKM = ifelse(
        threePseq,
        reads / (sum(reads) /  1000000),
        (reads / (exon_length / 1000)) / (sum(reads) /  1000000)
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
    ) %>%
    ### End of normalization
    dplyr::inner_join(
      RNA_that_varies,
      by = c("XF")
    ) %>%
    dplyr::inner_join(
      kdeg_corrs,
      by = "XF"
    ) %>%
    dplyr::group_by(XF) %>%
    dplyr::mutate(
      lkdeg_sd = sd(donorm_log_kdeg),
      lkdeg_range = max(donorm_log_kdeg) - min(donorm_log_kdeg)
    ) %>%
    dplyr::arrange(
      kdeg_corr
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      avg_sd = mean(lkdeg_sd),
      med_sd = median(lkdeg_sd),
      max_sd = quantile(lkdeg_sd)[4],
      max_range = max(lkdeg_range)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(XF) %>%
    dplyr::mutate(
      rownorm_log_kdeg = (donorm_log_kdeg - mean(donorm_log_kdeg)) / max_range
      #rownorm_log_kdeg = cell_line_avg - mean(cell_line_avg)
    ) %>%
    dplyr::mutate(
      RNA_order = factor(paste0("col", rank(log_RPKM)),
                         levels = paste0("col", 1:dplyr::n()))
    ) %>%
    dplyr::mutate(
      min_expression_cell = cell_line[RNA_order == "col1"],
      max_expression_cell = cell_line[RNA_order == paste0("col", dplyr::n()) ]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      XF, kdeg_corr, RNA_order, rownorm_log_kdeg, min_expression_cell, max_expression_cell
    ) %>%
    tidyr::pivot_wider(
      values_from = rownorm_log_kdeg,
      names_from = RNA_order
    )
  
  return(
    list(
      rownorm_cell = rownorm_heatmap_cell,
      rownorm_samp = rownorm_heatmap_samp
    )
  )
}



# Analysis pipeline ------------------------------------------------------------


analyze_nrseq_data <- function(
    dataset_str,
    metadf,
    use_arrow = TRUE,
    dropout_correct = TRUE,
    dropout_normalize = FALSE,
    normalize_across_tl = TRUE,
    pold_from_nolabel = TRUE,
    returnEZ = TRUE,
    kstrat = c("standard", "pulse-chase"),
    exclude_pulse_estimates = FALSE,
    suffix = "_refined_WTonly"
){
  
  kstrat <- match.arg(kstrat)
  
  if("tpulse" %in% colnames(metadf)){
    exclude_pulse_estimates <- exclude_pulse_estimates | !any(metadf$tpulse > 0 & metadf$tpulse < 8 & metadf$tchase == 0)
  }
  
  ### Run settings
  savedir <- paste0("C:/Users/isaac/Yale University/Simon Lab (General) - RNAdegDB_curation/EZbakR_output/", dataset_str, "/")
  
  
  ### Analysis
  
  if(use_arrow){
    
    cB <- arrow::open_dataset(paste0("C:/Users/isaac/Yale University/Simon Lab (General) - RNAdegDB_curation/fastq2EZbakR_output/", dataset_str,"/arrow_dataset"))
    
    ezbdo <- EZbakRArrowData(cB, metadf)
    
  }else{
    
    cB <- fread(paste0("C:/Users/isaac/Yale University/Simon Lab (General) - RNAdegDB_curation/fastq2EZbakR_output/", dataset_str,"/cB/cB.csv.gz"))
    
    ezbdo <- EZbakRData(cB, metadf)
    
  }
  
  
  pold_from_nolabel <- pold_from_nolabel & any(metadf$tl == 0)
  dropout_correct <- dropout_correct & any(metadf$tl == 0)
  
  ezbdo <- EstimateFractions(ezbdo,
                             features = "XF",
                             pold_from_nolabel = pold_from_nolabel)
  
  
  if(dropout_correct){
    
    ezbdo <- CorrectDropout(ezbdo)
    
  }
  
  if(dropout_normalize){
    
    ezbdo <- NormalizeForDropout(ezbdo,
                                 normalize_across_tls = normalize_across_tl)
    
  }
  
  
  
  QCchecks <- EZQC(ezbdo)
  
  
  ezbdo <- EstimateKinetics(ezbdo,
                            strategy = kstrat,
                            exclude_pulse_estimates = exclude_pulse_estimates)
  
  
  
  ### Save output
  dir.create(file.path(savedir))
  
  write_csv(
    ezbdo$kinetics$XF,
    paste0(savedir, dataset_str, "_sample_kdegs", suffix ,".csv")
  )
  

  saveRDS(
    QCchecks,
    paste0(savedir, dataset_str,"_EZbakR_QCchecks", suffix,".rds")
  )
  
  
  write_csv(
    ezbdo$mutation_rates$TC,
    paste0(savedir, dataset_str,"_Estimated_Mutrates", suffix, ".rds")
  )
  
  
  write_csv(
    metadf,
    paste0(savedir, dataset_str,"_metadf", suffix, ".rds")
  )
  
  
  saveRDS(
    ezbdo,
    paste0(savedir, dataset_str, "_EZbakRFit", suffix,".rds")
  )
  
  
  if(returnEZ){
    return(ezbdo)
  }
  
}



# Data -------------------------------------------------------------------------


# Datasets with questionable correlation with other dataset
questionable_datasets <- c(
  "Narain_etal_2021_U2OS", # Pretty high dropout it seems
  "Rawat_etal_2021_HeLa",
  "Zhou_etal_2024_HEK293T", # Bad replicate-to-replicate correlation especially in 1.5 hour label time data
  "Zhang_etal_2023_BE2C", # Weird anomaly with lots of genes having 0 mutations in thousands of reads
  "Mowery_etal_2018_Nalm6" 
)

# Samples with questionable correlation with other samples
questionable_samples <- c(
  "SRR131573505" # Only 1 hour Finkel (Calu3) data
)


##### Metadfs ######

# Williams et al. 2025 MDA-MB-453
metadf_Williams_etal_2025_MDAMB <- tibble(
  sample = c("SRR31699597", "SRR31699594", "SRR31699592", "SRR31699596",
             "SRR31699590", "SRR31699588",
             "SRR31699605", "SRR31699601", "SRR31699599"),
  tl = 1,
  compartment = c("total", "total", "nucleus", "nucleus",
                  "cytoplasm", "cytoplasm",
                  "total", "nucleus", "cytoplasm")
)



# Luo et al. 2020 HEK293T
metadf_Luo_etal_2020_HEK293T <- tibble(
  sample = c("SRR10887642", "SRR10887643",
             "SRR10887644", "SRR10887645",
             "SRR10887646", "SRR10887647",
             "SRR10887648", "SRR10887649", "SRR10887650",
             "SRR10887651", "SRR10887652", "SRR10887653"
  ),
  tl = c(2, 2,
         2, 2,
         0, 0,
         2, 2, 0,
         2, 2, 0),
  genotype = c("WT", "WT",
               "DCP2", "DCP2",
               "WT", "DCP2",
               "WT", "WT", "WT",
               "MSI2", "MSI2", "MSI2")
)


# Ietswaart et al. 2024 K562
metadf_Ietswaart_etal_2024_K562 <- tibble(
  sample = c("SRR20080264", "SRR20080263", "SRR20080262", "SRR20080261", "SRR20080260",
             "SRR20080287", "SRR20080286", "SRR20080285", "SRR20080284", "SRR20080283",
             "SRR20080253", "SRR20080252", "SRR20080251", "SRR20080250",
             "SRR20080243", "SRR20080242", "SRR20080241", "SRR20080240",
             "SRR20080236", "SRR20080232", "SRR20080231", "SRR20080230",
             "SRR24302461", "SRR24302460", "SRR24302459", "SRR24302458",
             "SRR24302454", "SRR24302457", "SRR24302453", "SRR24302452",
             "SRR24302448", "SRR24302447", "SRR24302446", "SRR24302449",
             "SRR24302443", "SRR24302442", "SRR24302440", "SRR24302439",
             "SRR24302437", "SRR24302436", "SRR24302435", "SRR24302434",
             paste0("SRR200802", 74:70),
             paste0("SRR200802", 97:93),
             paste0("SRR200802", 79:75),
             paste0("SRR20080", 302:298)),
  tl = c(0, 15, 30, 60, 120, # total
         0, 15, 30, 60, 120, # total
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 0, 60, 60,
         0, 15, 30, 60, 120, # cytoplasm
         0, 15, 30, 60, 120, # cytoplasm
         0, 15, 30, 60, 120, # nuclear
         0, 15, 30, 60, 120), # nuclear
  siRNA = c("none", "none", "none", "none", "none",
            "none", "none", "none", "none", "none",
            "scramble", "scramble", "scramble", "scramble",
            "DDX3X", "DDX3X", "DDX3X", "DDX3X",
            "PABPC4", "PABPC4", "PABPC4", "PABPC4",
            "scramble", "scramble", "scramble", "scramble",
            "DIS3", "DIS3", "DIS3", "DIS3",
            "EXOSC10", "EXOSC10", "EXOSC10", "EXOSC10",
            "PABPN1", "PABPN1", "PABPN1", "PABPN1",
            "ZFC3H1", "ZFC3H1", "ZFC3H1", "ZFC3H1",
            "none", "none", "none", "none", "none",
            "none", "none", "none", "none", "none",
            "none", "none", "none", "none", "none",
            "none", "none", "none", "none", "none"),
  compartment = c("total", "total", "total", "total", "total",
                  "total", "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  "total", "total", "total", "total",
                  rep("cytoplasm", times = 10),
                  rep("nucleus", times = 10))
) %>%
  dplyr::filter(
    siRNA == "none"
  )


# Swartzel et al. 2022 MOLM14
metadf_Swartzel_etal_2022_MOLM14 <- tibble(
  sample = paste0("SRR19981", 802:794),
  tl = c(0,
         2, 2,
         2, 2, 
         2, 2,
         2, 2),
  treatment = c("UT",
                "JCS1", "JCS1",
                "JCS2", "JCS2",
                "RG3039", "RG3039",
                "UT", "UT")
)


# Harada et al. 2022 MV411
metadf_Harada_etal_2022_MV411 <- tibble(
  sample = c('SRR15654897','SRR15654896','SRR15654895','SRR15655142',
             'SRR15655136','SRR15655135'),
  tl = c(1, 1, 1, 1,
         0, 0)
)



# Rawat et al. 2021 HeLa
metadf_Rawat_etal_2021_HeLa <- tibble(
  sample = c("SRR11067753", "SRR11067754", "SRR11067755",
             "SRR11067756", "SRR11067757", "SRR11067758",
             "SRR11067765", "SRR11067766", "SRR11067767",
             "SRR11067768", "SRR11067769", "SRR11067770"),
  tl = c(105, 105, 105,
         105, 105, 105,
         105, 105, 105,
         105, 105, 105),
  treatment = c("HS", "HS", "HS",
                "HS", "HS", "HS",
                "NHS", "NHS", "NHS",
                "NHS", "NHS", "NHS"),
  siRNA = c("Ctrl", "Ctrl", "Ctrl",
            "ZNF451", "ZNF451", "ZNF451",
            "Ctrl", "Ctrl", "Ctrl",
            "ZNF451", "ZNF451", "ZNF451")
)


# Narain et al. 2021 U2OS
metadf_Narain_etal_2021_U2OS <- tibble(
  sample = paste0("SRR131594", 26:49),
  tl = rep(c(0, 1, 2, 4), each = 6),
  treatment = rep(rep(c("Control", "Auxin"), each = 3), times = 4)
)


# Zuckerman et al. 2020 MCF7
metadf_Zuckerman_etal_2020_MCF7 <- tibble(
  sample = paste0("SRR103169", 40:57),
  tl = rep(c(0, 0, 2, 2, 4, 4, 8, 8, 8),
           times = 2),
  compartment = "cytoplasmic",
  siRNA = rep(c("NT", "NXF1"), each = 9)
)


# Finkel et al. 2021 Calu3
metadf_Finkel_etal_2021_Calu3 <- tibble(
  sample = c(paste0("SRR1357350", 4:9),
             paste0("SRR135735", 10:15)),
  tl = rep(c(0, 1, 2, 2, 3, 4), times = 2),
  hpi = c(0, 0, 0, 0, 0, 0,
          3, 4, 5, 5, 6, 7)
)


# Biancon et al. 2022 HEL
metadf_Biancon_etal_2022_HEL <- tibble(
  sample = c("SRR17803613", "SRR17803612", "SRR17803611", "SRR17803610",
             "SRR17803609", "SRR17803608",
             "SRR17803599", "SRR17803598", "SRR17803597"),
  tl = c(0, 2, 2, 2, 2, 2, 2, 2, 2),
  U2AF1 = c("WT", "WT", "WT", "WT", "WT",
            "S34F", "S34F", "Q157R", "Q157R")
)


# Mabin et al. 2025 HEK293T
metadf_Mabin_etal_2025_HEK293T <- metadf <- tibble(
  sample = c(paste0("DMSO_",1:3), "DMSO_nos4U"),
  tl = c(2, 2, 2, 0)
)


# Shimoda et al. 2021 (BC-1/BCBL-1)
metadf_Shimoda_etal_2021_BC1_BCBL1 <- tibble(
  sample = c(paste0("SRR1447530", 5:9),
             paste0("SRR144753", 10:16)),
  tl = 1,
  cell = rep(c("BC-1", "BCBL-1"), each = 6),
  treatment = rep(rep(c("Control", "WT_peptide", "Mutant_peptide"), each = 2),
                  times = 2)
)



# Leidendecker et al. 2020 MCC
metadf_Leidendecker_etal_2020_MCC <- tibble(
  sample = paste0("SRR114589", 17:28),
  tl = rep(rep(c(1, 6), each = 3), times = 2),
  treatment = rep(c("Ctrl", "GSK"), each = 6)
)


# Zhang et al. 2023 BE2C cells
metadf_Zhang_etal_2023_BE2C <- tibble(
  sample = paste0("SRR191650", 58:50),
  tl = 1,
  treatment = rep(c("DMSO", "ACY957", "dTAG13"), each = 3)
)


# Williams et al. 2025 MDA-MB-453
metadf_Williams_etal_2025_MDAMB453 <- tibble(
  sample = c("SRR31699597", "SRR31699594", "SRR31699592", "SRR31699596",
             "SRR31699590", "SRR31699588",
             "SRR31699605", "SRR31699601", "SRR31699599"),
  tl = 1,
  compartment = c("total", "total", "nucleus", "nucleus",
                  "cytoplasm", "cytoplasm",
                  "total", "nucleus", "cytoplasm")
)


# Sheppard et al. 2021 CH22 cells
metadf_Sheppard_etal_2021_CH22 <- tibble(
  sample = paste0("SRR1216863", 2:7),
  tl = 7,
  treatment = rep(c("DMSO", "THZ1"), each = 3)
)


# Mowery et al. 2018 Nalm6 cells
metadf_Mowery_etal_2018_Nalm6 <- tibble(
  sample = paste0("SRR799382", 1:8),
  tl = 5,
  treatment = c("DOX", "DOX", "Vehicle", "Vehicle", "Vehicle",
                "DOX", "DOX", "DOX"),
  genotype = c("S20_24E", "S20_24E", "WT", "WT", "WT",
               "WT", "WT", "WT")
)


# Bell et al. 2024 K562 cells
metadf_Bell_etal_2024_K562 <- tibble(
  sample = paste0("SRR280263", 41:36),
  tl = 1,
  treatment = rep(c("DMSO", "MED14-dTAG"), each = 3)
)


# Schofield et al. 2018 K562 cells
metadf_Schofield_etal_2018_K562 <- tibble(
  sample = c("SRR6262128", "SRR6262127", "SRR6262129"),
  tl = c(4, 0, 4) 
)


# Muhar et al. 2018 many cells
metadf_Muhar_etal_2018_many <- tibble(
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


# Thiecke et al. 2020 HeLa cells
metadf_Thiecke_etal_2020_HeLa <- tibble(
  sample = c("SRR11142916", "SRR11142918", "SRR11142920",
             "SRR11142922", "SRR11142924", "SRR11142926"),
  tl = 1,
  SCC1 = c("control", "control", "control",
           "depleted", "depleted", "depleted")
)


# Bauer et al. 2022 THP1 cells
metadf_Bauer_etal_2022_THP1 <- tibble(
  sample = c("SRR18583459", "SRR18583458", "SRR18583457", "SRR18583456", "SRR18583455",
             paste0("SRR185834", 44:48),
             paste0("SRR185834", 34:31),
             "SRR18583423", "SRR18583420", "SRR18583413"),
  tpulse = c(0, 8, 8, 8, 8, 
             0, 8, 8, 8, 8, 
             8, 8, 8, 8,
             1, 1, 1),
  tchase = c(0, 0, 1, 3, 6, 
             0, 0, 1, 3, 6, 
             0, 1, 3, 6,
             0, 0, 0)
)


# An et al. 2024 HEK293T cells
metadf_An_etal_2024_HEK293T <- tibble(
  sample = c("SRR30252747", "SRR30252744", "SRR30252741",
             "SRR30252738", "SRR24829286", "SRR24829283",
             "SRR24829281", "SRR24829277", "SRR24828554"),
  tpulse = c(24, 24, 24, 24, 24, 24, 24, 24, 1),
  tchase = c(0, 3, 6, 12, 0, 3, 6, 12, 0),
  siRNA = c("SCR", "SCR", "SCR", "SCR", "WT", "WT", "WT", "WT", "SCR")
)


# Zhou et al. 2024 HEK293T cells
metadf_Zhou_etal_2024_HEK293T <- tibble(
  sample = paste0("SRR300070", 27:10),
  tl = rep(c(0, 1.5, 3), each = 6),
  treatment = rep(rep(c("DMSO", "STM2457"),
                      each = 3),
                  times = 3)
)


# Rothamel et al. 2021 (THP1 cells)
  # KO is of ELAVL1
metadf_Rothamel_etal_2021_THP1 <- tibble(
  sample = paste0("SRR125398", 78:97),
  tpulse = 16,
  tchase = c(0, 0, 0, 0,
             1, 1, 3, 3,
             3, 3, 6, 6,
             6, 6, 8, 8,
             1, 1, 8, 8),
  genotype = c("KO", "KO", "WT", "WT",
               "KO", "WT", "KO", "KO",
               "WT", "WT", "KO", "KO",
               "WT", "WT", "KO", "KO",
               "KO", "WT", "WT", "KO")
)


# Wu et al. 2019 (K562)
metadf_Wu_etal_2019_K562 <- tibble(
  sample = c(paste0("SRR85712", 43:52),
             paste0("SRR85712", 41:42),
             paste0("SRR85712", 65:66)),
  tpulse = rep(c(24, 0), times = c(12, 2)),
  tchase = c(rep(c(0, 2, 4, 6), each = 3),
             0, 0)
)


# Zhu et al. 2024 HEK293T
metadf_Zhu_etal_2024_HEK293T <- tibble(
  sample = paste0("SRR291643", 36:13),
  tpulse = 24,
  tchase = rep(rep(c(0, 1, 2, 4, 8, 12),
                   each = 2), times = 2),
  gRNA = rep(
    c("CNOT3", "NT"),
    each = 12
  )
)



# Boulias et al. 2019 HEK293T 
  # KO = PCIF1 knockout
metadf_Boulias_etal_2019_HEK293T <- tibble(
  sample = paste0("SRR82444", 73:78),
  tpulse = 24,
  tchase = c(0, 12, 6,
             0, 12, 6),
  genotype = c("KO", "KO", "KO",
               "WT", "WT", "WT")
)



# Wang et al. 2020 THP-1 cells
metadf_Wang_etal_2020_THP1 <- dplyr::tibble(
  sample = paste0("SRR112946", 27:42),
  tpulse = 4,
  tchase = rep(c(0, 1.5, 4.5, 7.5), times = 4),
  treatment = rep(rep(c("control", "A5-KD"), each = 4),
                  times = 2)
)



##### Misc. ######

total_rna_datasets <- c(
  "Biancon_etal_2022_HEL",
  "Finkel_etal_2021_Calu3",
  "Harada_etal_2022_MV411",
  "Ietswaart_etal_2024_K562",
  "Luo_etal_2020_HEK293T",
  "Mowery_etal_2018_Nalm6",
  "Swartzel_etal_2022_MOLM14",
  "Zuckerman_etal_2020_MCF7"
)

# Different naming in another script
ietswaart_metadf <- metadf_Ietswaart_etal_2024_K562


# Functions --------------------------------------------------------------------


##### Identify variably stable genes #####

find_variable_kdeg <- function(
    kdegs,
    exon_file = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/GeneInfo/Exonic_sequence_stats.csv",
    use_ANOVA = TRUE,
    use_setunion = TRUE,
    use_adhoc_filters = TRUE,
    output_stats = TRUE,
    consider_dynamic_range = TRUE,
    nmin = 100,
    RPKmin = 100,
    range_sd_factor = 2,
    fn_lower_limit = 0.05,
    fn_upper_limit = 0.95
    ){
  


  ### Filter
  exonic_lengths <- read_csv(exon_file)
  
  kdegs <- kdegs %>%
    dplyr::inner_join(
    exonic_lengths,
    by = "XF"
  ) %>%
    dplyr::mutate(
      RPK = ifelse(
        threePseq,
        n,
        n / (length / 1000)
      )
    ) %>%
    dplyr::filter(
      (n > nmin) & 
       (RPK > RPKmin)
    )
  
  
  if(use_ANOVA){
    
    # lkdeg_actual_anova <- kdegs %>%
    #   dplyr::group_by(sample) %>%
    #   dplyr::mutate(
    #     log_kdeg_z = (log_kdeg - mean(log_kdeg)) / sd(log_kdeg)
    #   ) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::group_by(XF, cell_line) %>%
    #   # Only keep cell lines with multiple replicates of filtered data
    #   dplyr::filter(dplyr::n() > 1) %>%
    #   dplyr::group_by(XF) %>%
    #   # Only keep genes with multiple cell lines of filtered data
    #   dplyr::filter(length(unique(cell_line)) > 1) %>%
    #   dplyr::group_modify(
    #     ~ broom::tidy(
    #       aov(log_kdeg_z ~ cell_line, data = .x)
    #     )
    #   ) %>%
    #   ungroup() %>%
    #   filter(term == "cell_line")
    # 
    # 
    # anova_table <- lkdeg_actual_anova %>%
    #   transmute(
    #     XF,
    #     Fstat = statistic,
    #     pval = p.value,
    #     FDR = p.adjust(p.value, method = "BH")
    #   )
    
    
    # Conservative, home-spun ANOVA
    # Replaces mean within-group variance with
    # max within-group variance.
    anova_table <- kdegs %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        log_kdeg_z = (log_kdeg - mean(log_kdeg)) / sd(log_kdeg)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(XF, cell_line) %>%
      # Only keep cell lines with multiple replicates of filtered data
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::group_by(XF) %>%
      # Only keep genes with multiple cell lines of filtered data
      dplyr::filter(length(unique(cell_line)) > 1) %>%
      dplyr::group_by(XF, cell_line) %>%
      dplyr::mutate(
        within_group_var = var(log_kdeg_z)
      ) %>%
      dplyr::group_by(XF) %>%
      dplyr::mutate(
        within_group_var = max(within_group_var),
        between_group_var = var(log_kdeg_z)
      ) %>%
      dplyr::summarise(
        Fstat = mean(between_group_var) / mean(within_group_var),
        ntot = dplyr::n(),
        ncell = length(unique(cell_line))
      ) %>%
      dplyr::mutate(
        pval = pf(Fstat, ncell - 1, ntot - ncell, lower.tail = FALSE)
      ) %>%
      dplyr::mutate(
        FDR = p.adjust(pval, method = "BH")
      ) %>%
      dplyr::select(
        XF, Fstat, pval, FDR
      )
    
  }
  
  
  if(use_setunion){
    
    if(consider_dynamic_range){
      
      kdegs_su <- kdegs %>%
        dplyr::mutate(
          log_kdeg = case_when(
            log_kdeg < log(-log(1 - fn_lower_limit) / tl) ~ -Inf, # below DR
            log_kdeg > log(-log(1 - fn_upper_limit) / tl) ~ Inf, # above DR
            .default = log_kdeg
          )
        )
      
      
    }else{
      kdegs_su <- kdegs
    }
    
    
    
    lkdegrange_df <- kdegs_su %>%
      dplyr::group_by(
        XF, cell_line
      ) %>%
      dplyr::summarise(
        lo = min(log_kdeg) - range_sd_factor*sd(log_kdeg),
        hi = max(log_kdeg) + range_sd_factor*sd(log_kdeg),
        nsamp = dplyr::n(),
        avg_n = mean(n),
        avg_RPK = mean(RPK)
      ) %>%
      dplyr::filter(
        nsamp > 1
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(XF) %>%
      dplyr::filter(
        dplyr::n() > 1 &
          (max(lo) - min(hi) > 0) # Don't analyze things that have no chance
      )
    
    
    XF_in_df <- unique(
      lkdegrange_df$XF
    )
    
    XF_hits <- c()
    
    for(i in seq_along(XF_in_df)){
      
      df_subset <- lkdegrange_df %>%
        dplyr::filter(
          XF == XF_in_df[i]
        )
      
      interval_gene <- Intervals(
        df_subset %>%
          dplyr::ungroup() %>%
          dplyr::select(lo, hi) %>%
          as.matrix()
      )
      
      union_gene <- interval_union(interval_gene)
      
      if(length(union_gene) > 3){
        XF_hits <- c(XF_hits, XF_in_df[i])
      }
      
    }
    
  }
  
  
  if(output_stats){
    
    genewise_stats <- kdegs %>%
      dplyr::group_by(
        XF, cell_line
      ) %>%
      dplyr::summarise(
        lo_log_kdeg = min(log_kdeg),
        hi_log_kdeg = max(log_kdeg),
        nsamp = dplyr::n(),
        avg_n = mean(n),
        avg_RPK = mean(RPK)
      ) %>%
      dplyr::group_by(XF) %>%
      dplyr::summarise(
        lowest_hi_log_kdeg = min(hi_log_kdeg),
        highest_lo_log_kdeg = max(lo_log_kdeg),
        lowest_hi_cell_line = cell_line[hi_log_kdeg == min(hi_log_kdeg)],
        highest_lo_cell_line = cell_line[lo_log_kdeg == max(lo_log_kdeg)],
        med_avg_n = median(avg_n),
        med_avg_RPK = median(avg_RPK),
        min_avg_n = min(avg_n),
        min_avg_RPK = min(avg_RPK),
        max_avg_n = max(avg_n),
        max_avg_RPK = max(avg_RPK)
      )
    
  }
  
  
  if(use_adhoc_filters){
    
    XF_to_consider <- kdegs %>%
      dplyr::group_by(
        XF
      ) %>%
      dplyr::filter(
        !all( (log(2) / exp(log_kdeg)) < tl / 1 ) &
          !all( (log(2) / exp(log_kdeg)) < tl * 2 )
      )
    
  }
  
  outlist <- list()
  
  if(use_ANOVA){
    
    outlist[["ANOVA_res"]] <- anova_table
    
  }
  
  if(use_setunion){
    
    outlist[["SetOverlap_hits"]] <- XF_hits
    
  }
  
  if(output_stats){
    
    outlist[["Genewise_stats"]] <- genewise_stats
  }
      
  return(outlist)
  
}



##### Filter database #####

filter_database <- function(database,
                            filter_chrY = TRUE,
                            filter_rRNA = TRUE,
                            filter_chrM = TRUE,
                            filter_short = TRUE,
                            min_length = 100,
                            gtf_file = "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Genomes/Human_strict/hg38_Strict.gtf",
                            exon_file = "C:/Users/isaac/Documents/Simon_Lab/kdeg_curation/Data/GeneInfo/Exonic_sequence_stats.csv"){
  

  
  ### Load necessary input
  
  if(any(c(filter_chrY, filter_rRNA, filter_chrM))){
    gtf <- rtracklayer::import(gtf_file)
  }


  if(filter_short){
    exonic_lengths <- read_csv(exon_file)
  }
  
  
  if(filter_chrY){
    
    # Filter out Y-chr RNA (since not all cell types will have these)
    Y_XF <- gtf %>%
      as_tibble() %>%
      dplyr::filter(
        seqnames == "chrY"
      ) %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct() %>%
      unlist() %>%
      unname()
    
    database_sub <- database_sub %>%
      dplyr::filter(
        !(XF %in% Y_XF)
      )
    
  }
  

  
  
  if(filter_rRNA){
    
    # Filter out rRNA
    rRNA_XF <- gtf %>%
      as_tibble() %>%
      dplyr::filter(gene_biotype == "rRNA") %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct() %>%
      unlist() %>% unname()
    
    
    database_sub <- database_sub %>%
      dplyr::filter(
        !(XF %in% rRNA_XF)
      )
    
  }
  
  
  if(filter_chrM){
    
    M_XF <- gtf %>%
      as_tibble() %>%
      dplyr::filter(
        seqnames == "chrM"
      ) %>%
      dplyr::select(gene_id) %>%
      dplyr::distinct() %>%
      unlist() %>%
      unname()
    
    database_sub <- database_sub %>%
      dplyr::filter(
        !(XF %in% M_XF)
      )
    
    
  }
  

  
  if(filter_short){
    
    # Filter out short RNA
    short_XF <- exonic_lengths %>%
      dplyr::filter(length < min_length) %>%
      dplyr::select(
        XF
      ) %>%
      dplyr::distinct() %>%
      unlist() %>%
      unname()
    
    
    database_sub <- database_sub %>%
      dplyr::filter(
        !(XF %in% short_XF)
      )
    
  }


  
  return(database_sub)
  
}


##### Normalize for dropout #####

dropout_normalization <- function(
    kdegs,
    metadf,
    normalize_across_tls = TRUE,
    read_cutoff = 50,
    grouping_factors = NULL
    ){

  ### GENERAL STEPS:
  # 1) Get -s4U RPMs
  # 2) Get +s4U RPMs
  # 3) Calculate dropout (+s4U RPM)/(-s4U RPM)
  # 4) Fit dropout vs. estimated fraction new trend
  # 5) Correct fraction news and read counts accordingly.
  # 6) Return corrected EZbakRFractions object
  
  
  ### Figure out which fraction new estimates to use
  
  # Get fractions
  fractions <- kdegs %>%
    dplyr::rename(tl = tchase) %>%
    dplyr::mutate(fraction_highTC = 1 - exp(-kdeg*tl),
                  logit_fraction_highTC = EZbakR:::logit(fraction_highTC),
                  se_logit_fraction_highTC = se_log_kdeg * 
                    ((1 / fraction_highTC) + (1 / (1-fraction_highTC))) *
                    (1 - fraction_highTC)*
                    tl*kdeg) %>%
    dplyr::select(
      sample,
      XF,
      fraction_highTC,
      logit_fraction_highTC,
      se_logit_fraction_highTC,
      n
    )
  
  features_to_analyze <- "XF"
  
  ### Figure out column names to operate on
  
  fraction_cols <- colnames(fractions)
  
  fraction_of_interest <- fraction_cols[grepl("^fraction_high", fraction_cols)]
  
  # Other columns I will have to reference
  logit_fraction <- paste0("logit_", fraction_of_interest)
  logit_se <- paste0("se_", logit_fraction)
  
  
  ### Which columns should -s4U samples be grouped by?
  
  
  # Which column is label time column?
  # Could generalize to pulse-chases at some point
  tl_col <- "tl"
  
  
  if(!normalize_across_tls){
    grouping_factors <- unique(c(grouping_factors, tl_col))
  }
  
  
  
  ### Figure out which sample likely has lowest dropout in each group
  
  # Reference sample in each group is that which has the lowest
  # dropout and thus that other samples in that group will be
  # normalized with respect to
  reference_sample_info <- fractions %>%
    dplyr::inner_join(
      metadf,
      by = "sample"
    ) %>%
    dplyr::filter(tl > 0 & n > read_cutoff) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", grouping_factors)))) %>%
    dplyr::summarise(
      med_fn = ifelse(
        normalize_across_tls,
        stats::median(-log(1 - !!dplyr::sym(fraction_of_interest))/tl ), # Inferred kdeg
        stats::median(!!dplyr::sym(logit_fraction))
      )
    ) %>%
    dplyr::ungroup() %>%
    {if (length(grouping_factors) > 0) dplyr::group_by(., dplyr::across(dplyr::all_of(grouping_factors))) else .} %>%
    dplyr::mutate(
      normalization_reference = dplyr::case_when(
        med_fn == max(med_fn) ~ TRUE,
        .default = FALSE
      ),
      reference_samp = sample[med_fn == max(med_fn)],
      nsamps_in_group = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      sample, !!grouping_factors, normalization_reference, reference_samp, nsamps_in_group
    )
  
  
  ### Quantify dropout with respect to reference in each group
  
  
  fxn_wide <- fractions %>%
    dplyr::filter(n > read_cutoff & !!dplyr::sym(logit_fraction) > -Inf) %>%
    dplyr::select(sample, !!features_to_analyze, !!fraction_of_interest) %>%
    pivot_wider(
      names_from = "sample",
      values_from = all_of(fraction_of_interest)
    )
  
  pdos <- metadf %>%
    dplyr::filter(tl > 0) %>%
    dplyr::mutate(
      pdo = 0
    )
  
  for(i in seq_along(pdos$sample)){
    
    samp <- pdos$sample[i]
    
    if(reference_sample_info$normalization_reference[reference_sample_info$sample == samp]) next
    
    reference <- reference_sample_info$reference_samp[reference_sample_info$sample == samp]
    
    
    fxn_subset <- fxn_wide %>%
      dplyr::select(
        !!features_to_analyze,
        !!samp,
        !!reference
      )
    
    if(normalize_across_tls){
      
      ref_tl <- metadf[[tl_col]][metadf$sample == reference]
      samp_tl <- metadf[[tl_col]][metadf$sample == samp]
      
      # Infer reference fraction new assuming exponential decay kinetics
      fxn_subset <- fxn_subset %>%
        dplyr::mutate(
          # kdeg = -log(1 - fn)/tl
          # fn = 1 - exp(-kdeg*tl)
          !!reference := 1 - exp((log(1 - !!dplyr::sym(reference))/ref_tl) * samp_tl)
        )
      
    }
    
    pdos$pdo[i] <- fxn_subset %>%
      na.omit() %>%
      summarise(
        pdo = stats::optim(
          par = c(-2),
          fn = EZbakR:::do_norm_ll,
          reffn = !!dplyr::sym(reference),
          fndo = !!dplyr::sym(samp),
          sig = 0.2, # Should be some function of the uncertainties
          method = "L-BFGS-B",
          upper = 6,
          lower = -6
        )$par[1]
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(pdo) %>%
      unlist() %>%
      unname() %>%
      EZbakR:::inv_logit()
    
    
  }
  
  ### Normalize for dropout
  fractions_normalized <- fractions %>%
    dplyr::left_join(
      pdos %>%
        dplyr::select(sample, pdo),
      by = "sample"
    ) %>%
    dplyr::mutate(
      pdo = ifelse(is.na(pdo), 0, pdo)
    ) %>%
    dplyr::mutate(
      !!fraction_of_interest := (!!dplyr::sym(fraction_of_interest))/((1 - pdo) + !!dplyr::sym(fraction_of_interest) * pdo),
    ) %>%
    dplyr::mutate(
      !!logit_fraction := EZbakR:::logit(!!dplyr::sym(fraction_of_interest))
    ) %>%
    dplyr::group_by(
      sample
    ) %>%
    dplyr::mutate(
      num = mean(!!dplyr::sym(fraction_of_interest))*(1-pdo) + (1-mean(!!dplyr::sym(fraction_of_interest))),
      den = !!dplyr::sym(fraction_of_interest)*(1-pdo) + (1 - !!dplyr::sym(fraction_of_interest)),
      n = ceiling(n*(num/den))
    ) %>%
    dplyr::select(
      -num, -den
    )
  
  return(fractions_normalized)
  
}