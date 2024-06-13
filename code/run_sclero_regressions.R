# Test each taxon relative abundance with Sclerotinia incidence
run_sclero_regressions <- function(tax_sum, input) {
  
  # Prepare data for a particular taxonomic level
  # Get only taxa that pass prevalence and abundance cutoff
  prev_p25 <- data.frame("ASV_ID" = rownames(tax_sum),
                             "Absent" = rowSums(tax_sum==0)) %>%
    mutate(Present = ncol(tax_sum) - Absent) %>%
    mutate(Present_Perc = Present/ncol(tax_sum)*100) %>%
    filter(Present_Perc >= 25) # 285
  abund_gen <- data.frame("ASV_ID" = rownames(tax_sum),
                              "Mean_Abund" = rowMeans(tax_sum)) %>%
    filter(Mean_Abund > 0.0001)
  ts <- tax_sum %>%
    filter(rownames(.) %in% prev_p25$ASV_ID) %>%
    filter(rownames(.) %in% abund_gen$ASV_ID) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  df_sclero <- input$map_loaded %>%
    dplyr::select(sampleID, rep, pedigree, DiseaseIncidence) %>%
    left_join(., ts, by = "sampleID")
  
  # Prepare dataframe to store results
  results_sclero <- as.data.frame(matrix(NA, nrow = ncol(df_sclero), ncol = 7)) %>%
    set_names(c("OTU", "ShapiroData", "ShapiroResid", "ScleroEst", "ScleroR2", "ScleroP", "ScleroPFDR"))
  
  # Regress all taxa against Sclerotinia resistance
  for (i in 5:ncol(df_sclero)) {
    # OTU name
    results_sclero$OTU[i] <- names(df_sclero)[i]
    
    # Make data frame with the predictors and each OTU
    df <- df_sclero %>%
      dplyr::select(rep, pedigree, DiseaseIncidence, i) %>%
      set_names(c("rep", "pedigree", "sclero", "OTU"))
    
    # Shapiro Test (data)
    s1 <- shapiro.test(df$OTU)
    results_sclero$ShapiroData[i] <- s1$p.value
    
    # Model
    m <- lm(sclero ~ OTU, data = df)
    msum <- summary(m)
    results_sclero$ScleroEst[i] <- msum$coefficients[2,1]
    results_sclero$ScleroR2[i] <- msum$r.squared
    results_sclero$ScleroP[i] <- msum$coefficients[2,4]
    
    # Shapiro Test (residuals)
    s2 <- shapiro.test(m$residuals)
    results_sclero$ShapiroResid[i] <- s2$p.value
    
  }
  
  # Format output, correct p-values
  results_sclero <- results_sclero %>%
    drop_na(OTU) %>%
    mutate(ScleroPFDR = p.adjust(ScleroP, method = "fdr")) %>%
    mutate(SymbolGeno = ifelse(ScleroPFDR <= 0.001,
                               "***",
                               ifelse(ScleroPFDR > 0.001 & ScleroPFDR <= 0.01,
                                      "**",
                                      ifelse(ScleroPFDR > 0.01 & ScleroPFDR <= 0.05,
                                             "*",
                                             "N.S."))))
  
  # Print key results for Table
  sigSclero <- results_sclero %>%
    filter(ScleroP < 0.05)
  sigSclero_pos <- sigSclero %>%
    filter(ScleroEst > 0)
  sigSclero_neg <- sigSclero %>%
    filter(ScleroEst < 0)
  print("# Sig.")
  print(nrow(sigSclero))
  print("# Pos.")
  print(nrow(sigSclero_pos))
  print("# Neg.")
  print(nrow(sigSclero_neg))
  print("Max R2")
  print(max(results_sclero$ScleroR2))
  sigScleroPfdr <- results_sclero %>%
    filter(ScleroPFDR < 0.05)
  print("# Sig. Pfdr")
  print(nrow(sigScleroPfdr))
  
  # Save results df
  return(results_sclero)
  
}