# Run null model for taxa heritability
# Do 100 random bootstraps

# Loop through the NB model and etasq calc for each ASV
run_taxa_null_model <- function(ts_df, df_df, input){
for (i in 1:100) {
  # Set results dataframe
  results_ITS_null[[i]] <- as.data.frame(matrix(NA, nrow = ncol(df_df), ncol = 7)) %>%
    set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
                "GenotypePaov", "GenotypePnb"))
  
  # Randomly assign the samples to 95 groups by replicate field, 1 to 2 samples each
  rand_ITS_r1 <- input$map_loaded %>%
    filter(rep == "1MI") %>%
    mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                             rep(seq(from = 6, to = 95, by = 1), 2))) %>%
    mutate(pedigree_samp = sample(x = pedigree_rand,
                                  size = 185,
                                  replace = F))
  rand_ITS_r2 <- input$map_loaded %>%
    filter(rep == "2MI") %>%
    mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                             rep(seq(from = 6, to = 12, by = 1), 1), 
                             rep(seq(from = 13, to = 95, by = 1), 2))) %>%
    mutate(pedigree_samp = sample(x = pedigree_rand,
                                  size = 183,
                                  replace = F))
  rand_ITS <- rbind(rand_ITS_r1, rand_ITS_r2)
  
  # Combined dataframe
  df_ITS <- rand_ITS %>%
    dplyr::select(sampleID, rep, pedigree_samp) %>%
    mutate(pedigree_samp = as.factor(pedigree_samp)) %>%
    left_join(., ts_df, by = "sampleID")
  
  # Run the model for loop as before
  for (j in 4:ncol(df_ITS)) {
    
    # OTU name
    results_ITS_null[[i]]$OTU[j] <- names(df_ITS)[j]
    
    # Levene Test
    l2 <- leveneTest(df_ITS[,j] ~ df_ITS$pedigree_samp)
    results_ITS_null[[i]]$LeveneG[j] <- l2$`Pr(>F)`[1]
    
    # Make data frame with the predictors and each OTU
    df <- df_ITS %>%
      dplyr::select(rep, pedigree_samp, j) %>%
      mutate(pedigree_samp = as.factor(pedigree_samp)) %>%
      set_names(c("rep", "pedigree_samp", "OTU"))
    
    # Models
    m <- aov(OTU ~ rep + pedigree_samp, data = df)
    nb3 <- MASS::glm.nb(OTU ~ rep + pedigree_samp, data = df) # Not working but at least estimates theta
    nb4 <- glm(OTU ~ rep + pedigree_samp, data = df, family = negative.binomial(nb3$theta))
    c1 <- Anova(m, type = "II", singular.ok = TRUE)
    c4 <- Anova(nb4, test.statistic = "F")
    
    results_ITS_null[[i]]$GenotypePaov[j] <- c1$`Pr(>F)`[2]
    results_ITS_null[[i]]$GenotypePnb[j] <- c4$`Pr(>F)`[2]
    
    eta <- eta_sq_glm(nb4) %>%
      as.data.frame() %>%
      rownames_to_column(var = "variable") %>%
      set_names(c("variable", "eta"))
    
    results_ITS_null[[i]]$Heritability[j] <- eta$eta[2]
    
    # Shapiro Test
    s1 <- shapiro.test(m$residuals)
    s3 <- shapiro.test(nb4$residuals)
    results_ITS_null[[i]]$ShapiroAOV[j] <- s1$p.value
    results_ITS_null[[i]]$ShapiroNB4[j] <- s3$p.value
  }
  
  # Get the mean for that bootstrap
  results_ITS_meanh$Run[i] <- i
  results_ITS_meanh$MeanH[i] <- mean(results_ITS_null[[i]]$Heritability, na.rm = T)
}
  results_ITS_meanh <<- results_ITS_meanh
  results_ITS_null <<- results_ITS_null
}