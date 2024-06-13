# Run null model for diversity heritability
# Use aov not NB
# Do 100 random bootstraps

# Loop through the NB model and etasq calc for each ASV
run_div_null_model <- function(metric, df_df, input){
  for (i in 1:100) {
    # Set results dataframe
    results_null[[i]] <- as.data.frame(matrix(NA, nrow = ncol(df_df), ncol = 7)) %>%
      set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
                  "GenotypePaov", "GenotypePnb"))
    
    # Randomly assign the samples to 95 groups by replicate field, 1 to 2 samples each
    rand_r1 <- input$map_loaded %>%
      filter(rep == "1MI") %>%
      mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                               rep(seq(from = 6, to = 95, by = 1), 2))) %>%
      mutate(pedigree_samp = sample(x = pedigree_rand,
                                    size = 185,
                                    replace = F))
    rand_r2 <- input$map_loaded %>%
      filter(rep == "2MI") %>%
      mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                               rep(seq(from = 6, to = 12, by = 1), 1), 
                               rep(seq(from = 13, to = 95, by = 1), 2))) %>%
      mutate(pedigree_samp = sample(x = pedigree_rand,
                                    size = 183,
                                    replace = F))
    rand <- rbind(rand_r1, rand_r2)
    
    # Combined dataframe
    df <- rand %>%
      dplyr::select(sampleID, rep, pedigree_samp, all_of(metric)) %>%
      mutate(pedigree_samp = as.factor(pedigree_samp))
    
    # Run the model for loop as before
    for (j in 4:ncol(df)) {
      
      # OTU name
      results_null[[i]]$OTU[j] <- names(df)[j]
      
      # Levene Test
      l2 <- leveneTest(df[,j] ~ df$pedigree_samp)
      results_null[[i]]$LeveneG[j] <- l2$`Pr(>F)`[1]
      
      # Make data frame with the predictors and each OTU
      df <- df %>%
        dplyr::select(rep, pedigree_samp, j) %>%
        mutate(pedigree_samp = as.factor(pedigree_samp)) %>%
        set_names(c("rep", "pedigree_samp", "OTU"))
      
      # Models
      m <- aov(OTU ~ rep + pedigree_samp, data = df)
      results_null[[i]]$GenotypePaov[j] <- c1$`Pr(>F)`[2]
      results_null[[i]]$GenotypePnb[j] <- c4$`Pr(>F)`[2]
      
      eta <- eta_sq(m) %>%
        as.data.frame() %>%
        rownames_to_column(var = "variable") %>%
        set_names(c("variable", "eta"))
      
      results_null[[i]]$Heritability[j] <- eta$eta[2]
      
      # Shapiro Test
      s1 <- shapiro.test(m$residuals)
      results_null[[i]]$ShapiroAOV[j] <- s1$p.value
    }
    
    # Get the mean for that bootstrap
    results_meanh$Run[i] <- i
    results_meanh$MeanH[i] <- mean(results_null[[i]]$Heritability, na.rm = T)
  }
  results_meanh <<- results_meanh
  results_null <<- results_null
}

# For Bray-Curtis dissimilarity, use PERMANOVA
# Loop through the NB model and etasq calc for each ASV
run_com_null_model <- function(dissim, input, nruns){
  
  # Track time
  t1=Sys.time()
  
  # Set results dataframe
  results_null <- as.data.frame(matrix(NA, nrow = nruns, ncol = 4)) %>%
    set_names(c("Variable", "Heritability", "PERMDISP", "GenotypePpermanova")) %>%
    mutate("Variable" = "BC")
  
  for (i in 1:nruns) {
    # Progress printout
    print(paste('Beginning Run ', i,' of ', nruns,sep=''))
    
    # Randomly assign the samples to 95 groups by replicate field, 1 to 2 samples each
    rand_r1 <- input$map_loaded %>%
      filter(rep == "1MI") %>%
      mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                               rep(seq(from = 6, to = 95, by = 1), 2))) %>%
      mutate(pedigree_samp = sample(x = pedigree_rand,
                                    size = 185,
                                    replace = F))
    rand_r2 <- input$map_loaded %>%
      filter(rep == "2MI") %>%
      mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                               rep(seq(from = 6, to = 12, by = 1), 1), 
                               rep(seq(from = 13, to = 95, by = 1), 2))) %>%
      mutate(pedigree_samp = sample(x = pedigree_rand,
                                    size = 183,
                                    replace = F))
    rand <- rbind(rand_r1, rand_r2)
    
    # Combined dataframe
    df <- rand %>%
      dplyr::select(sampleID, rep, pedigree_samp) %>%
      mutate(pedigree_samp = as.factor(pedigree_samp))
    
    # PERMDISP Test
    l2 <- anova(betadisper(dissim, df$pedigree_samp))
    results_null$PERMDISP[i] <- l2$`Pr(>F)`[1]
    
    # Models
    m <- adonis2(dissim ~ df$rep + df$pedigree_samp)
    results_null$GenotypePpermanova[i] <- m$`Pr(>F)`[2]
    
    eta <- eta_sq_adonis(m) %>%
      as.data.frame() %>%
      rownames_to_column(var = "variable") %>%
      set_names(c("variable", "eta"))
    
    results_null$Heritability[i] <- eta$eta[2]*2
    
    # Time check	
    t2 = Sys.time()
    print(t2-t1)
  }
  results_null <<- results_null
}