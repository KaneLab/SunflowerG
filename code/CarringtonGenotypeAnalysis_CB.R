# Sunflower Rhizosphere Soil Microbes - Carrington 2017 - Genotype Effects
# 16S and ITS from one site (Carrington)
# 95 genotypes, 4 replicates
# Follow up on Pogoda et al. 2024 to compare archaea/bacteria and fungi
# by Cliff Bueno de Mesquita, May 2024

# Outline (can navigate with the outline in the top right):
# 1. Setup
# 2. Alpha
# 3. Beta
# 4. Heritability
# 5. Sclerotinia
# 6. Networks
# 7. SNPs/Microbes



##### 1. Setup ####
# Libraries
library(plyr) # For data manipulation
library(tidyverse) # For data manipulation
library(mctoolsr) # For microbial analyses
library(vegan) # For multivariate stats
library(RColorBrewer) # For colors
library(microseq) # For fastas
library(car) # For stats
library(MASS) # For stats
library(FSA) # For SE
library(lme4) # For LMER
library(emmeans) # For TukeyHSD
library(multcomp) # For cld
library(lmerTest) # For Sattherwaite df
library(afex) # For alternative mixed model
library(ggh4x) # For plots
library(ggrepel) # For plot text
library(cowplot) # For multipanel
library(phyloseq) # For networks
library(SpiecEasi) # For networks
library(igraph) # For networks
library(rnetcarto) # For networks
library(writexl) # Export
library(naniar) # For NA management
library(pheatmap) # For heatmaps
library(iCAMP) # NTI
library(ape) # Phylogenetics
library(phyloseq) # microbial analyses, can handle trees
library(picante) # Trees
library(conditionz)
library(taxize) # version 0.9.100.1 from archive

# Repo Directory
setwd("~/Documents/GitHub/SunflowerG/")

# Functions
`%notin%` <- Negate(`%in%`)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
source("code/effectSize.R")
source("code/run_taxa_null_model.R")
source("code/run_div_null_model.R")
source("code/run_sclero_regressions.R")
MyDiamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("diamond", clip=shapes("circle")$clip,
          plot=MyDiamond, parameters=list(vertex.frame.color="white",
                                          vertex.frame.width=1))

mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)



#### _Prepare Data ####
# Prepare data. Do once. After, can just load in the _Start Here section
# Use data with updated taxonomy by Cliff
# Sclerotinia resistance
resist <- read.csv("data/ScleroResistance_CB.csv") %>%
  dplyr::select(pedigree, DiseaseIncidence)



#### __16S ####
input_16S <- load_taxa_table("data/seqtab_wTax_mctoolsr_16S.txt", 
                             "data/Ha_16S_MappingFile.txt", 
                             filter_cat = "pedigree", 
                             filter_vals = c("Blank", "LibBlank")) # n = 383

# Filter out chloroplast, mitochondria, eukaryotes, domain NA
input_filt_16S <- input_16S
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                         taxa_to_remove = "Chloroplast") # None found
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                         taxa_to_remove = "Mitochondria") # 38 taxa removed
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                         taxa_to_remove = "Eukaryota") # None found
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                         taxa_to_remove = "NA",
                                         at_spec_level = 1) # 6 taxa removed

# Good. Now check for singletons and doubletons and remove if present
singdoub_16S <- data.frame("count" = rowSums(input_filt_16S$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 636 observations
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                         taxa_IDs_to_remove = singdoub_16S$ASV) # 636 taxa removed

sort(colSums(input_filt_16S$data_loaded))
mean(colSums(input_filt_16S$data_loaded)) # 13767.36
sd(colSums(input_filt_16S$data_loaded)) # 6730.455
se(colSums(input_filt_16S$data_loaded)) # 343.9102

# Check rarefaction curve
rarecurve(t(input_filt_16S$data_loaded), step = 1000)
rarecurve(t(input_filt_16S$data_loaded), step = 1000, label = FALSE)
rc_16S <- rarecurve(t(input_filt_16S$data_loaded), step = 1000, label = FALSE, tidy = TRUE)
ggplot(rc_16S, aes(Sample, Species, group = Site)) +
  geom_line() +
  geom_vline(xintercept = 5099, colour = "red") +
  theme_bw()

# Add unrarefied ASV richness and Shannon diversity
input_filt_16S$map_loaded$rich <- specnumber(input_filt_16S$data_loaded, 
                                             MARGIN = 2)
input_filt_16S$map_loaded$shannon <- vegan::diversity(input_filt_16S$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)

# Add resistance info
sum(resist$pedigree %in% input_filt_16S$map_loaded$pedigree) # 95, good
input_filt_16S$map_loaded <- input_filt_16S$map_loaded %>%
  dplyr::select(-7, -8, -9, -sample) %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., resist, by = "pedigree")
rownames(input_filt_16S$map_loaded) <- input_filt_16S$map_loaded$sampleID

# Rarefy at 5099, drop 10 samples
set.seed(110)
input_filt_rare_16S <- single_rarefy(input_filt_16S, 5099) # n = 373
sort(colSums(input_filt_rare_16S$data_loaded))
table(input_filt_rare_16S$map_loaded$pedigree)
length(table(input_filt_rare_16S$map_loaded$pedigree)) # 97 pedigrees

# Add rarefied richness and shannon
input_filt_rare_16S$map_loaded$rich <- specnumber(input_filt_rare_16S$data_loaded, 
                                                  MARGIN = 2)
input_filt_rare_16S$map_loaded$shannon <- vegan::diversity(input_filt_rare_16S$data_loaded, 
                                                           index = "shannon", 
                                                           MARGIN = 2)

# Get pedigrees with 4 replicates (better for analysis) (n = 85)
p4_16S <- as.data.frame(table(input_filt_rare_16S$map_loaded$pedigree)) %>%
  filter(Freq == 4)
input_n4_16S <- filter_data(input_filt_rare_16S,
                            filter_cat = "pedigree",
                            keep_vals = p4_16S$Var1) # n = 340
# Get pedigrees matching ITS dataset (circle back here after running ITS part)
input_n4_16S <- filter_data(input_n4_16S,
                            filter_cat = "pedigree",
                            keep_vals = p4_ITS$Var1) # n = 332

# Get pedigrees with >= 3 replicates (lose power for those with 4, but include more pedigrees)
p3_16S <- as.data.frame(table(input_filt_rare_16S$map_loaded$pedigree)) %>%
  filter(Freq >= 3)
input_n3_16S <- filter_data(input_filt_rare_16S,
                            filter_cat = "pedigree",
                            keep_vals = p3_16S$Var1) # n = 370
# Get pedigrees matching ITS dataset (circle back here after running ITS part)
input_n3_16S <- filter_data(input_n3_16S,
                            filter_cat = "pedigree",
                            keep_vals = p3_ITS$Var1) # n = 370
input_n3_16S <- filter_data(input_n3_16S,
                            filter_cat = "sampleID",
                            keep_vals = input_n3_ITS$map_loaded$sampleID) # 368
length(unique(input_n3_16S$map_loaded$pedigree)) # 95

# Save data
saveRDS(input_filt_16S, "data/input_filt_16S.rds")
saveRDS(input_filt_rare_16S, "data/input_filt_rare_16S.rds")
saveRDS(input_n4_16S, "data/input_n4_16S.rds")
saveRDS(input_n3_16S, "data/input_n3_16S.rds")



#### __ITS ####
input_ITS <- load_taxa_table("data/seqtab_wTax_mctoolsr_ITS.txt", 
                             "data/Ha_ITS_MappingFile.txt", 
                             filter_cat = "pedigree", 
                             filter_vals = c("Blank", "LibBlank")) # n = 382

# Check that Bacteria, Archaea, Domain NA are all filtered.
input_filt_ITS <- input_ITS
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "Bacteria") # None found
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "Archaea") # None found
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "NA",
                                         at_spec_level = 1) # 1111 taxa removed

# Important for ITS: filter out all non-fungi! keep k__Fungi
# (ITS includes DNA from plants and animals)
# Check other "kingdoms" here:
tax_sum_kingdom <- summarize_taxonomy(input = input_filt_ITS, level = 1, report_higher_tax = F)
# k__fungi is >99% of reads but still, for accuracy, we need to filter.
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_keep = "k__Fungi",
                                         at_spec_level = 1) # 349 taxa removed

# Good. Now check for singletons and doubletons and remove if present
singdoub_ITS <- data.frame("count" = rowSums(input_filt_ITS$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 5 observations
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_IDs_to_remove = singdoub_ITS$ASV) # 5 taxa removed

sort(colSums(input_filt_ITS$data_loaded))
mean(colSums(input_filt_ITS$data_loaded)) # 18042.02
sd(colSums(input_filt_ITS$data_loaded)) # 7083.585
se(colSums(input_filt_ITS$data_loaded)) # 362.4278

# Check rarefaction curve
rarecurve(t(input_filt_ITS$data_loaded), step = 1000)
rarecurve(t(input_filt_ITS$data_loaded), step = 1000, label = FALSE)
rc_ITS <- rarecurve(t(input_filt_ITS$data_loaded), step = 1000, label = FALSE, tidy = TRUE)
ggplot(rc_ITS, aes(Sample, Species, group = Site)) +
  geom_line() +
  geom_vline(xintercept = 2070, colour = "red") +
  geom_vline(xintercept = 4820, colour = "red") +
  theme_bw()
# 2070 is good for some samples but not for the most rich samples
# 4820 is better and would only drop 1 sample

# Add unrarefied ASV richness and Shannon diversity
input_filt_ITS$map_loaded$rich <- specnumber(input_filt_ITS$data_loaded, 
                                             MARGIN = 2)
input_filt_ITS$map_loaded$shannon <- vegan::diversity(input_filt_ITS$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)

# Add resistance info
sum(resist$pedigree %in% input_filt_ITS$map_loaded$pedigree) # 95, good
input_filt_ITS$map_loaded <- input_filt_ITS$map_loaded %>%
  dplyr::select(-7, -8, -9, -sample) %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., resist, by = "pedigree")
rownames(input_filt_ITS$map_loaded) <- input_filt_ITS$map_loaded$sampleID

# Rarefy at 4820.
set.seed(110)
input_filt_rare_ITS <- single_rarefy(input_filt_ITS, 4820) # n = 381
sort(colSums(input_filt_rare_ITS$data_loaded))
table(input_filt_rare_ITS$map_loaded$pedigree)

# Add rarefied richness and shannon
input_filt_rare_ITS$map_loaded$rich <- specnumber(input_filt_rare_ITS$data_loaded, 
                                                  MARGIN = 2)
input_filt_rare_ITS$map_loaded$shannon <- vegan::diversity(input_filt_rare_ITS$data_loaded, 
                                                           index = "shannon", 
                                                           MARGIN = 2)

# Get pedigrees with 4 replicates and in both datasets (better for analysis) (n = 83)
p4_ITS <- as.data.frame(table(input_filt_rare_ITS$map_loaded$pedigree)) %>%
  filter(Freq == 4)
input_n4_ITS <- filter_data(input_filt_rare_ITS,
                            filter_cat = "pedigree",
                            keep_vals = p4_ITS$Var1) # n = 372
# Get pedigrees matching 16S dataset
input_n4_ITS <- filter_data(input_n4_ITS,
                            filter_cat = "pedigree",
                            keep_vals = p4_16S$Var1) # n = 332
# Rerun after running 16S part
input_n4_ITS <- filter_data(input_n4_ITS,
                            filter_cat = "pedigree",
                            keep_vals = unique(input_n4_16S$map_loaded$pedigree)) # n = 332
sum(unique(input_n4_16S$map_loaded$pedigree) %notin% unique(input_n4_ITS$map_loaded$pedigree)) # 0
sum(unique(input_n4_16S$map_loaded$pedigree) %in% unique(input_n4_ITS$map_loaded$pedigree)) # 83

# Get pedigrees with >= 3 replicates and in both datasets (n = 95)
p3_ITS <- as.data.frame(table(input_filt_rare_ITS$map_loaded$pedigree)) %>%
  filter(Freq >= 3)
input_n3_ITS <- filter_data(input_filt_rare_ITS,
                            filter_cat = "pedigree",
                            keep_vals = p3_ITS$Var1) # n = 378
# Get pedigrees matching 16S dataset
input_n3_ITS <- filter_data(input_n3_ITS,
                            filter_cat = "pedigree",
                            keep_vals = p3_16S$Var1) # n = 378
# Rerun after running 16S part
input_n3_ITS <- filter_data(input_n3_ITS,
                            filter_cat = "pedigree",
                            keep_vals = unique(input_n3_16S$map_loaded$pedigree)) # n = 378
# Match sample IDs
input_n3_ITS <- filter_data(input_n3_ITS,
                            filter_cat = "sampleID",
                            keep_vals = input_n3_16S$map_loaded$sampleID) # 368
length(unique(input_n3_ITS$map_loaded$pedigree)) # 95
sum(unique(input_n3_16S$map_loaded$pedigree) %notin% unique(input_n3_ITS$map_loaded$pedigree)) # 0
sum(unique(input_n3_16S$map_loaded$pedigree) %in% unique(input_n3_ITS$map_loaded$pedigree)) # 95



# Save data
saveRDS(input_filt_ITS, "data/input_filt_ITS.rds")
saveRDS(input_filt_rare_ITS, "data/input_filt_rare_ITS.rds")
saveRDS(input_n4_ITS, "data/input_n4_ITS.rds")
saveRDS(input_n3_ITS, "data/input_n3_ITS.rds")



# Old Start Here
input_filt_16S <- readRDS("data/input_filt_16S.rds")
input_filt_rare_16S <- readRDS("data/input_filt_rare_16S.rds")
input_n4_16S <- readRDS("data/input_n4_16S.rds")
input_n3_16S <- readRDS("data/input_n3_16S.rds")

input_filt_ITS <- readRDS("data/input_filt_ITS.rds")
input_filt_rare_ITS <- readRDS("data/input_filt_rare_ITS.rds")
input_n4_ITS <- readRDS("data/input_n4_ITS.rds")
input_n3_ITS <- readRDS("data/input_n3_ITS.rds")

input_n3_ITS$taxonomy_loaded <- input_n3_ITS$taxonomy_loaded %>%
  mutate(taxonomy1 = gsub("k__", "", taxonomy1),
         taxonomy2 = gsub("p__", "", taxonomy2),
         taxonomy3 = gsub("c__", "", taxonomy3),
         taxonomy4 = gsub("o__", "", taxonomy4),
         taxonomy5 = gsub("f__", "", taxonomy5),
         taxonomy6 = gsub("g__", "", taxonomy6),
         taxonomy7 = gsub("s__", "", taxonomy7))

# Per discussion with Nolan, use the n3 dataset, want to use all 95 genotypes!

# Also filter out ASVs present in only 1 sample
# By definition these would only be present in 1 genotype which could skew results
tax_sum_asv_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
dim(tax_sum_asv_16S)
prev_16S <- data.frame("ASV_ID" = rownames(tax_sum_asv_16S),
                       "Absent" = rowSums(tax_sum_asv_16S==0)) %>%
  mutate(Present = ncol(tax_sum_asv_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_16S)*100)
# 15 taxa present in 100%
prev_16S_n1 <- prev_16S %>%
  filter(Present == 1)
input_n3_16S <- filter_taxa_from_input(input_n3_16S,
                                       taxa_IDs_to_remove = prev_16S_n1$ASV_ID) # 346 removed

tax_sum_asv_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
dim(tax_sum_asv_ITS)
prev_ITS <- data.frame("ASV_ID" = rownames(tax_sum_asv_ITS),
                       "Absent" = rowSums(tax_sum_asv_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_asv_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_ITS)*100)
# 21 taxa present in 100%
prev_ITS_n1 <- prev_ITS %>%
  filter(Present == 1)
input_n3_ITS <- filter_taxa_from_input(input_n3_ITS,
                                       taxa_IDs_to_remove = prev_ITS_n1$ASV_ID) # 70 removed

# Save those
#saveRDS(input_n3_16S, "data/input_filt_16S_noSing.rds")
#saveRDS(input_n3_ITS, "data/input_filt_ITS_noSing.rds")

#### _Start here ####
input_n3_16S <- readRDS("data/input_filt_16S_noSing.rds")
input_n3_ITS <- readRDS("data/input_filt_ITS_noSing.rds")
input_n3_16S$map_loaded$rich <- specnumber(input_n3_16S$data_loaded, 
                                           MARGIN = 2)
input_n3_16S$map_loaded$shannon <- vegan::diversity(input_n3_16S$data_loaded, 
                                                    index = "shannon", 
                                                    MARGIN = 2)
input_n3_ITS$map_loaded$rich <- specnumber(input_n3_ITS$data_loaded, 
                                                  MARGIN = 2)
input_n3_ITS$map_loaded$shannon <- vegan::diversity(input_n3_ITS$data_loaded, 
                                                           index = "shannon", 
                                                           MARGIN = 2)

# Export taxonomy
#write.csv(input_n3_16S$taxonomy_loaded, "data/16S_taxonomy.csv", row.names = F)
#write.csv(input_n3_ITS$taxonomy_loaded, "data/ITS_taxonomy.csv", row.names = F)



#### 2. Alpha ####
#### _16S ####
# Get descriptive info
min(input_n3_16S$map_loaded$rich) # 1676
max(input_n3_16S$map_loaded$rich) # 2463
mean(input_n3_16S$map_loaded$rich) # 2099.595
se(input_n3_16S$map_loaded$rich) # 7.217861
sd(input_n3_16S$map_loaded$rich) # 139

# Sclerotinia
m <- lm(DiseaseIncidence ~ rich, data = input_n3_16S$map_loaded)
summary(m) # NS
ggplot(input_n3_16S$map_loaded, aes(rich, DiseaseIncidence*100)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  labs(x = "ASV Richness",
       y = "Sclerotinia % incidence") +
  theme_bw()

m <- lm(DiseaseIncidence ~ shannon, data = input_n3_16S$map_loaded)
summary(m) # NS
ggplot(input_n3_16S$map_loaded, aes(shannon, DiseaseIncidence*100)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  labs(x = "Shannon diversity",
       y = "Sclerotinia % incidence") +
  theme_bw()

# Test and plot
leveneTest(input_n3_16S$map_loaded$rich ~ input_n3_16S$map_loaded$pedigree) # Not homogeneous
leveneTest(input_n3_16S$map_loaded$rich ~ input_n3_16S$map_loaded$rep) # Homogeneous
m <- lmer(rich ~ pedigree + (1|rep), data = input_n3_16S$map_loaded)
summary(m)
Anova(m)
anova(m)
aovtab <- anova(m)
mnull <- lmer(rich ~ 1 + (1|rep), data = input_n3_16S$map_loaded)
anova(mnull, m)
m2 <- mixed(rich ~ pedigree + (1|rep), data = input_n3_16S$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two way ANOVA
m <- aov(rich ~ rep + pedigree, data = input_n3_16S$map_loaded)
summary(m)
Anova(m, type = "II")
h <- eta_sq(m)[2]
h
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
TukeyHSD(m)
t <- emmeans(object = m, specs = "pedigree") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(input_n3_16S$map_loaded$rich)+
           (max(input_n3_16S$map_loaded$rich)-
              min(input_n3_16S$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "16S")

leveneTest(input_n3_16S$map_loaded$shannon ~ input_n3_16S$map_loaded$pedigree) # Homogeneous
leveneTest(input_n3_16S$map_loaded$shannon ~ input_n3_16S$map_loaded$rep) # Homogeneous
m <- lmer(shannon ~ pedigree + (1|rep), data = input_n3_16S$map_loaded)
summary(m)
Anova(m)
anova(m)
mnull <- lmer(shannon ~ 1 + (1|rep), data = input_n3_16S$map_loaded)
anova(mnull, m)
m2 <- mixed(shannon ~ pedigree + (1|rep), data = input_n3_16S$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two va
m <- aov(shannon ~ rep + pedigree, data = input_n3_16S$map_loaded)
summary(m)
Anova(m, type = "II")
h1 <- eta_sq(m)[2]
h1
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
TukeyHSD(m)
t1 <- emmeans(object = m, specs = "pedigree") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(input_n3_16S$map_loaded$shannon)+
           (max(input_n3_16S$map_loaded$shannon)-
              min(input_n3_16S$map_loaded$shannon))/10,
         #y = 5650,
         Dataset = "16S")

label_df_16S <- rbind(t, t1)
alpha_long_16S <- input_n3_16S$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_16S.pdf", width = 8, height = 6)
g1 <- ggplot(alpha_long_16S, aes(reorder(pedigree, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_point(size = 2, alpha = 0.75, pch = 16) +
  # geom_text(data = label_df_16S, aes(pedigree, y, label = str_trim(.group)), 
  #           size = 3, color = "black") +
  labs(x = NULL, y = "ZOTU richness") +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
g1
dev.off()



#### _ITS ####
# Get descriptive info
min(input_n3_ITS$map_loaded$rich) # 144
max(input_n3_ITS$map_loaded$rich) # 284
mean(input_n3_ITS$map_loaded$rich) # 219.5652
se(input_n3_ITS$map_loaded$rich) # 1.09572
sd(input_n3_ITS$map_loaded$rich) # 21

# Sclerotinia
m <- lm(DiseaseIncidence ~ rich, data = input_n3_ITS$map_loaded)
summary(m) # NS
ggplot(input_n3_ITS$map_loaded, aes(rich, DiseaseIncidence*100)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  labs(x = "ASV Richness",
       y = "Sclerotinia % incidence") +
  theme_bw()

m <- lm(DiseaseIncidence ~ shannon, data = input_n3_ITS$map_loaded)
summary(m) # NS
ggplot(input_n3_ITS$map_loaded, aes(shannon, DiseaseIncidence*100)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  labs(x = "Shannon diversity",
       y = "Sclerotinia % incidence") +
  theme_bw()

# Test and plot
leveneTest(input_n3_ITS$map_loaded$rich ~ input_n3_ITS$map_loaded$pedigree) # Not quite Homogeneous
leveneTest(input_n3_ITS$map_loaded$rich ~ input_n3_ITS$map_loaded$rep) # Homogeneous
m <- lmer(rich ~ pedigree + (1|rep), data = input_n3_ITS$map_loaded)
summary(m)
Anova(m)
anova(m)
aovtab <- anova(m)
mnull <- lmer(rich ~ 1 + (1|rep), data = input_n3_ITS$map_loaded)
anova(mnull, m)
m2 <- mixed(rich ~ pedigree + (1|rep), data = input_n3_ITS$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two-way ANOVA
m <- aov(rich ~ rep + pedigree, data = input_n3_ITS$map_loaded)
summary(m)
Anova(m, type = "II") # N.S.!
h2 <- eta_sq(m)[2]
h2
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
TukeyHSD(m)
t <- emmeans(object = m, specs = "pedigree") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(input_n3_ITS$map_loaded$rich)+
           (max(input_n3_ITS$map_loaded$rich)-
              min(input_n3_ITS$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "ITS")

leveneTest(input_n3_ITS$map_loaded$shannon ~ input_n3_ITS$map_loaded$pedigree) # Homogeneous
leveneTest(input_n3_ITS$map_loaded$shannon ~ input_n3_ITS$map_loaded$rep) # Homogeneous
m <- lmer(shannon ~ pedigree + (1|rep), data = input_n3_ITS$map_loaded)
summary(m)
Anova(m)
anova(m)
mnull <- lmer(shannon ~ 1 + (1|rep), data = input_n3_ITS$map_loaded)
anova(mnull, m)
m2 <- mixed(shannon ~ pedigree + (1|rep), data = input_n3_ITS$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two-way ANOVA
m <- aov(shannon ~ rep + pedigree, data = input_n3_ITS$map_loaded)
summary(m)
Anova(m, type = "II") # Pedigree marginal
h3 <- eta_sq(m)[2]
h3
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
TukeyHSD(m)
t1 <- emmeans(object = m, specs = "pedigree") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(input_n3_ITS$map_loaded$shannon)+
           (max(input_n3_ITS$map_loaded$shannon)-
              min(input_n3_ITS$map_loaded$shannon))/10,
         #y = 5650,
         Dataset = "ITS")

label_df_ITS <- rbind(t, t1)
alpha_long_ITS <- input_n3_ITS$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_ITS.pdf", width = 8, height = 6)
g2 <- ggplot(alpha_long_ITS, aes(reorder(pedigree, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_point(size = 2, alpha = 0.75, pch = 16) +
  # geom_text(data = label_df_ITS, aes(pedigree, y, label = str_trim(.group)), 
  #           size = 3, color = "black") +
  labs(x = NULL, y = "ZOTU richness") +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
g2
dev.off()



#### _Combined ####
alpha_long_16S <- alpha_long_16S %>%
  mutate(Dataset = "16S") %>%
  dplyr::select(pedigree, value, name, Dataset)
alpha_long_ITS <- alpha_long_ITS %>%
  mutate(Dataset = "ITS") %>%
  dplyr::select(pedigree, value, name, Dataset)
ped_order <- input_n3_16S$map_loaded %>%
  group_by(pedigree) %>%
  summarise(mean_asv_rich = mean(rich)) %>%
  ungroup() %>%
  arrange(mean_asv_rich)
alpha_long <- rbind(alpha_long_16S, alpha_long_ITS) %>%
  mutate(pedigree = as.factor(pedigree)) %>%
  droplevels() %>%
  mutate(pedigree = factor(pedigree,
                           levels = ped_order$pedigree))
label_df_long <- rbind(label_df_16S, label_df_ITS)
# Don't show sig. letters. Just state heritability
label_df_long <- data.frame(x = c(75, 75, 75, 75),
                            y = c(1700, 6.6, 150, 2.25),
                            label = c("H = 0.42", "H = 0.38", "H = 0.26", "H = 0.31"),
                            plabel = c("p < 0.001", "p < 0.001", "N.S.D", "p < 0.1"),
                            px = c(41.5, 41.5, 41.5, 41.5),
                            py = c(2400, 7.3, 275, 4.25),
                            Dataset = c("16S", "16S", "ITS", "ITS"),
                            name = c("rich", "shannon", "rich", "shannon"))
facet_names <- c("rich" = "ZOTU Richness",
                 "shannon" = "Shannon Diversity",
                 "16S" = "a) Archaea/Bacteria",
                 "ITS" = "b) Fungi")
pdf("InitialFigs/Alpha_Combined.pdf", width = 8, height = 6)
g3 <- ggplot(alpha_long, aes(pedigree, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, alpha = 1, pch = 16) +
  geom_text(data = label_df_long, aes(x, y, label = label),
            size = 4, color = "black") +
  geom_text(data = label_df_long, aes(px, py, label = plabel),
            size = 4, color = "black") +
  labs(x = NULL, y = "ASV richness") +
  facet_grid2(name ~ Dataset, scales = "free", independent = "y",
              labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        #axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
g3
dev.off()

pdf("FinalFigs/Figure1.pdf", width = 8, height = 6)
g3
dev.off()

fung_alph <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rich) %>%
  rename(fun_rich = rich)
alpha_fb <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rich) %>%
  rename(bac_rich = rich) %>%
  left_join(., fung_alph, by = "sampleID")
pdf("FinalFigs/FigureS1.pdf", width = 7, height = 5)
figS1 <- ggplot(alpha_fb, aes(bac_rich, fun_rich)) +
  geom_point(size = 2, alpha = 0.5, pch = 16) +
  geom_smooth(method = "lm") +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 2, label = "R^2 == 0.09"), 
            parse = TRUE, size = 3, check_overlap = TRUE) +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 5, label = "p < 0.001"), 
            size = 3, check_overlap = TRUE) +
  labs(x = "Prokaryotic ZOTU richness", 
       y = "Fungal ZOTU richness") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
figS1
dev.off()
m <- lm(fun_rich ~ bac_rich, data = alpha_fb)
summary(m) # R2 = 0.09, p < 0.001



#### 3. Beta ####
#### _16S ####
bc_16S <- calc_dm(input_n3_16S$data_loaded)
#saveRDS(bc_16S, "data/bc_16S.rds")
set.seed(1150)
ado1 <- adonis2(bc_16S ~ input_n3_16S$map_loaded$rep + input_n3_16S$map_loaded$pedigree)
ado1
# Both significant, 0.01***, 0.34***
h4 <- eta_sq_adonis(ado1)[2]*2
h4

set.seed(1150)
ado2 <- adonis2(bc_16S ~ input_n3_16S$map_loaded$pedigree + input_n3_16S$map_loaded$rep)
ado2
# Both significant, 0.34***, 0.01***, no effect of order.
anova(betadisper(bc_16S, input_n3_16S$map_loaded$rep)) # Dispersion homogeneous
anova(betadisper(bc_16S, input_n3_16S$map_loaded$pedigree)) # Dispersion not homogeneous

# PCoA
pcoa_16S <- cmdscale(bc_16S, k = nrow(input_n3_16S$map_loaded) - 1, eig = T)
env_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(`Chlorophyll concentration`, rich, DiseaseIncidence)
set.seed(100)
ef_16S <- envfit(pcoa_16S, env_16S, permutations = 999, na.rm = TRUE)
ef_16S # Chloro and rich sig. Sclero marginal.
ordiplot(pcoa_16S)
plot(ef_16S, p.max = 0.05, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_16S)
vec.df_16S <- as.data.frame(ef_16S$vectors$arrows*sqrt(ef_16S$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_16S,
         Dim2 = Dim2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Chlorophyll", "Richness"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_n3_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_n3_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls_16S <- ddply(input_n3_16S$map_loaded, c("rep"), find_hull) # No effect, don't show
g4 <- ggplot(input_n3_16S$map_loaded, aes(Axis01, Axis02)) +
  # geom_polygon(data = micro.hulls_16S, 
  #              aes(colour = rep),
  #              alpha = 0.1, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.25) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_16S,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", fontface = "bold") +
  geom_text(aes(x = 0.35, y = -0.175, label = "p = 0.001\nH = 0.34"), check_overlap = T) +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g4
pdf("InitialFigs/Beta_Bray_16S.pdf", width = 7, height = 5)
g4
dev.off()

# wascores
au_16S <- tax_sum_asv_16S %>%
  filter(rownames(.) %in% input_au_16S$taxonomy_loaded$taxonomy8)
vare.points_16S <- data.frame("Axis01" = input_n3_16S$map_loaded$Axis01,
                              "Axis02" = input_n3_16S$map_loaded$Axis02)
vare.wa_16S <- wascores(vare.points_16S, t(input_n3_16S$data_loaded))
vare.wa.top_16S <- as.data.frame(vare.wa_16S) %>%
  filter(rownames(.) %in% rownames(au_16S)) %>%
  mutate(PC1_abs = abs(Axis01),
         PC2_abs = abs(Axis02)) %>%
  arrange(desc(PC1_abs)) %>%
  mutate(PC1_top = c(rep("top", 10), rep("rest", nrow(.)-10))) %>%
  arrange(desc(PC2_abs)) %>%
  mutate(PC2_top = c(rep("top", 10), rep("rest", nrow(.)-10))) %>%
  filter(PC1_top == "top" | PC2_top == "top") %>%
  rownames_to_column(var = "OTU") %>%
  left_join(., input_au_16S$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  mutate(tax = paste(taxonomy4, OTU, sep = "_"))
g4_wa <- ggplot(input_n3_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.25) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_16S,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", fontface = "bold") +
  geom_text_repel(data = vare.wa.top_16S,
            aes(x = Axis01, y = Axis02, label = tax),
            size = 2, color = "blue", fontface = "bold", max.overlaps = 15) +
  geom_text(aes(x = 0.35, y = -0.175, label = "p = 0.001\nH = 0.34"), check_overlap = T) +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g4_wa

# PC Scores. Genotype and Sclerotinia
m <- aov(Axis01 ~ rep + pedigree, data = input_n3_16S$map_loaded)
summary(m)
Anova(m, type = "II")
h5 <- eta_sq(m)[2]
h5
m <- lm(DiseaseIncidence ~ Axis01, data = input_n3_16S$map_loaded)
summary(m)

m <- aov(Axis02 ~ rep + pedigree, data = input_n3_16S$map_loaded)
summary(m)
Anova(m, type = "II")
h6 <- eta_sq(m)[2]
h6
m <- lm(DiseaseIncidence ~ Axis02, data = input_n3_16S$map_loaded)
summary(m)



#### _ITS ####
bc_ITS <- calc_dm(input_n3_ITS$data_loaded)
#saveRDS(bc_ITS, "data/bc_ITS.rds")
set.seed(1150)
ado3 <- adonis2(bc_ITS ~ input_n3_ITS$map_loaded$rep + input_n3_ITS$map_loaded$pedigree)
ado3
h7 <- eta_sq_adonis(ado3)[2]*2
h7
# Both significant, 0.04***, 0.34***
set.seed(1150)
ado4 <- adonis2(bc_ITS ~ input_n3_ITS$map_loaded$pedigree + input_n3_ITS$map_loaded$rep)
ado4
# Both significant, 0.34***, 0.04***, no effect of order.
anova(betadisper(bc_ITS, input_n3_ITS$map_loaded$rep)) # Dispersion homogeneous
anova(betadisper(bc_ITS, input_n3_ITS$map_loaded$pedigree)) # Dispersion not homogeneous

# PCoA
pcoa_ITS <- cmdscale(bc_ITS, k = nrow(input_n3_ITS$map_loaded) - 1, eig = T)
env_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(`Chlorophyll concentration`, rich, DiseaseIncidence)
set.seed(100)
ef_ITS <- envfit(pcoa_ITS, env_ITS, permutations = 999, na.rm = TRUE)
ef_ITS # Chloro and rich sig. Sclero not.
ordiplot(pcoa_ITS)
plot(ef_ITS, p.max = 0.05, cex = 0.5)
multiplier_ITS <- ordiArrowMul(ef_ITS)
vec.df_ITS <- as.data.frame(ef_ITS$vectors$arrows*sqrt(ef_ITS$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_ITS,
         Dim2 = Dim2 * multiplier_ITS) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_ITS$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Chlorophyll", "Richness"))
pcoaA1_ITS <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2_ITS <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_n3_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_n3_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls_ITS <- ddply(input_n3_ITS$map_loaded, c("rep"), find_hull) # No effect, don't show
g5 <- ggplot(input_n3_ITS$map_loaded, aes(Axis01, Axis02)) +
  # geom_polygon(data = micro.hulls_ITS,
  #              aes(colour = rep),
  #              alpha = 0.1, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.25) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_ITS,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", fontface = "bold") +
  geom_text(aes(x = 0.15, y = -0.21, label = "p = 0.001\nH = 0.34"), check_overlap = T) +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g5
pdf("InitialFigs/Beta_Bray_ITS.pdf", width = 7, height = 5)
g5
dev.off()

# wascores
au_ITS <- tax_sum_asv_ITS %>%
  filter(rownames(.) %in% input_au_ITS$taxonomy_loaded$taxonomy8)
vare.points_ITS <- data.frame("Axis01" = input_n3_ITS$map_loaded$Axis01,
                              "Axis02" = input_n3_ITS$map_loaded$Axis02)
vare.wa_ITS <- wascores(vare.points_ITS, t(input_n3_ITS$data_loaded))
vare.wa.top_ITS <- as.data.frame(vare.wa_ITS) %>%
  filter(rownames(.) %in% rownames(au_ITS)) %>%
  mutate(PC1_abs = abs(Axis01),
         PC2_abs = abs(Axis02)) %>%
  arrange(desc(PC1_abs)) %>%
  mutate(PC1_top = c(rep("top", 10), rep("rest", nrow(.)-10))) %>%
  arrange(desc(PC2_abs)) %>%
  mutate(PC2_top = c(rep("top", 10), rep("rest", nrow(.)-10))) %>%
  filter(PC1_top == "top" | PC2_top == "top") %>%
  rownames_to_column(var = "OTU") %>%
  left_join(., input_au_ITS$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  mutate(tax = paste(taxonomy4, OTU, sep = "_"))
g5_wa <- ggplot(input_n3_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.25) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_ITS,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", fontface = "bold") +
  geom_text_repel(data = vare.wa.top_ITS,
                  aes(x = Axis01, y = Axis02, label = tax),
                  size = 2, color = "blue", fontface = "bold", max.overlaps = 15) +
  geom_text(aes(x = 0.15, y = -0.21, label = "p = 0.001\nH = 0.34"), check_overlap = T) +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g5_wa

# PC Scores. Genotype and Sclerotinia
m <- aov(Axis01 ~ rep + pedigree, data = input_n3_ITS$map_loaded)
summary(m)
Anova(m, type = "II")
h8 <- eta_sq(m)[2]
h8
m <- lm(DiseaseIncidence ~ Axis01, data = input_n3_ITS$map_loaded)
summary(m)

m <- aov(Axis02 ~ rep + pedigree, data = input_n3_ITS$map_loaded)
summary(m)
Anova(m, type = "II")
h9 <- eta_sq(m)[2]
h9
m <- lm(DiseaseIncidence ~ Axis02, data = input_n3_ITS$map_loaded)
summary(m)



#### _Combined ####
g4 <- g4 + 
  ggtitle("a) Archaea/Bacteria") + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = -1))
g5 <- g5 + 
  ggtitle("b) Fungi") + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = -1))
pdf("InitialFigs/Beta_Bray_Combined.pdf", width = 8, height = 4)
plot_grid(g4, g5, ncol = 2)
dev.off()

g4_wa <- g4_wa + 
  ggtitle("a) Archaea/Bacteria") + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = -1))
g5_wa <- g5_wa + 
  ggtitle("b) Fungi") + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = -1))
pdf("InitialFigs/Beta_Bray_Combined_wascores.pdf", width = 10, height = 6)
plot_grid(g4_wa, g5_wa, ncol = 2)
dev.off()
# ^ This is Figure 2.

set.seed(1156)
mantel(bc_16S, bc_ITS) # NS
pdf("FinalFigs/FigureS2.pdf", width = 8, height = 6)
qplot(bc_16S, bc_ITS, geom = "point", alpha = 0.001) +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 2, label = "r == -0.03"), 
            parse = TRUE, size = 3, check_overlap = TRUE) +
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 5, label = "p = 0.95"), 
            size = 3, check_overlap = TRUE) +
  labs(x = "Prokaryotic Bray-Curtis dissimilarity",
       y = "Fungal Bray-Curtis dissimilarity") +
  #xlim(0, 1) +
  #ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()



#### 4. Heritability ####
# Calculate heritability of ASVs and high order taxonomies
# Carefully consider filtering cutoffs
# Let's use 25% of samples and 0.01% mean rel abundance
# 16S has 5099. * 0.0001 
# ITS has 4820. * 0.0001
# Make new input object labeled "au" for abundant/ubiquitious



#### _16S ####
nrow(input_n3_16S$data_loaded) # 19274
prev_16S_p25 <- prev_16S %>%
  filter(Present_Perc >= 25) # 2344
input_au_16S <- filter_taxa_from_input(input_n3_16S,
                                       filter_thresh = 5099*0.0001)
input_au_16S <- filter_taxa_from_input(input_au_16S,
                                       taxa_IDs_to_keep = prev_16S_p25$ASV_ID)
nrow(input_au_16S$data_loaded) # 1726

# Loop through the NB model and etasq calc for each ASV
tax_sum_asv_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
ts_16S <- tax_sum_asv_16S %>%
  filter(rownames(.) %in% input_au_16S$taxonomy_loaded$taxonomy8) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_au_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S$ShapiroAOV[i] <- s1$p.value
  results_16S$ShapiroNB4[i] <- s3$p.value
}

results_16S <- results_16S %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(Level = "ASV",
         Dataset = "16S")
sigGeno_16S <- results_16S %>%
  filter(GenotypePFDRaov < 0.05) # 976 p, 864 pFDR
sigGenoNB_16S <- results_16S %>%
  filter(GenotypePFDRnb < 0.05) # 1410 p, 1385 pFDR

hist(results_16S$Heritability)
ggplot(results_16S, aes(Heritability)) +
  geom_histogram()
ggplot(results_16S, aes(Heritability)) +
  geom_density() +
  theme_bw()
mean(results_16S$Heritability) # 0.3906376
se(results_16S$Heritability) # 0.002089092
mean(sigGenoNB_16S$Heritability) # 0.4182238
se(sigGenoNB_16S$Heritability) # 0.001970485
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/nrow(input_au_16S$data_loaded)*100, digits = 2) # 80.24%



#### _ITS ####
nrow(input_n3_ITS$data_loaded) # 1425
prev_ITS_p25 <- prev_ITS %>%
  filter(Present_Perc >= 25) # 255
input_au_ITS <- filter_taxa_from_input(input_n3_ITS,
                                       filter_thresh = 4820*0.0001)
input_au_ITS <- filter_taxa_from_input(input_au_ITS,
                                       taxa_IDs_to_keep = prev_ITS_p25$ASV_ID)
nrow(input_au_ITS$data_loaded) # 240

# Loop through the NB model and etasq calc for each ASV
tax_sum_asv_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
ts_ITS <- tax_sum_asv_ITS %>%
  filter(rownames(.) %in% input_au_ITS$taxonomy_loaded$taxonomy8) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_au_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS$ShapiroAOV[i] <- s1$p.value
  results_ITS$ShapiroNB4[i] <- s3$p.value
}

results_ITS <- results_ITS %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(Level = "ASV",
         Dataset = "ITS")
sigGeno_ITS <- results_ITS %>%
  filter(GenotypePFDRaov < 0.05) # 98 p, 75 pFDR
sigGenoNB_ITS <- results_ITS %>%
  filter(GenotypePFDRnb < 0.05) # 222 p, 222 pFDR

hist(results_ITS$Heritability)
mean(results_ITS$Heritability) # 0.4066793
se(results_ITS$Heritability) # 0.005194114
mean(sigGenoNB_ITS$Heritability) # 0.4171681
se(sigGenoNB_ITS$Heritability) # 0.00496342
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/nrow(input_au_ITS$data_loaded)*100, digits = 2) # 92.5%



#### _Combined ####
results_16S <- results_16S %>%
  mutate(Dataset = "16S")
results_ITS <- results_ITS %>%
  mutate(Dataset = "ITS")
results_comb <- rbind(results_16S, results_ITS)
results_sum <- results_comb %>%
  group_by(Dataset) %>%
  summarise(mean = mean(Heritability),
            se = se(Heritability))
pdf("InitialFigs/Heritability_Combined.pdf", width = 7, height = 5)
ggplot(results_comb, aes(Dataset, Heritability)) +
  geom_violin(colour = "red", draw_quantiles = T, scale = "area") +
  geom_jitter(size = 1, alpha = 0.5, width = 0.1, pch = 16) +
  geom_point(data = results_sum, aes(Dataset, mean), size = 4, colour = "blue") +
  geom_errorbar(data = results_sum, aes(x = Dataset, ymax = mean + se, ymin = mean - se),
                inherit.aes = F, width = 0.1, colour = "blue") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()



#### _Check Taxa ####
# Merge taxonomy to results data frame
# Want to see which taxa were significant
# Any AMF heritable?
results_16S_wTax <- results_16S %>%
  left_join(., input_au_16S$taxonomy_loaded, by = c("OTU" = "taxonomy8"))
num_h_16S <- results_16S_wTax %>%
  filter(GenotypePFDRnb < 0.05) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 1385) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S") %>%
  arrange(prop)
top_taxa_16S <- num_h_16S %>%
  arrange(desc(prop)) %>%
  slice_head(n = 12)
num_h_16S$taxonomy2[num_h_16S$taxonomy2 %notin% top_taxa_16S$taxonomy2] <- "Other"
num_h_16S <- num_h_16S %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = sum(num_sig)) %>%
  mutate(tot = 1385) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S",
         Subset = "Sig.") %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$taxonomy2))))
g6 <- ggplot(num_h_16S, aes(x = Dataset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(12, "Paired")[11:1])) +
  labs(x = NULL,
       y = "% Heritable ZOTUs",
       fill = "Phylum") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"))
g6

View(input_au_ITS$taxonomy_loaded)
results_ITS_wTax <- results_ITS %>%
  left_join(., input_au_ITS$taxonomy_loaded, by = c("OTU" = "taxonomy8"))
num_h_ITS <- results_ITS_wTax %>%
  filter(GenotypePFDRnb < 0.05) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 222) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "ITS") %>%
  arrange(prop) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = taxonomy2)) %>%
  mutate(Subset = "Sig.")
g7 <- ggplot(num_h_ITS, aes(x = Dataset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(7, "Paired")) +
  labs(x = NULL,
       y = "% Heritable ZOTUs",
       fill = "Phylum") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"))
g7

# Check if this is just related to the starting number of taxa in each group
tot_count_16S <- input_n3_16S$taxonomy_loaded %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% top_taxa_16S$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_16S$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = nrow(input_n3_16S$taxonomy_loaded)) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S",
         Subset = "All") %>%
  arrange(desc(prop))
tot_count_ITS <- input_n3_ITS$taxonomy_loaded %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% as.character(num_h_ITS$taxonomy2),
                            "Other",
                            taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = sum(num_sig)) %>%
  mutate(tot = nrow(input_n3_ITS$taxonomy_loaded)) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "ITS",
         Subset = "All") %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = c("Other", as.character(num_h_ITS$taxonomy2)))) %>%
  arrange(desc(prop))
tested_count_16S <- input_au_16S$taxonomy_loaded %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% top_taxa_16S$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_16S$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = nrow(input_au_16S$taxonomy_loaded)) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S",
         Subset = "Tested")
tested_count_ITS <- input_au_ITS$taxonomy_loaded %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% num_h_ITS$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_ITS$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = nrow(input_au_ITS$taxonomy_loaded)) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "ITS",
         Subset = "Tested")

pdf("InitialFigs/Heritability_SigTaxa.pdf", width = 7, height = 5)
plot_grid(g6, g7, ncol = 1, align = "v")
dev.off()

# Plot with tot count and au count and sig count
levels_16S <- data.frame("tax" = as.character(tot_count_16S$taxonomy2)) %>%
  filter(tax != "NA") %>%
  filter(tax != "Other")
num_tax_16S <- rbind(tot_count_16S, tested_count_16S, num_h_16S) %>%
  mutate(Subset = factor(Subset, levels = c("All", "Tested", "Sig."))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(levels_16S$tax))))
nrow(input_n3_16S$taxonomy_loaded)
nrow(input_au_16S$taxonomy_loaded)
nrow(input_auh_16S$taxonomy_loaded)
n_text_16S <- data.frame(x = c("All", "Tested", "Sig."),
                         y = c(105, 105, 105),
                         label = c("n = 19274", "n = 1726", "n = 1385"))
g8 <- ggplot(num_tax_16S, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_16S, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(12, "Paired")[11:1])) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "a) Archaea/Bacteria") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g8
ITS_levels <- data.frame("tax" = as.character(tot_count_ITS$taxonomy2)) %>%
  filter(tax != "NA") %>%
  filter(tax != "Other")
num_tax_ITS <- rbind(tot_count_ITS, tested_count_ITS, num_h_ITS) %>%
  mutate(Subset = factor(Subset, levels = c("All", "Tested", "Sig."))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(ITS_levels$tax))))
nrow(input_n3_ITS$taxonomy_loaded)
nrow(input_au_ITS$taxonomy_loaded)
nrow(input_auh_ITS$taxonomy_loaded)
n_text_ITS <- data.frame(x = c("All", "Tested", "Sig."),
                         y = c(105, 105, 105),
                         label = c("n = 1425", "n = 240", "n = 222"))
g9 <- ggplot(num_tax_ITS, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_ITS, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(6, "Paired"))) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "b) Fungi") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g9

pdf("FinalFigs/Figure4.pdf", width = 7, height = 5)
plot_grid(g8, g9, ncol = 1, align = "v")
dev.off()

# Run section 7 to combine this info with SNP info for Figure 4



#### _Export ####
# Need to send Kyle data
# Send ASV by sampleID abundances table for ASVs that meet prev/abund/heritability cutoff
# Also add alpha and beta diversity metrics to the table
# 16S and ITS
input_auh_16S <- filter_taxa_from_input(input_au_16S,
                                        taxa_IDs_to_keep = sigGenoNB_16S$OTU)
nrow(input_auh_16S$data_loaded) # 1385
View(input_auh_16S$data_loaded)
map_subset_16S <- input_auh_16S$map_loaded %>%
  dplyr::select(sampleID, rich, shannon, Axis01, Axis02) %>%
  rename(PC1 = Axis01) %>%
  rename(PC2 = Axis02)
export_16S_copy <- as.data.frame(t(input_auh_16S$data_loaded)) %>%
  rownames_to_column(var = "sampleID") %>%
  arrange(sampleID) %>%
  separate(sampleID, into = c("Junk1", "Junk2", "Pedigree", "Rep"), remove = F) %>%
  mutate(newID = paste(Pedigree, Rep, sep = "_")) %>%
  dplyr::select(-Junk1, -Junk2, -Pedigree, -Rep) %>%
  arrange(newID) %>%
  left_join(., map_subset_16S, by = "sampleID") %>%
  dplyr::select(sampleID, newID, rich, shannon, PC1, PC2, everything())
export_16S_mean <- as.data.frame(t(input_auh_16S$data_loaded)) %>%
  rownames_to_column(var = "sampleID") %>%
  arrange(sampleID) %>%
  separate(sampleID, into = c("Junk1", "Junk2", "Pedigree", "Rep"), remove = F) %>%
  left_join(., map_subset_16S, by = "sampleID") %>%
  dplyr::select(Pedigree, rich, shannon, PC1, PC2, everything()) %>%
  dplyr::select(-Junk1, -Junk2, -Rep, -sampleID) %>%
  group_by(Pedigree) %>%
  summarise_if(is.numeric, mean) %>%
  arrange(Pedigree)
export_16S_wRep <- export_16S_copy %>%
  dplyr::select(-sampleID)
  
input_auh_ITS <- filter_taxa_from_input(input_au_ITS,
                                        taxa_IDs_to_keep = sigGenoNB_ITS$OTU)
nrow(input_auh_ITS$data_loaded) # 222
View(input_auh_ITS$data_loaded)
map_subset_ITS <- input_auh_ITS$map_loaded %>%
  dplyr::select(sampleID, rich, shannon, Axis01, Axis02) %>%
  rename(PC1 = Axis01) %>%
  rename(PC2 = Axis02)
export_ITS_copy <- as.data.frame(t(input_auh_ITS$data_loaded)) %>%
  rownames_to_column(var = "sampleID") %>%
  arrange(sampleID) %>%
  separate(sampleID, into = c("Junk1", "Junk2", "Pedigree", "Rep"), remove = F) %>%
  mutate(newID = paste(Pedigree, Rep, sep = "_")) %>%
  dplyr::select(-Junk1, -Junk2, -Pedigree, -Rep) %>%
  arrange(newID) %>%
  left_join(., map_subset_ITS, by = "sampleID") %>%
  dplyr::select(sampleID, newID, rich, shannon, PC1, PC2, everything())
export_ITS_mean <- as.data.frame(t(input_auh_ITS$data_loaded)) %>%
  rownames_to_column(var = "sampleID") %>%
  arrange(sampleID) %>%
  separate(sampleID, into = c("Junk1", "Junk2", "Pedigree", "Rep"), remove = F) %>%
  left_join(., map_subset_ITS, by = "sampleID") %>%
  dplyr::select(Pedigree, rich, shannon, PC1, PC2, everything()) %>%
  dplyr::select(-Junk1, -Junk2, -Rep, -sampleID) %>%
  group_by(Pedigree) %>%
  summarise_if(is.numeric, mean) %>%
  arrange(Pedigree)
export_ITS_wRep <- export_ITS_copy %>%
  dplyr::select(-sampleID)

# Export to .txt
#write.table(export_16S_mean, "data/div_asvs_mean_16S.txt", sep = "\t", row.names = F)
#write.table(export_16S_wRep, "data/div_asvs_wReps_16S.txt", sep = "\t", row.names = F)
#write.table(export_ITS_mean, "data/div_asvs_mean_ITS.txt", sep = "\t", row.names = F)
#write.table(export_ITS_wRep, "data/div_asvs_wReps_ITS.txt", sep = "\t", row.names = F)

# Reload
export_16S_mean <- read.delim("data/div_asvs_mean_16S.txt")
input_auh_16S <- filter_taxa_from_input(input_au_16S,
                                        taxa_IDs_to_keep = names(export_16S_mean))
nrow(input_auh_16S$taxonomy_loaded)
export_ITS_mean <- read.delim("data/div_asvs_mean_ITS.txt")
input_auh_ITS <- filter_taxa_from_input(input_au_ITS,
                                        taxa_IDs_to_keep = names(export_ITS_mean))
nrow(input_auh_ITS$taxonomy_loaded)

# How many of each taxonomic level?
input_auh_16S_bac <- filter_taxa_from_input(input_auh_16S,
                                            taxa_to_keep = "Bacteria",
                                            at_spec_level = 1)
input_auh_16S_arc <- filter_taxa_from_input(input_auh_16S,
                                            taxa_to_keep = "Archaea",
                                            at_spec_level = 1)
length(unique(input_auh_16S$taxonomy_loaded$taxonomy2))
length(unique(input_auh_16S$taxonomy_loaded$taxonomy3))
length(unique(input_auh_16S$taxonomy_loaded$taxonomy4))
length(unique(input_auh_16S$taxonomy_loaded$taxonomy5))
length(unique(input_auh_16S$taxonomy_loaded$taxonomy6))

length(unique(input_auh_16S_bac$taxonomy_loaded$taxonomy2))
length(unique(input_auh_16S_bac$taxonomy_loaded$taxonomy3))
length(unique(input_auh_16S_bac$taxonomy_loaded$taxonomy4))
length(unique(input_auh_16S_bac$taxonomy_loaded$taxonomy5))
length(unique(input_auh_16S_bac$taxonomy_loaded$taxonomy6))

length(unique(input_auh_16S_arc$taxonomy_loaded$taxonomy2))
length(unique(input_auh_16S_arc$taxonomy_loaded$taxonomy3))
length(unique(input_auh_16S_arc$taxonomy_loaded$taxonomy4)) # NA
length(unique(input_auh_16S_arc$taxonomy_loaded$taxonomy5)) # NA
length(unique(input_auh_16S_arc$taxonomy_loaded$taxonomy6)) # NA
View(input_auh_16S_arc$taxonomy_loaded)

length(unique(input_auh_ITS$taxonomy_loaded$taxonomy2))
length(unique(input_auh_ITS$taxonomy_loaded$taxonomy3))
length(unique(input_auh_ITS$taxonomy_loaded$taxonomy4))
length(unique(input_auh_ITS$taxonomy_loaded$taxonomy5))
length(unique(input_auh_ITS$taxonomy_loaded$taxonomy6))


#### _Higher Level ####
# Run for genus through phylum level
#### __16S ####
tax_sum_gen_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 6, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_gen_16S) # 535
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_gen_16S),
                           "Absent" = rowSums(tax_sum_gen_16S==0)) %>%
  mutate(Present = ncol(tax_sum_gen_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_gen_16S)*100) %>%
  filter(Present_Perc >= 25) # 285
abund_16S_gen <- data.frame("ASV_ID" = rownames(tax_sum_gen_16S),
                            "Mean_Abund" = rowMeans(tax_sum_gen_16S)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_16S <- tax_sum_gen_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_gen$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_gen <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S_gen$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S_gen$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S_gen$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S_gen$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S_gen$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S_gen$ShapiroAOV[i] <- s1$p.value
  results_16S_gen$ShapiroNB4[i] <- s3$p.value
}

results_16S_gen <- results_16S_gen %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_16S <- results_16S_gen %>%
  filter(GenotypePFDRaov < 0.05) # 136 p, 117 pFDR
sigGenoNB_16S <- results_16S_gen %>%
  filter(GenotypePFDRnb < 0.05) # 175 p, 169 pFDR

hist(results_16S_gen$Heritability)
mean(results_16S_gen$Heritability) # 0.361643
se(results_16S_gen$Heritability) # 0.005233573
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/(ncol(ts_16S)-1)*100, digits = 2) # 66.02%



tax_sum_fam_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 5, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_fam_16S) # 247
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_fam_16S),
                           "Absent" = rowSums(tax_sum_fam_16S==0)) %>%
  mutate(Present = ncol(tax_sum_fam_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_fam_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_fam <- data.frame("ASV_ID" = rownames(tax_sum_fam_16S),
                            "Mean_Abund" = rowMeans(tax_sum_fam_16S)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_16S <- tax_sum_fam_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_fam$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_fam <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S_fam$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S_fam$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S_fam$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S_fam$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S_fam$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S_fam$ShapiroAOV[i] <- s1$p.value
  results_16S_fam$ShapiroNB4[i] <- s3$p.value
}

results_16S_fam <- results_16S_fam %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_16S <- results_16S_fam %>%
  filter(GenotypePFDRaov < 0.05) # 86 p, 73 pFDR
sigGenoNB_16S <- results_16S_fam %>%
  filter(GenotypePFDRnb < 0.05) # 109 p, 104 pFDR

hist(results_16S_fam$Heritability)
mean(results_16S_fam$Heritability) # 0.3601051
se(results_16S_fam$Heritability) # 0.007226718
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/(ncol(ts_16S)-1)*100, digits = 2) # 67.1%



tax_sum_ord_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 4, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_ord_16S) # 194
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_ord_16S),
                           "Absent" = rowSums(tax_sum_ord_16S==0)) %>%
  mutate(Present = ncol(tax_sum_ord_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_ord_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_ord <- data.frame("ASV_ID" = rownames(tax_sum_ord_16S),
                            "Mean_Abund" = rowMeans(tax_sum_ord_16S)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_16S <- tax_sum_ord_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_ord$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_ord <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S_ord$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S_ord$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S_ord$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S_ord$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S_ord$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S_ord$ShapiroAOV[i] <- s1$p.value
  results_16S_ord$ShapiroNB4[i] <- s3$p.value
}

results_16S_ord <- results_16S_ord %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_16S <- results_16S_ord %>%
  filter(GenotypePFDRaov < 0.05) # 70 p, 63 pFDR
sigGenoNB_16S <- results_16S_ord %>%
  filter(GenotypePFDRnb < 0.05) # 90 p, 85 pFDR

hist(results_16S_ord$Heritability)
mean(results_16S_ord$Heritability) # 0.3620144
se(results_16S_ord$Heritability) # 0.008208965
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/(ncol(ts_16S)-1)*100, digits = 2) # 66.93%



tax_sum_cla_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 3, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_cla_16S) # 94
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_cla_16S),
                           "Absent" = rowSums(tax_sum_cla_16S==0)) %>%
  mutate(Present = ncol(tax_sum_cla_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_cla_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_cla <- data.frame("ASV_ID" = rownames(tax_sum_cla_16S),
                            "Mean_Abund" = rowMeans(tax_sum_cla_16S)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_16S <- tax_sum_cla_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_cla$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_cla <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S_cla$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S_cla$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S_cla$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S_cla$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S_cla$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S_cla$ShapiroAOV[i] <- s1$p.value
  results_16S_cla$ShapiroNB4[i] <- s3$p.value
}

results_16S_cla <- results_16S_cla %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_16S <- results_16S_cla %>%
  filter(GenotypePFDRaov < 0.05) # 38 p, 30 pFDR
sigGenoNB_16S <- results_16S_cla %>%
  filter(GenotypePFDRnb < 0.05) # 47 p, 45 pFDR

hist(results_16S_cla$Heritability)
mean(results_16S_cla$Heritability) # 0.3532372
se(results_16S_cla$Heritability) # 0.01008838
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/(ncol(ts_16S)-1)*100, digits = 2) # 62.5%



tax_sum_phy_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 2, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_phy_16S) # 40
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_phy_16S),
                           "Absent" = rowSums(tax_sum_phy_16S==0)) %>%
  mutate(Present = ncol(tax_sum_phy_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_phy_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_phy <- data.frame("ASV_ID" = rownames(tax_sum_phy_16S),
                            "Mean_Abund" = rowMeans(tax_sum_phy_16S)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_16S <- tax_sum_phy_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_phy$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_phy <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_16S)) {
  # OTU name
  results_16S_phy$OTU[i] <- names(df_16S)[i]
  
  # Levene Test
  l2 <- leveneTest(df_16S[,i] ~ df_16S$pedigree)
  results_16S_phy$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_16S %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_16S_phy$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_16S_phy$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_16S_phy$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_16S_phy$ShapiroAOV[i] <- s1$p.value
  results_16S_phy$ShapiroNB4[i] <- s3$p.value
}

results_16S_phy <- results_16S_phy %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_16S <- results_16S_phy %>%
  filter(GenotypePFDRaov < 0.05) # 15 p, 11 pFDR
sigGenoNB_16S <- results_16S_phy %>%
  filter(GenotypePFDRnb < 0.05) # 19 p, 16 pFDR

hist(results_16S_phy$Heritability)
mean(results_16S_phy$Heritability) # 0.3293332
se(results_16S_phy$Heritability) # 0.01646441
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_16S)/(ncol(ts_16S)-1)*100, digits = 2) # 51.61%



#### __ITS ####
tax_sum_gen_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 6, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_gen_ITS) # 289
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_gen_ITS),
                           "Absent" = rowSums(tax_sum_gen_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_gen_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_gen_ITS)*100) %>%
  filter(Present_Perc >= 25) # 285
abund_ITS_gen <- data.frame("ASV_ID" = rownames(tax_sum_gen_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_gen_ITS)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_ITS <- tax_sum_gen_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_gen$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_gen <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS_gen$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS_gen$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS_gen$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS_gen$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS_gen$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS_gen$ShapiroAOV[i] <- s1$p.value
  results_ITS_gen$ShapiroNB4[i] <- s3$p.value
}

results_ITS_gen <- results_ITS_gen %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_ITS <- results_ITS_gen %>%
  filter(GenotypePFDRaov < 0.05) # 29 p, 20 pFDR
sigGenoNB_ITS <- results_ITS_gen %>%
  filter(GenotypePFDRnb < 0.05) # 84 p, 84 pFDR

hist(results_ITS_gen$Heritability)
mean(results_ITS_gen$Heritability) # 0.3870184
se(results_ITS_gen$Heritability) # 0.008245881
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/(ncol(ts_ITS)-1)*100, digits = 2) # 87.5%



tax_sum_fam_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 5, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_fam_ITS) # 179
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_fam_ITS),
                           "Absent" = rowSums(tax_sum_fam_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_fam_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_fam_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_fam <- data.frame("ASV_ID" = rownames(tax_sum_fam_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_fam_ITS)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_ITS <- tax_sum_fam_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_fam$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_fam <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS_fam$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS_fam$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS_fam$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS_fam$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS_fam$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS_fam$ShapiroAOV[i] <- s1$p.value
  results_ITS_fam$ShapiroNB4[i] <- s3$p.value
}

results_ITS_fam <- results_ITS_fam %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_ITS <- results_ITS_fam %>%
  filter(GenotypePFDRaov < 0.05) # 28 p, 17 pFDR
sigGenoNB_ITS <- results_ITS_fam %>%
  filter(GenotypePFDRnb < 0.05) # 77 p, 77 pFDR

hist(results_ITS_fam$Heritability)
mean(results_ITS_fam$Heritability) # 0.3944223
se(results_ITS_fam$Heritability) # 0.008877199
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/(ncol(ts_ITS)-1)*100, digits = 2) # 88.51%



tax_sum_ord_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 4, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_ord_ITS) # 80
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_ord_ITS),
                           "Absent" = rowSums(tax_sum_ord_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_ord_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_ord_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_ord <- data.frame("ASV_ID" = rownames(tax_sum_ord_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_ord_ITS)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_ITS <- tax_sum_ord_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_ord$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_ord <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS_ord$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS_ord$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS_ord$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS_ord$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS_ord$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS_ord$ShapiroAOV[i] <- s1$p.value
  results_ITS_ord$ShapiroNB4[i] <- s3$p.value
}

results_ITS_ord <- results_ITS_ord %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_ITS <- results_ITS_ord %>%
  filter(GenotypePFDRaov < 0.05) # 14 p, 4 pFDR
sigGenoNB_ITS <- results_ITS_ord %>%
  filter(GenotypePFDRnb < 0.05) # 36 p, 36 pFDR

hist(results_ITS_ord$Heritability)
mean(results_ITS_ord$Heritability) # 0.3675886
se(results_ITS_ord$Heritability) # 0.01100514
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/(ncol(ts_ITS)-1)*100, digits = 2) # 83.72%



tax_sum_cla_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 3, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_cla_ITS) # 31
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_cla_ITS),
                           "Absent" = rowSums(tax_sum_cla_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_cla_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_cla_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_cla <- data.frame("ASV_ID" = rownames(tax_sum_cla_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_cla_ITS)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_ITS <- tax_sum_cla_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_cla$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_cla <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS_cla$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS_cla$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS_cla$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS_cla$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS_cla$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS_cla$ShapiroAOV[i] <- s1$p.value
  results_ITS_cla$ShapiroNB4[i] <- s3$p.value
}

results_ITS_cla <- results_ITS_cla %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_ITS <- results_ITS_cla %>%
  filter(GenotypePFDRaov < 0.05) # 5 p, 1 pFDR
sigGenoNB_ITS <- results_ITS_cla %>%
  filter(GenotypePFDRnb < 0.05) # 14 p, 14 pFDR

hist(results_ITS_cla$Heritability)
mean(results_ITS_cla$Heritability) # 0.3559972
se(results_ITS_cla$Heritability) # 0.01848159
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/(ncol(ts_ITS)-1)*100, digits = 2) # 70%



tax_sum_phy_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 2, 
                                      report_higher_tax = F,
                                      relative = T)
nrow(tax_sum_phy_ITS) # 14
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_phy_ITS),
                           "Absent" = rowSums(tax_sum_phy_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_phy_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_phy_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_phy <- data.frame("ASV_ID" = rownames(tax_sum_phy_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_phy_ITS)) %>%
  filter(Mean_Abund > 0.0001)

# Loop through the NB model and etasq calc for each ASV
ts_ITS <- tax_sum_phy_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_phy$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_phy <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df_ITS)) {
  # OTU name
  results_ITS_phy$OTU[i] <- names(df_ITS)[i]
  
  # Levene Test
  l2 <- leveneTest(df_ITS[,i] ~ df_ITS$pedigree)
  results_ITS_phy$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df_ITS %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results_ITS_phy$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results_ITS_phy$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results_ITS_phy$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results_ITS_phy$ShapiroAOV[i] <- s1$p.value
  results_ITS_phy$ShapiroNB4[i] <- s3$p.value
}

results_ITS_phy <- results_ITS_phy %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePFDRnb <= 0.001,
                             "***",
                             ifelse(GenotypePFDRnb > 0.001 & GenotypePFDRnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePFDRnb > 0.01 & GenotypePFDRnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno_ITS <- results_ITS_phy %>%
  filter(GenotypePFDRaov < 0.05) # 1 p, 0 pFDR
sigGenoNB_ITS <- results_ITS_phy %>%
  filter(GenotypePFDRnb < 0.05) # 6 p, 6 pFDR

hist(results_ITS_phy$Heritability)
mean(results_ITS_phy$Heritability) # 0.3389862
se(results_ITS_phy$Heritability) # 0.03055428
# What % of the input abund/ubiq ASVs were heritable?
round(nrow(sigGenoNB_ITS)/(ncol(ts_ITS)-1)*100, digits = 2) # 75%



#### __Plot ####
# Make multipanel plot showing heritability of all levels, color by NB Pfdr sig.
# Could do as before with more panels
# Or 16S and ITS groups bars ± SE, and state % sig
# Forgot to label the level in the dataframe
results_16S_gen <- results_16S_gen %>%
  mutate(Level = "Genus",
         Dataset = "16S")
results_16S_fam <- results_16S_fam %>%
  mutate(Level = "Family",
         Dataset = "16S")
results_16S_ord <- results_16S_ord %>%
  mutate(Level = "Order",
         Dataset = "16S")
results_16S_cla <- results_16S_cla %>%
  mutate(Level = "Class",
         Dataset = "16S")
results_16S_phy <- results_16S_phy %>%
  mutate(Level = "Phylum",
         Dataset = "16S")
results_ITS_gen <- results_ITS_gen %>%
  mutate(Level = "Genus",
         Dataset = "ITS")
results_ITS_fam <- results_ITS_fam %>%
  mutate(Level = "Family",
         Dataset = "ITS")
results_ITS_ord <- results_ITS_ord %>%
  mutate(Level = "Order",
         Dataset = "ITS")
results_ITS_cla <- results_ITS_cla %>%
  mutate(Level = "Class",
         Dataset = "ITS")
results_ITS_phy <- results_ITS_phy %>%
  mutate(Level = "Phylum",
         Dataset = "ITS")
results_comb_all <- rbind(results_16S, results_16S_gen, results_16S_fam, 
                          results_16S_ord, results_16S_cla, results_16S_phy,
                          results_ITS, results_ITS_gen, results_ITS_fam,
                          results_ITS_ord, results_ITS_cla, results_ITS_phy) %>%
  mutate(Level = factor(Level,
                        levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))) %>%
  mutate(Sig = ifelse(GenotypePFDRnb < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))
results_sum_all <- results_comb_all %>%
  group_by(Dataset, Level) %>%
  summarise(mean = mean(Heritability),
            se = se(Heritability))
test_res <- as.data.frame(matrix(NA, nrow = length(levels(results_comb_all$Level)), ncol = 3)) %>%
  set_names(c("Level", "t", "P"))
for (i in 1:6) {
  test_res$Level[i] <- levels(results_comb_all$Level)[i]
  t_df <- results_comb_all %>%
    filter(Level == levels(results_comb_all$Level)[i])
  m <- t.test(t_df$Heritability ~ t_df$Dataset)
  test_res$t[i] <- m$statistic
  test_res$P[i] <- m$p.value
}
lab_df <- data.frame(Level = c("ASV", "ASV", "Genus", "Genus", "Family", "Family",
                               "Order", "Order", "Class", "Class", "Phylum", "Phylum"),
                     Dataset = c("16S", "ITS", "16S", "ITS", "16S", "ITS",
                                 "16S", "ITS", "16S", "ITS", "16S", "ITS"),
                     label = c("a", "b", "a", "b", "a", "b", "a", "a", "a", "a", "a", "a"),
                     Heritability = c(rep(0.69, 12))) %>%
  mutate(Level = factor(Level,
                        levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum")))
pdf("InitialFigs/Heritability_Combined_All.pdf", width = 8, height = 5)
ggplot(results_comb_all, aes(Dataset, Heritability)) +
  geom_violin(colour = "black", draw_quantiles = T, scale = "area") +
  geom_jitter(size = 1, alpha = 0.5, width = 0.1, pch = 16, aes(colour = Sig)) +
  geom_point(data = results_sum_all, aes(Dataset, mean), size = 4, colour = "black") +
  geom_errorbar(data = results_sum_all, aes(x = Dataset, ymax = mean + se, ymin = mean - se),
                inherit.aes = F, width = 0.1, colour = "black") +
  geom_text(data = lab_df, aes(label = label)) +
  facet_wrap(~ Level, ncol = 6) +
  labs(colour = "GLM") +
  scale_colour_manual(values = c("#619CFF", "#F8766D")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()

# Save. Have handy for figure tweaking.
#saveRDS(results_comb_all, "data/results_comb_all.rds")
results_comb_all <- readRDS("data/results_comb_all.rds")
results_comb_all <- results_comb_all %>%
  mutate(Level = gsub("ASV", "ZOTU", Level)) %>%
  mutate(Level = factor(Level,
                        levels = c("ZOTU", "Genus", "Family", "Order", "Class", "Phylum"))) %>%
  droplevels()
lab_df <- data.frame(Level = c("ZOTU", "ZOTU", "Genus", "Genus", "Family", "Family",
                               "Order", "Order", "Class", "Class", "Phylum", "Phylum"),
                     Dataset = c("16S", "ITS", "16S", "ITS", "16S", "ITS",
                                 "16S", "ITS", "16S", "ITS", "16S", "ITS"),
                     label = c("a", "b", "a", "b", "a", "b", "a", "a", "a", "a", "a", "a"),
                     Heritability = c(rep(0.69, 12))) %>%
  mutate(Level = factor(Level,
                        levels = c("ZOTU", "Genus", "Family", "Order", "Class", "Phylum")))
results_sum_all <- results_comb_all %>%
  group_by(Dataset, Level) %>%
  summarise(mean = mean(Heritability),
            se = se(Heritability))

# Add histograms or density plots
bot <- ggplot(results_comb_all, aes(Heritability, colour = Dataset)) +
  geom_density(linewidth = 1.5) +
  scale_colour_manual(values = c("#FFC107", "#004D40")) +
  labs(x = "Heritability",
       y = "Density") +
  facet_wrap(~ Level, ncol = 6) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 9))
bot

# Combined
top <- ggplot(results_comb_all, aes(Dataset, Heritability)) +
  geom_violin(colour = "black", draw_quantiles = T, scale = "area") +
  geom_jitter(size = 1, alpha = 0.5, width = 0.1, pch = 16, aes(colour = Sig)) +
  geom_point(data = results_sum_all, aes(Dataset, mean), size = 4, colour = "black") +
  geom_errorbar(data = results_sum_all, aes(x = Dataset, ymax = mean + se, ymin = mean - se),
                inherit.aes = F, width = 0.1, colour = "black") +
  geom_text(data = lab_df, aes(label = label)) +
  facet_wrap(~ Level, ncol = 6) +
  labs(colour = "GLM") +
  scale_colour_manual(values = c("#619CFF", "#F8766D")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 9))
top_leg <- get_legend(top)
top <- top + theme(legend.position = "none")
bot_leg <- get_legend(bot)
bot <- bot + theme(legend.position = "none")
leg <- plot_grid(top_leg, bot_leg, ncol = 1)
leg
plot <- plot_grid(top, bot, ncol = 1, align = "v", rel_heights = c(0.54, 0.46))
plot
pdf("InitialFigs/Heritability_Combined_All_wDens.pdf", width = 8, height = 6)
plot_grid(plot, leg, ncol = 2, rel_widths = c(0.88, 0.12))
dev.off()
pdf("FinalFigs/Figure3.pdf", width = 8, height = 6)
plot_grid(plot, leg, ncol = 2, rel_widths = c(0.88, 0.12))
dev.off()


#### _NULL ####
# Need to test a null model and bootstrap it
# What if we randomly assigned the samples into 95 different groups 100 times?
# Need to also maintain the replicate structure
# Check if this number is lower than ~0.36
# If not then we have a statistical artifact problem due to the large number of genotypes

#### __16S ####
ped_count <- as.data.frame(table(input_n3_16S$map_loaded$pedigree)) # 12 n = 3, 83 n = 4
r1 <- subset(input_n3_16S$map_loaded, rep == "1MI")
r2 <- subset(input_n3_16S$map_loaded, rep == "2MI") # 2 less than R1
ped_count_r1 <- as.data.frame(table(r1$pedigree)) # 5 n = 1, 90 n = 2
ped_count_r2 <- as.data.frame(table(r2$pedigree)) # 7 n = 1, 87 n = 2
rand_16S_r1 <- input_n3_16S$map_loaded %>%
  filter(rep == "1MI") %>%
  mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                           rep(seq(from = 6, to = 95, by = 1), 2)))
rand_16S_r2 <- input_n3_16S$map_loaded %>%
  filter(rep == "2MI") %>%
  mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                           rep(seq(from = 6, to = 12, by = 1), 1), 
                           rep(seq(from = 13, to = 95, by = 1), 2)))
rand_16S <- rbind(rand_16S_r1, rand_16S_r2)
table(rand_16S$pedigree_rand)
nrow(input_n3_16S$data_loaded) # 19274
prev_16S_p25 <- prev_16S %>%
  filter(Present_Perc >= 25) # 2344
input_au_16S <- filter_taxa_from_input(input_n3_16S,
                                       filter_thresh = 5099*0.0001)
input_au_16S <- filter_taxa_from_input(input_au_16S,
                                       taxa_IDs_to_keep = prev_16S_p25$ASV_ID)
nrow(input_au_16S$data_loaded) # 1726

# Loop through the NB model and etasq calc for each ASV
tax_sum_asv_16S <- summarize_taxonomy(input_n3_16S, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
ts_16S <- tax_sum_asv_16S %>%
  filter(rownames(.) %in% input_au_16S$taxonomy_loaded$taxonomy8) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- rand_16S %>%
  dplyr::select(sampleID, rep, pedigree_rand) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_null <- list() 
results_16S_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
for (i in 1:100) {
  # Set results dataframe
  results_16S_null[[i]] <- as.data.frame(matrix(NA, nrow = ncol(df_16S), ncol = 7)) %>%
    set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
                "GenotypePaov", "GenotypePnb"))
  
  # Randomly assign the samples to 95 groups by replicate field, 1 to 2 samples each
  rand_16S_r1 <- input_n3_16S$map_loaded %>%
    filter(rep == "1MI") %>%
    mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                             rep(seq(from = 6, to = 95, by = 1), 2))) %>%
    mutate(pedigree_samp = sample(x = pedigree_rand,
                                  size = 185,
                                  replace = F))
  rand_16S_r2 <- input_n3_16S$map_loaded %>%
    filter(rep == "2MI") %>%
    mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                             rep(seq(from = 6, to = 12, by = 1), 1), 
                             rep(seq(from = 13, to = 95, by = 1), 2))) %>%
    mutate(pedigree_samp = sample(x = pedigree_rand,
                                  size = 183,
                                  replace = F))
  rand_16S <- rbind(rand_16S_r1, rand_16S_r2)
  
  # Combined dataframe
  df_16S <- rand_16S %>%
    dplyr::select(sampleID, rep, pedigree_samp) %>%
    mutate(pedigree_samp = as.factor(pedigree_samp)) %>%
    left_join(., ts_16S, by = "sampleID")
  
  # Run the model for loop as before
  for (j in 4:ncol(df_16S)) {
    
    # OTU name
    results_16S_null[[i]]$OTU[j] <- names(df_16S)[j]
    
    # Levene Test
    l2 <- leveneTest(df_16S[,j] ~ df_16S$pedigree_samp)
    results_16S_null[[i]]$LeveneG[j] <- l2$`Pr(>F)`[1]
    
    # Make data frame with the predictors and each OTU
    df <- df_16S %>%
      dplyr::select(rep, pedigree_samp, j) %>%
      mutate(pedigree_samp = as.factor(pedigree_samp)) %>%
      set_names(c("rep", "pedigree_samp", "OTU"))
    
    # Models
    m <- aov(OTU ~ rep + pedigree_samp, data = df)
    nb3 <- MASS::glm.nb(OTU ~ rep + pedigree_samp, data = df) # Not working but at least estimates theta
    nb4 <- glm(OTU ~ rep + pedigree_samp, data = df, family = negative.binomial(nb3$theta))
    c1 <- Anova(m, type = "II", singular.ok = TRUE)
    c4 <- Anova(nb4, test.statistic = "F")
    
    results_16S_null[[i]]$GenotypePaov[j] <- c1$`Pr(>F)`[2]
    results_16S_null[[i]]$GenotypePnb[j] <- c4$`Pr(>F)`[2]
    
    eta <- eta_sq_glm(nb4) %>%
      as.data.frame() %>%
      rownames_to_column(var = "variable") %>%
      set_names(c("variable", "eta"))
    
    results_16S_null[[i]]$Heritability[j] <- eta$eta[2]
    
    # Shapiro Test
    s1 <- shapiro.test(m$residuals)
    s3 <- shapiro.test(nb4$residuals)
    results_16S_null[[i]]$ShapiroAOV[j] <- s1$p.value
    results_16S_null[[i]]$ShapiroNB4[j] <- s3$p.value
  }
  
  # Get the mean for that bootstrap
  results_16S_meanh$Run[i] <- i
  results_16S_meanh$MeanH[i] <- mean(results_16S_null[[i]]$Heritability, na.rm = T)
}
hist(results_16S_meanh$MeanH)
mean(results_16S_meanh$MeanH) # 0.304948
se(results_16S_meanh$MeanH) # 0.0008883271
results_16S_totest <- results_16S %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 100/100 different.


### Genus
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_gen_16S),
                           "Absent" = rowSums(tax_sum_gen_16S==0)) %>%
  mutate(Present = ncol(tax_sum_gen_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_gen_16S)*100) %>%
  filter(Present_Perc >= 25) # 285
abund_16S_gen <- data.frame("ASV_ID" = rownames(tax_sum_gen_16S),
                            "Mean_Abund" = rowMeans(tax_sum_gen_16S)) %>%
  filter(Mean_Abund > 0.0001)
ts_16S <- tax_sum_gen_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_gen$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_ITS_null <- list()
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_16S, df_df = df_16S, input = input_n3_16S)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.2877438
se(results_ITS_meanh$MeanH) # 0.0006778203
results_16S_totest <- results_16S_gen %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 100/100 different.


### Family
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_fam_16S),
                           "Absent" = rowSums(tax_sum_fam_16S==0)) %>%
  mutate(Present = ncol(tax_sum_fam_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_fam_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_fam <- data.frame("ASV_ID" = rownames(tax_sum_fam_16S),
                            "Mean_Abund" = rowMeans(tax_sum_fam_16S)) %>%
  filter(Mean_Abund > 0.0001)
ts_16S <- tax_sum_fam_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_fam$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_ITS_null <- list()
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_16S, df_df = df_16S, input = input_n3_16S)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.2807439
se(results_ITS_meanh$MeanH) # 0.0008631375
results_16S_totest <- results_16S_fam %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 100/100 different.


### Order
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_ord_16S),
                           "Absent" = rowSums(tax_sum_ord_16S==0)) %>%
  mutate(Present = ncol(tax_sum_ord_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_ord_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_ord <- data.frame("ASV_ID" = rownames(tax_sum_ord_16S),
                            "Mean_Abund" = rowMeans(tax_sum_ord_16S)) %>%
  filter(Mean_Abund > 0.0001)
ts_16S <- tax_sum_ord_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_ord$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_null <- list()
results_16S_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_16S, df_df = df_16S, input = input_n3_16S)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.2787361
se(results_ITS_meanh$MeanH) # 0.001017611
results_16S_totest <- results_16S_ord %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 100/100 different.


### Class
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_cla_16S),
                           "Absent" = rowSums(tax_sum_cla_16S==0)) %>%
  mutate(Present = ncol(tax_sum_cla_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_cla_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_cla <- data.frame("ASV_ID" = rownames(tax_sum_cla_16S),
                            "Mean_Abund" = rowMeans(tax_sum_cla_16S)) %>%
  filter(Mean_Abund > 0.0001)
ts_16S <- tax_sum_cla_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_cla$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_null <- list()
results_16S_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_16S, df_df = df_16S, input = input_n3_16S)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.2738785
se(results_ITS_meanh$MeanH) # 0.001037529
results_16S_totest <- results_16S_cla %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 100/100 different.


### Phylum
prev_16S_p25 <- data.frame("ASV_ID" = rownames(tax_sum_phy_16S),
                           "Absent" = rowSums(tax_sum_phy_16S==0)) %>%
  mutate(Present = ncol(tax_sum_phy_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_phy_16S)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_16S_phy <- data.frame("ASV_ID" = rownames(tax_sum_phy_16S),
                            "Mean_Abund" = rowMeans(tax_sum_phy_16S)) %>%
  filter(Mean_Abund > 0.0001)
ts_16S <- tax_sum_phy_16S %>%
  filter(rownames(.) %in% prev_16S_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_16S_phy$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_16S, by = "sampleID")
results_16S_null <- list()
results_16S_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_16S, df_df = df_16S, input = input_n3_16S)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.2728629
se(results_ITS_meanh$MeanH) # 0.0009575493
results_16S_totest <- results_16S_phy %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_16S <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_16S$Run[i] <- i
  results_16S_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_16S <- rbind(results_16S_totest, results_16S_null_totest)
  m <- t.test(test_null_16S$Heritability ~ test_null_16S$Data)
  nulltest_16S$P[i] <- m$p.value 
}
sum(nulltest_16S$P < 0.05) # Yes. 98/100 different.


### Diversity
# Get mean ± SE for 100 runs of richness, shannon, PC1, PC2
df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, rich)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "rich", df_df = df_16S, input = input_n3_16S)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2463671
se(results_meanh$MeanH) # 0.002858537

df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, shannon)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "shannon", df_df = df_16S, input = input_n3_16S)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2355776
se(results_meanh$MeanH) # 0.003215132

df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, Axis01)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "Axis01", df_df = df_16S, input = input_n3_16S)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2525228
se(results_meanh$MeanH) # 0.003358773

df_16S <- input_n3_16S$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, Axis02)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "Axis02", df_df = df_16S, input = input_n3_16S)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2369572
se(results_meanh$MeanH) # 0.00256954



#### __ITS ####
# Prepare input data
ped_count <- as.data.frame(table(input_n3_ITS$map_loaded$pedigree)) # 12 n = 3, 83 n = 4
r1 <- subset(input_n3_ITS$map_loaded, rep == "1MI")
r2 <- subset(input_n3_ITS$map_loaded, rep == "2MI") # 2 less than R1
ped_count_r1 <- as.data.frame(table(r1$pedigree)) # 5 n = 1, 90 n = 2
ped_count_r2 <- as.data.frame(table(r2$pedigree)) # 7 n = 1, 87 n = 2
rand_ITS_r1 <- input_n3_ITS$map_loaded %>%
  filter(rep == "1MI") %>%
  mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                           rep(seq(from = 6, to = 95, by = 1), 2)))
rand_ITS_r2 <- input_n3_ITS$map_loaded %>%
  filter(rep == "2MI") %>%
  mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 2),
                           rep(seq(from = 6, to = 12, by = 1), 1), 
                           rep(seq(from = 13, to = 95, by = 1), 2)))
rand_ITS <- rbind(rand_ITS_r1, rand_ITS_r2)
table(rand_ITS$pedigree_rand)

nrow(input_n3_ITS$data_loaded) # 1425
prev_ITS_p25 <- prev_ITS %>%
  filter(Present_Perc >= 25) # 255
input_au_ITS <- filter_taxa_from_input(input_n3_ITS,
                                       filter_thresh = 4820*0.0001)
input_au_ITS <- filter_taxa_from_input(input_au_ITS,
                                       taxa_IDs_to_keep = prev_ITS_p25$ASV_ID)
nrow(input_au_ITS$data_loaded) # 240

tax_sum_asv_ITS <- summarize_taxonomy(input_n3_ITS, 
                                      level = 8, 
                                      report_higher_tax = F,
                                      relative = T)
ts_ITS <- tax_sum_asv_ITS %>%
  filter(rownames(.) %in% input_au_ITS$taxonomy_loaded$taxonomy8) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- rand_ITS %>%
  dplyr::select(sampleID, rep, pedigree_rand) %>%
  left_join(., ts_16S, by = "sampleID")
results_ITS_null <- list() 
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))

# Loop through the NB model and etasq calc for each ASV
for (i in 1:100) {
  # Set results dataframe
  results_ITS_null[[i]] <- as.data.frame(matrix(NA, nrow = ncol(df_ITS), ncol = 7)) %>%
    set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
                "GenotypePaov", "GenotypePnb"))
  
  # Randomly assign the samples to 95 groups by replicate field, 1 to 2 samples each
  rand_ITS_r1 <- input_n3_ITS$map_loaded %>%
    filter(rep == "1MI") %>%
    mutate(pedigree_rand = c(rep(seq(from = 1, to = 5, by = 1), 1), 
                             rep(seq(from = 6, to = 95, by = 1), 2))) %>%
    mutate(pedigree_samp = sample(x = pedigree_rand,
                                  size = 185,
                                  replace = F))
  rand_ITS_r2 <- input_n3_ITS$map_loaded %>%
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
    left_join(., ts_ITS, by = "sampleID")
  
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
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.345137
se(results_ITS_meanh$MeanH) # 0.0003746296
results_ITS_totest <- results_ITS %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # Yes. 100/100 different.
# Now run for all taxonomic levels.


### Genus
# Run null models
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_gen_ITS),
                           "Absent" = rowSums(tax_sum_gen_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_gen_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_gen_ITS)*100) %>%
  filter(Present_Perc >= 25) # 285
abund_ITS_gen <- data.frame("ASV_ID" = rownames(tax_sum_gen_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_gen_ITS)) %>%
  filter(Mean_Abund > 0.0001)
ts_ITS <- tax_sum_gen_ITS %>%
  filter(rownames(tax_sum_gen_ITS) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_gen$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_null <- list() 
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_ITS, df_df = df_ITS, input = input_n3_ITS)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.3419219
se(results_ITS_meanh$MeanH) # 0.0004314655

# Test if real different from random for the 100 runs
results_ITS_totest <- results_ITS_gen %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # Yes. 100/100 different.

### Family
# Run null models
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_fam_ITS),
                           "Absent" = rowSums(tax_sum_fam_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_fam_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_fam_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_fam <- data.frame("ASV_ID" = rownames(tax_sum_fam_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_fam_ITS)) %>%
  filter(Mean_Abund > 0.0001)
ts_ITS <- tax_sum_fam_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_fam$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_null <- list() 
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_ITS, df_df = df_ITS, input = input_n3_ITS)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.3525847
se(results_ITS_meanh$MeanH) # 0.0004832921

# Test if real different from random for the 100 runs
results_ITS_totest <- results_ITS_fam %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # Yes. 100/100 different.

### Order
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_ord_ITS),
                           "Absent" = rowSums(tax_sum_ord_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_ord_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_ord_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_ord <- data.frame("ASV_ID" = rownames(tax_sum_ord_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_ord_ITS)) %>%
  filter(Mean_Abund > 0.0001)
ts_ITS <- tax_sum_ord_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_ord$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_null <- list()
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_ITS, df_df = df_ITS, input = input_n3_ITS)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.3244539
se(results_ITS_meanh$MeanH) # 0.0005838461

# Test if real different from random for the 100 runs
results_ITS_totest <- results_ITS_ord %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # Sig. 99/100 different.

### Class
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_cla_ITS),
                           "Absent" = rowSums(tax_sum_cla_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_cla_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_cla_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_cla <- data.frame("ASV_ID" = rownames(tax_sum_cla_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_cla_ITS)) %>%
  filter(Mean_Abund > 0.0001)
ts_ITS <- tax_sum_cla_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_cla$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_null <- list()
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_ITS, df_df = df_ITS, input = input_n3_ITS)
View(results_ITS_null[[1]])
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.3220796
se(results_ITS_meanh$MeanH) # 0.0007941715

# Test if real different from random for the 100 runs
results_ITS_totest <- results_ITS_cla %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # 1. 99/100 not a difference from NULL!


### Phylum
prev_ITS_p25 <- data.frame("ASV_ID" = rownames(tax_sum_phy_ITS),
                           "Absent" = rowSums(tax_sum_phy_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_phy_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_phy_ITS)*100) %>%
  filter(Present_Perc >= 25) # 165
abund_ITS_phy <- data.frame("ASV_ID" = rownames(tax_sum_phy_ITS),
                            "Mean_Abund" = rowMeans(tax_sum_phy_ITS)) %>%
  filter(Mean_Abund > 0.0001)
ts_ITS <- tax_sum_phy_ITS %>%
  filter(rownames(.) %in% prev_ITS_p25$ASV_ID) %>%
  filter(rownames(.) %in% abund_ITS_phy$ASV_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts_ITS, by = "sampleID")
results_ITS_null <- list()
results_ITS_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_taxa_null_model(ts_df = ts_ITS, df_df = df_ITS, input = input_n3_ITS)
hist(results_ITS_meanh$MeanH)
mean(results_ITS_meanh$MeanH) # 0.3175758
se(results_ITS_meanh$MeanH) # 0.001204963

# Test if real different from random for the 100 runs
results_ITS_totest <- results_ITS_phy %>%
  dplyr::select(-GenotypePFDRaov, -GenotypePFDRnb, -SymbolGeno, -Level, -Dataset) %>%
  mutate(Data = "Real")
nulltest_ITS <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "P"))
for (i in 1:100) {
  nulltest_ITS$Run[i] <- i
  results_ITS_null_totest <- results_ITS_null[[i]] %>%
    mutate(Data = "Randomized")
  test_null_ITS <- rbind(results_ITS_totest, results_ITS_null_totest)
  m <- t.test(test_null_ITS$Heritability ~ test_null_ITS$Data)
  nulltest_ITS$P[i] <- m$p.value 
}
sum(nulltest_ITS$P < 0.05) # None. Never a difference from NULL!



### Diversity
# Get mean ± SE for 100 runs of richness, shannon, PC1, PC2
df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, rich)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "rich", df_df = df_ITS, input = input_n3_ITS)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2550137
se(results_meanh$MeanH) # 0.00335851

df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, shannon)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "shannon", df_df = df_ITS, input = input_n3_ITS)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2535926
se(results_meanh$MeanH) # 0.003581935

df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, Axis01)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "Axis01", df_df = df_ITS, input = input_n3_ITS)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.2514883
se(results_meanh$MeanH) # 0.003047386

df_ITS <- input_n3_ITS$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree, Axis02)
results_null <- list()
results_meanh <- as.data.frame(matrix(NA, 100, 2)) %>%
  set_names(c("Run", "MeanH"))
run_div_null_model(metric = "Axis02", df_df = df_ITS, input = input_n3_ITS)
View(results_null[[1]])
hist(results_meanh$MeanH)
mean(results_meanh$MeanH) # 0.1176695
se(results_meanh$MeanH) # 0.001593909



#### 5. Sclerotinia ####
# Correlations between all taxa levels and Sclerotinia resistance
# Not really needed here, already done by Pogoda et al., but can compare results
# Also test fungi which was not done by Pogoda et al.
# Still not the best dataset for testing this as none of the plants actually had Sclero
# Just using Sclero incidence from previous surveys
# Also check how many sig. taxa are heritable
# Wrote a function to streamline this process, check the function code for details

#### _16S ####
results_sclero_16S_asv <- run_sclero_regressions(tax_sum = tax_sum_asv_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_asv %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S$OTU)
results_sclero_16S_gen <- run_sclero_regressions(tax_sum = tax_sum_gen_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_gen %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S_gen$OTU)
results_sclero_16S_fam <- run_sclero_regressions(tax_sum = tax_sum_fam_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_fam %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S_fam$OTU)
results_sclero_16S_ord <- run_sclero_regressions(tax_sum = tax_sum_ord_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_ord %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S_ord$OTU)
results_sclero_16S_cla <- run_sclero_regressions(tax_sum = tax_sum_cla_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_cla %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S_cla$OTU)
results_sclero_16S_phy <- run_sclero_regressions(tax_sum = tax_sum_phy_16S,
                                                 input = input_n3_16S)
sigher <- results_sclero_16S_phy %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_16S_phy$OTU)



#### _ITS ####
results_sclero_ITS_asv <- run_sclero_regressions(tax_sum = tax_sum_asv_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_asv %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS$OTU)
results_sclero_ITS_gen <- run_sclero_regressions(tax_sum = tax_sum_gen_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_gen %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS_gen$OTU)
results_sclero_ITS_fam <- run_sclero_regressions(tax_sum = tax_sum_fam_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_fam %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS_fam$OTU)
results_sclero_ITS_ord <- run_sclero_regressions(tax_sum = tax_sum_ord_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_ord %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS_ord$OTU)
results_sclero_ITS_cla <- run_sclero_regressions(tax_sum = tax_sum_cla_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_cla %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS_cla$OTU)
results_sclero_ITS_phy <- run_sclero_regressions(tax_sum = tax_sum_phy_ITS,
                                                 input = input_n3_ITS)
sigher <- results_sclero_ITS_phy %>%
  filter(ScleroP < 0.05) %>%
  filter(OTU %in% results_ITS_phy$OTU)



##### 6. Networks ####
# Use SPIEC-EASI to build a network for Archaea/Bacteria and Fungi
# Again, carefully consider filtering cutoffs
# Here use the top 100 most abundant taxa for each
# Also not really needed here, not really part of genotype/heritability/sclero story
# Can compare 16S and ITS though

# Basic network plot with phyloseq.
top_16S <- as.data.frame(rowMeans(input_n3_16S$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_abund_16S <- filter_taxa_from_input(input_n3_16S,
                                          taxa_IDs_to_keep = rownames(top_16S),
                                          at_spec_level = 7)
nrow(input_abund_16S$data_loaded) # 100
otu_16S <- otu_table(input_abund_16S$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_abund_16S$taxonomy_loaded))
map_16S <- sample_data(input_abund_16S$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
se_mb_16S <- spiec.easi(input_phy_16S, 
                        method = 'mb', 
                        lambda.min.ratio = 1e-2,
                        nlambda = 20, 
                        pulsar.params = list(rep.num = 50))
net_16S <- adj2igraph(getRefit(se_mb_16S),  
                      vertex.attr=list(name = taxa_names(input_phy_16S)))
pdf("InitialFigs/Network_16S.pdf", width = 7, height = 5)
set.seed(1)
plot_network(net_16S, 
             input_phy_16S, 
             type = 'taxa', 
             color = "taxonomy2",
             shape = "taxonomy1",
             point_size = 3,
             label = NULL) +
  labs(shape = "Domain",
       color = "Phylum") +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(17,17,17,17,
                                                            16,17,17,17,
                                                            17,17,17,17)))) +
  ggtitle("16S Network (top 100)") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

View(input_n3_ITS$taxonomy_loaded)
input_n3_ITS$taxonomy_loaded <- input_n3_ITS$taxonomy_loaded %>%
  mutate(taxonomy1 = gsub("k__", "", taxonomy1),
         taxonomy2 = gsub("p__", "", taxonomy2),
         taxonomy3 = gsub("c__", "", taxonomy3),
         taxonomy4 = gsub("o__", "", taxonomy4),
         taxonomy5 = gsub("f__", "", taxonomy5),
         taxonomy6 = gsub("g__", "", taxonomy6),
         taxonomy7 = gsub("s__", "", taxonomy7))
top_ITS <- as.data.frame(rowMeans(input_n3_ITS$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_abund_ITS <- filter_taxa_from_input(input_n3_ITS,
                                          taxa_IDs_to_keep = rownames(top_ITS),
                                          at_spec_level = 7)
nrow(input_abund_ITS$data_loaded) # 100
otu_ITS <- otu_table(input_abund_ITS$data_loaded, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(input_abund_ITS$taxonomy_loaded))
map_ITS <- sample_data(input_abund_ITS$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
se_mb_ITS <- spiec.easi(input_phy_ITS, 
                        method = 'mb', 
                        lambda.min.ratio = 1e-2,
                        nlambda = 20, 
                        pulsar.params = list(rep.num = 50))
net_ITS <- adj2igraph(getRefit(se_mb_ITS),  
                      vertex.attr=list(name = taxa_names(input_phy_ITS)))
pdf("InitialFigs/Network_ITS.pdf", width = 7, height = 5)
set.seed(1)
plot_network(net_ITS, 
             input_phy_ITS, 
             type = 'taxa', 
             color = "taxonomy3",
             point_size = 3,
             label = NULL) +
  labs(color = "Class") +
  ggtitle("ITS Network (top 100)") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

# Need to create unique prok and euk ESV (ASV) IDs
# Unique ASV IDs are needed for rownames of $data_loaded and $taxonomy_loaded
# Also filter 16S data to ITS data
input_abund_16S <- filter_data(input_abund_16S,
                               filter_cat = "sampleID",
                               keep_vals = input_n3_ITS$map_loaded$sampleID)
input_abund_16S_toMerge <- input_abund_16S
num_16S <- seq(1:nrow(input_abund_16S_toMerge$data_loaded))
ASVlabel_16S <- rep("ASV_Prok_", nrow(input_abund_16S_toMerge$data_loaded))
IDs_16S <- paste(ASVlabel_16S, num_16S, sep = "")
rownames(input_abund_16S_toMerge$data_loaded) <- IDs_16S
rownames(input_abund_16S_toMerge$taxonomy_loaded) <- IDs_16S
input_abund_16S_toMerge$taxonomy_loaded$taxonomy9 <- IDs_16S

input_abund_ITS_toMerge <- input_abund_ITS
num_ITS <- seq(1:nrow(input_abund_ITS_toMerge$data_loaded))
ASVlabel_ITS <- rep("ASV_Euk_", nrow(input_abund_ITS_toMerge$data_loaded))
IDs_ITS <- paste(ASVlabel_ITS, num_ITS, sep = "")
rownames(input_abund_ITS_toMerge$data_loaded) <- IDs_ITS
rownames(input_abund_ITS_toMerge$taxonomy_loaded) <- IDs_ITS
input_abund_ITS_toMerge$taxonomy_loaded$taxonomy9 <- IDs_ITS

input_filt_abund_combined <- input_abund_16S_toMerge
input_filt_abund_combined$data_loaded <- rbind(input_filt_abund_combined$data_loaded,
                                               input_abund_ITS_toMerge$data_loaded)
input_filt_abund_combined$taxonomy_loaded <- rbind(input_filt_abund_combined$taxonomy_loaded,
                                                   input_abund_ITS_toMerge$taxonomy_loaded)
nrow(input_filt_abund_combined$data_loaded) # 200 total taxa

# Convert mctoolsr to phyloseq
names(input_filt_abund_combined$taxonomy_loaded) <- c("Domain", "Phylum", "Class", "Order",
                                                      "Family", "Genus", "Species", "ASV_ID", "ASV_ID2")
otu <- otu_table(input_filt_abund_combined$data_loaded, taxa_are_rows = T)
tax <- tax_table(as.matrix(input_filt_abund_combined$taxonomy_loaded))
map <- sample_data(input_filt_abund_combined$map_loaded)
input.phy <- phyloseq(otu, tax, map)
se.mb2 <- spiec.easi(input.phy, 
                     method='mb', 
                     lambda.min.ratio=1e-2,
                     nlambda=20,
                     pulsar.params=list(rep.num=50))
net <- adj2igraph(getRefit(se.mb2),  
                  vertex.attr=list(name=taxa_names(input.phy)))

# Phyloseq plot
plot_network(net, 
             input.phy, 
             type = 'taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL)
plot_network(net, 
             input.phy, 
             type = 'taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL,
             layout.method = layout.circle)

# Stats
E(net) # 466
V(net) # 200
transitivity(net) # Average clustering coefficient. 0.225
deg <- degree(net, mode="all")
mean(deg) # 4.66

# Add taxonomic information and color
# Add phylum, which is in the taxonomy_loaded table
# Check phylum numbers
table(input_filt_abund_combined$taxonomy_loaded$Phylum) # 18

# Check match
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match

# Make factor
input_filt_abund_combined$taxonomy_loaded$Phylum <- as.factor(input_filt_abund_combined$taxonomy_loaded$Phylum)

# Confer to network
V(net)$phylum = input_filt_abund_combined$taxonomy_loaded$Phylum

# Check levels
levels(input_filt_abund_combined$taxonomy_loaded$Phylum) # There are 18 phyla

# Set n to number of levels
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$Phylum))

# Save taxonomy and colors in tax
tax <- input_filt_abund_combined$taxonomy_loaded

# Get colors for n levels
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$Phylum
levels(tax$Phylum)
tax$color <- recode_factor(tax$color,
                           "Acidobacteriota" = colrs[1],
                           "Actinobacteriota" = colrs[2],
                           "Ascomycota" = colrs[3],
                           "Bacteroidota" = colrs[4],
                           "Basidiomycota" = colrs[5],
                           "Chloroflexi" = colrs[6],
                           "Chytridiomycota" = colrs[7],
                           "Crenarchaeota" = colrs[8],
                           "Firmicutes" = colrs[9],
                           "Gemmatimonadota" = colrs[10],
                           "Mortierellomycota" = colrs[11],
                           "Mucoromycota" = colrs[12],
                           "Myxococcota" = colrs[13],
                           "NA" = colrs[14],
                           "Nitrospirota" = colrs[15],
                           "Planctomycetota" = colrs[16],
                           "Proteobacteria" = colrs[17],
                           "Verrucomicrobiota" = colrs[18])
V(net)$color <- as.character(tax$color)

# Layout in circle
par(mar = c(0,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color, 
     vertex.size = deg*2, 
     vertex.shape = "circle", 
     vertex.frame.color = "black",
     vertex.label = NA, 
     #edge.color = ifelse(cor.matrix$r > 0, "#619CFF","#F8766D"),
     edge.curved = 0.2,
     edge.width = 0.2,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -2, y = 0.7, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)

# Refine taxonomy, shapes, color edged, other stylistic things

# Check Modules
rng_adj <- igraph::get.adjacency(net, sparse = FALSE)
netcarto(rng_adj) # 9 modules
edge.betweenness.community(net) # 28 groups
fastgreedy.community(net) # 110 groups
walktrap.community(net) # 101 groups
spinglass.community(net) # Can't work with unconnected graph
leading.eigenvector.community(net) # 70 groups
label.propagation.community(net) # 63 groups
cluster_louvain(net) # 65 groups

# Color edge by weight
se.mb2$lambda
bm <- symBeta(getOptBeta(se.mb2), mode = "maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
E(net)$weight <- weights
E(net)[weight > 0]$color <-"blue" # positive blue
E(net)[weight < 0]$color <-"red"  # negative red
edge_attr(net)

# Vertex shape by role
adj <- as.matrix(as_adjacency_matrix(net))
role <- netcarto(adj)
role <- role[[1]]
v.names <- data.frame(name = V(net)$name)
v.names.ni <- v.names %>%
  filter(name %notin% role$name) %>%
  mutate(module = NA,
         connectivity = NA,
         participation = NA,
         role = NA)
role <- rbind(role, v.names.ni) %>%
  mutate(role = as.factor(role)) %>%
  droplevels() %>%
  mutate(shape = recode_factor(role,
                               #"Peripheral Hub" = "square",
                               "Connector Hub" = "square",
                               #"Kinless Hub" = "square",
                               "Connector" = "diamond",
                               #"Kinless" = "triangle",
                               "Peripheral" = "circle",
                               "Ultra peripheral" = "circle")) %>%
  mutate(shape = as.character(shape)) %>%
  mutate(shape = replace_na(shape, "circle"))
vertex_attr(net)
V(net)$shape <- role$shape
hubcon <- role %>%
  filter(role == "Peripheral Hub" | role == "Connector Hub") %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ASV_ID2"))

# Info for labels
length(V(net))
length(E(net))
round(mean(deg), 1)
round(transitivity(net), 3)

#### _Color, Shape
# Color edges, shape by role

# Plot Circle
par(mar = c(2,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -1.7, y = 1.1, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.4, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.8, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1.5, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -2, cex = 0.8)
mtext("Edges = 466", side = 1, line = -1, cex = 0.8)
mtext("Mean Degrees = 4.7", side = 1, line = 0, cex = 0.8)
mtext("Clustering Coefficient = 0.225", side = 1, line = 1, cex = 0.8)

# Update taxonomy
# Need to label Archaea, Bacteria, Fungi (A_, B_, F_)
View(input_filt_abund_combined$taxonomy_loaded)
input_filt_abund_combined$taxonomy_loaded <- input_filt_abund_combined$taxonomy_loaded %>%
  mutate(Phylum = as.character(Phylum)) %>%
  mutate(taxonomy = paste(Domain, Phylum, sep = "_")) %>%
  mutate(taxonomy = gsub("Archaea", "A", taxonomy)) %>%
  mutate(taxonomy = gsub("Bacteria", "B", taxonomy)) %>%
  mutate(taxonomy = gsub("Fungi", "F", taxonomy))

# Run as above but with taxonomy instead of Phylum
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
input_filt_abund_combined$taxonomy_loaded$taxonomy <- as.factor(input_filt_abund_combined$taxonomy_loaded$taxonomy)
V(net)$taxonomy = input_filt_abund_combined$taxonomy_loaded$taxonomy
levels(input_filt_abund_combined$taxonomy_loaded$taxonomy) # There are 17 phyla
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$taxonomy))
tax <- input_filt_abund_combined$taxonomy_loaded
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$taxonomy
levels(tax$taxonomy)
tax$color <- recode_factor(tax$color,
                           "A_Crenarchaeota" = colrs[1],
                           "B_Acidobacteriota" = colrs[2],
                           "B_Actinobacteriota" = colrs[3],
                           "B_Bacteroidota" = colrs[4],
                           "B_Chloroflexi" = colrs[5],
                           "B_Firmicutes" = colrs[6],
                           "B_Gemmatimonadota" = colrs[7],
                           "B_Myxococcota" = colrs[8],
                           "B_Nitrospirota" = colrs[9],
                           "B_Planctomycetota" = colrs[10],
                           "B_Proteobacteria" = colrs[11],
                           "B_Verrucomicrobiota" = colrs[12],
                           "F_Ascomycota" = colrs[13],
                           "F_Basidiomycota" = colrs[14],
                           "F_Chytridiomycota" = colrs[15],
                           "F_Mortierellomycota" = colrs[16],
                           "F_Mucoromycota" = colrs[17],
                           "F_NA" = colrs[18])
V(net)$color <- as.character(tax$color)

# And add ASV ID as vertex attribute
V(net)$ASV_ID2 <- input_filt_abund_combined$taxonomy_loaded$ASV_ID2

#### _Good Taxonomy
# Replot, Circle, shape by Role
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.8, y = 1.1, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.8, y = -0.35, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.8, y = -0.7, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 466", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 4.7", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.225", side = 1, line = 1.5, cex = 0.8)

# Replot, Circle, shape by Domain. Save this one.
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
V(net)$Domain = input_filt_abund_combined$taxonomy_loaded$Domain
V(net)$Domain <- gsub("Archaea", "square", V(net)$Domain)
V(net)$Domain <- gsub("Bacteria", "circle", V(net)$Domain)
V(net)$Domain <- gsub("Fungi", "triangle", V(net)$Domain)

pdf("InitialFigs/Network_Combined.pdf", width = 8, height = 6)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.8, y = 1.1, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(22, 2), rep(21, 10), rep(24, 15)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.8, y = -0.6, c("Archaea", "Bacteria", "Fungi"), pch = c(22, 21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.8, y = -1, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.2)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 466", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 4.7", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.225", side = 1, line = 1.5, cex = 0.8)
dev.off()



#### _Betweenness
# Plot degree versus betweenness for each network
se <- se.mb2
net <- adj2igraph(getRefit(se),  vertex.attr=list(name=taxa_names(input.phy)))
bw <- data.frame("Degree" = degree(net, mode="all"),
                 "Betweenness" = betweenness(net)) %>%
  mutate("ASV" = rownames(.)) %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("ASV" = "ASV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum))

nb.cols <- length(levels(bw$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.2, 35.2),
                     breaks = seq(from = 0, to = 35, by = 1)) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid.minor.x = element_blank())

# Check ASVs, add labels
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  geom_text(data = bw,
            aes(x = Degree, y = Betweenness, label = ASV), 
            size = 3, inherit.aes = F) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.2, 35.2),
                     breaks = seq(from = 0, to = 35, by = 1)) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid.minor.x = element_blank())

bw_asvs <- bw %>%
  filter(Betweenness > 500) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          ifelse(Class != "NA",
                                                 Class,
                                                 ifelse(Phylum != "NA",
                                                        Phylum,
                                                        Domain)))))) %>%
  mutate("HightTax_ASV" = paste(HighTax, ASV_ID, sep = "_"))

pdf("InitialFigs/Network_Betweenness.pdf", width = 7, height = 5)
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  geom_text_repel(data = bw_asvs,
                  #min.segment.length = 0,
                  aes(x = Degree, y = Betweenness, label = HightTax_ASV), 
                  size = 2, inherit.aes = F,
                  position = position_jitter(width = 0.15)) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.2, 35.2),
                     breaks = seq(from = 0, to = 35, by = 1)) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(16,17,17,17,
                                                            17,17,17,17,
                                                            17,17,17,17,
                                                            15,15,15,15,
                                                            15,15)))) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid.minor.x = element_blank())
dev.off()



#### _Participation
# Plot participation coefficient and within-module degree (Barnes et al.)
# Dashed lines at 0.61 and 2.2
# z-score (within module degree) is "connectivity"
# participation coefficient P is "participation
roles <- role %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ASV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum),
         role = as.factor(role)) %>%
  filter(is.na(module) == F) %>%
  droplevels()
nb.cols <- length(levels(roles$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(16,17,17,17,
                                                            17,17,17,17,
                                                            17,17,17,15,
                                                            15,15,15,15)))) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid = element_blank())

# Check ASVs, add labels
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text(data = roles,
            aes(x = participation, y = connectivity, label = name), 
            size = 3, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(16,17,17,17,
                                                            17,17,17,17,
                                                            17,17,17,15,
                                                            15,15,15,15)))) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid = element_blank())

pz_asvs <- roles %>%
  filter(connectivity > 2.5 | participation > 0.62) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          Class))))

pdf("InitialFigs/Network_Participation.pdf", width = 7, height = 5)
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(data = pz_asvs,
                  min.segment.length = 0,
                  max.overlaps = 20,
                  aes(x = participation, y = connectivity, label = HighTax), 
                  size = 2, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(16,17,17,17,
                                                            17,17,17,17,
                                                            17,17,17,15,
                                                            15,15,15,15)))) +
  theme_bw() +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(0,0,0,0),
        panel.grid = element_blank())
dev.off()

keystone <- roles %>%
  filter(participation > 0.62 | connectivity > 2.5)

#write_xlsx(keystone, format_headers = F, "data/keystone_taxa.xlsx")



#### 7. SNPs/Microbes ####
# This is a really cool aspect of having BOTH microbial data and plant genomic data
# Look for associations between SNPs and ASVs
# Again use the abund/ubiq ASVs
# Need help from others (Kyle, Cloe, Brian etc. to run this)
# Most is done outside of R on the super computer
# Kyle will run
# Import files from Kyle

# Old
# snp_16S <- read.delim("data/AllSignificantSites16SWithGeneFunction.txt")
# nrow(snp_16S) # 1873
# sum(snp_16S$GeneIDs == "NoGenes") # 403
# sum(snp_16S$GeneIDs == "NoGenes")/nrow(snp_16S)
# snp_16S_taxa <- snp_16S %>%
#   group_by(ASV) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count)) %>%
#   mutate("ASV_ID" = gsub("_16S", "", ASV)) %>%
#   left_join(., input_n3_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
# snp_16S_genes <- snp_16S %>%
#   group_by(GeneIDs) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count))
# 
# snp_ITS <- read.delim("data/AllSignificantSitesITSWithGeneFunction.txt")
# nrow(snp_ITS) # 25921
# sum(snp_ITS$GeneIDs == "NoGenes") # 5528
# sum(snp_ITS$GeneIDs == "NoGenes")/nrow(snp_ITS)
# snp_ITS_taxa <- snp_ITS %>%
#   group_by(ASV) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count)) %>%
#   mutate("ASV_ID" = gsub("_ITS", "", ASV)) %>%
#   left_join(., input_n3_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
# snp_ITS_genes <- snp_ITS %>%
#   group_by(GeneIDs) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count))
# 
# # Update with PC1 and PC2 included
# # Import files from Kyle
# snp_16S <- read.delim("data/Window16SProteinIdentitiesNoGenesRM.txt", header = F)
# nrow(snp_16S) # 1820
# sum(snp_16S$GeneIDs == "NoGenes") # 403
# sum(snp_16S$GeneIDs == "NoGenes")/nrow(snp_16S)
# snp_16S_taxa <- snp_16S %>%
#   group_by(ASV) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count)) %>%
#   mutate("ASV_ID" = gsub("_16S", "", ASV)) %>%
#   left_join(., input_n3_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
# snp_16S_genes <- snp_16S %>%
#   group_by(GeneIDs) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count))



# Update from Kyle September 2, 2024 (use this!)
# Update is that some ZOTUs were transformed if necessary, other were not
# Uses the same PC1/PC2 as last update
snp_16S <- read.delim("data/16SCandidateSites_wTaxonomy.txt", header = F) %>%
  set_names(c("Test", "SNP", "Position", "Pvalue", "Gene1", "Gene2", "Taxonomy")) %>%
  filter(is.na(SNP) == FALSE) %>%
  separate(Test, into = c("OTU_ID", "Junk"), remove = F, sep = "_16S") %>%
  dplyr::select(-Junk)
nrow(snp_16S) # 640 (of 1654 rows, more than 1385 because some with multiple hits)
length(unique(snp_16S$OTU_ID)) # 375 ZOTUs of 1385 tested (27%)
sum(snp_16S$Gene1 == "NoGenes") # 433
sum(snp_16S$Gene1 == "NoGenes")/nrow(snp_16S) # 68% of sig. SNPs with no associated genes
# Richness, Shannon, PC1, PC2 no significant SNPs
snp_16S_taxa <- snp_16S %>%
  left_join(., input_auh_16S$taxonomy_loaded, by = c("OTU_ID" = "taxonomy8")) %>%
  mutate(Taxonomy = paste(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                          taxonomy5, taxonomy6, taxonomy7, OTU_ID, sep = ";")) %>%
  group_by(Test, Taxonomy) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
snp_16S_genes <- snp_16S %>%
  filter(Gene1 != "NoGenes") %>%
  group_by(Gene1) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
snp_16S_genes_taxlist <- snp_16S %>%
  filter(Gene1 != "NoGenes")
length(unique(snp_16S_genes_taxlist$OTU_ID)) # 162 ZOTUs with SNPs with genes

snp_ITS <- read.delim("data/ITSCandidateSites_wTaxonomy.txt", header = F) %>%
  dplyr::select(-V8) %>%
  set_names(c("Test", "SNP", "Position", "Pvalue", "Gene1", "Gene2", "Taxonomy")) %>%
  filter(is.na(SNP) == FALSE) %>%
  separate(Test, into = c("OTU_ID", "Junk"), remove = F, sep = "_ITS") %>%
  dplyr::select(-Junk)
nrow(snp_ITS) # 403 (out of 537 rows, more than 222 because some with multiple hits)
length(unique(snp_ITS$OTU_ID)) # 92 ZOTUs of of 222 tested (41%)
sum(snp_ITS$Gene1 == "NoGenes") # 239
sum(snp_ITS$Gene1 == "NoGenes")/nrow(snp_ITS) # 59% of sig. SNPs with no associated genes
# Richness, Shannon, PC1, PC2 no significant SNPs
snp_ITS_taxa <- snp_ITS %>%
  left_join(., input_auh_ITS$taxonomy_loaded, by = c("OTU_ID" = "taxonomy8")) %>%
  mutate(Taxonomy = paste(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                          taxonomy5, taxonomy6, taxonomy7, OTU_ID, sep = ";")) %>%
  group_by(Test, Taxonomy) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
snp_ITS_genes <- snp_ITS %>%
  filter(Gene1 != "NoGenes") %>%
  group_by(Gene1) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
snp_ITS_genes_taxlist <- snp_ITS %>%
  filter(Gene1 != "NoGenes")
length(unique(snp_ITS_genes_taxlist$OTU_ID)) # 51 ZOTUs with SNPs with genes

# Make Figure 5 - 
# Now the left column was the Sig. column in figure 4
# Make the dataframes
SNP_tested_16S <- num_h_16S %>%
  arrange(desc(prop)) %>%
  mutate(Subset = "Heritable")
sum(SNP_tested_16S$prop)
SNP_sig_16S <- input_n3_16S$taxonomy_loaded %>%
  filter(taxonomy8 %in% snp_16S$OTU_ID) %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% top_taxa_16S$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_16S$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 375) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S",
         Subset = "Sig. SNP")
sum(SNP_sig_16S$prop)
SNP_sig_gene_16S <- input_n3_16S$taxonomy_loaded %>%
  filter(taxonomy8 %in% snp_16S_genes_taxlist$OTU_ID) %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% top_taxa_16S$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_16S$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 162) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "16S",
         Subset = "Sig. SNP w/Gene")
sum(SNP_sig_gene_16S$prop)

SNP_tested_ITS <- num_h_ITS %>%
  arrange(desc(prop)) %>%
  mutate(Subset = "Heritable")
sum(SNP_tested_ITS$prop)
SNP_sig_ITS <- input_n3_ITS$taxonomy_loaded %>%
  filter(taxonomy8 %in% snp_ITS$OTU_ID) %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% num_h_ITS$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_ITS$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 92) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "ITS",
         Subset = "Sig. SNP")
sum(SNP_sig_ITS$prop)
SNP_sig_gene_ITS <- input_n3_ITS$taxonomy_loaded %>%
  filter(taxonomy8 %in% snp_ITS_genes_taxlist$OTU_ID) %>%
  mutate(taxonomy2 = ifelse(taxonomy2 %notin% num_h_ITS$taxonomy2,
                            "Other",
                            taxonomy2)) %>%
  mutate(taxonomy2 = factor(taxonomy2, levels = num_h_ITS$taxonomy2)) %>%
  group_by(taxonomy2) %>%
  summarise(num_sig = n()) %>%
  mutate(tot = 51) %>%
  mutate(prop = num_sig/tot) %>%
  mutate(Dataset = "ITS",
         Subset = "Sig. SNP w/Gene")
sum(SNP_sig_gene_ITS$prop)

# Plot with heritable, sig SNP, sig SNP with gene
SNP_tax_16S <- rbind(SNP_tested_16S, SNP_sig_16S, SNP_sig_gene_16S) %>%
  mutate(Subset = factor(Subset, levels = c("Heritable", "Sig. SNP", "Sig. SNP w/Gene"))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(levels_16S$tax))))
n_text_16S <- data.frame(x = c("Heritable", "Sig. SNP", "Sig. SNP w/Gene"),
                         y = c(105, 105, 105),
                         label = c("n = 1385", "n = 375", "n = 162"))
g10 <- ggplot(SNP_tax_16S, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_16S, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(12, "Paired")[11:1])) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "a) Archaea/Bacteria") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g10
SNP_tax_ITS <- rbind(SNP_tested_ITS, SNP_sig_ITS, SNP_sig_gene_ITS) %>%
  mutate(Subset = factor(Subset, levels = c("Heritable", "Sig. SNP", "Sig. SNP w/Gene"))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(ITS_levels$tax))))
n_text_ITS <- data.frame(x = c("Heritable", "Sig. SNP", "Sig. SNP w/Gene"),
                         y = c(105, 105, 105),
                         label = c("n = 222", "n = 92", "n = 51"))
g11 <- ggplot(SNP_tax_ITS, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_ITS, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(6, "Paired"))) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "b) Fungi") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g11

pdf("InitialFigs/SNP_taxa.pdf", width = 7, height = 5)
plot_grid(g10, g11, ncol = 1, align = "v")
dev.off()



# Combine all counts from heritability testing and SNP testing into 1 figure!
# Otherwise the Heritability column is shown in both figures, which isn't good
# Use stacked bars for now but could also consider an alluvial plot
SNP_tax_16S <- rbind(tot_count_16S, tested_count_16S, SNP_tested_16S, SNP_sig_16S, SNP_sig_gene_16S) %>%
  mutate(Subset = factor(Subset, 
                         levels = c("All", "Tested", "Heritable", "Sig. SNP", "Sig. SNP w/Gene"))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(levels_16S$tax))))
n_text_16S <- data.frame(x = c("All", "Tested", "Heritable", "Sig. SNP", "Sig. SNP w/Gene"),
                         y = c(105, 105, 105, 105, 105),
                         label = c("n = 19274", "n = 1726", "n = 1385", "n = 375", "n = 162"))
g12 <- ggplot(SNP_tax_16S, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_16S, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(12, "Paired")[11:1])) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "a) Archaea/Bacteria") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g12
SNP_tax_ITS <- rbind(tot_count_ITS, tested_count_ITS, SNP_tested_ITS, SNP_sig_ITS, SNP_sig_gene_ITS) %>%
  mutate(Subset = factor(Subset, 
                         levels = c("All", "Tested", "Heritable", "Sig. SNP", "Sig. SNP w/Gene"))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(ITS_levels$tax))))
n_text_ITS <- data.frame(x = c("All", "Tested", "Heritable", "Sig. SNP", "Sig. SNP w/Gene"),
                         y = c(105, 105, 105, 105, 105),
                         label = c("n = 1425", "n = 240", "n = 222", "n = 92", "n = 51"))
g13 <- ggplot(SNP_tax_ITS, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_ITS, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(6, "Paired"))) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "b) Fungi") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g13

pdf("InitialFigs/Heritability_TaxaCountswGWAS.pdf", width = 7, height = 5)
plot_grid(g12, g13, ncol = 1, align = "v")
dev.off()

# Code copied from above
# Plot with tot count and au count and sig count
num_tax_ITS <- rbind(tot_count_ITS, tested_count_ITS, num_h_ITS) %>%
  mutate(Subset = factor(Subset, levels = c("All", "Tested", "Sig."))) %>%
  mutate(taxonomy2 = factor(taxonomy2, 
                            levels = c("Other", "NA", rev(ITS_levels$tax))))
nrow(input_n3_ITS$taxonomy_loaded)
nrow(input_au_ITS$taxonomy_loaded)
nrow(input_auh_ITS$taxonomy_loaded)
n_text_ITS <- data.frame(x = c("All", "Tested", "Sig."),
                         y = c(105, 105, 105),
                         label = c("n = 1425", "n = 240", "n = 222"))
g9 <- ggplot(num_tax_ITS, aes(x = Subset, y = prop*100, fill = taxonomy2)) +
  geom_bar(stat = "identity") +
  geom_text(data = n_text_ITS, aes(x = x, y = y, label = label), 
            inherit.aes = F, size = 3) +
  scale_fill_manual(values = c("grey90", "grey40", brewer.pal(6, "Paired"))) +
  labs(x = NULL,
       y = "% ZOTU Count",
       fill = "Phylum",
       title = "b) Fungi") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = -2))
g9

# Export for sharing
#write.csv(snp_16S_taxa, "data/snp_16S_taxa.csv", row.names = F)
#write.csv(snp_ITS_taxa, "data/snp_ITS_taxa.csv", row.names = F)
#write.csv(snp_16S_genes, "data/snp_16S_genes.csv", row.names = F)
#write.csv(snp_ITS_genes, "data/snp_ITS_genes.csv", row.names = F)



# Other figures
# Conceptual figure with diagram and some microbe shapes and gene text
# Plot with the significant red dots from Manhattan plots
# Heatmap with p-value info?
# Heatmap with total number of archaea, bacteria, fungal hits by chromosome
# Heatmap with number of unique archaea, bacteria, fungal ZOTUs by chromosome

# Total counts
snp_arc <- snp_16S %>%
  filter(grepl("Archaea", Taxonomy))
snp_bac <- snp_16S %>%
  filter(grepl("Bacteria", Taxonomy))
snp_fun <- snp_ITS
chr_arc <- as.data.frame(table(snp_arc$SNP)) %>%
  set_names(c("ChromosomeID", "Archaea"))
chr_bac <- as.data.frame(table(snp_bac$SNP)) %>%
  set_names(c("ChromosomeID", "Bacteria"))
chr_fun <- as.data.frame(table(snp_fun$SNP)) %>%
  set_names(c("ChromosomeID", "Fungi"))
chr_counts <- chr_arc %>%
  full_join(., chr_bac, by = "ChromosomeID") %>%
  full_join(., chr_fun, by = "ChromosomeID") %>%
  replace_na_with(value = 0) %>%
  mutate(Chrome = "Chromosome",
         Number = as.integer(substr(ChromosomeID, start = 11, stop = 12))) %>%
  mutate(Chromosome = paste(Chrome, Number, sep = "")) %>%
  arrange(Number) %>%
  column_to_rownames(var = "Chromosome") %>%
  dplyr::select(Archaea, Bacteria, Fungi)
pheatmap(mat = chr_counts,
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         angle_col = 315,
         display_numbers = chr_counts,
         filename = "InitialFigs/Kyle_ChromosomeCountsTotal.png",
         width = 4,
         height = 5)

# ZOTU counts
chr_arc <- snp_arc %>%
  group_by(SNP) %>%
  summarise(Archaea = length(unique(OTU_ID))) %>%
  set_names(c("ChromosomeID", "Archaea"))
chr_bac <- snp_bac %>%
  group_by(SNP) %>%
  summarise(Bacteria = length(unique(OTU_ID))) %>%
  set_names(c("ChromosomeID", "Bacteria"))
chr_fun <- snp_fun %>%
  group_by(SNP) %>%
  summarise(Fungi = length(unique(OTU_ID))) %>%
  set_names(c("ChromosomeID", "Fungi"))
chr_zotus <- chr_arc %>%
  full_join(., chr_bac, by = "ChromosomeID") %>%
  full_join(., chr_fun, by = "ChromosomeID") %>%
  replace_na_with(value = 0) %>%
  mutate(Chrome = "Chromosome",
         Number = as.integer(substr(ChromosomeID, start = 11, stop = 12))) %>%
  mutate(Chromosome = paste(Chrome, Number, sep = "")) %>%
  arrange(Number) %>%
  column_to_rownames(var = "Chromosome") %>%
  dplyr::select(Archaea, Bacteria, Fungi)
pheatmap(mat = chr_zotus,
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         angle_col = 315,
         display_numbers = chr_zotus,
         filename = "InitialFigs/Kyle_ChromosomeCountsZOTUs.png",
         width = 4,
         height = 5)



#### 8. Inhibitors ####
# Check heritability/genotype effects of 6 inhibitor OTUs 
# For paper with Lara Vimercati/Alisha Quandt
inhibitors <- c("OTU_10441", "OTU_12372", "OTU_978", "OTU_4850", "OTU_275", "OTU_1672")
tax_sum_OTU_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                      level = 8, 
                                      report_higher_tax = F)
input_6 <- filter_taxa_from_input(input_filt_rare_16S,
                                  taxa_IDs_to_keep = inhibitors)
View(input_6$taxonomy_loaded)
facet_names_6 <- c("OTU_275" = "Pseudomonas\nsp.",
                   "OTU_978" = "Rhodococcus\nsp.",
                   "OTU_12372" = "Bacillus\nsp.",
                   "OTU_1672" = "Variovorax\nparadoxus",
                   "OTU_4850" = "Paenibacillus\npolymyxa",
                   "OTU_10441" = "Sphingobacterium\nkitahiroshimense")
ts6 <- tax_sum_OTU_16S %>%
  filter(rownames(.) %in% input_6$taxonomy_loaded$taxonomy8) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df6 <- input_6$map_loaded %>%
  dplyr::select(sampleID, rep, pedigree) %>%
  left_join(., ts6, by = "sampleID")
prev_all <- df6 %>%
  pivot_longer(cols = c(4:ncol(df6)), names_to = "OTU") %>%
  group_by(OTU) %>%
  summarise(otu_n = sum(value > 0)) %>%
  mutate(otu_per = otu_n/368*100)
prev_6 <- df6 %>%
  pivot_longer(cols = c(4:ncol(df6)), names_to = "OTU") %>%
  group_by(pedigree, OTU) %>%
  summarise(pedigree_n = n(),
            otu_n = sum(value > 0)) %>%
  mutate(otu_per = otu_n/pedigree_n*100) %>%
  mutate(OTU = factor(OTU,
                      levels = c("OTU_275", "OTU_978", "OTU_12372", "OTU_1672",
                                 "OTU_4850", "OTU_10441")))
pdf("InitialFigs/Inhibitors_GenotypePrevalence95.pdf", width = 7, height = 4)
ggplot(prev_6, aes(pedigree, otu_per)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ OTU, ncol = 6, labeller = as_labeller(facet_names_6)) +
  labs(x = "Genotype (n = 95)",
       y = "% Prevalence") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 8, face = "italic"),
        strip.background = element_rect(linewidth = 0.2, fill = "white"))
dev.off()

# This is why we need to use NB
hist(df6$OTU_10441)
hist(df6$OTU_12372)
hist(df6$OTU_1672)
hist(df6$OTU_275)
hist(df6$OTU_4850)
hist(df6$OTU_978)

results6 <- as.data.frame(matrix(NA, nrow = ncol(df6), ncol = 7)) %>%
  set_names(c("OTU", "Heritability", "LeveneG", "ShapiroAOV", "ShapiroNB4",
              "GenotypePaov", "GenotypePnb"))
for (i in 4:ncol(df6)) {
  # OTU name
  results6$OTU[i] <- names(df6)[i]
  
  # Levene Test
  l2 <- leveneTest(df6[,i] ~ df6$pedigree)
  results6$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df6 %>%
    dplyr::select(rep, pedigree, i) %>%
    set_names(c("rep", "pedigree", "OTU"))
  
  # Models
  m <- aov(OTU ~ rep + pedigree, data = df)
  nb3 <- MASS::glm.nb(OTU ~ rep + pedigree, data = df) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ rep + pedigree, data = df, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "II", singular.ok = TRUE)
  c4 <- Anova(nb4, test.statistic = "F")
  
  results6$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results6$GenotypePnb[i] <- c4$`Pr(>F)`[2]
  
  eta <- eta_sq_glm(nb4) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results6$Heritability[i] <- eta$eta[2]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s3 <- shapiro.test(nb4$residuals)
  results6$ShapiroAOV[i] <- s1$p.value
  results6$ShapiroNB4[i] <- s3$p.value
}

results6 <- results6 %>%
  drop_na(OTU) %>%
  mutate(GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr")) %>%
  mutate(SymbolGeno = ifelse(GenotypePnb <= 0.001,
                             "***",
                             ifelse(GenotypePnb > 0.001 & GenotypePnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePnb > 0.01 & GenotypePnb <= 0.05,
                                           "*",
                                           "N.S."))))
sigGeno <- results6 %>%
  filter(GenotypePaov < 0.05) # 1 p, 0 pFDR
sigGenoNB <- results6 %>%
  filter(GenotypePnb < 0.05) # 6 p, 6 pFDR

otu_tax_16S <- input_6$taxonomy_loaded %>%
  group_by(taxonomy8) %>%
  slice_head(n = 1) %>%
  as.data.frame()
sigGeno_otus <- sigGenoNB$OTU
otherGeno_otus <- results6$OTU[results6$OTU %notin% sigGenoNB$OTU]
df6_long_geno <- df6 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  mutate(name = factor(name,
                       levels = c(sigGeno_otus, otherGeno_otus))) %>%
  mutate(name = as.character(name)) %>%
  mutate(Genotype = "Genotype") %>%
  mutate(name = factor(name,
                       levels = c("OTU_275", "OTU_978", "OTU_12372", "OTU_1672",
                                  "OTU_4850", "OTU_10441")))
sigLab_geno <- df6_long_geno %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/10) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  mutate(name = factor(name,
                       levels = c(sigGeno_otus, otherGeno_otus))) %>%
  arrange(name = desc(levels(df6_long_geno$name))) %>%
  mutate(name = as.character(name)) %>%
  mutate(name = factor(name,
                       levels = c("OTU_275", "OTU_978", "OTU_12372", "OTU_1672",
                                  "OTU_4850", "OTU_10441")))
facet_names_6 <- c("OTU_275" = "Pseudomonas\nsp.",
                   "OTU_978" = "Rhodococcus\nsp.",
                   "OTU_12372" = "Bacillus\nsp.",
                   "OTU_1672" = "Variovorax\nparadoxus",
                   "OTU_4850" = "Paenibacillus\npolymyxa",
                   "OTU_10441" = "Sphingobacterium\nkitahiroshimense",
                   "Genotype" = "Genotype (n = 95)")
pdf("InitialFigs/Inhibitors_GenotypeAbundance95.pdf", width = 7, height = 5)
ggplot(df6_long_geno, aes(pedigree, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0, height = 0) +
  geom_text(data = subset(sigLab_geno, SymbolGeno != "N.S."),
            aes(x = 48, y = y*100, label = SymbolGeno),
            inherit.aes = F, size = 8) +
  facet_grid(name ~ Genotype, scales = "free_y", labeller = as_labeller(facet_names_6)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.2), 
                                        add = c(0, 0))) +
  labs(x = "Genotype",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),        
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 6, face = "italic"),
        strip.background = element_rect(linewidth = 0.2, fill = "white"))
dev.off()



#### 9. NTI ####
# Look at Nearest Taxon Index to infer assembly mechanisms
# Just do for 16S as ITS trees are questionable
# Hypothesis that genotype effects lead to deterministic assembly
# Read in and filter repset fasta
nrow(input_n3_16S$taxonomy_loaded)
cs15 <- as.data.frame(rowSums(input_n3_16S$data_loaded)) %>%
  set_names("Reads") %>%
  filter(Reads >= 15)
input_15 <- filter_taxa_from_input(input_n3_16S,
                                   taxa_IDs_to_keep = rownames(cs15))
nrow(input_15$taxonomy_loaded) # 10522
View(input_15$taxonomy_loaded)
repset <- readFasta("~/Desktop/Sunflower/Cloe/16S/16S_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% input_15$taxonomy_loaded$taxonomy8)
#writeFasta(repset, "~/Desktop/Sunflower/Cloe/16S/16S_rep_set_zotus_filt15.fa")

# Transfer to microbe, run qiime commands for alignment and tree
# screen align_seqs.py -i  16S_rep_set_zotus_filt15.fa -m muscle -o ./
# screen make_phylogeny.py -i aligned.fasta -o ./rep_phylo.tre
# Transfer tree back here

# The NTI values between −2 and 2 indicate stochastic community assembly, whereas NTI values less than −2 or higher than 2 indicate that deterministic processes play a more important role in structuring the community (deterministic environmental filtering)

# Run on microbe
#saveRDS(input_15, "data/input_15.rds")
# Use dataframe with taxa with >= 15 reads in whole dataset (10522)
bac_asv <- as.data.frame(t(input_15$data_loaded))
tree <- read.tree("data/rep_phylo.tre")
colnames(bac_asv)
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$tip.label
sum(colnames(bac_asv) %in% tree$tip.label)
sum(tree$tip.label %in% colnames(bac_asv))
#tree <- prune.sample(bac_asv, tree)
#bac_asv <- bac_asv[,tree$tip.label]
phy.dist <- cophenetic(tree)
NTI <- NTI.p(bac_asv, 
             phy.dist, 
             nworker = 4, 
             memo.size.GB = 50,
             weighted = TRUE, 
             rand = 1000,
             check.name = TRUE, 
             output.MNTD = FALSE,
             sig.index = "NTI",
             silent = FALSE)
#saveRDS(NTI, "NTI.rds")
NTI <- readRDS("data/NTI.rds") %>%
  rownames_to_column(var = "sampleID")
hist(NTI$NTI)
t.test(NTI$NTI, mu = 0) # p < 0.05, different than 0!
nrow(NTI) # 368
sum(NTI$NTI > 2) # 127 of 368 above the 2 cutoff
input_15$map_loaded <- input_15$map_loaded %>%
  left_join(., NTI, by = "sampleID")
rownames(input_n3_16S$map_loaded)
rownames(input_15$map_loaded) <- input_15$map_loaded$sampleID

# Get descriptive info
min(input_15$map_loaded$NTI) # -4.186063
max(input_15$map_loaded$NTI) # 5.074002
mean(input_15$map_loaded$NTI) # 1.421827
se(input_15$map_loaded$NTI) # 0.07257843
sd(input_15$map_loaded$NTI) # 1.392296

# Test and plot versus Sclerotinia, chlorophyll, richness, shannon
summary(lm(NTI ~ DiseaseIncidence, data = input_15$map_loaded))
summary(lm(NTI ~ `Chlorophyll concentration`, data = input_15$map_loaded))
summary(lm(NTI ~ rich, data = input_15$map_loaded))
summary(lm(NTI ~ shannon, data = input_15$map_loaded))
plot(input_15$map_loaded$DiseaseIncidence, input_15$map_loaded$NTI)
plot(input_15$map_loaded$`Chlorophyll concentration`, input_15$map_loaded$NTI)
plot(input_15$map_loaded$rich, input_15$map_loaded$NTI)
plot(input_15$map_loaded$shannon, input_15$map_loaded$NTI)
pdf("InitialFigs/NTI_Sclero.pdf", width = 7 , height = 5)
ggplot(input_15$map_loaded, aes(DiseaseIncidence, NTI)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
dev.off()
pdf("InitialFigs/NTI_Chloro.pdf", width = 7 , height = 5)
ggplot(input_15$map_loaded, aes(`Chlorophyll concentration`, NTI)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
dev.off()

# Test and plot
leveneTest(input_15$map_loaded$NTI ~ input_15$map_loaded$pedigree) # Homogeneous
leveneTest(input_15$map_loaded$NTI ~ input_15$map_loaded$rep) # Not Homogeneous
m <- lmer(NTI ~ pedigree + (1|rep), data = input_15$map_loaded)
summary(m)
Anova(m)
anova(m)
aovtab <- anova(m)
mnull <- lmer(NTI ~ 1 + (1|rep), data = input_15$map_loaded)
anova(mnull, m)
m2 <- mixed(NTI ~ pedigree + (1|rep), data = input_15$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two way ANOVA
m <- aov(NTI ~ rep + pedigree, data = input_15$map_loaded)
summary(m)
Anova(m, type = "II")
h <- eta_sq(m)[2]
h
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
TukeyHSD(m)
t <- emmeans(object = m, specs = "pedigree") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "NTI",
         y = max(input_15$map_loaded$NTI)+
           (max(input_15$map_loaded$NTI)-
              min(input_15$map_loaded$NTI))/10)

pdf("InitialFigs/NTI_15reads.pdf", width = 7, height = 5)
ggplot(input_15$map_loaded, aes(reorder(pedigree, NTI, mean), NTI)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = -2, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.75, pch = 16, aes(colour = rep)) +
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6)) +
  #geom_text(data = t, aes(pedigree, y, label = str_trim(.group)), 
  #          size = 1, color = "black") +
  labs(x = NULL, y = "NTI") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank())
dev.off()

ggplot(input_15$map_loaded, aes(reorder(rep, NTI, mean), NTI)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = -2, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  #geom_text(data = t, aes(pedigree, y, label = str_trim(.group)), 
  #          size = 1, color = "black") +
  labs(x = NULL, y = "NTI") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank())

# Use the same style as the alpha diversity plot. Figure 5
pdf("FinalFigs/Figure5.pdf", width = 8, height = 4)
ggplot(input_15$map_loaded, aes(reorder(pedigree, NTI, mean), NTI)) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -2, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, alpha = 1, pch = 16) +
  labs(x = "Genotype", y = "NTI") +
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
dev.off()



#### 10. Sloan ####
# Run Sloane Neutral Model
# Run with pctax package, use github version for plot editing
#devtools::install_github("Asa12138/pctax", force = TRUE)
library(pctax)

# 16S
ncm_res <- ncm(input_15$data_loaded, model = "nls")
table(ncm_res$ncm_data$group) # 9884 above, 364 below, 274 in
9884/10522*100
364/10522*100
274/10522*100

# ITS
nrow(input_n3_ITS$taxonomy_loaded) # 1425
cs15_ITS <- as.data.frame(rowSums(input_n3_ITS$data_loaded)) %>%
  set_names("Reads") %>%
  filter(Reads >= 15)
input_15_ITS <- filter_taxa_from_input(input_n3_ITS,
                                       taxa_IDs_to_keep = rownames(cs15_ITS))
nrow(input_15_ITS$taxonomy_loaded) # 981

ncm_res_ITS <- ncm(input_15_ITS$data_loaded, model = "nls")
table(ncm_res_ITS$ncm_data$group) # 558 above, 166 below, 257 in
558/981*100
166/981*100
257/981*100

# Combined plot (Figure 6)
fig6a <- plot(ncm_res,
     mycols = c(Above = "#069870", Below = "#e29e02", In = "#1e353a"),
     text_position = NULL,
     pie_text_params = list(size = 0)) +
  ggtitle("a) Archaea/Bacteria") +
  theme(plot.title = element_text(hjust = 0, vjust = -1))
fig6a

fig6b <- plot(ncm_res_ITS,
              mycols = c(Above = "#069870", Below = "#e29e02", In = "#1e353a"),
              text_position = NULL,
              pie_text_params = list(size = 0)) +
  ggtitle("b) Fungi") +
  theme(plot.title = element_text(hjust = 0, vjust = -1))
fig6b

# Custom plot, use this
lib_ps("patchwork", library = FALSE)
mycols = c("Above" = "#069870", "Below" = "#e29e02", "In" = "#1e353a",
           lincol <- "#4a8fb8")

out <- ncm_res[[2]]
text_position <- c(log(max(out$p)) - 2, 0.05)
out$group %>%
  table() %>%
  as.data.frame() -> ad
colnames(ad) <- c("type", "n")
p1 <- ggplot() +
  geom_line(data = out, aes(x = log(p), y = freq.pred), 
            linewidth = 1.2, linetype = 1, col = lincol) +
  geom_line(data = out, aes(x = log(p), y = Lower), 
            linewidth = 1.2, linetype = 2, col = lincol) +
  geom_line(data = out, aes(x = log(p), y = Upper), 
            linewidth = 1.2, linetype = 2, col = lincol) +
  geom_point(data = out, aes(x = log(p), y = freq, color = group), size = 1) +
  labs(x = "log10(mean relative abundance)",
       y ="Occurrence frequency") +
  scale_colour_manual(values = mycols) +
  annotate("text", 
           x = text_position[1], 
           y = text_position[2], 
           label = paste("Nm = ", sprintf("%.0f", ncm_res[[1]][1] * ncm_res[[1]][3]), 
                         sep = ""), size = 4) +
  annotate("text", 
           x = text_position[1], 
           y = text_position[2] + 0.1,
           label = paste0("R^2 ==", round(ncm_res[[1]][2], 3)), size = 4, parse = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggtitle("a) Archaea/Bacteria") +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(), legend.background = element_rect(I(0)),
        plot.title = element_text(hjust = 0.5, vjust = -1))
p1
pie <- pcutils::gghuan(ad, name = FALSE, text_params = list(size = 0)) +
  xlim(0.2, 3.3) +
  scale_fill_manual(values = mycols) +
  theme(plot.background = element_rect(I(0), linetype = 0),
        panel.background = element_rect(I(0))) # 去除图片背景白色
pie
fig6a <- p1 + 
  patchwork::inset_element(pie, left = -0.1, bottom = 0.5, right = 0.3, top = 1)
fig6a

out <- ncm_res_ITS[[2]]
text_position <- c(log(max(out$p)) - 2.5, 0.05)
out$group %>%
  table() %>%
  as.data.frame() -> ad
colnames(ad) <- c("type", "n")
p1 <- ggplot() +
  geom_line(data = out, aes(x = log(p), y = freq.pred), 
            linewidth = 1.2, linetype = 1, col = lincol) +
  geom_line(data = out, aes(x = log(p), y = Lower), 
            linewidth = 1.2, linetype = 2, col = lincol) +
  geom_line(data = out, aes(x = log(p), y = Upper), 
            linewidth = 1.2, linetype = 2, col = lincol) +
  geom_point(data = out, aes(x = log(p), y = freq, color = group), size = 1) +
  labs(x = "log10(mean relative abundance)",
       y ="Occurrence frequency") +
  scale_colour_manual(values = mycols) +
  annotate("text", 
           x = text_position[1], 
           y = text_position[2], 
           label = paste("Nm = ", sprintf("%.0f", ncm_res_ITS[[1]][1] * ncm_res_ITS[[1]][3]), 
                         sep = ""), size = 4) +
  annotate("text", 
           x = text_position[1], 
           y = text_position[2] + 0.1,
           label = paste0("R^2 ==", round(ncm_res_ITS[[1]][2], 3)), size = 4, parse = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggtitle("b) Fungi") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.4),
        legend.title = element_blank(), legend.background = element_rect(I(0)),
        plot.title = element_text(hjust = 0.5, vjust = -1))
p1
pie <- pcutils::gghuan(ad, name = FALSE, text_params = list(size = 0)) +
  xlim(0.2, 3.3) +
  scale_fill_manual(values = mycols) +
  theme(plot.background = element_rect(I(0), linetype = 0),
        panel.background = element_rect(I(0))) # 去除图片背景白色
pie
fig6b <- p1 + 
  patchwork::inset_element(pie, left = -0.1, bottom = 0.5, right = 0.3, top = 1)
fig6b

pdf("FinalFigs/Figure6.pdf", width = 8, height = 4)
plot_grid(fig6a, fig6b)
dev.off()

#### 11. IAA ####
# Check which taxa have genes for IAA
# From Kyle: Check these:
# 1. Acidobacteriota Holophagae (no)
# 2. Actinobacteriota Solirubrobacterales (no)
# 3. Proteobacteria Hyphomonadaceae (no)
# 4. Acidobacteriota Holophagae (no)
# 5. Myxococcota Polyangiales (no)
# Are there genomes of those taxa in IMG? Yes.
# Holophagae: 20
# Solirubrobacterales: 28
# Hyphomonadaceae: 63
# Polyangiales: 37

# Check for more detailed taxonomy. BLAST on NCBI
# OTU_88, OTU_1044, OTU_1124, OTU_1156, and OTU_1514
iaa_repset <- readFasta("~/Desktop/Kane/Cloe/16S/16S_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% c("OTU_88", "OTU_1044", "OTU_1124", "OTU_1156", "OTU_1514"))
# OTU_88: 100% with many uncultured bacteria
# OTU_1044: 100% with many uncultured bacteria
# OTU_1124: 100% with many uncultured bacteria
# OTU_1156: 100% with many uncultured bacteria
# OTU_1514: 100% with many uncultured bacteria

# Searched MetaCyc for IAA biosynthesis pathways
# There are 6 pathways on there, 4 of which are performed by bacteria
# III: K00466, K01426, K21801
# IV: K01721, K20807, (K01426, K21801)
# V: K01501
# VI: K14265, K04103, K11817, K22417

# Searched IMG for these KOs on Dec 5 2024
# Downloaded list of genomes with those KOs

# K00466 (iaaM) - those taxa don't have
iaa <- read.delim("data/K00466.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus != "Uncultured") # removes 2
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))

# K01426 (amiE)
iaa <- read.delim("data/K01426.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus != "Uncultured") %>% # removes 114
  filter(Domain %in% c("Archaea", "Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



# K21801 (iaaH)
iaa <- read.delim("data/K21801.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus != "Uncultured") %>% # removes 114
  filter(Domain %in% c("Archaea", "Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
taxonomy_info <- Filter(function(x) !is.logical(x), taxonomy_info)
taxonomy_info_K21801 <- taxonomy_info
#saveRDS(taxonomy_info_K21801, "data/taxonomy_info_K21801.rds")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = length(taxonomy_info), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:length(taxonomy_info)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_taxonomy_K21801 <- iaa_taxonomy
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



# K01721
iaa <- read.delim("data/K01721.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus %notin% c("Uncultured", "unclassified", "co-assembly", "putative")) %>%
  filter(Domain %in% c("Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
iaa_genus_toclassify <- iaa_genus %>%
  filter(Genus %notin% iaa_taxonomy_K21801$Genus)
taxonomy_info <- classification(iaa_genus_toclassify$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



# K20807
iaa <- read.delim("data/K20807.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus %notin% c("Uncultured", "unclassified", "co-assembly", "putative")) %>%
  filter(Domain %in% c("Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



# K01501
iaa <- read.delim("data/K01501.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus %notin% c("Uncultured", "unclassified", "co-assembly", "putative")) %>%
  filter(Domain %in% c("Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



# K04103
iaa <- read.delim("data/K04103.txt") %>%
  mutate(Genome.Name = gsub("Candidatus ", "", Genome.Name)) %>% # remove Candidatus
  separate(Genome.Name, into = c("Genus", "Species", "Strain1", "Strain2"), 
           sep = " ", remove = F) %>%
  filter(Genus %notin% c("Uncultured", "unclassified", "co-assembly", "putative")) %>%
  filter(Domain %in% c("Bacteria"))
iaa_genus <- iaa %>%
  group_by(Genus) %>%
  summarise(count = n())
taxonomy_info <- classification(iaa_genus$Genus, db = "ncbi")
iaa_taxonomy <- as.data.frame(matrix(data = NA, nrow = nrow(iaa_genus), ncol = 5)) %>%
  set_names(c("Phylum", "Class", "Order", "Family", "Genus"))
for (i in 1:nrow(iaa_genus)) {
  resP <- taxonomy_info[[i]] %>%
    filter(rank == "phylum") %>%
    pull(name)
  iaa_taxonomy$Phylum[i] <- if (length(resP) == 0) NA else resP
  resC <- taxonomy_info[[i]] %>%
    filter(rank == "class") %>%
    pull(name)
  iaa_taxonomy$Class[i] <- if (length(resC) == 0) NA else resC
  resO <- taxonomy_info[[i]] %>%
    filter(rank == "order") %>%
    pull(name)
  iaa_taxonomy$Order[i] <- if (length(resO) == 0) NA else resO
  resF <- taxonomy_info[[i]] %>%
    filter(rank == "family") %>%
    pull(name)
  iaa_taxonomy$Family[i] <- if (length(resF) == 0) NA else resF
  iaa_taxonomy$Genus[i] <- names(taxonomy_info)[i]
}
iaa_taxonomy <- iaa_taxonomy %>%
  left_join(., iaa_genus, by = "Genus")
iaa_family <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(count = sum(count))
hyph <- iaa_family %>%
  filter(Family == "Hyphomonadaceae")
iaa_order <- iaa_taxonomy %>%
  group_by(Phylum, Class, Order) %>%
  summarise(count = sum(count))
soli <- iaa_order %>%
  filter(Order == "Solirubrobacterales")
poly <- iaa_order %>%
  filter(Order == "Polyangiales")
iaa_class <- iaa_taxonomy %>%
  group_by(Phylum, Class) %>%
  summarise(count = sum(count))
holo <- iaa_class %>%
  filter(Class == "Holophagae")
iaa_phylum <- iaa_taxonomy %>%
  group_by(Phylum) %>%
  summarise(count = sum(count))



#### 12. NCBI ####
# Submit ZOTU sequences to NCBI
# First need to take the repset fasta from Jason and subset it to the used ZOTUs
# Used ZOTU IDs are from the rarefied dataset with n = 368 samples
repset_16S <- microseq::readFasta("~/Desktop/Kane/Cloe/16S/16S_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% input_n3_16S$taxonomy_loaded$taxonomy8)
nrow(repset_16S) == nrow(input_n3_16S$taxonomy_loaded)
microseq::writeFasta(repset_16S, "data/repset_16S.fasta")

repset_ITS <- microseq::readFasta("~/Desktop/Kane/Cloe/ITS/ITS_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% input_n3_ITS$taxonomy_loaded$taxonomy8)
nrow(repset_ITS) == nrow(input_n3_ITS$taxonomy_loaded)
microseq::writeFasta(repset_ITS, "data/repset_ITS.fasta")

# Now assign ASVs to samples (not comprehensive, just first sample)
info <- input_n3_16S$data_loaded
names(info) <- input_n3_16S$map_loaded$sampleID
for (i in 1:ncol(info)) {
  for (j in 1:nrow(info)) {
    ifelse(info[j, i] > 0, info[j, i] <- names(info)[i], info[j, i] <- "")
  }
}
info_cat <- info
info_cat <- info_cat %>%
  mutate_all(na_if, "") %>%
  mutate(unite(., "sample_name", c(names(info)), sep = ", ")) %>%
  mutate(sample_name = gsub("NA, ", "", sample_name)) %>%
  mutate(sample_name = gsub(", NA", "", sample_name)) %>%
  rownames_to_column(var = "Sequence_ID") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_cat %>%
  separate(sample_name, into = c("sample_name", "Junk"), sep = ", ") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_first %>%
  filter(Sequence_ID %in% repset_16S$Header)
write_tsv(info_first, file = "data/biosample_assignment_16S.tsv")

info <- input_n3_ITS$data_loaded
names(info) <- input_n3_ITS$map_loaded$sampleID
for (i in 1:ncol(info)) {
  for (j in 1:nrow(info)) {
    ifelse(info[j, i] > 0, info[j, i] <- names(info)[i], info[j, i] <- "")
  }
}
info_cat <- info
info_cat <- info_cat %>%
  mutate_all(na_if, "") %>%
  mutate(unite(., "sample_name", c(names(info)), sep = ", ")) %>%
  mutate(sample_name = gsub("NA, ", "", sample_name)) %>%
  mutate(sample_name = gsub(", NA", "", sample_name)) %>%
  rownames_to_column(var = "Sequence_ID") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_cat %>%
  separate(sample_name, into = c("sample_name", "Junk"), sep = ", ") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_first %>%
  filter(Sequence_ID %in% repset_ITS$Header)
write_tsv(info_first, file = "data/biosample_assignment_ITS.tsv")


# End Script