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
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(mctoolsr) # Microbial analyses
library(vegan) # Multivariate stats
library(RColorBrewer) # Colors
library(microseq) # for fastas
library(car) # stats
library(MASS) # stats
library(FSA) # se
library(lme4) # LMER
library(emmeans) # TukeyHSD
library(multcomp) # cld
library(lmerTest) # For Sattherwaite df
library(afex) # For alternative mixed model
library(ggh4x) # For plots


# Repo Directory
setwd("~/Documents/GitHub/SunflowerG/")

# Functions
`%notin%` <- Negate(`%in%`)
source("code/effectSize.R")



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
length(table(input_filt_rare_16S$map_loaded$pedigree)) # 98 pedigrees

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

# Save data
saveRDS(input_filt_16S, "data/input_filt_16S.rds")
saveRDS(input_filt_rare_16S, "data/input_filt_rare_16S.rds")
saveRDS(input_n4_16S, "data/input_n4_16S.rds")



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

# Save data
saveRDS(input_filt_ITS, "data/input_filt_ITS.rds")
saveRDS(input_filt_rare_ITS, "data/input_filt_rare_ITS.rds")
saveRDS(input_n4_ITS, "data/input_n4_ITS.rds")



#### _Start Here ####
input_filt_16S <- readRDS("data/input_filt_16S.rds")
input_filt_rare_16S <- readRDS("data/input_filt_rare_16S.rds")
input_n4_16S <- readRDS("data/input_n4_16S.rds")

input_filt_ITS <- readRDS("data/input_filt_ITS.rds")
input_filt_rare_ITS <- readRDS("data/input_filt_rare_ITS.rds")
input_n4_ITS <- readRDS("data/input_n4_ITS.rds")



#### 2. Alpha ####
#### _16S ####
# Get descriptive info
min(input_n4_16S$map_loaded$rich) # 1676
max(input_n4_16S$map_loaded$rich) # 2465
mean(input_n4_16S$map_loaded$rich) # 2096.584
se(input_n4_16S$map_loaded$rich) # 7.755601
plot(as.factor(input_n4_16S$map_loaded$pedigree), input_n4_16S$map_loaded$rich)
plot(input_n4_16S$map_loaded$DiseaseIncidence, input_n4_16S$map_loaded$rich)

# Test and plot
leveneTest(input_n4_16S$map_loaded$rich ~ input_n4_16S$map_loaded$pedigree) # Not homogeneous
leveneTest(input_n4_16S$map_loaded$rich ~ input_n4_16S$map_loaded$rep) # Homogeneous
m <- lmer(rich ~ pedigree + (1|rep), data = input_n4_16S$map_loaded)
summary(m)
Anova(m)
aovtab <- anova(m)
mnull <- lmer(rich ~ 1 + (1|rep), data = input_n4_16S$map_loaded)
anova(mnull, m)
m2 <- mixed(rich ~ pedigree + (1|rep), data = input_n4_16S$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two va
m <- aov(rich ~ rep + pedigree, data = input_n4_16S$map_loaded)
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
         y = max(input_n4_16S$map_loaded$rich)+
           (max(input_n4_16S$map_loaded$rich)-
              min(input_n4_16S$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "16S")

leveneTest(input_n4_16S$map_loaded$shannon ~ input_n4_16S$map_loaded$pedigree) # Homogeneous
leveneTest(input_n4_16S$map_loaded$shannon ~ input_n4_16S$map_loaded$rep) # Homogeneous
m <- lmer(shannon ~ pedigree + (1|rep), data = input_n4_16S$map_loaded)
summary(m)
Anova(m)
anova(m)
mnull <- lmer(shannon ~ 1 + (1|rep), data = input_n4_16S$map_loaded)
anova(mnull, m)
m2 <- mixed(shannon ~ pedigree + (1|rep), data = input_n4_16S$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two va
m <- aov(shannon ~ rep + pedigree, data = input_n4_16S$map_loaded)
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
         y = max(input_n4_16S$map_loaded$shannon)+
           (max(input_n4_16S$map_loaded$shannon)-
              min(input_n4_16S$map_loaded$shannon))/10,
         #y = 5650,
         Dataset = "16S")

label_df_16S <- rbind(t, t1)
alpha_long_16S <- input_n4_16S$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_16S.pdf", width = 8, height = 6)
g1 <- ggplot(alpha_long_16S, aes(reorder(pedigree, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  # geom_text(data = label_df_16S, aes(pedigree, y, label = str_trim(.group)), 
  #           size = 3, color = "black") +
  labs(x = NULL, y = "ASV richness") +
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
min(input_n4_ITS$map_loaded$rich) # 145
max(input_n4_ITS$map_loaded$rich) # 285
mean(input_n4_ITS$map_loaded$rich) # 219.6717
se(input_n4_ITS$map_loaded$rich) # 1.195132
plot(as.factor(input_n4_ITS$map_loaded$pedigree), input_n4_ITS$map_loaded$rich)
plot(input_n4_ITS$map_loaded$DiseaseIncidence, input_n4_ITS$map_loaded$rich)

# Test and plot
leveneTest(input_n4_ITS$map_loaded$rich ~ input_n4_ITS$map_loaded$pedigree) # Homogeneous
leveneTest(input_n4_ITS$map_loaded$rich ~ input_n4_ITS$map_loaded$rep) # Homogeneous
m <- lmer(rich ~ pedigree + (1|rep), data = input_n4_ITS$map_loaded)
summary(m)
Anova(m)
aovtab <- anova(m)
mnull <- lmer(rich ~ 1 + (1|rep), data = input_n4_ITS$map_loaded)
anova(mnull, m)
m2 <- mixed(rich ~ pedigree + (1|rep), data = input_n4_ITS$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two va
m <- aov(rich ~ rep + pedigree, data = input_n4_ITS$map_loaded)
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
         y = max(input_n4_ITS$map_loaded$rich)+
           (max(input_n4_ITS$map_loaded$rich)-
              min(input_n4_ITS$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "ITS")

leveneTest(input_n4_ITS$map_loaded$shannon ~ input_n4_ITS$map_loaded$pedigree) # Homogeneous
leveneTest(input_n4_ITS$map_loaded$shannon ~ input_n4_ITS$map_loaded$rep) # Homogeneous
m <- lmer(shannon ~ pedigree + (1|rep), data = input_n4_ITS$map_loaded)
summary(m)
Anova(m)
anova(m)
mnull <- lmer(shannon ~ 1 + (1|rep), data = input_n4_ITS$map_loaded)
anova(mnull, m)
m2 <- mixed(shannon ~ pedigree + (1|rep), data = input_n4_ITS$map_loaded, method = "PB")
m2
anova(m2)
# No way to get all Sum of Squares for mixed model, so just use two variables
m <- aov(shannon ~ rep + pedigree, data = input_n4_ITS$map_loaded)
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
         y = max(input_n4_ITS$map_loaded$shannon)+
           (max(input_n4_ITS$map_loaded$shannon)-
              min(input_n4_ITS$map_loaded$shannon))/10,
         #y = 5650,
         Dataset = "ITS")

label_df_ITS <- rbind(t, t1)
alpha_long_ITS <- input_n4_ITS$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_ITS.pdf", width = 8, height = 6)
g2 <- ggplot(alpha_long_ITS, aes(reorder(pedigree, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  # geom_text(data = label_df_ITS, aes(pedigree, y, label = str_trim(.group)), 
  #           size = 3, color = "black") +
  labs(x = NULL, y = "ASV richness") +
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
ped_order <- input_n4_16S$map_loaded %>%
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
                            label = c("H = 0.41", "H = 0.38", "H = 0.26", "H = 0.30"),
                            plabel = c("p < 0.05", "p < 0.05", "N.S.D", "N.S.D"),
                            px = c(41.5, 41.5, 41.5, 41.5),
                            py = c(2400, 7.3, 275, 4.25),
                            Dataset = c("16S", "16S", "ITS", "ITS"),
                            name = c("rich", "shannon", "rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
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
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
g3
dev.off()



#### 3. Beta ####



#### 4. Heritability ####
# Calculate heritability of ASVs and high order taxonomies
# Carefully consider filtering cutoffs
# Let's use 25% of samples and 0.1% rel abundance



#### 5. Sclerotinia ####
# Correlations between ASVs and Sclerotinia resistance



##### 6. Networks ####
# Use SPIEC-EASI to build a network for Archaea/Bacteria and Fungi
# Again, carefully consider filterting cutoffs
# Here use the top 100 most abundant taxa for each



#### 7. SNPs/Microbes ####
# This is a really cool aspect of having BOTH microbial data and plant genomic data
# Look for associations between SNPs and ASVs
# Again use the abund/ubiq ASVs
# Need help from others (Cloe, Brian etc. to run this)
# Most is done outside of R