# Reassign Taxonomy for Carrington 2017 Study
# by Cliff Bueno de Mesquita, May 2024
# Note: This was performed on the Innes Server.
# Including this in the GitHub repo to show what was done but not including all inputs.

# For 16S use SILVA v138.1 (released August 27, 2020)
# For ITS use UNITE v10.0 (released April 4, 2024)
# This will be a substantial update since Jason et al. used Greengenes 13.5 and UNITE 8.
# Lot's of things changed since then! 
# Need to get repset and ASV table, merge them, then assign taxonomy with DADA2
# Need to get sequences as column names and samples as rownames for assignTaxonomy()

#### Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(dada2)
library(microseq)

# Directory
setwd("~/prev_SunflowerRhizosphereStudy/")



#### 16S ####
# Import data
seqtab_16S_orig <- read.delim("OTU_Tables/16S/zotutab_wTax_noChloroMito.txt", skip = 1)
seqtab_16S <- read.delim("OTU_Tables/16S/zotutab_wTax_noChloroMito.txt", skip = 1) %>%
  dplyr::select(-taxonomy)
repset_16S <- readFasta("OTU_Tables/16S/16S_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% seqtab_16S$X.OTU.ID)
seqtab_16S <- seqtab_16S %>%
  left_join(., repset_16S, by = c("X.OTU.ID" = "Header")) %>%
  dplyr::select(-X.OTU.ID) %>%
  column_to_rownames(var = "Sequence") %>%
  t()

# Assign Taxonomy
tax_16S <- assignTaxonomy(seqtab_16S, "~/Sunflower_GxE/OTU_DBs/silva_nr99_v138.1_wSpecies_train_set.fa.gz", 
                          tryRC = TRUE,
                          multithread = TRUE)

# Rest of the pipeline
# Check match
nrow(tax_16S)
ncol(seqtab_16S)
sum(rownames(tax_16S) != colnames(seqtab_16S)) # Order is same so can easily reassign OTU ID

# Save file and then transfer back to personal computer
saveRDS(tax_16S, "taxonomy/tax_final_SILVA138.rds")

# Flip table
seqtab.t <- as.data.frame(t(seqtab_16S))

# Add ASV numbers to table (use Jason's ids!)
rownames(seqtab.t) <- seqtab_16S_orig$X.OTU.ID

# Add ASV numbers to taxonomy
taxonomy_16S <- as.data.frame(tax_16S) %>%
  mutate(ASV_ID = seqtab_16S_orig$X.OTU.ID)
taxonomy_for_mctoolsr_16S <- unite(taxonomy_16S, "taxonomy", 
                                   c("Kingdom", "Phylum", "Class", "Order",
                                     "Family", "Genus", "Species", "ASV_ID"),
                                   sep = ";")

# Merge taxonomy and table
seqtab_wTax <- seqtab.t %>%
  mutate(taxonomy = taxonomy_for_mctoolsr_16S$taxonomy) %>%
  rownames_to_column(var = "#ASV_ID") %>%
  dplyr::select("#ASV_ID", everything())

# Set name of table in mctoolsr format and save
out_fp <- "taxonomy/seqtab_wTax_mctoolsr_16S.txt"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Compare taxonomies
tax_comp_16S <- seqtab_16S_orig %>%
  dplyr::select(X.OTU.ID, taxonomy) %>%
  mutate(OldTax = taxonomy) %>%
  dplyr::select(-taxonomy) %>%
  mutate(NewTax = taxonomy_for_mctoolsr_16S$taxonomy)



#### ITS ####
# Import data
seqtab_ITS_orig <- read.delim("OTU_Tables/ITS/zotutab_wTax_noChloroMito.txt", skip = 1)
seqtab_ITS <- read.delim("OTU_Tables/ITS/zotutab_wTax_noChloroMito.txt", skip = 1) %>%
  dplyr::select(-taxonomy)
repset_ITS <- readFasta("OTU_Tables/ITS/ITS_rep_set_zotus_filt_relabeled.fa") %>%
  filter(Header %in% seqtab_ITS$X.OTU.ID)
seqtab_ITS <- seqtab_ITS %>%
  left_join(., repset_ITS, by = c("X.OTU.ID" = "Header")) %>%
  dplyr::select(-X.OTU.ID) %>%
  column_to_rownames(var = "Sequence") %>%
  t()

# Assign Taxonomy
tax_ITS <- assignTaxonomy(seqtab_ITS, "~/Sunflower_GxE/OTU_DBs/sh_general_release_dynamic_s_all_04.04.2024.fasta", 
                          tryRC = TRUE,
                          multithread = TRUE)

# Rest of the pipeline
# Check match
nrow(tax_ITS)
ncol(seqtab_ITS)
sum(rownames(tax_ITS) != colnames(seqtab_ITS)) # Order is same so can easily reassign OTU ID

# Save file and then transfer back to personal computer
saveRDS(tax_ITS, "taxonomy/tax_final_UNITE10.rds")

# Flip table
seqtab.t <- as.data.frame(t(seqtab_ITS))

# Add ASV numbers to table (use Jason's ids!)
rownames(seqtab.t) <- seqtab_ITS_orig$X.OTU.ID

# Add ASV numbers to taxonomy
taxonomy_ITS <- as.data.frame(tax_ITS) %>%
  mutate(ASV_ID = seqtab_ITS_orig$X.OTU.ID)
taxonomy_for_mctoolsr_ITS <- unite(taxonomy_ITS, "taxonomy", 
                                   c("Kingdom", "Phylum", "Class", "Order",
                                      "Family", "Genus", "Species", "ASV_ID"),
                                    sep = ";")

# Merge taxonomy and table
seqtab_wTax <- seqtab.t %>%
  mutate(taxonomy = taxonomy_for_mctoolsr_ITS$taxonomy) %>%
  rownames_to_column(var = "#ASV_ID") %>%
  dplyr::select("#ASV_ID", everything())
  
# Set name of table in mctoolsr format and save
out_fp <- "taxonomy/seqtab_wTax_mctoolsr_ITS.txt"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Compare taxonomies
tax_comp_ITS <- seqtab_ITS_orig %>%
  dplyr::select(X.OTU.ID, taxonomy) %>%
  mutate(OldTax = taxonomy) %>%
  dplyr::select(-taxonomy) %>%
  mutate(NewTax = taxonomy_for_mctoolsr_ITS$taxonomy)



# End Script