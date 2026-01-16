library(tidyverse)

setwd("~/8.heterozygosity/")

# 1. Filter to only include hybrids
# Read in admixture information
admixture <- read_csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  filter(!is.na(K1),
         K1 >= 0.25 & K1 <= 0.75) %>%
  select(newID, K1) %>%
  mutate(newID = paste("BWS_", newID, sep = ""))

write_delim(admixture[,1], "data/hybridsonly-25-75.txt", delim = "\t", col_names = F)

# 2. Filter to only include fixed SNPs with --positions
# These must be obtained from supplemental data repo
afdsnps <- read_delim("AFD-table-parents.txt") %>%
  filter(AFD == 1) %>%
  select(chr, ps) %>%
  mutate(chr = paste("scaffold_", chr, sep = ""))

write_delim(afdsnps, "data/insnps_AFD1.txt", delim = "\t", col_names = F)
