setwd("6.pigments/")
library(tidyverse)

# Read in gemma snps
gemmasnps <- read_csv("data/significant-gemmasnps.csv")

# Read in AFD table
# These must be obtained from supplemental data repo
afds <- read_delim("AFD-table-parents.txt")

# Merge
outtable <- gemmasnps %>%
  left_join(., afds, by = c("chr", "ps"))
write_csv(outtable, "data/merged-afds-gemmasnps.csv")


# Calculate % fixed differenced for each trait
sumtable <- outtable %>%
  group_by(phenotype, chr) %>%
  summarise(n_total = n(),
            n_fixed = sum(AFD == 1, na.rm = TRUE),
            pct_fixed = 100 * n_fixed / n_total) %>%
  filter(n_total > 5)
write_csv(sumtable, "summary-afds-gemmasnps.csv")  
