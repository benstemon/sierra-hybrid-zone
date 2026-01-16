# General use
library(tidyverse)

setwd("~/6.pigments/")


# PARENT INFORMATION FOR BACKGROUND SNPS
# Parent 1 (dav)
snptable_p1 <- read_delim("dav_allele_freq.frq", delim = "\t") %>%
 rename(chr = "CHROM", ps = "POS") %>%
 select(c(chr, ps, allele_1)) %>%
 mutate(allele_1 = as.numeric(str_extract(allele_1, "[0-9.]+")))

# Parent 2 (new)
snptable_p2 <- read_delim("new_allele_freq.frq", delim = "\t") %>%
 rename(chr = "CHROM", ps = "POS") %>%
 select(c(chr, ps, allele_1)) %>%
 mutate(allele_1 = as.numeric(str_extract(allele_1, "[0-9.]+")))

# Combine these dfs and find AFD for each site
afd_parent_table <- full_join(snptable_p1, snptable_p2, by = c("chr", "ps")) %>%
 mutate(AFD = abs(allele_1.x - allele_1.y)) %>%
 mutate(chr = str_remove(chr, "scaffold_"))
write_delim(afd_parent_table, file = "AFD-table-parents.txt")