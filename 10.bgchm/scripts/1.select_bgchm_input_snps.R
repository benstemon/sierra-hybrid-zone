# This script generates the SNPs to be used in genomic clines analyses
# 1. Collect SNPs in trait loci with MAF ≥ 0.1. Filter as necessary
# 2. Collect background SNPs with MAF ≥ 0.1 and > 20kb from any trait locus or GWAS SNP

library(tidyverse)

setwd("~/10.bgchm/")



# SETUP
################################################
# Generate MAF file with no parents
t <- read_delim("ulmm.NOPARENTS_normalized_LDpruned.modified_hue.assoc.txt",
                delim = "\t") %>%
  select(rs, af) %>%
  separate(rs, into = c("tmp1", "chr", "ps"), sep = ":", remove = F) %>%
  mutate(chr = as.numeric(str_remove(chr, "scaffold_")),
         ps = as.numeric(ps)) %>%
  select(-tmp1, -rs)

tt <- t %>%
  mutate(maf = case_when(af > 0.5 ~ 1-af,
                         af < 0.5 ~ af,
                         af == 0.5 ~ 0.5,
                         TRUE ~ NA)) %>%
  select(-af)

write_csv(tt, file = "data/MAF-snps-NOPARENTS.csv")

# Read in parent AFD file generated from vcftools
# Combine with hybrid MAF file generated from NOPARENTS GEMMA output (unfiltered)
# Get this from the online supplemental data repo
allsnps <- read_delim("AFD-table-parents.txt",
                      delim = " ") %>%
  select(chr, ps, AFD) %>%
  left_join(read_csv("MAF-snps-NOPARENTS.csv"),
            by = c("chr", "ps"))

# Read in candidate trait loci
# file must be filtered for duplicate loci
# Get this from the online supplemental data repo
traitloci <- read_csv("candidate_loci_20kb-sierras-NOPARENTS-normalized-LDpruned-DUPFILTER-perm0.0001.csv")

# Read in list of SNPs for GWAS hits
# Then merge it to the allsnps object
sigsnps <- read_delim("all-significant-snps-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned-WITHPVE.tsv",
                      delim = "\t") %>%
  mutate(phenotype = str_remove(phenotype, "PC2covariate.")) %>%
  select(c(chr, ps, phenotype)) %>%
  filter(!duplicated(across(c(chr, ps)))) %>%
  left_join(., allsnps, by = c("chr", "ps"))

################################################

# Choosing SNPs, filtering, etc.
################################################
# 1. Collect SNPs in trait loci with MAF ≥ 0.1. Filter as necessary
# Convert to data.tables
setDT(allsnps)
setDT(traitloci)

# Add required columns to allsnps to represent it as a 1bp interval
allsnps[, `:=`(bpstart = ps, bpend = ps)]

# Set keys for both data.tables
setkey(allsnps, chr, bpstart, bpend)
setkey(traitloci, chr, bpstart, bpend)

# Find SNPs in trait loci, filtering by MAF
TRAITsnps <- foverlaps(allsnps, traitloci, nomatch = 0L) %>%
  select(-c(i.bpstart, i.bpend, modelblast)) %>%
  filter(maf >= 0.1)

# Then subsample 1086 and 2532 so as not to completely dominate trait loci
# Of 60,031 SNPs, 56483 (94%) are from those two scaffolds
# We subsample to roughly 60% total
TRAITsnps_subsampled <- TRAITsnps %>%
  group_by(chr) %>%
  group_modify(~ {frac <- case_when(.y$chr %in% c(1086, 2532) ~ .10,
                                    TRUE ~ 1)
  sample_frac(.x, size = frac)}) %>%
  ungroup()



# 2. Collect background SNPs with MAF ≥ 0.1 and > 20kb from any GWAS SNP
# Convert to data.tables
setDT(sigsnps)

# Make 40kb windows around significant SNPs
sigsnps[, `:=`(bpstart = ps - 20000, bpend = ps + 20000)]

# Set keys for sigsnps. Already done for allsnps
setkey(sigsnps, chr, bpstart, bpend)

# Find overlap -- all SNPs within 20kb of GWAS hits
OVERLAPsnps <- foverlaps(allsnps, sigsnps, nomatch = 0L) %>%
  distinct(chr, i.ps, .keep_all = T) %>%
  select(chr, i.ps, i.AFD, i.maf, i.bpstart, i.bpend) %>%
  rename("ps" = i.ps)

# Anti-join -- all SNPs NOT within the 20kb GWAS window
# ALSO, no SNPs within trait loci
# Filter by MAF and subsample by chromosome
BACKGROUNDsnps_subsampled <- anti_join(allsnps, OVERLAPsnps, by = c("chr", "ps")) %>%
  anti_join(., TRAITsnps, by = c("chr", "ps")) %>%
  filter(maf >= 0.1) %>%
  group_by(chr) %>%
  sample_frac(size = 0.02)
################################################

# Writing
################################################
# Write full trait snps
write_csv(TRAITsnps_subsampled, file = "data/data-fullinfo-traitsnps-maf0.1-2scaff0.1-bgchm.csv")

# Write full background snps
write_csv(BACKGROUNDsnps_subsampled, file = "data/data-fullinfo-backgroundsnps-maf0.1-sub.02-bgchm.csv")


# Merge objects and write chr, bp to disk for bgchm
finalsnps <- bind_rows(TRAITsnps_subsampled %>% select(chr, ps),
                       BACKGROUNDsnps_subsampled %>% select(chr, ps)) %>%
  arrange(chr, ps) %>%
  mutate(chr = paste0("scaffold_", chr, ""))

# Ensure no funny business -- no duplicated SNPs  
finalsnps[which(duplicated(finalsnps)),]

# Write final snps to output 
write_delim(finalsnps,
            file = "data/subsampled-maf0.1-bgcinput.txt",
            delim = "\t", col_names = F)
################################################

