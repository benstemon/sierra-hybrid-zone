# This script first extracts significant SNPs for each phenotype
# Then it identifies candidate genes based on proximity to those SNPs

library(tidyverse)

####
# ULMM LOCI
####
# Set directory
setwd("~/4.gemma/")

#Phase 1. Extract significant SNPs
################################################################################
# List of all the output files
# These must be obtained from supplemental data repo
filelist = list.files(getwd(), pattern = "*.assoc.txt")

# Path to directory with permutation results:
# These must be obtained from supplemental data repo
permfiles = list.files("quantiles_pvals/", full.names = T)


# For loop to filter to only SNPs above the significance threshold
for (i in 1:length(filelist)){
  
  # Set up the list to hold significant SNPs, if the first iteration
  if(i == 1){
    sigsnplist <- list()
  }
  
  # Specify the infile
  infile <- filelist[i]
  
  # Specify the phenotype being examined
  phenotype <- str_extract(infile, "(?<=LDpruned\\.).*?(?=\\.assoc)") %>%
    str_remove(., "avg_")
  
  # Read in the permutation file
  perminfo <- read_delim(permfiles[which(str_detect(string = permfiles, pattern = phenotype))], delim = "\t")
  
  # Set p-values
  # First is bonferroni based on number of SNPs
  # Second is based on permutations test, top 0.0001% of top 1% of permutation p-values
  p_bonf = -log10(.05/3021696)
  p_perm_0.0001 = -log10(as.numeric(perminfo[1,2]))
  
  # Read in the GWAS results
  # # Clean results, and only keep snps with transformed p > cutoff
  # sigsnps <- read_delim(infile, delim = "\t") %>%
  #   na.exclude() %>%
  #   mutate(chr = str_extract_all(rs, "(?<=:)(.*?)(?=:)")) %>%
  #   mutate(chr = as.numeric(str_remove_all(chr, "scaffold_"))) %>%
  #   mutate(ps = as.numeric(str_extract_all(rs, "(?<=:)[^:]+$"))) %>%
  #   select(., -c(rs, n_miss)) %>%
  #   mutate(ptransform = -log10(p_wald)) %>%
  #   filter(ptransform > p_perm_0.0001) %>%
  #   mutate(phenotype = phenotype)
  
  # Clean results, and only keep snps with transformed p > cutoff
  sigsnps <- read_delim(infile, delim = "\t") %>%
    na.exclude() %>%
    mutate(chr = str_extract_all(rs, "(?<=:)(.*?)(?=:)")) %>%
    mutate(chr = as.numeric(str_remove_all(chr, "scaffold_"))) %>%
    mutate(ps = as.numeric(str_extract_all(rs, "(?<=:)[^:]+$"))) %>%
    select(., -c(rs, n_miss)) %>%
    mutate(ptransform = -log10(p_wald)) %>%
    filter(ptransform > p_bonf) %>%
    mutate(phenotype = phenotype)
  
  # Save this to the significant SNP list
  if(nrow(sigsnps) > 0){
    sigsnplist[[i]] <- sigsnps
  }
}

# Concatenate all the significant SNP data frames and write to file
all_sigsnps <- bind_rows(sigsnplist)
# write_delim(all_sigsnps, file = "all-significant-snps-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned.tsv", delim = "\t",
#             quote = "none")

write_delim(all_sigsnps, file = "all-significant-snps-pwald-bonferroni-sierras-NOPARENTS-normalized-LDpruned.tsv", delim = "\t",
            quote = "none")
################################################################################


# Phase 2. Find candidate loci
################################################################################
# Do the significant SNPs occur within any interesting gene models?

# Read in the significant SNPs
# These must be obtained from supplemental data repo
sigsnps <- read_delim("all-significant-snps-pwald-bonferroni-sierras-NOPARENTS-normalized-LDpruned.tsv", delim = "\t")

# Read in the P. davidsonii functional annotations
# Note that these are already prefiltered to single isoform, and exclude other things
# So that only gene models and some tRNAs remain
# This is a custom .gff file made from the P. davidsonii reference genome
# Must be obtained from supplemental data repo
annotations <- read_delim("single_isoform_davidsonii_FUNCTIONAL-INCLUDED.gff",
                          delim = "\t", col_names = F) %>%
  select(-c(2,6,7,8)) %>%
  rename(chr=X1, model=X3, bpstart=X4, bpend=X5, info=X9) %>%
  mutate(chr = str_remove(chr, "scaffold_"))

# Keep just gene models
genes <- annotations %>%
  filter(model == "mRNA")
  
# Keep just exons
exons <- annotations %>%
  filter(model == "exon")


# Find significant SNPs that fall within these gene models and exons
# Also extract the info column for those significant hits
sites_in_genic <- sigsnps %>%
  mutate(in_genemodel = map2_lgl(chr, ps, ~ any(.x == genes$chr & .y >=genes$bpstart & .y <= genes$bpend))) %>%
  filter(in_genemodel == TRUE) %>%
  mutate(in_exon = map2_lgl(chr, ps, ~ any(.x == exons$chr & .y >= exons$bpstart & .y <= exons$bpend))) %>%
  mutate(matched_row_index = map2_int(chr, ps, ~ {
    matching_rows <- which(.x == genes$chr & .y >= genes$bpstart & .y <= genes$bpend)
    if (length(matching_rows) > 0) matching_rows[1] else NA
  })) %>%
  mutate(info = case_when(!is.na(matched_row_index) ~ genes$info[matched_row_index],
                          TRUE ~ NA)) %>%
  select(-matched_row_index) %>%
  mutate(modelID = str_extract(info, "^[^;]+"),
         modelID = str_remove(modelID, "ID="),
         modelblast = str_extract(info, "Note=[^;]+"),
         modelblast = str_remove(modelblast, "Note=")) %>%
  select(-info)
  
  
# Write this to output for further analysis
# write_delim(sites_in_genic,
#             file = "all-significant-snps-in-genes-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned.tsv", delim = "\t")
write_delim(sites_in_genic,
            file = "all-significant-snps-in-genes-pwald-bonferroni-sierras-NOPARENTS-normalized-LDpruned.tsv", delim = "\t")



# Find all gene models within 20kb of significant SNPs
candidate_loci <- sigsnps %>%
  mutate(matched_rows = map2(chr, ps, ~ {
    which(.x == genes$chr & .y >= genes$bpstart - 20000 & .y <= genes$bpend + 20000)
    })) %>%
  unnest_longer(matched_rows, values_to = "gene_row") %>%
  filter(!is.na(gene_row)) %>%
  mutate(gene_chr = genes$chr[gene_row],
         bpstart = genes$bpstart[gene_row],
         bpend = genes$bpend[gene_row],
         geneID = genes$info[gene_row]) %>%
  select(chr, bpstart, bpend, phenotype, geneID) %>%
  distinct(phenotype, geneID, .keep_all = T) %>%
  mutate(modelID = str_extract(geneID, "^[^;]+"),
         modelID = str_remove(modelID, "ID="),
         modelblast = str_extract(geneID, "Note=[^;]+"),
         modelblast = str_remove(modelblast, "Note=")) %>%
  select(-geneID)

# Write list of candidate loci to disk -- this keeps all hits within phenotypes
# Even if there are duplicate gene hits across phenotypes (colocalization)
# write_csv(candidate_loci, file = "candidate_loci_20kb-sierras-NOPARENTS-normalized-LDpruned-NODUPFILTER-perm0.0001.csv")
write_csv(candidate_loci, file = "candidate_loci_20kb-sierras-NOPARENTS-normalized-LDpruned-NODUPFILTER-bonferroni.csv")

# Filter again, keeping only one example of each possible gene (no colocalization)
t <- candidate_loci %>%
  distinct(modelID, .keep_all = T)
# write_csv(t, file = "candidate_loci_20kb-sierras-NOPARENTS-normalized-LDpruned-DUPFILTER-perm0.0001.csv")

################################################################################


# Phase 3: Write a summary table to broadly summarize GWAS hits for each trait
################################################################################
indata <- read_delim("all-significant-snps-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned.tsv",
                     delim = "\t")

# For each trait, find min-max bp for each chromosome with significant hits
summary_hits <- indata %>%
  group_by(phenotype, chr) %>%
  summarize(min_ps = min(ps, na.rm = TRUE),
            max_ps = max(ps, na.rm = TRUE),
            snp_count = n(),
            .groups = "drop") %>%
  arrange(phenotype, chr)

write.csv(summary_hits, file = "RESULT-GWAS-sighits-summary-permutation.csv")
################################################################################


# Calculate PVE for top SNPs in association peaks

# Read in top SNPs
indata <- read_delim("all-significant-snps-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned.tsv",
                     delim = "\t")

# Code to calculate PVE (from Shim et al. 2015)
indata$pve = (2 * (indata$beta^2) * (indata$af) * (1 - indata$af)) / 
  (2 * (indata$beta^2) * (indata$af) * (1 - indata$af) + (indata$se^2) * 2 * 283 * (indata$af) * (1 - indata$af))

# Write to file
write_delim(indata, file = "all-significant-snps-pwald-perm0.0001-sierras-NOPARENTS-normalized-LDpruned-WITHPVE.tsv",
          delim = "\t")

