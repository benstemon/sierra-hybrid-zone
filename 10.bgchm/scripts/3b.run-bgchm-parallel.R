# This script estimates cline SDs and hybrid indexes with a random
# subset of the full data

library(bgchm)

# Load in each of the data sets
hybrids <- as.matrix(read.csv("1000snp-recoded-hybrids.tsv", sep = "\t", row.names = 1))
genp0 <- as.matrix(read.csv("1000snp-recoded-p0.tsv", sep = "\t", row.names = 1))
genp1 <- as.matrix(read.csv("1000snp-recoded-p1.tsv", sep = "\t", row.names = 1))


# Estimate hybrid indexes with default HMC settings
# Concurrently estimate allele frequencies
h_out <- est_hi(Gx = hybrids, G0 = genp0, G1 = genp1,
                model = "genotype", ploidy = "diploid")

# Write hybrid index estimate to file to easily access later
write.table(file = "h_est.txt", h_out$hi, row.names = F, quote = F)


# Fit hierarchical genomic cline model for subset of loci with the estimated HIs
# Also estimate parental allele frequencies concurrently
# Uses 4000 iterations and 2000 burnin to ensure appropriate ESS
# Main goal is to estimate the cline SDs
gc_out <- est_genocl(Gx = hybrids, G0 = genp0, G1 = genp1, H = h_out$hi[,1],
                     model = "genotype", ploidy = "diploid", hier = T, n_iters = 4000)


# Write cline SDs to file, as they are used as point estimates later
sink(file = "SD-point-estimates.txt")
print("SDc")
print(gc_out$SDc)
print("SDv")
print(gc_out$SDv)
sink()


# Save HI convergence results to quickly scan n_eff and Rhat
options(max.print = 3000)
sink(file = "hi_convergence_results.txt")
h_out$hi_hmc
sink()

# Same for genomic clines
options(max.print = 300000)
sink(file = "genocline_convergence_results.txt")
gc_out$gencline_hmc
sink()

# Save HI and clines files to further assess convergence if necessary
save(h_out, file = "forconvergence_HI.RData")
save(gc_out, file = "forconvergence_GC.RData")


