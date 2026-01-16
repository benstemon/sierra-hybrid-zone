library(bgchm)
library(tidyverse)

setwd("~/Desktop/testing_bgchm/clineresults/")

# Combining .rda plots from parallel runs
################################################################################
# List .rda outfiles
rdafiles <- list.files(pattern="clinesOut")

# Load the first .rda file and create storage objects
load(rdafiles[1])

# Number of loci (hybrid matrix from .rda)
nloci <- ncol(hybrids)

# Set up the output matrices
finalgradient <- matrix(nrow = nloci, ncol = 3)
finalcenter <- matrix(nrow = nloci, ncol = 3)

# For loop over all .rda files to fill in estimates for gradient and center
# Note that the .rda objects may not be loaded in the correct numeric order
# so SNP locations must be obtained from the .rda object itself
for(i in 1:length(rdafiles)){
        # Load each .rda file
        load(rdafiles[i])
        
        # Place gradient and center in corresponding rows (SNP index from .rda)
        finalgradient[snpstart:snpend,] <- gradientout
        finalcenter[snpstart:snpend,] <- centerout
}

# Save the combined clines object
save(list = ls(), file = "combinedClines.rda")

# Write the hybrid index results (slightly diff format from h_est)
 write_delim(h, file = "hybrid_index-sierras.txt", delim = "\t")

# Impose sum-to-zero constraints
sz_out <- sum2zero(center = finalcenter, v = finalgradient, transform = T)

# Write out CIs for center and gradient for regular and sum-to-zero constrained clines
write.csv(finalcenter, "result-cline-centers.csv")
write.csv(finalgradient, "result-cline-gradients.csv")
write.csv(sz_out$center, "result-cline-centers-s2z.csv")
write.csv(sz_out$gradient, "result-cline-gradients-s2z.csv")
################################################################################





