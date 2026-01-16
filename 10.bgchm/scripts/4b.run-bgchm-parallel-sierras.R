# This script is parallelized to run bgchm on chunks of SNPs specified in batch script

library(bgchm)

# Load in each of the data sets
hybrids <- as.matrix(read.csv("genotypes-recoded-hybrids.tsv", sep = "\t", row.names = 1))
genp0 <- as.matrix(read.csv("genotypes-recoded-p0.tsv", sep = "\t", row.names = 1))
genp1 <- as.matrix(read.csv("genotypes-recoded-p1.tsv", sep = "\t", row.names = 1))


# Read in prior estimates of hybrid indexes, and specify cline SD point estimates
h <- read.table("h_est.txt", header = T)
sdc <- 0.6315415
sdv <- 0.1289992


# The parallelization on the github is a little clunky, so I'm doing it this way
# Parameters specifying the range of SNPs are done on the job submission side
myargs <- commandArgs(trailingOnly = T)
snpstart <- as.numeric(myargs[1])
snpend <- as.numeric(myargs[2])
batchnumber <- myargs[3]

# Subset genotype objects
sGhyb <- hybrids[,snpstart:snpend]
sGP0 <- genp0[,snpstart:snpend]
sGP1 <- genp1[,snpstart:snpend]

# Make table to write with individual SDs and convergence metrics
tableout <- data.frame(center_mean = numeric(),
                       center_SD = numeric(),
                       center_ESS = numeric(),
                       center_Rhat = numeric(),
                       gradient_mean = numeric(),
                       gradient_SD = numeric(),
                       gradient_ESS = numeric(),
                       gradient_Rhat = numeric())

# For loop to estimate cline parameters one SNP at a time
for(i in 1:ncol(sGhyb)){
  # Estimate genomic cline, using previous estimates for hybrid index and cline SDs
  tmpgc <- est_genocl(Gx = sGhyb[,i], G0 = sGP0[,i], G1 = sGP1[,i], H = h[,1],
                      model = "genotype", ploidy = "diploid", hier = F,
                      SDc = sdc, SDv = sdv)
  
  # Add to the outfiles
  # If on the first iteration, generate new data.frames
  # Otherwise, append
  if(i == 1){
    gradientout <- tmpgc$gradient
    centerout <- tmpgc$center
  }else{
    gradientout <- rbind(gradientout, tmpgc$gradient)
    centerout <- rbind(centerout, tmpgc$center)
  }
  
  # Now summarize the object
  invisible(show(tmpgc$gencline_hmc))
  
  
  # Add mean, sd, ESS and Rhat for center and gradient
  tableout[i,1:2] <- tmpgc$gencline_hmc@.MISC$summary$msd[1] # center mean and sd
  tableout[i,3] <- tmpgc$gencline_hmc@.MISC$summary$ess[1] # center ESS
  tableout[i,4] <- tmpgc$gencline_hmc@.MISC$summary$rhat[1] # center Rhat
  tableout[i,5:6] <- tmpgc$gencline_hmc@.MISC$summary$msd[2] # slope mean and sd
  tableout[i,7] <- tmpgc$gencline_hmc@.MISC$summary$ess[2] # slope ESS
  tableout[i,8] <- tmpgc$gencline_hmc@.MISC$summary$rhat[2] # slope Rhat
}

# Add batch number and write the table to output
tableout$batchid <- batchnumber
write.csv(tableout, file = paste0("convergencemetrics_batch", batchnumber, ".csv", sep = ""))

# Save the R object with batch ID included
save(list = ls(), file = paste("clinesOut", batchnumber, ".rda", sep = ""))


