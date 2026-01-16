# Pigment analysis

Additionally uses data from:
[1.phenotype_curation](~/1.phenotype_curation), [7.ancestryHMM](~/7.ancestryHMM), and [9.visualmodels](~/9.visualmodels), 
particularly for making some of the plots (step 5)


#### 1. Generate parent allele frequency files
* see [`1.generate-parent-freqs.sh`](scripts/1.generate-parent-freqs.sh)
* Uses vcftools and the imputed .vcf

#### 2. Generate parent AFD table, for all SNPs
* see [`2.create-parent-AFD-tables.R`](scripts/2.create-parent-AFD-tables.R)
* uses files generated from step 1. These must be either generated on your own or downloaded from the supplemental data repo

#### 3. Combine AFD values for all significant GEMMA SNPs
* see [`3.find-afds-gemmasnps.R`](scripts/3.find-afds-gemmasnps.R)

#### 4. Grab genotype data at the DFR-like locus
* see [`4.grab_DFR_genotypes.sh`](scripts/4.grab_DFRlike_genotypes.sh)
* Uses [`ANR_GEMMAhits.txt`](data/DFRlike_GEMMAhits.txt) to generate [`final_ANR_genotype_matrix.tsv`](data/final_DFRlike_genotype_matrix.tsv)

#### 5. Genotype at DFR-like and F3'5'H, and ANOVAs
* see [`5.plot-pigment-f3p5ph-genotype.R`](scripts/5.plot-pigment-f3p5ph-genotype.R)
* uses [`pigment_intensity_values.txt`](data/pigment_intensity_values.txt), [`pcr-f3p5ph-genotypes.csv`](data/pcr-f3p5ph-genotypes.csv), [`ANR_GEMMAhits.txt`](data/DFRlike_GEMMAhits.txt), [`final_ANR_genotype_matrix.tsv`](data/final_DFRlike_genotype_matrix.tsv), and the full AFD table, either generated from step 2 or from supplemental data repo

