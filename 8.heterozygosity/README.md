# Heterozygosity inference

Additionally uses data from:
[1.phenotype_curation](~/1.phenotype_curation/)
[10.bgchm](~/10.bgchm/)
AFD table from supplemental data repo


#### 1. Generate file to only analyze fixed opposite parent SNPs in hybrids
* see [`1.prepare-admixbins-hybrids-heterozygosity.R`](scripts/1.prepare-admixbins-hybrids-heterozygosity.R)

#### 2. Estimate heterozygosity metrics
* see [`2.ARRAY_excess_heterozygosity-sierras-25-75.sh`](scripts/2.ARRAY_excess_heterozygosity-sierras-25-75.sh)
* uses vcftools, the imputed .vcf, and output from step 1

#### 3. Plot heterozygosity across the genome
* see [`3.plot_heterozygosity_excess-allsnps.R`](scripts/3.plot_heterozygosity_excess-allsnps.R)

