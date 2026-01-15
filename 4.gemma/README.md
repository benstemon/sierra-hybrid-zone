# GWAS

Uses:
2.QC_MAPPING_VARIANTS
3.admixture-pca


#### 1. Clean up samples included and transform traits as necessary
* see [`1.prepare_gemma_phenotypes-sierras.R`](scripts/1.prepare_gemma_phenotypes-sierras.R)
* Reorders phenotypes to match sample order in vcf
* Filters out parent taxa and the lone sample with genetic data but no phenotypic data
* Tests each phenotype for normality and, if significantly different from normality, attempts several different transformations. If no transformation makes it better, keep default.
* Generates a list of excluded samples [`exluded_samples_GWAS.txt`](data/exluded_samples_GWAS.txt) which is used downstream to generate the relatedness matrix without including samples not included in the GWAS (parent species and the one missing phenotype hybrid sample)

#### 2. Phase and impute data
* see [`2.phase-impute-beagle-gemma-sierras.sh`](scripts/2.phase-impute-beagle-gemma-sierras.sh)
* Uses `finalfilter.vcf.gz` -- see supplemental data repo
* Uses BEAGLE and bcftools

#### 3. Convert .vcf to bimbam and generate GEMMA relatedness matrix
* see [`3.vcf2bimbam-generate-GTmatrix-sierras-NOPARENTS-LDpruned.sh`](scripts/3.vcf2bimbam-generate-GTmatrix-sierras-NOPARENTS-LDpruned.sh)
* Converts the .vcf into bimbam format, and then generates the genotype matrix, excluding parents and the one missing phenotype hybrid sample
* Uses bcftools, plink2, gemma, and qctool
* Generated relatedness matrix is found here: [`GTMATRIX-NOPARENT-LDpruned.cXX.txt`](data/GTMATRIX-NOPARENT-LDpruned.cXX.txt)

#### 4. Run the ULMM in GEMMA
* see [`4.run-gemma-ULMM-sierras-NOPARENTS-normalized-LDpruned.sh`](scripts/4.run-gemma-ULMM-sierras-NOPARENTS-normalized-LDpruned.sh)
* Uses GEMMA, the phenotype files, and the relatedness matrix

#### 5. Prefilter GEMMA output
* see [`5.prefilter-gemma-output.py`](scripts/5.prefilter-gemma-output.py)
* The output files are gigantic, to the point it's nearly impossible to plot anything with R. This script removes markers that are exceedingly not significant to reduce the computational load of plotting. Default is p > 0.3.

#### 6. GWAS permutations
* see [`6a.permutations-wrapper.sh`](scripts/6a.permutations-wrapper.sh) and [`6b.permutations-alltraits.sh`](scripts/6b.permutations-alltraits.sh)
* Main code `6b.permutations-alltraits.sh` shuffles phenotype values and reruns GEMMA for a phenotype. Built-in default to do 100 iterations as an array. It then extracts the top 1% of values from each iteration.
* The wrapper code `6a.permutations-wrapper.sh` just submits the main code for each phenotype designated

#### 7. Generate top quantiles from permutation results
* see [`7.generate-quantiles.py`](scripts/7.generate-quantiles.py)

#### 8. Plot ULMM GWAS results
* see [`8.plot_gwas_ULMM-sierras.R`](scripts/8.plot_gwas_ULMM-sierras.R)

#### 9. Extract significant GEMMA results
* see [`9.extract_significant_snps-ULMM.R`](scripts/9.extract_significant_snps-ULMM.R)


