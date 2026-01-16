# Bayesian Genomic Clines

Additionally uses data from:
Gemma results, AFD table, 

#### 1. Generate list of SNPs to use for genomic clines analysis
* see [`1.select_bgchm_input_snps.R`](scripts/1.select_bgchm_input_snps.R)
* Only considers SNPs with MAF ≥ 0.1
* Finds all SNPs in trait loci, and after subsampling the overrepresented pseudochromosomes, collect SNPs from trait loci and random sample from across the genome, assuming MAF ≥ 0.1 and > 20kb from any GWAS SNP.

#### 2. Prepare data for input into bgchm
* see [`2a.prepare-bgchm-input-sierras.sh`](scripts/2a.prepare-bgchm-input-sierras.sh) and companion script [`2b.vcf2bgchm-gzipupdate_v2.py`](scripts/2b.vcf2bgchm-gzipupdate_v2.py)
* The first script uses bcftools to filter the imputed .vcf to only include the designated SNPs.
* The second script is the main workhorse that converts the .vcf to bgchm format. Works on zipped, phased .vcf files.

#### 3. Get cline SD estimates from subset of data
* see [`3a.run-bgchm-parallel.sh`](scripts/3a.run-bgchm-parallel.sh) and [`3b.run-bgchm-parallel.R`](scripts/3b.run-bgchm-parallel.R)
* The first script is the wrapper, and randomly selects 1000 SNPs from the input set. The second script is running the code for bgchm.
* Input data need to be generated from step 2, or downloaded from supplemental data repo.

#### 4. Run bgchm
* see [`4a.ARRAY-run-bgchm-parallel-sierras.sh`](scripts/4a.ARRAY-run-bgchm-parallel-sierras.sh) and [`4b.run-bgchm-parallel-sierras.R`](scripts/4b.run-bgchm-parallel-sierras.R)
* Like before, the first script is a wrapper, and the second actually runs bgchm.
The wrapper script by default cuts up data into batches of 5000 SNPs and submits jobs. It must be run as an array, where the array is the number of jobs you need to submit, contingent on the batch size in SNPs. As an example, if you have 25,000 SNPs and have set batch_size to 10,000, you would need to submit three jobs. This data set has 49019 SNPs, with batches of 5000, so the array is set to ```--array=0-9%10```. THE ARRAY MUST START WITH 0.
* The second script has cline SD estimates from step 2 hard-coded in. It is a revamped version of the script available on the bgchm developer repo [here](https://github.com/zgompert/bgc-hm). It writes out each batch to an R object file, which will later be combined into one file for analysis.

#### 5. Combine bgchm subsets and write results for cline center and slope estimates
* see [`5.combine-bgchm-parallel.R`](scripts/5.combine-bgchm-parallel.R)

#### 6. Plotting cline results, permutation test, and a few other things
* see [`6.post-cluster-plots-bgchm.R`](scripts/6.post-cluster-plots-bgchm.R)
* There's a lot going on in here -- some of the plots produced are not present in the manuscript, but could be useful to someone
