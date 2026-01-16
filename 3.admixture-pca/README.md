## Admixture and Genomic PCA

Note: these scripts make use of data from [`1.phenotype_curation`](~/1.phenotype_curation), and [`11.elevation_plots`]((../11.elevation_plots/))

#### 1. Convert vcf to admixture data format and run admixture
* see [`1.vcf2admixture.sh`](scripts/1.vcf2admixture.sh)

#### 2. Run PCA in plink2
* see [`2.run_genomic_pca-sierras.sh`](scripts/2.run_genomic_pca-sierras.sh)

#### 3. Plotting genomic PCA
* see [`3.plot-pca-sierras.R`](scripts/3.plot-pca-sierras.R)