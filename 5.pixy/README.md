# Genomic differentiation and divergence

Additionally uses data from:
[4.gemma](~/4.gemma)
Reference genome annotations (get from supplemental data repo)

#### 1. Generate all-sites .vcf file
* see [`1.generate-parent-allsites-file.sh`](scripts/1.generate-parent-allsites-file.sh)
* Uses the P. davidsonii genome, the mapped reads, and two ID files containing the names of parent individuals, [`davidsonii-fst-pops.txt`](data/davidsonii-fst-pops.txt) and [`newberryi-fst-pops.txt`](data/newberryi-fst-pops.txt)
* Uses bcftools, and vcftools

#### 2. Run pixy genome-wide and on f3'5'h region
* see [`2a.run_pixy_nomask-100kb.sh`](scripts/2a.run_pixy_nomask-100kb.sh) and [`2b.run_pixy_nomask_f3p5ph.sh`](scripts/2b.run_pixy_nomask_f3p5ph.sh)
* Genome-wide runs on entire genome, in 100kb windows. The f3'5'h region runs just on the region encompassing the floral hue association peak, in 100bp windows.
* Uses a file indicating all parent taxa, [`pixy_parentpops_sierras.txt`](data/pixy_parentpops_sierras.txt)

#### 3. Plotting results
* see [`3.plot_pixy_sierras_parents.R`](scripts/3.plot_pixy_sierras_parents.R)

