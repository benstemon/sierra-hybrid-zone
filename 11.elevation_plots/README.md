# Fix elevation data and make some plots

Additionally uses data from: [1.phenotype_curation](~/1.phenotype_curation/)

#### 1. Estimate elevation for each sample
* see [`1.fix_elevation_elevatr.R`](scripts/1.fix_elevation_elevatr.R)
* The original data were clearly off, from the gps. So, this was necessary to get better estimates of elevation for each sample. Uses elevatr package.

#### 2. GLM for elevation and admixture proportion with new elevation estimates
* see [`2.elevation-investigation.R`](scripts/2.elevation-investigation.R)