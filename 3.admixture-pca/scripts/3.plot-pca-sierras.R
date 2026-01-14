setwd("data/")
library(tidyverse)
library(ggpubr)

# File setup and organization
################################################################################
# Read in the subpopulations phenotype file
# This will help with plotting
subpopulations <- read_csv("~/1.phenotype_curation/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  filter(region == "sierras")

# Read in eigenvectors and eigenvalues
# Then, bind to the subpopulations data.frame
pca <- read_delim("PCA.eigenvec", col_names = T, delim = "\t") %>%
  select(-`#FID`) %>%
  rename(newID = IID) %>%
  mutate(newID = str_remove(newID, "BWS_"),
         pop = str_extract(newID, "^.{3}")) %>%
  left_join(., subpopulations, by = c('newID')) %>% 
  filter(!newID == "133-013") %>%
  mutate(subpopulation = case_when(subpopulation == "davidsonii_virginialakes" ~ "dav (VL)",
                                    subpopulation == "davidsonii_gemlakes" ~ "dav (GL)",
                                    subpopulation == "newberryi_virginialakes" ~ "new (VL)",
                                    subpopulation == "newberryi_gemlakes" ~ "new (GL)",
                                    subpopulation == "hybrids_virginialakes" ~ "hybrid (VL)",
                                    subpopulation == "hybrids_gemlakes" ~ "hybrid (GL)"))

# Write these PC results to a new file (helps with downstream analyses)
# write.csv(pca, file = "PCA_values_phenotypes_admixture-sierras.csv", row.names = F)

# Check which sample is missing in subpopulations (admixture)
# Turns out its 133-013 -- filter that out before plotting
#anti_join(pca, subpopulations, by = "newID")

# Reset levels so plotting order makes more sense
mylevels <- c("dav (VL)", "hybrid (VL)", "new (VL)",
              "dav (GL)", "hybrid (GL)", "new (GL)")
pca$subpopulation <- factor(pca$subpopulation, levels=mylevels)



# Read Eigenvalues
eigenval <- scan("PCA.eigenval")

# Convert eigenvalues to % variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
# write.csv(pve, file = "PVE-eigenvalues.csv", row.names = F)
################################################################################


# Basic plot of the PC space
################################################################################
mycols <- c("dav (VL)" = "mediumorchid3",
              "hybrid (VL)" = "springgreen4",
              "new (VL)" = "palevioletred2",
              "dav (GL)" = "mediumpurple1",
              "hybrid (GL)" = "darkseagreen3",
              "new (GL)" = "firebrick1")

# PCs 1 and 2
ggplot(data = pca %>% arrange(subpopulation), aes(x = PC1, y = PC2, col = subpopulation)) +
  geom_point(aes(col = subpopulation)) +
  scale_color_manual(values = mycols) +
  labs(x = paste0("PC1 (", round(pve[1,2], digits = 2), "%)", sep = ""),
       y = paste0("PC2 (", round(pve[2,2], digits = 2), "%)", sep = ""),
       color =NULL) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.5),
        legend.background = element_rect(fill = NA, color = "black"),
        legend.key=element_blank())
ggsave(filename = "PLOT-genomic-PCA-sierras.svg", dpi = 400, width = 4.75, height = 4.5, units = "in")

# PCs 2 and 3
ggplot(data = pca %>% arrange(subpopulation), aes(x = PC2, y = PC3, col = subpopulation)) +
  geom_point(aes(col = subpopulation)) +
  scale_color_manual(values = mycols) +
  labs(x = paste0("PC2 (", round(pve[2,2], digits = 2), "%)", sep = ""),
       y = paste0("PC3 (", round(pve[3,2], digits = 2), "%)", sep = ""),
       title = "Genomic PCA - Sierras") +
  theme_bw()
# ggsave(filename = "PLOT-genomic-PCA-sierras-PC2+PC3.png", dpi = 400, width = 4.7, height = 3.5, units = "in")

# PCs 1 and 3
ggplot(data = pca %>% arrange(subpopulation), aes(x = PC1, y = PC3, col = subpopulation)) +
  geom_point(aes(col = subpopulation)) +
  scale_color_manual(values = mycols) +
  labs(x = paste0("PC1 (", round(pve[1,2], digits = 2), "%)", sep = ""),
       y = paste0("PC3 (", round(pve[3,2], digits = 2), "%)", sep = ""),
       title = "Genomic PCA - Sierras") +
  theme_bw()
# ggsave(filename = "PLOT-genomic-PCA-sierras-PC1+PC3.png", dpi = 400)
################################################################################
