# This script runs visual models on flower vs. leaf contrast for honeybees and hummingbirds

library(tidyverse)
library(ggbeeswarm)
library(zoo)


setwd("~9.visualmodels/")

# DATA I/O
################################################################################
# Import lists of parents and hybrids for which we have genetic data
p1list <- read_lines("~/10.bgchm/data/p0file_newberryi-sierras.txt") %>%
  str_remove("BWS_")

p2list <- read_lines("~/10.bgchm/data/p1file_davidsonii-sierras.txt") %>%
  str_remove("BWS_")

hybridlist <- read_lines("~/10.bgchm/data/samplelist-bgchm-sierras.txt") %>%
  str_remove("BWS_") %>%
  discard(~ .x %in% c(p1list, p2list))


# Read in matrix with information about admixture proportion
admixmat <- read_csv("~1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  filter(!is.na(K1)) %>%
  rename(sample = "newID")

# Read in matrix with F3'5'H genotypes
f3p5phmat <- read_csv("~/6.pigments/data/pcr-f3p5ph-genotypes.csv") %>%
  select(-class)

# Read in the spectral data previously generated
specdata <- read_csv("data/final_averaged_specdata_300-700-fixneg-only.csv") %>%
  rename("wavelength" = rowname)

################################################################################

# PLOTTING VISUAL MODELS
################################################################################
# Generated from step 1
setwd("data/")
infiles <- list.files(pattern = "result-vismodel*")

# PLOTTING CONTRAST IN HYBRIDS VS PARENTS BOXPLOT
for(i in infiles){
  
  # Get the basename for naming plot
  namepattern = str_remove(i, "result-vismodel-contrast-") %>% str_remove(., ".csv")
  
  # Read in matrix to plot, make new naming class
  plotmat <- read_csv(i) %>%
    mutate(species = case_when(sample %in% p1list ~ "newberryi",
                               sample %in% p2list ~ "davidsonii",
                               sample %in% hybridlist ~ "hybrid")) %>%
    left_join(f3p5phmat, by = "sample") %>%
    pivot_longer(cols = c(bg_dSbird, bg_dSbee, bg_dLbird, bg_dLbee),
                 names_to = c(".value", "visitor"),
                 names_pattern = "bg_(d[SL])(bird|bee)") %>%
    rename(chromatic_contrast = dS, achromatic_contrast = dL)
  
  # Get background type
  backgroundtype = unique(plotmat$background)
  
  # Plot contrast in hybrids vs. each parent for bees and birds
  ggplot(plotmat, aes(x = species, y = chromatic_contrast, col = species)) +
    geom_violin() +
    geom_beeswarm(size = 0.01) +
    facet_wrap(~visitor) +
    scale_color_manual(values = c("#6245ba", "goldenrod", "#c9368c")) +
    labs(x = "Species",
         y = "Chromatic contrast",
         # title = paste0("Background: ", backgroundtype, sep = "")) +
         title = NULL) +
    # scale_y_continuous(limits = c(0,16), n.breaks = 6) +
    theme_bw() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 7),
          title = element_text(size = 10))
  ggsave(paste0("PLOT-chromatic-contrast-allsamples-violin-", namepattern, ".pdf", sep = ""),
         height = 2.5, width = 3.5, units = "in", dpi = 900)
}


################################################################################