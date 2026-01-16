# This script runs visual models on flower vs. leaf contrast for honeybees and hummingbirds

library(pavo)
library(tidyverse)


setwd("~/9.visualmodels/")

# DATA I/O
################################################################################
# Import floral spectra, limiting from 300-700 nm wavelength
# Filter to move negatives to 0
# Obtain the spec data from the online repo
rawspec <- getspec(where = "specdata-flowers/",
                    ext = "txt", lim = c(300,700))
rawspec_filtered <- as.rspec(procspec(rawspec, fixneg = "zero"))
# Filter data. Original models only used fixneg
# Alternative filters for min and max (span 0-1 can be tried too?)



# Import leaf spectra, limiting from 300-800 nm wavelength
# Filter to move negatives to 0, and only up to 700 nm
# Obtain the spec data from the online repo
leafdata <- getspec(where = "specdata-leaf/", ext = "txt", lim = c(300, 700))
leafdata_filtered <- as.rspec(procspec(leafdata, fixneg = "zero"))
write_csv(leafdata_filtered, "specdata_leafdata_300-700-fixneg-only.csv")


# Import lists of parents and hybrids for which we have genetic data
p1list <- read_lines("~/10.bgchm/data/p0file_newberryi-sierras.txt") %>%
  str_remove("BWS_")

p2list <- read_lines("~/10.bgchm/data/p1file_davidsonii-sierras.txt") %>%
  str_remove("BWS_")

hybridlist <- read_lines("~/10.bgchm/data/samplelist-bgchm-sierras.txt") %>%
  str_remove("BWS_") %>%
  discard(~ .x %in% c(p1list, p2list))

################################################################################

# GENERATING AVERAGED FLORAL SPECTRA AND FILTERING TO SEQUENCED SAMPLES
################################################################################
# Find A and B names
a_names <- names(rawspec_filtered)[grep("A", names(rawspec_filtered))]
b_names <- names(rawspec_filtered)[grep("B", names(rawspec_filtered))]

# Make new df with only flowers to be averaged (also sort and add wl column to front)
avg_flowers <- rawspec_filtered %>%
  select(contains(sort(c(a_names, b_names)))) %>%
  mutate(wl = rawspec_filtered$wl, .before = names(.)[1])

# Verify that the lengths of the dfs are correct
length(a_names) + length(b_names) +1 == ncol(avg_flowers)

# Find flowers that do not need averaged (and do not keep the wl column)
noavg_flowers <- rawspec_filtered %>%
  select(-contains(names(avg_flowers))) %>%
  select(-contains("wl"))

# Aggregate spectra -- find average spectra for samples that need it
# Combine averaged and unaveraged spectra into new data set and clean up names
renamed_specdata <- cbind(aggspec(avg_flowers, by = 2, FUN = mean), noavg_flowers) %>%
  rename_all(~gsub("A", "", .)) %>%
  select(-wl) %>%
  `rownames<-`(rawspec_filtered$wl)
colnames(renamed_specdata)

# Rename samples so there are always three digits after the dash
# This finds a dash followed by one or more digits at end of string (-(\\d+)$)
# Then replaces it with - and same digits, left-padded with 0 to make it 3 digits
colnames(renamed_specdata) <- colnames(renamed_specdata) %>%
  str_replace("-(\\d+)$", ~ paste0("-", str_pad(str_match(.x, "-(\\d+)$")[,2], 3, pad = "0")))
colnames(renamed_specdata)


# Finally, filter out all samples not in parents or hybrid lists
finalspec_filtered <- renamed_specdata %>%
  select(any_of(c(p1list, p2list, hybridlist)))

write_csv(finalspec_filtered %>% rownames_to_column(var = "rowname"),
          file = "data/final_averaged_specdata_300-700-fixneg-only.csv")

################################################################################


# Questions/notes about I/O
# 1. I don't have to run separate models -- the result is labeled, so I can filter after
# 2. I assume it is OK to average spectra in the same way as for calculating hue?


# RUN VISUAL MODELS -- CHROMATIC AND ACHROMATIC CONTRASTS
################################################################################

# Read in prepared data
# The background spectra should have first col = "wl", second col = data
# For samples, no wl column
finalspec_filtered <- read_csv("data/final_averaged_specdata_300-700-fixneg-only.csv") %>%
  column_to_rownames(., var = "rowname")
leafdata_filtered <- read_csv("data/specdata_leafdata_300-700-fixneg-only.csv")


finalmat <- data.frame()
# Run visual modes -- background is sample 1 and input flower is sample 2
# Then, calculate color distances (contrasts) based on visual model
# Applies receptor-noise model of Vorobyev et al. 1998 to calculate color distances
# Noise based on relative photoreceptor densities
# Weber fraction (AKA receptor noise) -- this is default value

# Set parameters... -- can set multiple and iterate through if desired
birdvis_params <- c("avg.v")
illum_params <- c("ideal")
birdrecept_params <- list(c(1,1,1,2))
background_options <- list(leaf = leafdata_filtered)
# background_options <- list(leaf = leafdata_filtered, sand = sanddata_filtered, tree = treedata_filtered)


# MY ORIGINAL -- USES LEAF BG
final_output <- list()
# Loop through all combinations of parameter sets
for (birdvis_param in birdvis_params) {
  for (illum_param in illum_params) {
    for (birdrecept_param in birdrecept_params) {
      for (bgname in names(background_options)) {
        
        backgroundspectra <- background_options[[bgname]]
        finalmat <- data.frame()
        
        for (i in 1:ncol(finalspec_filtered)) {
          
          # Bird visual model
          birdmodel <- vismodel(
            rspecdata = cbind(backgroundspectra, finalspec_filtered[i]),
            visual = birdvis_param,
            illum = illum_param,
            qcatch = "Qi",
            achromatic = "l",
            vonkries = FALSE,
            relative = FALSE
          )
          
          birdcoldist <- coldist(
            modeldata = birdmodel,
            achromatic = T,
            n = birdrecept_param,
            weber = 0.1
          ) %>%
            setNames(c("background", "sample", "bg_dSbird", "bg_dLbird"))
          
          # Bee visual model
          beemodel <- vismodel(
            rspecdata = cbind(backgroundspectra, finalspec_filtered[i]),
            visual = "apis",
            # visual = pseudomasaris_sens,
            illum = illum_param,
            qcatch = "Qi",
            achromatic = "l",
            vonkries = FALSE,
            relative = FALSE
          )
          
          beecoldist <- coldist(
            modeldata = beemodel,
            achromatic = T,
            n = c(1, 2, 1),
            weber = 0.1
          ) %>%
            setNames(c("background", "sample", "bg_dSbee", "bg_dLbee"))
          
          # Combine
          finalmat <- rbind(finalmat, cbind(birdcoldist, beecoldist[, 3:4]))
        }
        
        # Build file name
        param_string <- paste0(
          birdvis_param, "-",
          illum_param, "-",
          paste(birdrecept_param, collapse = ""), "-",
          bgname)
        
        # Write CSV
        write_csv(finalmat, file = paste0("result-vismodel-contrast-fixneg-only-", param_string, ".csv"))
      }
    }
  }
}
################################################################################

# STATISTICS FOR JNDs
################################################################################

# Find average color distance (JND) per species per visitor
# Read in JNDs
finalmat <- read.csv("data/result-vismodel-contrast-fixneg-only-avg.v-ideal-1112-leaf.csv") %>%
  mutate(species = case_when(sample %in% p1list ~ "newberryi",
                             sample %in% p2list ~ "davidsonii",
                             sample %in% hybridlist ~ "hybrid"))

# Summarize
summarymat <- finalmat %>%
  group_by(species) %>%
  summarize(bird_chromatic_contrast = mean(bg_dSbird, na.rm = TRUE),
            bee_chromatic_contrast = mean(bg_dSbee, na.rm = TRUE))
write.csv(summarymat, file = "data/RESULT-meanvalues-per-species-avg.v-ideal-1112-leaf.csv")

# Run ANOVA and Tukey Test
anova_bird <- aov(bg_dSbird ~ species, data = finalmat)
tukey_bird <- TukeyHSD(anova_bird)
write.csv(tukey_bird$species, "data/RESULT-tukeytest-BIRD.v-ideal-1112-leaf.csv")

anova_bee <- aov(bg_dSbee ~ species, data = finalmat)
tukey_bee <- TukeyHSD(anova_bee)
write.csv(tukey_bee$species, "data/RESULT-tukeytest-BEE.v-ideal-1112-leaf.csv")
################################################################################


# RUN VISUAL MODELS -- COLOR SPACE
################################################################################
setwd("~/Desktop/sierras_hybridzone/14.visualmodels/")
finalspec_filtered_colorspace <- as.rspec(read_csv("data/final_averaged_specdata_300-700-fixneg-only.csv") %>%
  rename(wl = "rowname"))


# Run color space models
# Bee models
beemodel <- vismodel(rspecdata = finalspec_filtered_colorspace, visual = "apis",
                     illum = "ideal", qcatch = "Ei", achromatic = "l",
                     relative = F, vonkries = T,  bkg = "green")

beecolspace <- colspace(beemodel, space = "hexagon")

# Bird models
# birdmodel <- vismodel(rspecdata = finalspec_filtered_colorspace, visual = "avg.v",
#                      illum = "ideal", qcatch = "Ei", achromatic = "l",
#                      relative = T, vonkries = T,  bkg = "green")
# 
# birdcolspace <- colspace(birdmodel, space = "tcs")

################################################################################


admixmat <- read_csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  filter(!is.na(K1)) %>%
  rename(sample = "newID")

# BEE PLOT
# Add info for coloring
beeplot <- beecolspace %>%
  rownames_to_column(var = "sample") %>%
  left_join(., admixmat, by = "sample")

# Define hexagon outline
hex_outline <- tibble(angle = seq(0, 2 * pi, length.out = 7) + pi / 6,
                      x = cos(angle), y = sin(angle))

# Add points for drawing sectors
hex_midpoints <- tibble(x = (hex_outline$x[1:6] + hex_outline$x[2:7]) / 2,
                        y = (hex_outline$y[1:6] + hex_outline$y[2:7]) / 2,
                        group = 1:6)

# Add labels inside sectors
hex_labels <- tibble(x = cos(seq(0, 2 * pi, length.out = 7) + pi / 6)[1:6] * 0.7,
                     y = sin(seq(0, 2 * pi, length.out = 7) + pi / 6)[1:6] * 0.8,
                     label = c("Blue-Green", "Blue", "UV-Blue", "UV", "UV-Green", "Green"))
                     
  

# Plot
ggplot(beeplot, aes(x = x, y = y, color = rgb)) +
  geom_segment(data = hex_midpoints, aes(x = x, y = y, xend = 0, yend = 0, group = group), 
               color = "grey", linewidth = 0.5) +
  geom_path(data = hex_outline, aes(x, y), inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_text(data = hex_labels, aes(x, y, label = label),
            inherit.aes = FALSE, color = "black", size = 1.5) +
  geom_point(x = 0, y = 0, shape = 15, size = 2, color = "grey", inherit.aes = FALSE) +
  geom_point(size = 1, pch = 16, alpha = 0.8, stroke = 0) + # These are the data
  scale_color_identity() +
  labs(x = NULL,
       y = NULL,
       title = NULL) +
  coord_fixed() +
  facet_wrap(~class) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        # panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", linewidth = 0))
ggsave("plot-colorspace-beemodel-ideal-fixneg-only.pdf",
       height = 3, width = 5, units = "in", dpi = 900)
