library(tidyverse)
library(ggh4x)
library(ggtern)
library(ggridges)
library(patchwork)
setwd("~/7.ancestryHMM/data")

# POSTERIOR ESTIMATES -> making ternary plot
################################################################################
# Process the posterior files
# Need to obtain these from the supplemental data repo
posteriorlist <- list.files(pattern = "\\.posterior$")

for(i in 1:length(posteriorlist)){
  samplename <- posteriorlist[i] %>% str_remove(., "BWS_") %>% str_remove(., ".posterior")
  
  infile <- read_delim(posteriorlist[i], delim = "\t") %>%
    mutate(gt = case_when(
      `2,0` >= `1,1` & `2,0` >= `0,2` & `2,0` >= 0.8 ~ "0",
      `1,1` >= `2,0` & `1,1` >= `0,2` & `1,1` >= 0.8 ~ "1",
      `0,2` >= `2,0` & `0,2` >= `1,1` & `0,2` >= 0.8 ~ "2",
      TRUE ~ NA
    )) %>%
    select(c(chrom, position, gt)) %>%
    rename(!!samplename := gt)
  
  if(i == 1){
    posteriordf <- infile
  }else{
    posteriordf <- cbind(posteriordf, infile[,3])
  }
}

# Calculate average dav, het, and new proportion for each sample
post_proportions <- posteriordf %>%
  select(-c(chrom, position)) %>%
  pivot_longer(cols = starts_with("1"), names_to = "ID", values_to = "gt") %>%
  group_by(ID) %>%
  reframe(
    total_gt = sum(!is.na(gt)),
    new = mean(gt == 0, na.rm = TRUE),
    het = mean(gt == 1, na.rm = TRUE),
    dav = mean(gt == 2, na.rm = TRUE)
  ) %>%
  mutate(pop = str_extract(ID,  ".*(?=-)"))

write.csv(post_proportions, file = "result-posterior-proportions-sierras.csv", row.names = F)
post_proportions <- read_csv("result-posterior-proportions-sierras.csv")

# Plot these values in a ternary plot
# New plot with VL vs. GL
# This is new part of fig. 1
newdata <- post_proportions %>%
  mutate(Population = case_when(pop == "129" ~ "VL",
                            TRUE ~ "GL")) %>%
  arrange(Population == "VL")

ggtern(data = newdata, aes(x = dav, y = het, z = new, col = Population)) +
  geom_point(alpha = 0.90, pch = 16) +
  theme_custom(col.T = "black",
               tern.panel.background = "white",
               col.grid.minor = "grey") +
  theme_showarrows()
ggsave(filename = "PLOT-ternary-withpops-sierras-NEWORDER.svg",
       device = "svg", height = 3, units = "in", dpi= 600)

################################################################################


# Investigating excess het region
################################################################################
# Compare ancestry at Hue locus (not viterbi) compared to genome-wide proportion
# Reminder: 0 is new homozygous, 1 is het, 2 is dav homozygous
# Read in sample
# Need 
setwd("outfiles_0.75-sierras/")
posteriorlist <- list.files(pattern = "\\.posterior$")
comparisonmat <- data.frame()

for(i in 1:length(posteriorlist)){
  samplename <- posteriorlist[i] %>% str_remove(., "BWS_") %>% str_remove(., ".posterior")
  
  infile <- read_delim(posteriorlist[i], delim = "\t") %>%
    mutate(gt = case_when(`2,0` >= 0.8 ~ "0",
                          `1,1` >= 0.8 ~ "1",
                          `0,2` >= 0.8 ~ "2",
                          TRUE ~ NA)) %>%
    select(c(chrom, position, gt)) %>%
    setNames(c("chr", "ps", "gt")) %>%
    mutate(chr = str_remove(chr, "scaffold_")) %>%
    rename(!!samplename := gt)
  
  # Filter to hue locus
  huelocus <- infile %>%
    filter(chr == "1086" & ps >= 34150604 & ps <= 46470553)
  
  # Filter to rest of genome
  restofgenome <- infile %>%
    anti_join(., huelocus)
  
  # Calculate average dav, het, and new proportion -- Hue locus
  proportions_hue <- huelocus %>%
    select(-chr, -ps) %>%
    pivot_longer(cols = starts_with("1"), names_to = "ID", values_to = "gt") %>%
    group_by(ID) %>%
    reframe(
      total_gt = sum(!is.na(gt)),
      new_hue = mean(gt == 0, na.rm = TRUE),
      het_hue = mean(gt == 1, na.rm = TRUE),
      dav_hue = mean(gt == 2, na.rm = TRUE)
    ) %>%
    mutate(pop = str_extract(ID,  ".*(?=-)"))
  
  # Calculate average dav, het, and new proportion -- Rest of genome
  proportions_genome <- restofgenome %>%
    select(-chr, -ps) %>%
    pivot_longer(cols = starts_with("1"), names_to = "ID", values_to = "gt") %>%
    group_by(ID) %>%
    reframe(
      total_gt = sum(!is.na(gt)),
      new = mean(gt == 0, na.rm = TRUE),
      het = mean(gt == 1, na.rm = TRUE),
      dav = mean(gt == 2, na.rm = TRUE)
    ) %>%
    mutate(pop = str_extract(ID,  ".*(?=-)"))
  
  # Combine and add to cumulative data.frame
  t <- full_join(proportions_hue, proportions_genome, by = c("ID", "pop")) %>%
    select(ID, pop, new, het, dav, new_hue, het_hue, dav_hue)
  
  comparisonmat <- rbind(comparisonmat, t)
  
}

# Write combined matrix to disk
write_csv(comparisonmat, file = "result-hue-vs-genome-comparisonmat.csv")
comparisonmat <- read_csv("result-hue-vs-genome-comparisonmat.csv")

# Add admixture proportion and use that as X axis instead:
admixture <- read_csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  rename(ID = "newID")


# Make data.frame for density plots
hybrid_densityplot <- comparisonmat %>%
  mutate(dav_combhue = dav_hue + (het_hue/2),
         dav_combgenome = dav + (het/2)) %>%
  left_join(., admixture, by = "ID") %>%
  filter(!is.na(K1)) %>%
  mutate(densityID = case_when(dav_combhue >= 0.25 & dav_combhue <= 0.75 ~ "hybrid",
                               dav_combhue < 0.25 ~ "newberryi",
                               dav_combhue > 0.25 ~ "davidsonii"))

# Plot via Density map and scatterplot
# Total ancestry -- homozygous dav + 1/2 het
ggplot(hybrid_densityplot, aes(x = dav_combgenome, y = dav_combhue)) +
  geom_density_2d_filled(contour_var = "density", alpha = 0.9) +
  geom_point(alpha = 0.95, size = 1, aes(color = rgb), pch = 16) +
  scale_color_identity() +
  # geom_point(data = outliers_newberryi, inherit.aes = T, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white") +
  coord_fixed() +
  scale_x_continuous(limits= c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits= c(0,1), expand = c(0,0)) +
  labs(x = "Genome-wide Davidsonii ancestry (ancHMM)",
       y = "Hue Locus Davidsonii ancestry (ancHMM)",
       title = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), "cm"))
ggsave("PLOT-densityplot-genomewide-vs-hue-dav-ancestry.png",
       width = 3.5, height = 3.5, dpi = 600)

# Plot with admixture proportion
ggplot(hybrid_densityplot, aes(x = K2, y = dav_combhue)) +
  geom_density_2d_filled(contour_var = "density", alpha = 0.9) +
  geom_point(alpha = 0.95, size = 1, aes(color = rgb), pch = 16) +
  scale_color_identity() +
  #geom_point(data = outliers_newberryi, inherit.aes = T, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white") +
  coord_fixed() +
  scale_x_continuous(limits= c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits= c(0,1), expand = c(0,0)) +
  labs(x = "Genome-wide Davidsonii ancestry (Admixture)",
       y = "Hue Locus Davidsonii ancestry (ancHMM)",
       title = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), "cm"))
ggsave("PLOT-densityplot-genomewide-vs-hue-dav-ancestry-AMIXTURE.pdf",
       width = 3.5, height = 3.5, dpi = 900)


################################################################################

