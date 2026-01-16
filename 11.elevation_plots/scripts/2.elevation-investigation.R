library(tidyverse)
library(elevatr)
library(ggpmisc)
library(rcompanion)
library(pROC)
setwd("~/11.elevation_plots/")


# Read in new elevation data
elevation_data <- read_delim("data/sierras-pops-NEWelevation.txt", delim = "\t") %>%
  rename(newID = accession)

# Read in lat/lon data
latlon_data <- read_delim("data/sierras-pops-elevation.txt", delim = "\t") %>%
  select(-elevation) %>%
  rename(newID = accession)

# Read in phenotype data
phenodata <- read_csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  filter(region == "sierras")

# Combine together, add additional designations. Add factors to plot in the right order
alldata <- left_join(phenodata, elevation_data, by = "newID") %>%
  left_join(., latlon_data, by = "newID") %>%
  mutate(location = case_when(pop %in% c("131", "133", "135") & lat > 37.400 ~ "Gem Lakes North",
                              pop %in% c("131", "132") & lat < 37.400 ~ "Gem Lakes South",
                              TRUE ~ "Virginia Lakes")) %>%
  transform(location = factor(location, levels = c("Virginia Lakes", "Gem Lakes North", "Gem Lakes South")))



# USING GEM LAKES AS A SINGLE POPULATION
################################################################################
otherdata <- alldata %>%
  mutate(location2 = case_when(location %in% c("Gem Lakes North", "Gem Lakes South") ~ "Gem Lakes",
         TRUE ~ location)) %>%
  transform(location2 = factor(location2, levels = c("Virginia Lakes", "Gem Lakes")))

# Use a binomial GLM
gl_data <- otherdata %>% filter(location2 == "Gem Lakes")
vl_data <- otherdata %>% filter(location2 == "Virginia Lakes")

glm_GL <- glm(K1 ~ newelevation, data = gl_data, family = binomial(link = "logit"))
glm_VL <- glm(K1 ~ newelevation, data = vl_data, family = binomial(link = "logit"))


# Testing fit and significant interactions: nagelkerkes R2, anova (chi-squared)
r2_GL <- signif(nagelkerke(glm_GL)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anova_GL <- signif(anova(glm_GL, test = "LRT")$`Pr(>Chi)`[2], digits = 3)

r2_VL <- signif(nagelkerke(glm_VL)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anova_VL <- signif(anova(glm_VL, test = "LRT")$`Pr(>Chi)`[2], digits = 3)


# Make text objects of these so they can be added to the plot
plot_text <- data.frame(r2 = c(paste0("R^2 = ",r2_VL), paste0("R^2 = ",r2_GL)),
                        chi = c(paste0("p = ",anova_VL), paste0("p = ",anova_GL)),
                        K1 = c(0.35, 0.35), newelevation = c(11250, 10250),
                        location2 = c(levels(otherdata$location2))) %>%
  transform(location2 = factor(location2, levels = c("Virginia Lakes", "Gem Lakes")))
plot_text$r2[[1]] <- list(bquote(italic(R)^2 == .(r2_VL)))
plot_text$r2[[2]] <- list(bquote(italic(R)^2 == .(r2_GL)))

ggplot(otherdata %>% filter(!is.na(K1)), aes(x = newelevation, y = K1)) +
  geom_point() +
  facet_wrap(~location2, nrow = 1, ncol = 2, strip.position = "top") +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(y = "Admixture Proportion (newberryi)",
       x = "Elevation (ft)") +
  theme_bw() +
  geom_text(data = plot_text, mapping = aes(label = r2), parse = T) +
  geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.1))+
  theme(strip.text.y = element_text(size = 10)) 

ggsave(filename = "PLOT-elevation-by-admixture-proportion-sierras-GEMLAKESMERGED.pdf", dpi = 900,
       height = 3, width = 5, units = "in")
################################################################################


