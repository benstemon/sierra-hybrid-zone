# This is a collection of plots based on phenotypic data
setwd("~/1.phenotype_curation/")
library(tidyverse)

# correlations
library(corrplot)

# random forest, imputation
library(randomForest)
library(mice)

# ggridges and plotting
library(ggridges)
library(cowplot)

# regressions
library(ggpmisc)
library(rcompanion)

# additional PCA plotting options
library(ggbiplot)
library(geomtextpath)
library(scales)

# Data input and tidying
################################################################################
# Read in summarized phenotype information
sumdf <- read_delim("data/SUMMARIZED_sierras_phenotype_data.txt", delim = "\t")

# Read in table that includes rgb -- for some reason not in the summarized information
rgbdf <- read_delim("data/ALLSAMPLES_sierras_color_chroma_hue_averageflowers.csv", delim = ",") %>%
 select(c("newID", "rgb"))

# Read in admixture proportion data (from admixture analysis directory)
samples <- read.table("~/3.admixture-pca/data/extracted.admixture_prep_LD-50-10-0.1.fam", header = F)[,1] %>%
 str_remove(., "BWS_")
qmat <- read.table("~/3.admixture-pca/data/extracted.admixture_prep_LD-50-10-0.1.2.Q")
colnames(qmat) <- paste("K", seq(1:ncol(qmat)), sep="")
qmat$newID <- samples

# Merge all three data frames
merged_df <- left_join(sumdf, qmat, by = "newID") %>%
 left_join(., rgbdf, by = "newID")

# Now, for some additional labeling, based on the admixture proportions
`%nin%` <- negate(`%in%`)
admix_included <- merged_df %>%
 mutate(region = "sierras") %>%
 mutate(subpopulation = case_when(pop == "129" & K1 > 0.999 ~ "newberryi_virginialakes",
                                  pop == "131" & K1 > 0.999 ~ "newberryi_gemlakes",
                                  pop == "129" & K2 > 0.999 ~ "davidsonii_virginialakes",
                                  pop == "135" & K2 > 0.999 ~ "davidsonii_gemlakes")) %>%
 mutate(subpopulation = case_when(pop == "129" & K1 > 0 & subpopulation %nin% c("newberryi_virginialakes", "davidsonii_virginialakes") ~ "hybrids_virginialakes",
                                  pop %in% c("131", "132", "133", "135") & K1 > 0 & subpopulation %nin% c("newberryi_gemlakes", "davidsonii_gemlakes") ~ "hybrids_gemlakes",
                                  TRUE ~ subpopulation)) %>%
 mutate(class = case_when(subpopulation %in% c("newberryi_virginialakes", "newberryi_gemlakes") ~ "newberryi",
                          subpopulation %in% c("davidsonii_virginialakes", "davidsonii_gemlakes") ~ "davidsonii",
                          subpopulation %in% c("hybrids_virginialakes", "hybrids_gemlakes") ~ "hybrid"))

nrow(admix_included %>% filter(class == "newberryi"))
nrow(admix_included %>% filter(class == "davidsonii"))
nrow(admix_included %>% filter(class == "hybrid"))



# Write this combined matrix to a new file
write.csv(admix_included,
         file = "data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv",
         row.names = F)


# Read in the curated phenotype data set
setwd("1.phenotype_curation/")
admix_included <- read.csv("data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  select(-total_brightness)

################################################################################


# 1. RF to ID traits that differentiate parent species
################################################################################
# Impute missing data
imputed_matrix <- complete(mice(admix_included[,3:28], method = "pmm"))

# Generate tibble for RF input -- trying to differentiate parent species via traits
rf_input <- cbind(imputed_matrix, admix_included$class) %>%
  rename(class = `admix_included$class`) %>%
  filter(class != "hybrid") %>%
  mutate(class = as.factor(class))

# Generate RF model
rf_model <- randomForest(data = rf_input, class ~ ., importance = TRUE, ntree = 5000)

# SUMMARIZATION AND PERFORMANCE
# # #
# GINI, ACCURACY -- VARIABLE IMPORTANCE
importance_rf <- as_tibble(importance(rf_model), rownames = "Variable") %>%
  mutate(Variable = str_replace_all(Variable, "_", " ") %>%
           str_remove("avg") %>%
           str_to_title()) %>%
  mutate(MeanDecreaseAccuracy = ifelse(MeanDecreaseAccuracy < 0, 0, MeanDecreaseAccuracy)) %>%
  mutate(RelativeImportance = MeanDecreaseAccuracy / sum(MeanDecreaseAccuracy))

write_csv(importance_rf %>% arrange(desc(MeanDecreaseAccuracy)), "results-accuracy-gini.csv")


# Gini decrease
a <- ggplot(importance_rf, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = NULL,
       y = "Mean Decrease in Gini") +
  theme_bw()
ggsave("PLOT-Gini-decrease-RF.png", a, height = 4, width = 4, dpi = 600)

# Accuracy decrease
b <- ggplot(importance_rf, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = NULL,
       y = "Mean Decrease in Accuracy") +
  theme_bw()
ggsave("PLOT-Accuracy-decrease-RF.png", b, height = 4, width = 4, dpi = 600)

# Relative importance
c <- ggplot(importance_rf, aes(x = reorder(Variable, RelativeImportance), y = RelativeImportance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = NULL,
       y = "Relative Importance") +
  theme_bw()
ggsave("PLOT-RelativeImportance-RF.pdf", c, height = 4, width = 4, dpi = 600)


# # #

# # #
# CONFUSION MATRIX, PERFORMANCE (this is based on OOB error)
rf_confusion <- rf_model$confusion
rf_confusion[1,1:2] <- rf_confusion[1,1:2]/47
rf_confusion[2,1:2] <- rf_confusion[2,1:2]/27
rf_confusion
write.csv(file = "confusion-matrix-oob.csv", as.data.frame(rf_confusion))
# # #

# Together, these results show me which variables are most important
# And also that the model has an acceptably low error rate

################################################################################

# 2. Correlation matrices
################################################################################
# List important traits based on variables with > 2.5% relative importance
# using only complete pairwise observations when there is missing data
# I removed crease to top height because it is highly correlated (> 95%) with height at mouth
#avg_crease_to_top_height
important_traits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                       "avg_mouth_height_width_ratio")
reduced_matrix <- admix_included %>%
  select(all_of(important_traits), K1) %>%
  rename_with(~ .x %>%
                str_replace_all("_", " ") %>%
                str_remove_all("avg") %>%
                str_to_title())
  
pheno_corrplot <- cor(reduced_matrix, use = "pairwise.complete.obs")
write.csv(pheno_corrplot, file = "phenotype-correlation-matrix-sierras-RFIMPORTANT-characters.csv")
png("PLOT-correlation-RFIMPORTANT-sierras.png", units = "in", width = 5, height = 5, res = 400)
corrplot(pheno_corrplot, method="color", type = "upper", diag = FALSE,
         addCoef.col = "white", number.cex = 0.45, tl.col = "black",
         tl.cex = 0.6, cl.cex = 0.6)
dev.off()


# Make a correlation matrix between only the "select" phenotypic traits that we like
# using only complete pairwise observations when there is missing data
reduced_matrix <- admix_included %>%
  select(c("avg_long_stamen_length","avg_short_stamen_length","avg_pistil_length",
           "avg_floral_tube_length","avg_tube_width_at_mouth","avg_tube_height_at_mouth",
           "avg_platform_angle","avg_longstamen_pistil_ratio","avg_longstamen_tube_ratio",
           "avg_floral_crease_ratio","avg_throat_inflation_ratio","avg_nectary_area",
           "chroma","modified_hue","K1")) %>%
  rename_with(~ .x %>%
                str_replace_all("_", " ") %>%
                str_remove_all("avg") %>%
                str_to_title())
pheno_corrplot <- cor(reduced_matrix, use = "pairwise.complete.obs")
write.csv(pheno_corrplot, file = "phenotype-correlation-matrix-sierras-SELECT-characters.csv")
png("PLOT-correlation-SELECTtraits-sierras.png", units = "in", width = 5, height = 5, res = 400)
corrplot(pheno_corrplot, method="color", type = "upper", diag = FALSE,
         addCoef.col = "white", number.cex = 0.45,tl.col = "black",
         tl.cex = 0.6, cl.cex = 0.6)
dev.off()


# Make a matrix for the full set of traits, too.
pheno_corrplot <- cor(admix_included %>% select(-c("newID", "pop", "K2", "rgb", "region", "subpopulation", "class")),
                      use = "pairwise.complete.obs")
write.csv(pheno_corrplot, file = "phenotype-correlation-matrix-sierras-ALL-characters.csv")
png("PLOT-correlation-ALLtraits-sierras.png", units = "in", width = 7.5, height = 7.5, res = 400)
corrplot(pheno_corrplot, method="color", type = "upper", diag = FALSE,
         addCoef.col = "white", number.cex = 0.45,tl.col = "black",
         tl.cex = 0.6, cl.cex = 0.6)
dev.off()

################################################################################


# 3.1 GGridges plots of most important traits for differentiation according to RF
################################################################################
# re-list the important traits (remember text about what these came from)
important_traits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                      "avg_mouth_height_width_ratio")

# Filter to only samples with genetic data
ggridgesplot <- admix_included %>%
  filter(!is.na(K1))

# Set up collection of plots
plot_list <- list()
for(i in important_traits){
  #Condition for plotting text
  if(length(plot_list) %in% c(0,4,8)){
    textsize = 6
  }else(textsize = 0)
  
  tmpmat <- ggridgesplot %>% 
    drop_na(all_of(i)) %>% 
    select(all_of(i), class) %>%
    rename_with(~ .x %>% 
                  str_replace_all("_", " ") %>% 
                  str_remove_all("avg") %>% 
                  str_to_title())
  
  trait <- colnames(tmpmat)[[1]]
  
  # Compute means for trait for each species
  means <- tmpmat %>%
    group_by(Class) %>%
    summarize(mean_val = mean(.data[[trait]], na.rm = TRUE))
    
  # Decide gradient orientation based on mean
  if(mean(means$mean_val[means$Class == "newberryi"]) < 
     mean(means$mean_val[means$Class == "davidsonii"])){
    low_col <- "palevioletred1"
    high_col <- "mediumorchid4"}
  else{
    low_col <- "mediumorchid4"
    high_col <- "palevioletred1"}
  
  a <- ggplot(data = tmpmat,
         aes(y = Class, x = .data[[trait]], fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.5, quantile_lines = TRUE, quantiles = 2) +
    scale_fill_gradient(low = low_col, high = high_col) +
    scale_x_continuous(expand = c(0,0,0,0)) +
    scale_y_discrete(expand = c(0,0,0.15,1.5)) + #1.1 for 6 classes
    coord_cartesian(clip = "off") +
    labs(x = colnames(tmpmat)[[1]]) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_text(size = textsize), axis.text.x = element_text(size = 4),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          axis.ticks.y = element_blank(), axis.title.x = element_text(size = 5))
  
  # append to list
  plot_list[[length(plot_list)+1]] <- a
}

# stitch together in groups of 4
stitched_plots1 <- plot_grid(
  plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
  plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
  plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
  nrow = 3,
  ncol = 4,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1.32,1,1,1) # 2.6 for 6 classes.
)
ggsave("results - ggridges/PLOT-ggridges-RFimportant-traits-CLASS.png",
       stitched_plots1, height = 3.75, width = 5, dpi = 600)

################################################################################

# 3.2 GGridges plots of all other traits for supplemental figure
################################################################################
# re-list the important traits (remember text about what these came from)
removetraits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                      "avg_mouth_height_width_ratio")

# List columns that are actually trait data, but not in the set of "important" traits
important_traits <- setdiff(colnames(admix_included[3:28]), removetraits)

# Filter to only samples with genetic data
ggridgesplot <- admix_included %>%
  filter(!is.na(K1))

# Set up collection of plots
plot_list <- list()
for(i in important_traits){
  #Condition for plotting text
  if(length(plot_list) %in% c(0,4,8,12,16,20,24)){
    textsize = 6
  }else(textsize = 0)
  
  tmpmat <- ggridgesplot %>% 
    drop_na(all_of(i)) %>% 
    select(all_of(i), class) %>%
    rename_with(~ .x %>% 
                  str_replace_all("_", " ") %>% 
                  str_remove_all("avg") %>% 
                  str_to_title())
  
  trait <- colnames(tmpmat)[[1]]
  
  # Compute means for trait for each species
  means <- tmpmat %>%
    group_by(Class) %>%
    summarize(mean_val = mean(.data[[trait]], na.rm = TRUE))
  
  # Decide gradient orientation based on mean
  if(mean(means$mean_val[means$Class == "newberryi"]) < 
     mean(means$mean_val[means$Class == "davidsonii"])){
    low_col <- "palevioletred1"
    high_col <- "mediumorchid4"}
  else{
    low_col <- "mediumorchid4"
    high_col <- "palevioletred1"}
  
  a <- ggplot(data = tmpmat,
              aes(y = Class, x = .data[[trait]], fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.5, quantile_lines = TRUE, quantiles = 2) +
    scale_fill_gradient(low = low_col, high = high_col) +
    scale_x_continuous(expand = c(0,0,0,0)) +
    scale_y_discrete(expand = c(0,0,0.15,1.5)) + #1.1 for 6 classes
    coord_cartesian(clip = "off") +
    labs(x = colnames(tmpmat)[[1]]) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_text(size = textsize), axis.text.x = element_text(size = 4),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          axis.ticks.y = element_blank(), axis.title.x = element_text(size = 5))
  
  # append to list
  plot_list[[length(plot_list)+1]] <- a
}

# stitch together in groups of 4
stitched_plots1 <- plot_grid(
  plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
  plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
  plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
  plot_list[[13]],
  nrow = 4,
  ncol = 4,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1.32,1,1,1) # 2.6 for 6 classes.
)
ggsave("results - ggridges/PLOT-ggridges-UNIMPORTANTTRAITS-CLASS.pdf",
       stitched_plots1, height = 5, width = 5, dpi = 900)

################################################################################

# 4. Plotting Admixture proportion vs. trait variation in hybrids
################################################################################
# filter to hybrid df
hybrid_df <- admix_included %>%
  filter(class == "hybrid")

# re-list the important traits (remember text about what these came from)
important_traits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                      "avg_mouth_height_width_ratio")

# Set up collection of plots for K1 vs. important traits
# LINEAR REGRESSION
plot_list_LINEAR <- list()
for(i in important_traits){
  a <- ggplot(hybrid_df, aes(x = .data[[i]], y = K1)) +
    geom_point(size = 0.5) +
    coord_cartesian(y = c(0,1.06)) +
    stat_poly_line() +
    stat_poly_eq(use_label(c("adj.R2","p", "n")), label.y = 0.98, size = 3) +
    scale_y_continuous(expand = c(0.1,0,0,0.1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_bw()
  
  # append to list
  plot_list_LINEAR[[length(plot_list_LINEAR)+1]] <- a
}

# stitch together in groups of 5
stitched_plots_LINEAR <- plot_grid(
  plot_list_LINEAR[[1]],plot_list_LINEAR[[2]],plot_list_LINEAR[[3]],plot_list_LINEAR[[4]],
  plot_list_LINEAR[[5]],plot_list_LINEAR[[6]],plot_list_LINEAR[[7]],plot_list_LINEAR[[8]],
  plot_list_LINEAR[[9]],plot_list_LINEAR[[10]],plot_list_LINEAR[[11]],plot_list_LINEAR[[12]],
  nrow = 3,
  ncol = 4,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1,1,1,1) # 2.6 for 6 classes.
)
ggsave("results - regressions/PLOT-regressions-RFimportant-traits-LINEAR.png",
       stitched_plots_LINEAR, height = 7.5, width = 11, dpi = 600)



# LOGISTIC REGRESSION
plot_list_LOGISTIC <- list()
for(i in important_traits){
  #Condition for plotting text
  if(length(plot_list_LOGISTIC) %in% c(0,5)){
    textsize = 10
  }else(textsize = 0)
  
  # Generate logistic regression GLM
  myglm <- glm(data = hybrid_df, reformulate(i, response = "K1"), family = binomial(link = "logit"))
  
  # Testing fit and significant interactions: nagelkerkes R2, anova (chi-squared)
  r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
  anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
  
  # Text for annotating the plot
  x_coord <- min(hybrid_df[[i]], na.rm = TRUE)
  plot_text <- data.frame(r2 = paste0("R^2 = ", r2GLM),
                          chi = paste0("Chi = ", anovaGLM),
                          K1 = 1.11, Value = x_coord)
  plot_text <- setNames(plot_text, c("r2", "chi", "K1", i))
  
  b <- ggplot(hybrid_df, aes(x = .data[[i]], y = K1)) +
    geom_point(size = 0.5) +
    coord_cartesian(y = c(0,1.06)) +
    stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
    scale_y_continuous(expand = c(0.1,0,0,0.1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    geom_text(data = plot_text, mapping = aes(label = r2),hjust = 0, size = 3) +
    geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.06), hjust = 0, size = 3) +
    theme_bw()
  
  # append to list
  plot_list_LOGISTIC[[length(plot_list_LOGISTIC)+1]] <- b
}

# stitch together in groups of 5
stitched_plots_LOGISTIC <- plot_grid(
  plot_list_LOGISTIC[[1]],plot_list_LOGISTIC[[2]],plot_list_LOGISTIC[[3]],plot_list_LOGISTIC[[4]],
  plot_list_LOGISTIC[[5]],plot_list_LOGISTIC[[6]],plot_list_LOGISTIC[[7]],plot_list_LOGISTIC[[8]],
  plot_list_LOGISTIC[[9]],plot_list_LOGISTIC[[10]],plot_list_LOGISTIC[[11]],plot_list_LOGISTIC[[12]],
  nrow = 3,
  ncol = 4,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1,1,1,1) # 2.6 for 6 classes.
)
ggsave("results - regressions/PLOT-regressions-RFimportant-traits-LOGISTIC.png",
       stitched_plots_LOGISTIC, height = 7.5, width = 11, dpi = 600)
################################################################################


# 5. Floral PCA -- all samples (with genotype data)
################################################################################
# re-list the important traits (remember text about what these came from)
important_traits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                      "avg_mouth_height_width_ratio")

# Filter to only samples with genetic data and to most important traits
pca_input <- admix_included %>%
  filter(!is.na(class))

impute_input <- pca_input %>%
  select(all_of(important_traits))

# Impute missing data, subset to most important traits 
imputed_matrix <- complete(mice(impute_input, method = "pmm"))

# Scale the data
scaled_data <- scale(imputed_matrix)

# Perform the PCA
pca_result <- prcomp(scaled_data)

# Extract PC scores, loadings, and eigenvalues
pc_scores <- as.data.frame(pca_result$x)
loadings <- as.data.frame(pca_result$rotation)
eigenvalues <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Write the loadings and eigenvalues to a new file
write.csv(loadings, "loadings-allsamples-floral-pca-IMPORTANTTRAITS.csv")
write.csv(eigenvalues, "eigenvalues-allsamples-floral-pca-IMPORTANTTRAITS.csv", row.names = F)

# Generate new matrix that includes PC scores
pc_included <- pc_scores %>% 
  mutate(newID = pca_input$newID) %>% 
  full_join(., pca_input, by = 'newID')

# Write this to disk
write_csv(pc_included, file = "PCA-results-IMPORTANTTRAITS.csv")


# Read back in if desired
pc_included <- read_csv("PCA-results-IMPORTANTTRAITS.csv")
loadings <- read.csv("loadings-allsamples-floral-pca-IMPORTANTTRAITS.csv", row.names = 1)
eigenvalues <- read.csv("eigenvalues-allsamples-floral-pca-IMPORTANTTRAITS.csv")[,1]

# Plot the PCA, coloring by rgb.
ggplot(pc_included, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = rgb)) +
  scale_color_identity() +
  labs(title = "Floral trait PCA -- all samples -- color point",
       x = paste0("PC1 (", round(eigenvalues[1],digits = 2), "% Variance)"),
       y = paste0("PC2 (", round(eigenvalues[2],digits = 2), "% Variance)"),
       col = "Parent species") +
  theme_bw()
ggsave("PLOT-floral-PCA-allsamples-RGBcolor-IMPORTANTTRAITS.png", width = 4.75, height = 4.5, dpi = 600)


# Plot the PCA, coloring by class
mycols <- c("newberryi"="palevioletred1", "davidsonii"="mediumorchid4", "hybrid"="mediumorchid2")
ggplot(pc_included, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = class)) +
  scale_color_manual(values = mycols) +
  labs(title = "Floral trait PCA -- all samples -- class point",
       x = paste0("PC1 (", round(eigenvalues[1],digits = 2), "% Variance)"),
       y = paste0("PC2 (", round(eigenvalues[2],digits = 2), "% Variance)"),
       col = "Parent species") +
  theme_bw()
ggsave("PLOT-floral-PCA-allsamples-CLASScolor-IMPORTANTTRAITS.png", width = 6, height = 4.5, dpi = 600)



# Looking more closely at loadings
# ggbiplot multiplies endpoint by 2.5, with scale = 0
loadingsplot <- loadings %>%
  rownames_to_column(var = "var")

ggplot(pc_included, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = K1)) +
  scale_color_gradient(low="steelblue3", high="palevioletred3") +
  labs(x = paste0("PC1 (", round(eigenvalues[1],digits = 2), "% Variance)"),
       y = paste0("PC2 (", round(eigenvalues[2],digits = 2), "% Variance)"),
       col = "P. newberryi\nancestry") +
  geom_segment(data = loadingsplot, aes(x = 0, y = 0, xend = PC1 * 2.5, yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "orange") +
  #geom_text(data = loadingsplot, aes(x = PC1 * 5.6, y = PC2 * 5.6, label = var), size = 3) +
  theme_bw()
  # theme(legend.position.inside = "left",
  #       legend.position = c(0.15,0.25),
  #       legend.background = element_blank())
ggsave("PLOT-floral-PCA-allsamples-ADMIXTUREcolor-IMPORTANTTRAITS-withloadings1.pdf",
       width = 5.5, height = 4.5, dpi = 600)




# Correlation between PC axes and Admixture Proportion
# PC1
myglm <- glm(data = pc_included, K1 ~ PC1, family = binomial(link = "logit"))
r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- 2.5
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        K1 = 0.92, Value = x_coord,
                        PC1 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

a <- ggplot(data = pc_included, aes(x = PC1, y = K1)) +
  geom_point(size = 0.5) +
  labs(x = "Floral PC1",
       y= "Proportion P. newberryi ancestry") +
  coord_cartesian(y = c(0,1)) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()


# PC2
myglm <- glm(data = pc_included, K1 ~ PC2, family = binomial(link = "logit"))
r2GLM <- formatC(signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3),format = "e", digits = 2)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- min(pc_included$PC2, na.rm = TRUE)
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        K1 = 0.92,
                        PC2 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

b <- ggplot(data = pc_included, aes(x = PC2, y = K1)) +
  geom_point(size = 0.5) +
  labs(x = "Floral PC2",
       y= "Proportion P. newberryi ancestry") +
  coord_cartesian(y = c(0,1)) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()



# stitch together the two regressions
stitched_PC_vs_K1 <- plot_grid(
  a, b,
  nrow = 1,
  ncol = 2,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1,1)
)
ggsave("PLOT-regressions-PCA-vs-K1-LOGISTIC.png",
       stitched_PC_vs_K1, height = 3, width = 6, dpi = 600)



################################################################################



# 6. Floral PCA plotting, vs. admixture, vs. elevation
# Separated by hybrid zone
################################################################################
# Read in elevation data
elevation_data <- read_delim("~/11.elevation_plots/sierras-pops-NEWelevation.txt", delim = "\t") %>%
  rename(newID = accession)

# Read back PCA results
# pc_included <- read_csv("PCA-results-IMPORTANTTRAITS.csv")
# loadings <- read.csv("loadings-allsamples-floral-pca-IMPORTANTTRAITS.csv", row.names = 1)
# eigenvalues <- read.csv("eigenvalues-allsamples-floral-pca-IMPORTANTTRAITS.csv")[,1]

# Merge elevation data into pc_included
pc_included <- left_join(pc_included, elevation_data, by = "newID")

# Make VL and GL pc_includeds, and rescale elevation per pop
pc_included_VL <- pc_included %>%
  filter(pop == 129) %>%
  mutate(elevation_scaled = rescale(newelevation, to = c(0,1)))

pc_included_GL <- pc_included %>%
  filter(pop != 129) %>%
  mutate(elevation_scaled = rescale(newelevation, to = c(0,1)))


# A. 
# Correlation between PC1 and Admixture Proportion -- (Virginia Lakes)
myglm <- glm(data = pc_included_VL, K1 ~ PC1, family = binomial(link = "logit"))
r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- 2.5
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        K1 = 0.88, Value = x_coord,
                        PC1 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

a <- ggplot(data = pc_included_VL, aes(x = PC1, y = K1)) +
  geom_point(size = 0.5, col = "darkblue") +
  labs(x = "Floral PC1",
       y= "Proportion P. newberryi ancestry") +
  coord_cartesian(y = c(0,1), x = c(min(pc_included$PC1), max(pc_included$PC1))) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()

a
###


# B. 
# Correlation between PC1 and elevation -- (Virginia Lakes)
myglm <- glm(data = pc_included_VL, elevation_scaled ~ PC1, family = binomial(link = "logit"))
r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- 2.5
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        elevation_scaled = 0.88, Value = x_coord,
                        PC1 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

b <- ggplot(data = pc_included_VL, aes(x = PC1, y = elevation_scaled)) +
  geom_point(size = 0.5, col = "darkblue") +
  labs(x = "Floral PC1",
       y= "Elevation (scaled)") +
  coord_cartesian(y = c(0,1), x = c(min(pc_included$PC1), max(pc_included$PC1))) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     position = "right") +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = elevation_scaled-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()

b
###

# C.
# Correlation between PC1 and Admixture Proportion -- (Gem Lakes)
myglm <- glm(data = pc_included_GL, K1 ~ PC1, family = binomial(link = "logit"))
r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- 2.5
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        K1 = 0.88, Value = x_coord,
                        PC1 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

c <- ggplot(data = pc_included_GL, aes(x = PC1, y = K1)) +
  geom_point(size = 0.5, col = "skyblue") +
  labs(x = "Floral PC1",
       y= "Proportion P. newberryi ancestry") +
  coord_cartesian(y = c(0,1), x = c(min(pc_included$PC1), max(pc_included$PC1))) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = K1-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()

c
###

# Correlation between PC1 and elevation -- (Gem Lakes)
myglm <- glm(data = pc_included_GL, elevation_scaled ~ PC1, family = binomial(link = "logit"))
r2GLM <- signif(nagelkerke(myglm)$Pseudo.R.squared.for.model.vs.null[3], digits = 3)
anovaGLM <- signif(anova(myglm, test = "LRT")$`Pr(>Chi)`[2], digits = 3)
x_coord <- 2.5
plot_text <- data.frame(chi = paste0("p = ", anovaGLM),
                        elevation_scaled = 0.88, Value = x_coord,
                        PC1 = x_coord)
plot_text$r2 <- list(bquote(italic(R)^2 == .(r2GLM)))

d <- ggplot(data = pc_included_GL, aes(x = PC1, y = elevation_scaled)) +
  geom_point(size = 0.5, col = "skyblue") +
  labs(x = "Floral PC1",
       y= "Elevation (scaled)") +
  coord_cartesian(y = c(0,1), x = c(min(pc_included$PC1), max(pc_included$PC1))) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  scale_y_continuous(expand = c(0.01,0,0,0.01),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     position = "right") +
  geom_text(data = plot_text, aes(label = r2), parse = TRUE, hjust = 0, size = 3, col = "red") +
  geom_text(data = plot_text, mapping = aes(label = chi, y = elevation_scaled-0.06), hjust = 0, size = 3, col = "red") +
  theme_bw()

d
###


# stitch together the regressions
stitched_PC_regressions_split <- plot_grid(
  a, b, c, d,
  nrow = 2,
  ncol = 2,
  #align = "hv",
  #axis = "none",
  greedy = TRUE,
  rel_widths = c(1,1,1,1)
)
stitched_PC_regressions_split

ggsave("PLOT-regressions-PC-regressions-GLM.pdf",
       stitched_PC_regressions_split, height = 6, width = 6, dpi = 600)

################################################################################
