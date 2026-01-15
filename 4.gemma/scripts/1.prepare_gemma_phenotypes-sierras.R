# This script prepared data for GEMMA GWAS input, including transforming traits
# To be normally distributed when necessary
library(tidyverse)

setwd("~4.gemma/")

# Read in the curated phenotype data set with PC axes
phenotypes <- read_delim("~/3.admixture-pca/data/PCA_values_phenotypes_admixture-sierras.csv",
                         delim = ",")

# Read in the sample order file for the .vcf (.vcf header)
sample_order <- read_delim("~/2.QC_MAPPING_VARIANTS/data/samplename_changer.txt", delim = "\t", col_names = F) %>%
  select(-X1) %>%
  rename("newID" = X2) %>%
  mutate(newID = str_remove(newID, "BWS_"))


# We have a sample (133-013) with genetic data but no phenotype data, due to a mistake
# Let's create a row for this sample filled with missing data
phenotypes[nrow(phenotypes)+1, 1:ncol(phenotypes)] <- NA
phenotypes[nrow(phenotypes), 1] <- "133-013"


# Reorder the phenotypes data.frame so it matches the sample order
reordered_phenotypes <- phenotypes %>%
  inner_join(sample_order, by = "newID") %>%
  arrange(match(newID, sample_order$newID)) %>%
  select(-c(pop.x, pop.y, K2, rgb, total_brightness, region, subpopulation))

# Write PC2 for the covariate matrix
#covariates <- reordered_phenotypes %>%
#  select(newID, PC1, PC2) %>%
#  mutate(newID = 1)
#write_delim(covariates, col_names = F, delim = "\t",
#            file = "GEMMA-sierras-PC1-2-covariate.txt")


# Optional: Filter out parent taxa and the sample missing phenotypic data
filtered_reordered_phenotypes <- reordered_phenotypes %>%
  filter(class == "hybrid") %>%
  select(-c(class, K1))


# Optional: Write the filtered out taxa to a new file
excludethese <- reordered_phenotypes %>%
  anti_join(filtered_reordered_phenotypes) %>%
  select(newID) %>%
  mutate(newID = paste0("BWS_", newID, sep = ""))
  
write_delim(excludethese, file = "excluded_samples_GWAS.txt", delim = "\t", col_names = F)
  

# This is the final filtered phenotype matrix -- prior to testing for normality
finaldf <- filtered_reordered_phenotypes[,c(1,12:ncol(filtered_reordered_phenotypes))]


# Assess each trait for normality
################################################################################
# Copy of original data
df_normalized <- finaldf

# Initialize a data frame to track transformations
transform_log <- tibble(trait = character(), transformed = logical(), method = character(),
                        reason = character(), original_p = numeric(), transformed_p = numeric())

# Identify phenotype columns (all but sample ID)
phenotypes <- setdiff(names(finaldf), "newID")

# Iterate over each phenotype
for (trait in phenotypes) {
  
  # Split out missing data
  vec <- finaldf[[trait]]
  vec_non_na <- vec[!is.na(vec)]
  
  # Test for normality of original
  shap_p_orig <- shapiro.test(vec_non_na)$p.value
  
  # If original data are approximately normal, leave it
  if (shap_p_orig > 0.05) {
    transform_log <- transform_log %>%
      add_row(trait = trait, transformed = FALSE, method = "none", reason = "Approximately normal",
              original_p = shap_p_orig, transformed_p = NA)
    next
  }
  
  
  # Try multiple transformations
  transformations <- list(
    none = vec_non_na,
    log = if (all(vec_non_na > 0)) log(vec_non_na) else NULL,
    sqrt = if (all(vec_non_na >= 0)) sqrt(vec_non_na) else NULL,
    boxcox = tryCatch({
      bc <- boxCox(lm(vec_non_na ~ 1), lambda = seq(-2, 2, 0.1), plotit = FALSE)
      best_lambda <- bc$x[which.max(bc$y)]
      if (best_lambda == 0) log(vec_non_na) else (vec_non_na ^ best_lambda - 1) / best_lambda
    }, error = function(e) NULL)
  )
  
  # Evaluate Shapiro-Wilk p-values for each available transformation
  shap_results <- map_dbl(transformations, ~{
    if (is.null(.x)) return(NA_real_)
    tryCatch(shapiro.test(.x)$p.value, error = function(e) NA_real_)
  })
  
  # Identify best transformation
  best_method <- names(which.max(shap_results))
  best_p <- max(shap_results, na.rm = TRUE)
  
  if (best_method == "none" || is.na(best_p) || best_p <= shap_p_orig) {
    # Keep original data
    transform_log <- transform_log %>%
      add_row(
        trait = trait, transformed = FALSE, method = "none",
        reason = "No transformation improved normality",
        original_p = shap_p_orig, transformed_p = best_p
      )
  } else {
    # Use transformed data
    vec_transformed <- vec
    vec_transformed[!is.na(vec)] <- transformations[[best_method]]
    df_normalized[[trait]] <- vec_transformed
    
    reason_text <- case_when(
      best_method == "log" ~ "All values > 0",
      best_method == "sqrt" ~ "All values ≥ 0",
      best_method == "boxcox" ~ "Negative or mixed values, used Box-Cox",
      TRUE ~ "Unknown"
    )
    
    transform_log <- transform_log %>%
      add_row(
        trait = trait, transformed = TRUE, method = best_method,
        reason = reason_text,
        original_p = shap_p_orig, transformed_p = best_p
      )
  }
}

# This creates
# 1. df_normalized, which has all the final data after transformation
# 2. transform_log, which tracks for each trait what was done (if anything) and why
View(transform_log)


ggplot(df_normalized, aes(x = avg_platform_angle))+
  geom_density()
ggplot(finaldf, aes(x = avg_platform_angle))+
  geom_density()

################################################################################

# Write the transformation log to disk
write_csv(transform_log, file = "phenotype_transformation_log_noparents.csv")



# Recursively write matrices for each phenotype to be used in GEMMA
for (i in 2:ncol(df_normalized)){
  tmpdf <- df_normalized[,c(1,i)]
  write_delim(tmpdf, delim = "\t",
              file = paste("data/phenotypes_noparents_normalized/GEMMA_INPUTFILE.", colnames(tmpdf)[2], ".txt", sep = ""),
              col_names = F)
}

