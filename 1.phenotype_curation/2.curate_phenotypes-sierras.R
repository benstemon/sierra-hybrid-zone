library(tidyverse)

setwd("1.phenotype_curation")

##########################
#### DATA PREPARATION ####
##########################
##############################################################################
# Read in spec data and prepare for later join
# These data were curated previously in a different script
spec <- read_delim("data/ALLSAMPLES_sierras_color_chroma_hue_averageflowers.csv", delim = ",") %>%
  select(c(newID, total_brightness, chroma, modified_hue))



# Read in raw internal and external floral measurements
internal <- read_delim("data/RAW_sierras_phenotype_INTERNAL.txt", delim = "\t")
external <- read_delim("data/RAW_sierras_phenotype_EXTERNAL.txt", delim = "\t")


# Combine floral measurements into a single spreadsheet, and calculate additional metrics
combined <- full_join(internal, external, by = c("ID", "pop"))

# Check that all the pop+ID columns are the same in combined, external, and internal
nrow(anti_join(combined, external, by = c("pop", "ID")))
nrow(anti_join(combined, internal, by = c("pop", "ID")))


combined <- combined %>%
  mutate(stamen_length_ratio = long_stamen_length / short_stamen_length) %>% #long stamen: short stamen
  mutate(longstamen_pistil_ratio = long_stamen_length / pistil_length) %>% #pistil length: long stamen
  mutate(longstamen_tube_ratio = long_stamen_length / floral_tube_length) %>% #long stamen: tube length
  mutate(pistil_tube_ratio = pistil_length / floral_tube_length) %>% #pistil length: tube length
  mutate(floral_crease_ratio = crease_to_top_height / crease_to_bottom_height) %>% #top crease height: bottom crease height
  mutate(mouth_height_width_ratio = tube_height_at_mouth / tube_width_at_mouth) %>%
  mutate(throat_inflation_ratio = tube_height_at_inflection / tube_height_at_base)

write.csv(combined, "data/ALLSAMPLES_RAW_sierras_phenotypes.csv", row.names = F)

# Now that we have per-flower metrics, summarize across flowers for each sample
summarized <- combined %>%
  mutate(newID = str_replace(ID, "[AB]", "")) %>% #generate new ID without letters
  group_by(newID, pop) %>%
  mutate(nectary_join = paste(nectary_area_1, nectary_area_2, sep = ",")) %>%
  summarize(avg_long_stamen_length = mean(long_stamen_length, na.rm = T),
            avg_short_stamen_length = mean(short_stamen_length, na.rm = T),
            avg_pistil_length = mean(pistil_length, na.rm = T),
            nectary_join = paste(nectary_join, collapse = ","),
            avg_floral_tube_length = mean(floral_tube_length, na.rm = T),
            avg_crease_to_top_height = mean(crease_to_top_height, na.rm = T),
            avg_crease_to_bottom_height = mean(crease_to_bottom_height, na.rm = T),
            avg_tube_height_at_base = mean(tube_height_at_base, na.rm = T),
            avg_tube_height_at_inflection = mean(tube_height_at_inflection, na.rm = T),
            avg_mouth_angle = mean(mouth_angle, na.rm = T),
            avg_platform_angle = mean(platform_angle, na.rm = T),
            avg_upper_petal_angle = mean(upper_petal_angle, na.rm = T),
            avg_tube_width_at_mouth = mean(tube_width_at_mouth, na.rm = T),
            avg_tube_width_at_base = mean(tube_width_at_base, na.rm = T),
            avg_lower_petal_angle = mean(lower_petal_angle, na.rm = T),
            avg_tube_height_at_mouth = mean(tube_height_at_mouth, na.rm = T),
            avg_mouth_constriction_angle = mean(mouth_constriction_angle, na.rm = T),
            avg_stamen_length_ratio = mean(stamen_length_ratio, na.rm = T),
            avg_longstamen_pistil_ratio = mean(longstamen_pistil_ratio, na.rm = T),
            avg_longstamen_tube_ratio = mean(longstamen_tube_ratio, na.rm = T),
            avg_pistil_tube_ratio = mean(pistil_tube_ratio, na.rm = T),
            avg_floral_crease_ratio = mean(floral_crease_ratio, na.rm = T),
            avg_mouth_height_width_ratio = mean(mouth_height_width_ratio, na.rm = T),
            avg_throat_inflation_ratio = mean(throat_inflation_ratio, na.rm = T),) %>%
  separate_wider_delim(nectary_join, delim = ",",
                       names = c("nec1", "nec2", "nec3", "nec4"),
                       too_few = "align_start") %>%
  mutate_all(~ifelse(. %in% c("NaN", "NA"), NA, .)) %>%
  mutate_at(c("nec1", "nec2", "nec3", "nec4"), as.numeric) %>%
  mutate(avg_nectary_area = mean(c(nec1, nec2, nec3, nec4), na.rm = T)) %>%
  mutate_all(~ifelse(. %in% c("NaN"), NA, .)) %>%
  select(-c(nec1, nec2, nec3, nec4)) %>%
  #full_join(., nectar, by = "newID") %>%
  left_join(., spec, by = "newID")


# Check that all the pop+ID columns are the same in summarized and spec
nrow(anti_join(summarized, spec, by = c("newID")))
nrow(anti_join(spec, summarized, by = c("newID")))
# There are 4 extra samples in spec, but they are pop 130 and not useful


# Write this summarized data set into a new file to be used for downstream analyses
write_delim(summarized, file = "data/SUMMARIZED_sierras_phenotype_data.txt", delim = "\t")

##############################################################################

  