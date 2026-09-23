#################################################################################################
#################################################################################################
       ##### BEGINNING OF STRUCTURAL DIVERSITY AND ASSOCIATIONAL RESISTANCE R SCRIPT #####
#################################################################################################
#################################################################################################
### Load in required packages ###
install.packages("spatialreg")
library(spatialreg)
library(tidyverse)
library(lme4)
library(broom)
library(performance)
library(spdep)
library(sf)

### Load in data ###
mergeddata_climate_joinedALL4 <- read_csv("mergeddata_climate_joinedALL4.csv")
mergeddata_climate_joinedGEN4 <- read_csv("mergeddata_climate_joinedGEN4.csv")
mergeddata_climate_joinedSPEC4 <- read_csv("mergeddata_climate_joinedSPEC4.csv")

#################################################################################################
#################################################################################################
                    ##### PIECEWISE STRUCTURAL EQUATION MODELING #####
#################################################################################################
#################################################################################################

### Load in shapefiles for calculating neighborhood list ###
merged_with_geometry_sfALL <- st_read("path/to/location/merged_with_geometry_sfALL.shp")

# Check the CRS
st_crs(merged_with_geometry_sfALL)

# Check geometry is intact
sum(is.na(merged_with_geometry_sfALL$geometry))  # Should be 0

# View structure
head(merged_with_geometry_sfALL)

#### creating neighborhood list for overall insect dataset ###
county_centroids_points <- st_centroid(merged_with_geometry_sfALL)
county_nearest_neighbors <- knn2nb(knearneigh(county_centroids_points, k = 5)) 
nbw <- nb2listw(county_nearest_neighbors, style = "W")

############################################
### overall insects, s_diameter_cv_diff ###
############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(Total_pests ~ s_prop_hosts + s_pop_density +  s_diameter_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMALL)

# --- reg2: spatial regression for diameter_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  s_diameter_cv_diff %~~% s_prop_hosts,
  Total_pests %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

##########################################
### overall insects, s_height_cv_diff ###
##########################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(Total_pests ~ s_prop_hosts + s_pop_density +  s_height_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMALL)

# --- reg2: spatial regression for height_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  Total_pests %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

############################################
### overall insects, s_height_mean_diff ###
############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(Total_pests ~ s_prop_hosts + s_pop_density +  s_height_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMALL)

# --- reg2: spatial regression for height_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  Total_pests %~~% s_lat,
  s_height_mean_diff %~~% s_pop_density)

# --- Inspect SEM ---
summary(modelList)

##############################################
### overall insects, s_diameter_mean_diff ###
##############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(Total_pests ~ s_prop_hosts + s_pop_density +  s_diameter_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMALL)

# --- reg2: spatial regression for diameter_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMALL,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  s_diameter_mean_diff %~~% s_prop_hosts,
  Total_pests %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

#########################################################################################################################
### generalist insects ###
#########################################################################################################################
### Load in shapefiles for calculating neighborhood list ###
merged_with_geometry_sfGEN <- st_read("path/to/location/merged_with_geometry_sfGEN.shp")

# Check the CRS
st_crs(merged_with_geometry_sfGEN)

# Check geometry is intact
sum(is.na(merged_with_geometry_sfGEN$geometry))  # Should be 0

# View structure
head(merged_with_geometry_sfGEN)

# Creating neighborhood list
county_centroids_points <- st_centroid(merged_with_geometry_sfGEN)
county_nearest_neighbors <- knn2nb(knearneigh(county_centroids_points, k = 5)) 
nbw <- nb2listw(county_nearest_neighbors, style = "W")

###############################################
### generalist insects, s_diameter_cv_diff ###
###############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(generalist_count ~ s_prop_gen_hosts + s_pop_density +  s_diameter_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMGEN)

# --- reg2: spatial regression for diameter_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_gen_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  generalist_count %~~% s_lat,
  s_diameter_cv_diff %~~% s_prop_gen_hosts)

# --- Inspect SEM ---
summary(modelList)

#############################################
### generalist insects, s_height_cv_diff ###
#############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(generalist_count ~ s_prop_gen_hosts + s_pop_density +  s_height_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMGEN)

# --- reg2: spatial regression for height_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_gen_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  generalist_count %~~% s_lat,
  s_height_cv_diff %~~% s_prop_gen_hosts)

# --- Inspect SEM ---
summary(modelList)

###############################################
### generalist insects, s_height_mean_diff ###
###############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(generalist_count ~ s_prop_gen_hosts + s_pop_density +  s_height_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMGEN)

# --- reg2: spatial regression for height_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_gen_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  generalist_count %~~% s_lat,
  s_height_mean_diff %~~% s_pop_density
)

# --- Inspect SEM ---
summary(modelList)

#################################################
### generalist insects, s_diameter_mean_diff ###
#################################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(generalist_count ~ s_prop_gen_hosts + s_pop_density + s_diameter_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMGEN)

# --- reg2: spatial regression for diameter_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_gen_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_lat,
  data = SCALED_CLEAN_SEMGEN,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  s_diameter_mean_diff %~~% s_prop_gen_hosts, 
  generalist_count %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

#########################################################################################################################
### specialist insects ###
#########################################################################################################################
### Load in shapefiles for calculating neighborhood list ###
merged_with_geometry_sfSPEC <- st_read("path/to/location/merged_with_geometry_sfSPEC.shp")

# Check the CRS
st_crs(merged_with_geometry_sfSPEC)

# Check geometry is intact
sum(is.na(merged_with_geometry_sfSPEC$geometry))  # Should be 0

# View structure
head(merged_with_geometry_sfSPEC)

# Creating neighborhood list
county_centroids_points <- st_centroid(merged_with_geometry_sfSPEC)
county_nearest_neighbors <- knn2nb(knearneigh(county_centroids_points, k = 5)) 
nbw <- nb2listw(county_nearest_neighbors, style = "W")

###############################################
### specialist insects, s_diameter_cv_diff ###
###############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(specialist_count ~ s_prop_spec_hosts + s_pop_density +  s_diameter_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMSPEC)

# --- reg2: spatial regression for diameter_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_spec_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  s_diameter_cv_diff %~~% s_pop_density,
  specialist_count %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

#############################################
### specialist insects, s_height_cv_diff ###
#############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(specialist_count ~ s_prop_spec_hosts + s_pop_density +  s_height_cv_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMSPEC)

# --- reg2: spatial regression for height_cv_t_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_cv_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_spec_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  s_height_cv_diff %~~% s_pop_density,
  specialist_count %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

###############################################
### specialist insects, s_height_mean_diff ###
###############################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(specialist_count ~ s_prop_spec_hosts + s_pop_density +  s_height_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMSPEC)

# --- reg2: spatial regression for height_difference_NHminH ---
reg2 <- lagsarlm(
  s_height_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_spec_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  specialist_count %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

#################################################
### specialist insects, s_diameter_mean_diff ###
#################################################

# --- reg1: mixed model with random effect ---
reg1 <- glmer(specialist_count ~ s_prop_spec_hosts + s_pop_density +  s_diameter_mean_diff +
                s_MEAN_MAT + s_MEAN_MAP + s_tree_div + (1|US_L3NAME), family = poisson(link = "log"),
              data = SCALED_CLEAN_SEMSPEC)

# --- reg2: spatial regression for diameter_difference_NHminH ---
reg2 <- lagsarlm(
  s_diameter_mean_diff ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- reg3: spatial regression for prop_hosts ---
reg3 <- lagsarlm(
  s_prop_spec_hosts ~ s_tree_div + s_MEAN_MAT + s_MEAN_MAP + s_lat,
  data = SCALED_CLEAN_SEMSPEC,
  listw = nbw,
  type = "lag"
)

# --- Piecewise SEM ---
modelList <- psem(
  reg1,
  reg2,
  reg3,
  specialist_count %~~% s_lat)

# --- Inspect SEM ---
summary(modelList)

#################################################################################################
#################################################################################################
##### INTERACTION EFFECTS #####
### Investigating the effect of structural metrics on insect richness ###
### across varying forest types and temperature bins (low, moderate, high) ### 
#################################################################################################
#################################################################################################

### Load in required packages ###
library(dplyr)
library(ggplot2)
library(emmeans)
install.packages("ggrepel")
library(ggrepel)
library(car)
library(sjPlot)
library(sjmisc)
library(readr)
library(lme4)
library(lmerTest)
library(patchwork)

### Load in datasets ###
SCALED_CLEAN_SEMALL <- read_csv("C:/Users/sromy2/Downloads/IIP project\\SCALED_CLEAN_SEMALL.csv")
SCALED_CLEAN_SEMGEN <- read_csv("C:/Users/sromy2/Downloads/IIP project\\SCALED_CLEAN_SEMGEN.csv")
SCALED_CLEAN_SEMSPEC <- read_csv("C:/Users/sromy2/Downloads/IIP project\\SCALED_CLEAN_SEMSPEC.csv")

#changing datatype for US_L3NAME
SCALED_CLEAN_SEMALL$US_L3NAME <- factor(SCALED_CLEAN_SEMALL$US_L3NAME)
SCALED_CLEAN_SEMGEN$US_L3NAME <- factor(SCALED_CLEAN_SEMGEN$US_L3NAME)
SCALED_CLEAN_SEMSPEC$US_L3NAME <- factor(SCALED_CLEAN_SEMSPEC$US_L3NAME)

#################################
### CALCULATING SAMPLE SIZES ###
#################################

sample_sizes <- SCALED_CLEAN_SEMALL %>%
  count(broad_group, name = "n_counties")

sample_sizes <- SCALED_CLEAN_SEMGEN %>%
  count(broad_group, name = "n_counties")

sample_sizes <- SCALED_CLEAN_SEMSPEC %>%
  count(broad_group, name = "n_counties")

#########################################################################################################################
### total insects ###
#########################################################################################################################

###################################################
### interaction ran with SCALED_CLEAN_SEMALL and s_diameter_cv_diff ###
###################################################
# --- fit model ---
fit <- lmer(
  Total_pests ~ s_prop_hosts + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_tree_div +
    s_diameter_cv_diff * broad_group * temp_bin +
    (1 | US_L3NAME),
  data = SCALED_CLEAN_SEMALL
)

# --- get emtrends and tidy results (with p-values) ---
emm_slope <- emtrends(fit, ~ broad_group * temp_bin, var = "s_diameter_cv_diff")
emm_slope_df <- summary(emm_slope, infer = TRUE) %>% as.data.frame()

# quick check: what columns exist
print(names(emm_slope_df))

# --- ensure factors and initial sig flag ---
emm_slope_df <- emm_slope_df %>%
  mutate(
    broad_group = factor(broad_group),
    temp_bin = factor(temp_bin, levels = c("Low", "Moderate", "High")),
    sig = ifelse("p.value" %in% names(.), p.value <= 0.05, NA)
  )

emm_slope_df <- emm_slope_df %>%
  left_join(sample_sizes, by = "broad_group") %>%
  mutate(
    broad_group_label = paste0(broad_group, " (n=", n_counties, ")")
  )

# --- compute min, max, and range per forest group across bins ---
range_df <- emm_slope_df %>%
  group_by(broad_group, temp_bin) %>%
  summarize(
    min_val = min(lower.CL, na.rm = TRUE),
    max_val = max(upper.CL, na.rm = TRUE),
    range_val = max_val - min_val,
    .groups = "drop"
  )

# identify groups with extreme ranges (min < -5 or max > 5)
removed_groups_df <- range_df %>%
  filter(max_val > 5 | min_val < -5) %>%
  distinct(broad_group)

# show which groups are removed
if (nrow(removed_groups_df) == 0) {
  message("No groups meet the removal threshold (max>5 or min < -5).")
} else {
  message("Groups that will be removed (extreme slopes):")
  print(removed_groups_df)
}

removed_groups <- as.character(removed_groups_df$broad_group)

# --- filter the dataset to drop extreme groups ---
emm_slope_df_filtered <- emm_slope_df %>%
  filter(!(broad_group %in% removed_groups))

# --- recalculate significance cleanly after filtering ---
emm_slope_df_filtered <- emm_slope_df_filtered %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "NA",
      p.value <= 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Significant", "Not Significant", "NA"))
  )

# --- create annotation text for removed groups ---
if (length(removed_groups) > 0) {
  removed_text <- paste0(
    "Removed groups (max > 5 or min < -5): ",
    paste0(removed_groups_df$broad_group, collapse = "; ")
  )
} else {
  removed_text <- "No forest groups removed."
}

# --- plot ---
pos <- position_dodge(width = 0.6)

p <- ggplot(emm_slope_df_filtered %>% filter(!is.na(s_diameter_cv_diff.trend)),
            aes(x = s_diameter_cv_diff.trend,
                y = broad_group_label,    # <<<<< NEW LABEL HERE
                alpha = sig)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_point(size = 3, position = pos) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, position = pos) +
  
  # highlight significant points
  geom_point(
    data = subset(emm_slope_df_filtered, p.value <= 0.05),
    aes(x = s_diameter_cv_diff.trend, y = broad_group_label),
    shape = 21, fill = "blue", size = 2, stroke = 0.6, position = pos
  ) +
  
  facet_wrap(~ temp_bin, scales = "fixed") +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.3, "NA" = 0.3),
                     guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Effect of Diameter Variance on Invasive Insect Richness",
    subtitle = removed_text,
    x = "Effect (slope)",
    y = "Forest Type (n = sample size)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    plot.subtitle = element_text(color = "darkgray", size = 10)
  )

print(p)

############################################################################
##### interaction ran with SCALED_CLEAN_SEMALL and s_height_cv_diff #####
############################################################################
# --- fit model ---
fit <- lmer(
  Total_pests ~ s_prop_hosts + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_tree_div +
    s_height_cv_diff * broad_group * temp_bin +
    (1 | US_L3NAME),
  data = SCALED_CLEAN_SEMALL
)

# --- get emtrends and tidy results (with p-values) ---
emm_slope <- emtrends(fit, ~ broad_group * temp_bin, var = "s_height_cv_diff")
emm_slope_df <- summary(emm_slope, infer = TRUE) %>% as.data.frame()

# --- ensure factors and initial sig flag ---
emm_slope_df <- emm_slope_df %>%
  mutate(
    broad_group = factor(broad_group),
    temp_bin = factor(temp_bin, levels = c("Low", "Moderate", "High")),
    sig = ifelse("p.value" %in% names(.), p.value <= 0.05, NA)
  )

emm_slope_df <- emm_slope_df %>% left_join(sample_sizes, by = "broad_group") %>%
  mutate(broad_group_label = paste0(broad_group, " (n=", n_counties, ")")
  )

# --- compute min, max, and range per forest group across bins ---
range_df <- emm_slope_df %>%
  group_by(broad_group, temp_bin) %>%
  summarize(
    min_val = min(lower.CL, na.rm = TRUE),
    max_val = max(upper.CL, na.rm = TRUE),
    range_val = max_val - min_val,
    .groups = "drop"
  )

# identify groups with extreme ranges (min < -5 or max > 5)
removed_groups_df <- range_df %>%
  filter(max_val > 5 | min_val < -5) %>%
  distinct(broad_group)

# show which groups are removed
if (nrow(removed_groups_df) == 0) {
  message("No groups meet the removal threshold (max>5 or min < -5).")
} else {
  message("Groups that will be removed (extreme slopes):")
  print(removed_groups_df)
}

removed_groups <- as.character(removed_groups_df$broad_group)

# --- filter the dataset to drop extreme groups ---
emm_slope_df_filtered <- emm_slope_df %>%
  filter(!(broad_group %in% removed_groups))

# --- recalculate significance cleanly after filtering ---
emm_slope_df_filtered <- emm_slope_df_filtered %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "NA",
      p.value <= 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Significant", "Not Significant", "NA"))
  )

# --- create annotation text for removed groups ---
if (length(removed_groups) > 0) {
  removed_text <- paste0(
    "Removed groups (max > 5 or min < -5): ",
    paste0(removed_groups_df$broad_group, collapse = "; ")
  )
} else {
  removed_text <- "No forest groups removed."
}

# --- plot ---
pos <- position_dodge(width = 0.6)

p <- ggplot(emm_slope_df_filtered %>% filter(!is.na(s_height_cv_diff.trend)),
            aes(x = s_height_cv_diff.trend,
                y = broad_group_label, #removed "color = broad_group,"
                alpha = sig)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_point(size = 3, position = pos) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, position = pos) +
  
  # emphasize significant points (black outline)
  geom_point(
    data = subset(emm_slope_df_filtered, p.value <= 0.05),
    aes(x = s_height_cv_diff.trend, y = broad_group_label),
    shape = 21, fill = "blue", size = 2, stroke = 0.6, position = pos
  ) +
  
  facet_wrap(~ temp_bin, scales = "fixed") +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.3, "NA" = 0.3),
                     guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Effect of Height Variance on Invasive Insect Richness",
    subtitle = removed_text,
    x = "Effect (slope)",
    y = "Forest Type (n = sample size)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none", #changed from "bottom"
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    plot.subtitle = element_text(color = "darkgray", size = 10)
  )

print(p)

############################################################################
##### interaction ran with SCALED_CLEAN_SEMALL and s_height_mean_diff #####
############################################################################
# --- fit model ---
fit <- lmer(
  Total_pests ~ s_prop_hosts + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_tree_div +
    s_height_mean_diff * broad_group * temp_bin +
    (1 | US_L3NAME),
  data = SCALED_CLEAN_SEMALL
)

# --- get emtrends and tidy results (with p-values) ---
emm_slope <- emtrends(fit, ~ broad_group * temp_bin, var = "s_height_mean_diff")
emm_slope_df <- summary(emm_slope, infer = TRUE) %>% as.data.frame()

# quick check: what columns exist
#print(names(emm_slope_df))

# --- ensure factors and initial sig flag ---
emm_slope_df <- emm_slope_df %>%
  mutate(
    broad_group = factor(broad_group),
    temp_bin = factor(temp_bin, levels = c("Low", "Moderate", "High")),
    sig = ifelse("p.value" %in% names(.), p.value <= 0.05, NA)
  )

emm_slope_df <- emm_slope_df %>% left_join(sample_sizes, by = "broad_group") %>%
  mutate(broad_group_label = paste0(broad_group, " (n=", n_counties, ")")
  )

# --- compute min, max, and range per forest group across bins ---
range_df <- emm_slope_df %>%
  group_by(broad_group, temp_bin) %>%
  summarize(
    min_val = min(lower.CL, na.rm = TRUE),
    max_val = max(upper.CL, na.rm = TRUE),
    range_val = max_val - min_val,
    .groups = "drop"
  )

# identify groups with extreme ranges (min < -5 or max > 5)
removed_groups_df <- range_df %>%
  filter(max_val > 5 | min_val < -5) %>%
  distinct(broad_group)

# show which groups are removed
if (nrow(removed_groups_df) == 0) {
  message("No groups meet the removal threshold (max>5 or min < -5).")
} else {
  message("Groups that will be removed (extreme slopes):")
  print(removed_groups_df)
}

removed_groups <- as.character(removed_groups_df$broad_group)

# --- filter the dataset to drop extreme groups ---
emm_slope_df_filtered <- emm_slope_df %>%
  filter(!(broad_group %in% removed_groups))

# --- recalculate significance cleanly after filtering ---
emm_slope_df_filtered <- emm_slope_df_filtered %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "NA",
      p.value <= 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Significant", "Not Significant", "NA"))
  )

# --- create annotation text for removed groups ---
if (length(removed_groups) > 0) {
  removed_text <- paste0(
    "Removed groups (max > 5 or min < -5): ",
    paste0(removed_groups_df$broad_group, collapse = "; ")
  )
} else {
  removed_text <- "No forest groups removed."
}

# --- plot ---
pos <- position_dodge(width = 0.6)

p <- ggplot(emm_slope_df_filtered %>% filter(!is.na(s_height_mean_diff.trend)),
            aes(x = s_height_mean_diff.trend,
                y = broad_group_label, #removed "+ color = broad_group,"
                alpha = sig)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_point(size = 3, position = pos) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, position = pos) +
  
  # emphasize significant points (black outline)
  geom_point(
    data = subset(emm_slope_df_filtered, p.value <= 0.05),
    aes(x = s_height_mean_diff.trend, y = broad_group_label),
    shape = 21, fill = "blue", size = 2, stroke = 0.6, position = pos
  ) +
  
  facet_wrap(~ temp_bin, scales = "fixed") +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.3, "NA" = 0.3),
                     guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Effect of Height on Invasive Insect Richness",
    subtitle = removed_text,
    x = "Effect (slope)",
    y = "Forest Type (n = sample size)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    plot.subtitle = element_text(color = "darkgray", size = 10)
  )

print(p)

############################################################################
##### interaction ran with SCALED_CLEAN_SEMALL and s_diameter_mean_diff #####
############################################################################
# --- fit model ---
fit <- lmer(
  Total_pests ~ s_prop_hosts + s_MEAN_MAT + s_MEAN_MAP + s_pop_density + s_tree_div +
    s_diameter_mean_diff * broad_group * temp_bin +
    (1 | US_L3NAME),
  data = SCALED_CLEAN_SEMALL
)

# --- get emtrends and tidy results (with p-values) ---
emm_slope <- emtrends(fit, ~ broad_group * temp_bin, var = "s_diameter_mean_diff")
emm_slope_df <- summary(emm_slope, infer = TRUE) %>% as.data.frame()

# --- ensure factors and initial sig flag ---
emm_slope_df <- emm_slope_df %>%
  mutate(
    broad_group = factor(broad_group),
    temp_bin = factor(temp_bin, levels = c("Low", "Moderate", "High")),
    sig = ifelse("p.value" %in% names(.), p.value <= 0.05, NA)
  )

emm_slope_df <- emm_slope_df %>% left_join(sample_sizes, by = "broad_group") %>%
  mutate(broad_group_label = paste0(broad_group, " (n=", n_counties, ")")
  )

# --- compute min, max, and range per forest group across bins ---
range_df <- emm_slope_df %>%
  group_by(broad_group, temp_bin) %>%
  summarize(
    min_val = min(lower.CL, na.rm = TRUE),
    max_val = max(upper.CL, na.rm = TRUE),
    range_val = max_val - min_val,
    .groups = "drop"
  )

# identify groups with extreme ranges (min < -5 or max > 5)
removed_groups_df <- range_df %>%
  filter(max_val > 5 | min_val < -5) %>%
  distinct(broad_group)

# show which groups are removed
if (nrow(removed_groups_df) == 0) {
  message("No groups meet the removal threshold (max>5 or min < -5).")
} else {
  message("Groups that will be removed (extreme slopes):")
  print(removed_groups_df)
}

removed_groups <- as.character(removed_groups_df$broad_group)

# --- filter the dataset to drop extreme groups ---
emm_slope_df_filtered <- emm_slope_df %>%
  filter(!(broad_group %in% removed_groups))

# --- recalculate significance cleanly after filtering ---
emm_slope_df_filtered <- emm_slope_df_filtered %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "NA",
      p.value <= 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Significant", "Not Significant", "NA"))
  )

# --- create annotation text for removed groups ---
if (length(removed_groups) > 0) {
  removed_text <- paste0(
    "Removed groups (max > 5 or min < -5): ",
    paste0(removed_groups_df$broad_group, collapse = "; ")
  )
} else {
  removed_text <- "No forest groups removed."
}

# --- plot ---
pos <- position_dodge(width = 0.6)

p <- ggplot(emm_slope_df_filtered %>% filter(!is.na(s_diameter_mean_diff.trend)),
            aes(x = s_diameter_mean_diff.trend,
                y = broad_group_label, #removed "+ color = broad_group,"
                alpha = sig)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_point(size = 3, position = pos) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, position = pos) +
  
  # emphasize significant points (black outline)
  geom_point(
    data = subset(emm_slope_df_filtered, p.value <= 0.05),
    aes(x = s_diameter_mean_diff.trend, y = broad_group_label),
    shape = 21, fill = "blue", size = 2, stroke = 0.6, position = pos
  ) +
  
  facet_wrap(~ temp_bin, scales = "fixed") +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.3, "NA" = 0.3),
                     guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Effect of Diameter on Invasive Insect Richness",
    subtitle = removed_text,
    x = "Effect (slope)",
    y = "Forest Type (n = sample size)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none", #changed from "bottom"
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    plot.subtitle = element_text(color = "darkgray", size = 10)
  )

print(p)

#################################################################################################
#################################################################################################
##### END OF STRUCTURAL DIVERSITY AND ASSOCIATIONAL RESISTANCE R SCRIPT #####
#################################################################################################
#################################################################################################
