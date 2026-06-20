#Title: Master's Project modelling
#Author: Maxine Chapman
#Date: 03/11/26

# ============================================================
# 1. SETUP
# ============================================================

# ------------------------------
# 1.1 Libraries
# ------------------------------
library(dplyr)
library(terra)
library(hydrographr)
library(tmap)
library(viridis)
library(sf)
library(ggplot2)
library(ggsci)
library(patchwork)
library(tidyr)
library(purrr)
library(stringr)

# ------------------------------
# 1.2 Directories
# ------------------------------
# Define working directory (server or local environment)
#wdir <- "/home/mchapman/masters_project" #server
wdir <- "C:/Users/mchapman/Documents/masters_project" #local
setwd(wdir)

# Define main directories for model outputs, figures, and data
figure_dir <- file.path(wdir, "figures")
data_dir <- file.path(wdir, "data")
results_dir <- file.path(wdir, "results")

pdp_dir <- file.path(results_dir, "partial_dependence")
plot_dir <- file.path(pdp_dir, "plots")
#dir.create(plot_dir)
predict_dir <- file.path(data_dir, "raw", "predict_tables")
change_fig_dir <- file.path(figure_dir, "change_figures")
#dir.create(change_fig_dir)
change_dir <- file.path(data_dir, "intermediate", "change_tifs")
layer_dir <- file.path(data_dir, "intermediate", "layer_dir")


# ============================================================
# 2. KEY DATASETS
# ============================================================

# Load occurrence dataset and modelling results
combined_df <- read.csv(file.path(data_dir, "raw", "occurrence_data", "global_df.csv"))
results_df <- read.csv(file.path(results_dir, "working_df2.csv"))

# Filter occurrence data to include only successfully modelled species
combined_df <- combined_df %>% filter(species %in% results_df$species) 

# Extract list of species to process
inverts_list <- unique(results_df$species)
length(inverts_list)


# ============================================================
# 3. CREATING HABITAT SUITABILITY MAPS
# ============================================================

# ---------------------------------
# 3a. Overall map
# ---------------------------------

#check files are there
#list.files(model_dir)

#set aesthetics
pal <- colorRampPalette(c(
  "saddlebrown",
  "tan3",
  "white",
  "#bdb76b",
  "#6b8e23"
))(100)

scale_settings <- list(
  col.scale = tm_scale_continuous(
    limits = c(-100, 100),
    midpoint = 0,
    values = pal
  ),
  col.legend = tm_legend(
    title = "Change in Number of Species Present",
    orientation = "landscape"   # horizontal legend
  )
)

tmap_mode("plot") #static mode of plotting

# ----------------------------
# Build file paths
# ----------------------------
ssp1_path <- file.path(layer_dir, "ssp1_layered_predictions.tif")
ssp3_path <- file.path(layer_dir, "ssp3_layered_predictions.tif")
ssp5_path <- file.path(layer_dir, "ssp5_layered_predictions.tif")

# ----------------------------
# Read rasters safely
# ----------------------------
ssp1 <- try(rast(ssp1_path), silent = TRUE)
ssp3 <- try(rast(ssp3_path), silent = TRUE)
ssp5 <- try(rast(ssp5_path), silent = TRUE)


# ----------------------------
# Aggregate for plotting
# ----------------------------
ssp1_plot <- aggregate(ssp1, fact = 5, fun = mean)
ssp3_plot <- aggregate(ssp3, fact = 5, fun = mean)
ssp5_plot <- aggregate(ssp5, fact = 5, fun = mean)

# ----------------------------
# FUTURE MAPS
# ----------------------------
map1 <- tm_shape(ssp1_plot) +
  do.call(tm_raster, scale_settings) +
  tm_layout(
    legend.show = FALSE,
    # title = paste0(species_title, " (2041 - 2070)"),
    # title.position = c("right", "top"),
    # title.size = 1,
    # inner.margins = c(0.08, 0.02, 0.02, 0.02),
    bg.color = "lightsteelblue1"
  ) +
  tm_scalebar(position = "bottom") +
  tm_credits("SSP126.", position = c("LEFT", "BOTTOM")) 

map2 <- tm_shape(ssp3_plot) +
  do.call(tm_raster, scale_settings) +
  tm_layout(legend.show = FALSE,
            bg.color = "lightsteelblue1") +
  tm_scalebar(position = "bottom") +
  tm_credits("SSP370.", position = c("LEFT", "BOTTOM"))

map3 <- tm_shape(ssp5_plot) +
  do.call(tm_raster, scale_settings) +
  tm_layout(legend.show = FALSE,
            bg.color = "lightsteelblue1") +
  tm_scalebar(position = "bottom") +
  tm_credits("SSP585.", position = c("LEFT", "BOTTOM"))

map4 <- tm_shape(ssp5_plot) +
  do.call(tm_raster, scale_settings) +
  tm_layout(legend.show = TRUE,
            legend.only = TRUE,
            legend.outside = TRUE,
            legend.outside.position = "bottom",
            legend.just = "center")

# ----------------------------
# Combine maps
# ----------------------------
multiple <- tmap_arrange(map1, map2, map3, map4, nrow = 4)

# ----------------------------
# Save output
# ----------------------------
output_file <- file.path(figure_dir,
                         "overall_map.png")

try(
  tmap_save(multiple,
            filename = output_file,
            width = 6,
            height = 8,
            units = "in",
            dpi = 300),
  silent = TRUE
)

# ----------------------------
# 3b. Habitat map per Species
# ----------------------------
`for (each_species in seq_along(inverts_list)) {
  
  each_species <- 1 #TESTING ONLY
  
  species_name <- inverts_list[each_species]
  species_title <- gsub("_", " ", species_name) #change format for figure title
  
  message("Processing: ", species_name)
  
  # ----------------------------
  # Build file paths
  # ----------------------------
  ssp1_path <- file.path(change_dir, paste0(species_name, "_ssp1_change.tif"))
  ssp3_path <- file.path(change_dir, paste0(species_name, "_ssp3_change.tif"))
  ssp5_path <- file.path(change_dir, paste0(species_name, "_ssp5_change.tif"))
  
  # ----------------------------
  # Check if all raster files exist
  # ----------------------------
  if (!all(file.exists(ssp1_path, ssp3_path, ssp5_path))) {
    message("  → Skipping ", species_name, ": Missing raster file(s)")
    next
  }
  
  # ----------------------------
  # Read rasters safely
  # ----------------------------
  ssp1 <- try(rast(ssp1_path), silent = TRUE)
  ssp3 <- try(rast(ssp3_path), silent = TRUE)
  ssp5 <- try(rast(ssp5_path), silent = TRUE)
  
  if (inherits(ssp1, "try-error") |
      inherits(ssp3, "try-error") |
      inherits(ssp5, "try-error") 
      ) {
    message("  → Skipping ", species_name, ": Error reading raster")
    next
  }
  
  # ----------------------------
  # Aggregate for plotting
  # ----------------------------
  ssp1_plot <- aggregate(ssp1, fact = 5, fun = mean)
  ssp3_plot <- aggregate(ssp3, fact = 5, fun = mean)
  ssp5_plot <- aggregate(ssp5, fact = 5, fun = mean)
  
  # ----------------------------
  # FUTURE MAPS
  # ----------------------------
  map1 <- tm_shape(ssp1_plot) +
    do.call(tm_raster, scale_settings) +
    tm_layout(
      legend.show = FALSE,
      title = paste0(species_title, " (2041 - 2070)"),
      title.position = c("right", "top"),
      title.size = 1,
      inner.margins = c(0.08, 0.02, 0.02, 0.02),
      bg.color = "lightsteelblue1"
    ) +
    tm_scalebar(position = "bottom") +
    tm_credits("SSP126.", position = c("LEFT", "BOTTOM")) 
  
  map2 <- tm_shape(ssp3_plot) +
    do.call(tm_raster, scale_settings) +
    tm_layout(legend.show = FALSE,
              bg.color = "lightsteelblue1") +
    tm_scalebar(position = "bottom") +
    tm_credits("SSP370.", position = c("LEFT", "BOTTOM"))
  
  map3 <- tm_shape(ssp5_plot) +
    do.call(tm_raster, scale_settings) +
    tm_layout(legend.show = FALSE,
              bg.color = "lightsteelblue1") +
    tm_scalebar(position = "bottom") +
    tm_credits("SSP585.", position = c("LEFT", "BOTTOM"))
  
  map4 <- tm_shape(ssp5_plot) +
    do.call(tm_raster, scale_settings) +
    tm_layout(legend.show = TRUE,
              legend.only = TRUE,
              legend.outside = TRUE,
              legend.outside.position = "bottom",
              legend.just = "center")

  # ----------------------------
  # Combine maps
  # ----------------------------
  multiple <- tmap_arrange(map1, map2, map3, map4, nrow = 4)
  
  # ----------------------------
  # Save output
  # ----------------------------
  output_file <- file.path(change_fig_dir,
                           paste0(species_name, "_suitability_change.pdf"))
  
  try(
    tmap_save(multiple,
              filename = output_file,
              width = 6,
              height = 8,
              units = "in",
              dpi = 300),
    silent = TRUE
  )
  
  message("  ✓ Finished ", species_name)
}

# ============================================================
# 4. PLOTTING EXPOSURE
# ============================================================

#read in df
exposure_df <- read.csv(file.path(results_dir, "exposure_df2.csv"))

#plot exposure scores for each variable split by scenario
p_exposure <- ggplot(exposure_df, aes(x = variable, y = exposure, fill = variable)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_fill_npg() +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    x = "Climate variable",
    y = "Exposure",
    title = "Exposure distribution by climate variable and scenario"
  )

#save plot
ggsave(
  filename = file.path(figure_dir, "exposure_boxplot.png"),
  plot = p_exposure,
  width = 8,
  height = 6,
  dpi = 600
)

#plot exposure values as proportion vulnerable in each scenario

#change format of data first
prop_data <- exposure_df %>%
  group_by(scenario, species_name) %>%
  summarise(vulnerable = max(check), .groups = "drop") %>%
  group_by(scenario) %>%
  summarise(
    Vulnerable = sum(vulnerable),
    Safe       = 96 - sum(vulnerable),
    .groups    = "drop"
  ) %>%
  pivot_longer(cols = c(Safe, Vulnerable),
               names_to  = "Status",
               values_to = "Count") %>%
  mutate(Status = factor(Status, levels = c("Vulnerable", "Safe")))


#create plot
p_vulnerable <- ggplot(prop_data, aes(x = scenario, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.6) +
  # geom_text(
  #   aes(label = ifelse(Proportion > 0.03, scales::percent(Proportion, accuracy = 1), "")),
  #   position = position_stack(vjust = 0.5),
  #   size = 3.5, colour = "white", fontface = "bold"
  # ) +
  scale_fill_npg() +
  theme_bw(base_size = 12) +
  labs(
    x        = "Scenario",
    y        = "Species counts",
    fill     = NULL
  )

ggsave(
  filename = file.path(figure_dir, "exposure_vulnerable.png"),
  plot = p_vulnerable,
  width = 8,
  height = 6,
  dpi = 600
)

  # ============================================================
  # 5. PLOTTING SENSITIVITY
  # ============================================================

#plot sensitivity scores for each variable split by scenario
p_sensitivity <- ggplot(exposure_df, aes(x = scenario, y = sensitivity, fill = scenario)) +
  geom_boxplot(width = 0.6, outlier.size = 1.5) +
  scale_fill_npg() +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    x = "Scenario",
    y = "Sensitivity",
    title = "Sensitivity distribution across climate scenarios"
  )

ggsave(
  filename = file.path(figure_dir, "sensitivity_boxplot.png"),
  plot = p_sensitivity,
  width = 6,
  height = 4,
  dpi = 600
)

#overall pattern
results_df <- read.csv(file.path(results_dir, "results.csv"))

results_df %>% count(pattern)

# ============================================================
# 6. PARTIAL DEPENDENCY PLOTS
# ============================================================

# ----------------------------
# Overall plot
# ----------------------------

#read in data
files <- list.files(pdp_dir, pattern = "_variable_importance_percent.csv$", full.names = TRUE)

#get species list for successful species
species <- results_df$species

# get species name from filename
file_species <- basename(files) %>%
  str_remove("_variable_importance_percent.csv$")

#filter
filtered_files <- files[file_species %in% species]

#combine individual species into overall df
pre_df <- map_dfr(filtered_files, function(f) {
  read.csv(f) %>%
    mutate(
      file = basename(f),
      species = str_remove(file, "_variable_importance_percent.csv$")
    )
}) %>%
  select(species, variable, importance_percent)

#calculate mean and sd
plot_df <- pre_df %>%
  group_by(variable) %>%
  summarise(
    mean_importance = mean(importance_percent, na.rm = TRUE),
    sd_importance = sd(importance_percent, na.rm = TRUE)
  )

#update variable names
v_lookup <- c(
  "bio01_1981.2010_observed_mean" = "Annual Mean Temperature",
  "bio05_1981.2010_observed_mean" = "Maximum Temperature of Warmest Month",
  "bio06_1981.2010_observed_mean" = "Minimum Temperature of Coldest Month",
  "bio07_1981.2010_observed_mean" = "Annual Temperature Range",
  "bio12_1981.2010_observed_mean" = "Annual Precipitation",
  "bio13_1981.2010_observed_mean" = "Precipitation of Wettest Month",
  "bio14_1981.2010_observed_mean" = "Precipitation of Driest Month",
  "bio15_1981.2010_observed_mean" = "Precipitation Seasonality",
  "bio18_1981.2010_observed_mean" = "Precipitation of Warmest Quarter",
  "bio17_1981.2010_observed_mean" = "Precipitation of Driest Quarter",
  "c100_y2020" = "Tree and Shrub Cover",
  "accumulation_mean" = "Accumulation Mean"
)

plot_df$variable_name <- v_lookup[plot_df$variable]

#make plot
pdp_plot <- ggplot(plot_df, aes(x = reorder(variable_name, mean_importance),
                   y = mean_importance)) +
  stat_summary(fun = mean, geom = "col", fill = "#4DBBD5") +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Mean % Importance"
  ) +
  theme_bw(base_size = 11) +
  geom_errorbar(aes(ymin = mean_importance - sd_importance,
                    ymax = mean_importance + sd_importance),
                width = 0.2) 

#save
ggsave(file.path(figure_dir, "1_pdp_plot.png"),
       pdp_plot,
       width = 10,
       height = 4
)

# ----------------------------
# Individual pdp plots
# ----------------------------

for (each_species in 1:length(inverts_list)) {
  
  # each_species <- 1
  
  sp <- inverts_list[each_species]
  
  #select 5 files
  species_files <- list.files(
    pdp_dir,
    pattern = paste0("^", sp, "_pdp_.*\\.csv$"),
    full.names = TRUE
  )
  
  if (length(species_files) == 0) next
  
  species_files <- species_files[1:min(5, length(species_files))]
  
  #read in files for plotting
  pd_list <- lapply(species_files, read.csv)
  
  #extract variable names from pdp_long with positions var1, var2, var3
  vars <- pdp_long %>%
    dplyr::filter(species_name == sp) %>%
    dplyr::arrange(position) %>%
    dplyr::pull(variable)
  
  vars <- vars[1:min(5, length(vars))]
  
  plots <- list()
  
  for (i in 1:length(pd_list)) {
    
    pd <- pd_list[[i]]
    variable_name <- vars[i]
    
    #make plots like this v1 to v5
    plots[[i]] <- ggplot(pd, aes(x = predictor_value, y = yhat_mean)) +
      geom_line(size = 1.2) +
      theme_classic() +
      labs(
        x = variable_name,
        y = "Presence probability",
        title = paste(i, "-", variable_name)
      )
    
  }
  
  #save a png of all plots together in pdp_dir
  combined_plot <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = gsub("_", " ", sp)
    )
  
  ggsave(
    file.path(plot_dir, paste0(sp, "_pdp_plots.png")),
    combined_plot,
    width = 10,
    height = 8
  )
  
  #informative message
  message("Saved PDP plots for ", sp)
  
}

# ============================================================
# 7. ENVIRONMENTAL VARIABLES
# ============================================================

#read in data
present_predict <- read.csv(file.path(predict_dir, "present_predict_gbif.csv"))
ssp1 <- read.csv(file.path(predict_dir, "ssp126_predict_table.csv"))
ssp3 <- read.csv(file.path(predict_dir, "ssp370_predict_table.csv"))
ssp5 <- read.csv(file.path(predict_dir, "ssp585_predict_table.csv"))

# Add scenario labels
present_predict$scenario <- "Present"
ssp1$scenario <- "SSP1"
ssp3$scenario <- "SSP3"
ssp5$scenario <- "SSP5"

prepare_scenario <- function(df, scenario_name) {
  
  colnames(df) <- gsub("^(bio[0-9]+).*", "\\1", colnames(df))
  
  df$scenario <- scenario_name
  
  return(df)
}


present_predict <- prepare_scenario(present_predict, "present")
ssp1 <- prepare_scenario(ssp1, "ssp1")
ssp3 <- prepare_scenario(ssp3, "ssp3")
ssp5 <- prepare_scenario(ssp5, "ssp5")

# Combine datasets
all_data <- bind_rows(
  present_predict,
  ssp1,
  ssp3,
  ssp5
)

# Identify BIO columns
bio_cols <- grep("^bio", colnames(all_data), value = TRUE)

# Keep only scenario + BIO variables and sample rows for speed
bio_data <- all_data %>%
  select(scenario, all_of(bio_cols)) %>%
  slice_sample(n = 10000)

bio_data <- all_data %>%
  select(scenario, all_of(bio_cols))

# Convert to long format
bio_long <- bio_data %>%
  pivot_longer(
    cols = -scenario,
    names_to = "bio_variable",
    values_to = "value"
  )

# BIO variable labels
bio_lookup <- c(
  "bio01" = "Annual Mean Temperature",
  "bio05" = "Maximum Temperature of Warmest Month",
  "bio06" = "Minimum Temperature of Coldest Month",
  "bio07" = "Annual Temperature Range",
  "bio12" = "Annual Precipitation",
  "bio13" = "Precipitation of Wettest Month",
  "bio14" = "Precipitation of Driest Month",
  "bio15" = "Precipitation Seasonality",
  "bio18" = "Precipitation of Warmest Quarter",
  "bio17" = "Precipitation of Driest Quarter"
)

bio_long$bio_label <- bio_lookup[bio_long$bio_variable]

# Plot boxplots
bioclim <- ggplot(bio_long, aes(x = scenario, y = value, fill = scenario)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~bio_label, scales = "free_y") +
  scale_fill_npg() +
  theme_minimal() +
  labs(
    x = "Scenario",
    y = "Value",
    title = "Comparison of Bioclimatic Variables Across Climate Scenarios") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

ggsave(file.path(figure_dir, "bioclim_variables.png"),
  bioclim,
  width = 10,
  height = 8
)

# ============================================================
# 8. RANGE BOXPLOT
# ============================================================

range_df <- read.csv(file.path(results_dir, "range_df.csv"))

scenario_lookup <- c(
"1" = "ssp1",
"2" = "ssp3",
"3" = "ssp5"
)

range_df$scenario <- scenario_lookup[as.character(range_df$scenario)]

range_plot <- ggplot(range_df,
       aes(x = scenario,
           y = range_change,
           fill = scenario)) +
  scale_fill_npg() +
  geom_boxplot() +
  labs(
    x = "Scenario",
    y = "Range Change (%)"
  ) +
  scale_y_continuous(
    limits = c(-500, 2100)
  ) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

ggsave(file.path(figure_dir, "range_plot.png"))
