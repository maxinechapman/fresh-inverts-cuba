#Title: Master's Project analysis
#Author: Maxine Chapman
#Date: 03/11/26

# ==============================================================
# 01. SETUP: LIBRARIES AND DIRECTORIES
# ==============================================================

library(dplyr)
library(terra)
library(tidyr)
library(sf)
library(hydrographr)
library(stringr)
library(purrr)
library(readr)

wdir <- "~/masters_project"
setwd(wdir)

#set key directories
data_dir <- file.path(wdir, "data")
results_dir <- file.path(wdir, "results")
occ_dir <- file.path(data_dir, "raw", "occurrence_data")
pdp_dir <- file.path(results_dir, "partial_dependence")
model_dir1 <- file.path(data_dir, "intermediate","model_tables")
model_dir2 <- file.path(data_dir, "intermediate","model_tables_2")
predict_dir <- file.path(data_dir, "raw", "predict_tables")

# ==============================================================
# 02. KEY FILES
# ==============================================================

results_df <- read.csv(file.path(results_dir, "results.csv")) #modelling results
global_df <- read.csv(file.path(occ_dir, "global_df.csv")) #gbif observation data

# ==============================================================
# 03. CREATING WORKING_DF
# ==============================================================

# Combine species data
species_df <- global_df %>%
  select(-X, -decimalLatitude, -decimalLongitude)
species_df <- distinct(species_df) # species df

results_df <- results_df %>% mutate(species = species_name) %>%
  select(-species_name, -X, -species_index, -status) #edit columns

working_df <- inner_join(results_df, species_df, by = "species")  
#working_df <- inner_join(working_df, tss_log, by = "species")

combined_df <- global_df %>% filter(species %in% working_df$species)
length(unique(combined_df$species)) #occurrences only for species with good models

# Remove duplicates and filter unwanted species
working_df <- working_df[!duplicated(working_df[c("species", "order")]), ]
working_df <- working_df %>%
  filter(!((species == "Littoridinops_monroensis" & order == "Littorinimorpha") |
             (species == "Nereina_punctulata" & order == "Cycloneritida") |
             (species == "Vitta_virginea" & order == "Cycloneritida"))
  )
working_df <- working_df[!(working_df$order == ""), ]
working_df <- working_df %>% drop_na(species)

# Ensure only one of each species
check <- working_df %>% count(species)
check
# dupes <- check %>% filter(n >1)
# fix <- working_df %>% filter(working_df$species %in% dupes$species) %>% select(species, order)
write.csv(working_df, file.path(results_dir,  "working_df2.csv"))

working_df <- read.csv(file.path(results_dir, "working_df2.csv"))
head(working_df)

# Initial exploration
working_df %>% count(order) %>% arrange(n, sort = TRUE)
working_df %>% count(class)
working_df %>% count(cubanendemic)

# ==============================================================
# 04. GENERATE GPKG FILES
# ==============================================================

#create directory to store output
points_dir <- file.path(results_dir, "points_dir_2")
dir.create(points_dir)
# unlink(list.files(points_dir, full.names = TRUE, recursive = TRUE))

#create df with only cuba occurrences for results species
occ_df <- read.csv(file.path(occ_dir, "species_occurrences_cuba.csv"))
gpkg_df <- occ_df %>% filter(occ_df$species %in% working_df$species)
nrow(gpkg_df)

occ_sf <- st_as_sf(
  gpkg_df,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

species_list <- unique(occ_sf$species) #get list of species

#make files using loop
for (sp in species_list) {
  
  sp_clean <- gsub("[^A-Za-z0-9_]", "_", sp)
  file_name <- file.path(points_dir, paste0(sp_clean, ".gpkg"))
  
  st_write(
    occ_sf[occ_sf$species == sp, ],
    dsn = file_name,
    delete_dsn = TRUE,
    quiet = TRUE
  )
}

# ==============================================================
# 05. CHANGE IN HABITAT SUITABILITY
# ==============================================================

#calculate changes in habitat suitability between present and future

model_dir <- file.path(data_dir, "intermediate", "model_tables") #directory
working_df <- read.csv(file.path(results_dir, "working_df.csv")) #key df
change_dir <- file.path(data_dir, "intermediate", "change_tifs") #directory #output dir
head(working_df)
dir.create(change_dir)

future_scenarios <- c("ssp1", "ssp3", "ssp5")

for (each_species in 1:nrow(working_df)){
  
  #each_species <- 1 #TESTING ONLY
  
  species_name <- working_df[each_species, "species"] #get species name
  mean_tss <- working_df[each_species, "threshold_tss"] #select threshold
  
  message("Processing: ", species_name)
  
  #read in model output for present
  filename <- paste0(
    species_name, "_present_ensemble_predictions.csv"
  )
  pres_df <- read.csv(file.path(model_dir, filename))
  
  # Loop through scenarios
  for (scenario in 1:length(future_scenarios)) {
    
    #scenario <- 1 #TESTING ONLY
    
    current <- future_scenarios[scenario] #scenario name
    
    #read in model output for future
    future_model_path <- file.path(model_dir, 
                                   paste0(species_name, "_", current, "_ensemble_predictions.csv"))
    fut_df <- read.csv(future_model_path)
    
    #join dfs and calculate changes
    new_df <- pres_df %>%
      select(subcID, pres_prob = pres_prob) %>%
      inner_join(
        fut_df %>% select(subcID, fut_prob = pres_prob),
        by = "subcID"
      ) %>%
      mutate(change = fut_prob - pres_prob)
    
    #save changes
    new_filename <- paste0(species_name, "_", current, "_change.csv")
    write.csv(new_df, file.path(model_dir, new_filename), row.names = FALSE)
    
    message("Creating raster: ", species_name, " ", current)
    
    #reclass rasters
    reclassTB <- new_df %>% select(subcID, prob = change)
    reclassTB <- reclassTB %>% mutate(new = as.integer(prob*100))
    
    #reclassTB %>% count(new)
    
    reclass_raster(
      data = reclassTB,
      rast_val = "subcID",
      new_val = "new",
      raster_layer = file.path(data_dir, "raw", "subcatchments.tif"),
      recl_layer = file.path(change_dir, paste0(
        species_name, "_", current, "_change.tif"
      )),
      bigtiff = TRUE
    )
    
  }}

# ==============================================================
# 06. RECLASSIFY RASTER TIFS
# ==============================================================

# ----------------------------
# 6a. Individual Plots
# ----------------------------

#create directories for input and output
model_dir <- model_dir2 #input dir
suitability_dir <- file.path(data_dir, "intermediate", "suitability_tifs_2") #output dir
dir.create(suitability_dir)

#list of scenarios for looping and naming outputs
scenario_names <- c("ssp1", "ssp3", "ssp5")

#go through each species and reclassify raster by pres/abs
for (each_species in 1:nrow(working_df)){
  
  species_name <- working_df[each_species, "species"] #select species
  mean_tss <- working_df[each_species, "threshold_tss"] #select threshold
  
  message("Processing: ", species_name)
  
  for (scenario in 1:length(scenario_names)) {
    
    current <- scenario_names[scenario] #select scenario
    
    filename <- paste0(
      species_name, "_", current, "_ensemble_predictions.csv"
    )
    pred_df <- read.csv(file.path(model_dir, filename)) #read in file for correct sp and scenario
    
    # Create raster
    reclassTB <- pred_df %>% select(subcID, prob = pres_prob)
    reclassTB <- reclassTB %>% mutate(new = as.integer(ifelse(prob > mean_tss, 1, 0)))
    
    #reclassTB %>% count(new)
    
    reclass_raster(
      data = reclassTB,
      rast_val = "subcID",
      new_val = "new",
      raster_layer = file.path(data_dir, "raw", "subcatchments.tif"),
      recl_layer = file.path(suitability_dir, paste0(
        species_name, "_", current, "_predictions.tif"
      )),
      bigtiff = TRUE
    )
  }
}

# ----------------------------
# 6b. Overall Plot
# ----------------------------

#list of scenarios for looping and naming outputs
scenario_names <- c("ssp1", "ssp3", "ssp5")
change_dir <- file.path(data_dir, "intermediate", "layer_dir") #directory #output dir

#for present only
# model_dir <- file.path(data_dir, "intermediate","model_tables") #input dir
# current <- "present"

#loop over scenarios
for (scenario in seq_along(scenario_names)) {
  
  current <- scenario_names[scenario]
  message("Processing: ", current, " scenario")
  
  long_list <- vector("list", nrow(working_df))
  
  for (i in seq_len(nrow(working_df))) {
    
    species_name <- working_df[i, "species"]
    mean_tss <- working_df[i, "threshold_tss"]
    
    message("Processing: ", species_name)
    
    filename <- paste0(species_name, "_", current, "_ensemble_predictions.csv")
    pred_df <- read.csv(file.path(model_dir, filename))
    
    long_list[[i]] <- pred_df %>%
      transmute(
        subcID,
        presence = as.integer(pres_prob >= mean_tss)
      )
  }
  
  # Combine all species
  long_df <- bind_rows(long_list)
  
  # Sum presences per subcatchment
  summary_df <- long_df %>%
    group_by(subcID) %>%
    summarise(sum = sum(presence), .groups = "drop")
  
  #save file
  write.csv(summary_df, file.path(change_dir, paste0(current, "_sum.csv")) )
  
  # Sanity check
  # message("Range: ", paste(range(summary_df$sum), collapse = " - "))
  # 
  # # Raster
  # reclass_raster(
  #   data = summary_df,
  #   rast_val = "subcID",
  #   new_val = "sum",
  #   raster_layer = file.path(data_dir, "subcatchments.tif"),
  #   recl_layer = file.path(
  #     suitability_dir, 
  #     paste0(current, "_layered_predictions.tif")
  #   ),
  #   bigtiff = TRUE
  # )
}

min(summary_df$sum)
max(summary_df$sum)
sd(summary_df$sum)

#write loop to sort through scenarios and make files with present-future

present_df <- read.csv(file.path(change_dir, "present_sum.csv"))

for (scenario in seq_along(scenario_names)) {
  
  current <- scenario_names[scenario]
  message("Processing: ", current, " scenario")
  
  filename <- paste0(current, "_sum.csv")
  sum_df <- read.csv(file.path(change_dir, filename))
  
  # Join by subcatchment ID
  change_df <- present_df %>%
    select(subcID, present_sum = sum) %>%
    left_join(
      sum_df %>% select(subcID, future_sum = sum),
      by = "subcID"
    ) %>%
    mutate(change = future_sum - present_sum)
  
  reclass_raster(
    data = change_df,
    rast_val = "subcID",
    new_val = "change",
    raster_layer = file.path(data_dir, "raw","subcatchments.tif"),
    recl_layer = file.path(
      change_dir,
      paste0(current, "_layered_predictions.tif")
    ),
    bigtiff = TRUE
  )
  
  
}

# ==============================================================
# 06. CALCULATE SENSITIVITY AND EXPOSURE
# ==============================================================

#read in dfs
global_predict <- read.csv(file.path(predict_dir,"present_predict_gbif.csv"))
global_df <- read.csv(file.path(occ_dir, "global_df.csv"))

ssp1 <- read.csv(file.path(predict_dir, "ssp126_predict_table_2.csv"))
ssp3 <- read.csv(file.path(predict_dir, "ssp370_predict_table_2.csv"))
ssp5 <- read.csv(file.path(predict_dir, "ssp585_predict_table_2.csv"))

#list future scenarios
scenarios <- list(
  ssp1 = ssp1,
  ssp3 = ssp3,
  ssp5 = ssp5
)

scenario_names <- names(scenarios)
scenario_names1 <- c("ssp126", "ssp370", "ssp585")

#create df to store exposure and sensitivity values
exposure_df <- data.frame(
  species_name = character(),
  scenario = character(),
  variable = character(),
  mean_pres = numeric(),
  sd_pres = numeric(),
  mean_fut = numeric(),
  sd_fut = numeric(),
  exposure = numeric(),
  sensitivity = numeric(),
  stringsAsFactors = FALSE
)

#list climate variables
climate_vars <- c(
  "bio05",  # Max temperature of warmest month (BIO5)
  "bio06",  # Min temperature of coldest month (BIO6)
  "bio13",  # Precipitation of wettest month (BIO13)
  "bio14",  # Precipitation of driest month (BIO14)
  "bio15",  # Precipitation seasonality (BIO15)
  "bio17",  # Precipitation of driest quarter
  "bio18"   # Precipitation of warmest quarter
)

#correct naming structure for variable data
variables_list <- paste0(
  climate_vars,
  "_1981.2010_observed_mean"
)

#main loop to calculate exposure and sensitivity for each species
for (each_species in 1:nrow(working_df)){
  
  #each_species <- 1 #TESTING ONLY
  
  species_name <- working_df[each_species, "species"] #get species name
  mean_tss <- working_df[each_species, "threshold_tss"] #get threshold to set pres/abs
  
  message("Processing: ", species_name)
  
  # Select present data
  sp_data <- global_df %>% 
    dplyr::filter(species == species_name) #global occurrence data for species
  
  filename <- paste0(
    species_name, "_present_ensemble_predictions.csv"
  )
  pred_df <- read.csv(file.path(model_dir1, filename)) #read in model output for present
  
  #reclassify raster
  reclassTB <- pred_df %>% select(subcID, prob = pres_prob)
  reclassTB <- reclassTB %>% mutate(new = as.integer(ifelse(prob > mean_tss, 1, 0)))
  reclassTB %>% count(new)
  
  #get subcIDs classified as present
  subc_ids <- reclassTB %>% filter(new == 1) %>% select(subcID) 
  occ_pts_ids <- reclassTB %>%
    filter(new == 1) %>%
    pull(subcID)
  
  #get environmental data for present subcatchments in the present
  df_pres <- global_predict %>%
    filter(subcID %in% occ_pts_ids) 

  # Loop through scenarios
  for (scenario in 1:length(scenario_names)) {
    
    future_predict <- scenarios[[scenario]] #select correct environmental predictor table for scenario
    current <- scenario_names[scenario] #scenario name
    
    df_fut <- future_predict %>%
      filter(subcID %in% occ_pts_ids) #get environmental data for the same subcatchments in the future
    
    # Calculate sensitivity (different for each scenario)
    future_model_path <- file.path(model_dir2, 
                                   paste0(species_name, "_", current, "_ensemble_predictions.csv"))
    
    df_present_model <- pred_df #present model predictions
    df_future_model  <- read.csv(future_model_path) #future model predictions
    
    sensitivity_val <- 
      (sum(df_present_model$pres_prob, na.rm = TRUE) - 
         sum(df_future_model$pres_prob, na.rm = TRUE)) /
      sum(df_present_model$pres_prob, na.rm = TRUE)
    
    # Loop through climate variables to calculate exposure (different for scenario/variable combination)
    for (var in 1:length(variables_list)) {
      
      var_name <- variables_list[var] #select variable name
      vars_future <- paste0(
        climate_vars,
        "_2071.2100_mpi.esm1.2.hr_", scenario_names1[scenario], "_v2_1_mean"
      ) #get variable column name
      var_fut <- vars_future[var]
      
      mean_pres <- mean(df_pres[[var_name]], na.rm = TRUE) #present mean of var
      sd_pres   <- sd(df_pres[[var_name]], na.rm = TRUE) #present sd of var
      mean_fut  <- mean(df_fut[[var_fut]], na.rm = TRUE) #future mean of var
      
      exposure_val <- (mean_fut - mean_pres)/(2*sd_pres) #calculate exposure
      
      #update df
      exposure_df <- rbind(
        exposure_df,
        data.frame(
          species_name = species_name,
          scenario = current,
          variable = climate_vars[var],
          mean_pres = mean_pres,
          sd_pres = sd_pres,
          mean_fut = mean_fut,
          exposure = exposure_val,
          sensitivity = sensitivity_val,
          stringsAsFactors = FALSE
        ))
      
      message("Completed variable ", climate_vars[var], " for species ", species_name, " in scenario ", current)
    }
  }
}
 
write.csv(exposure_df, file.path(results_dir, "exposure_df2.csv")) #save results

# ==============================================================
# 07. EXPOSURE ANALYSIS
# ==============================================================

#read in resutls
exposure_df <- read.csv(file.path(results_dir, "exposure_df.csv"))

#should only be these columns but lets make sure
exposure_df <- exposure_df %>% select(species_name, scenario, variable, mean_pres, sd_pres, mean_fut, exposure, sensitivity)

#check overall
exposure_df %>% group_by(variable) %>% 
  summarise(mean_exposure = mean(exposure), mean_sensitivity = mean(sensitivity))

#for one species
exposure_df %>% filter(species_name == "Xiphocentron_cubanum")

#add columns to easily see what species are classed as vulnerable (multiple over 1SD, or any over 2SD)
exposure_df <- exposure_df %>%
  mutate(
    moderate = abs(exposure) > 0.5,  # > 1 SD
    extreme  = abs(exposure) > 1     # > 2 SD
  )

#count number of extreme values per set of scenario and species
vulnerability_summary <- exposure_df %>%
  group_by(species_name, scenario) %>%
  summarise(
    n_moderate = sum(moderate, na.rm = TRUE),
    n_extreme  = sum(extreme, na.rm = TRUE),
    .groups = "drop"
  )

#create overall vulnerability column 
vulnerability_summary <- vulnerability_summary %>%
  mutate(
    check = case_when(
      n_extreme >= 1 ~ 1,         # 1 variable greater than 2SD
      n_moderate >= 2 ~ 1,        # 2 variables greater than 1SD
      TRUE ~ 0
    )
  )

vulnerability_summary %>% filter(check == 1)

#add to original df
exposure_df <- exposure_df %>%
  left_join(
    vulnerability_summary %>% select(species_name, scenario, check),
    by = c("species_name", "scenario")
  )

#sensitivity
exposure_df <- exposure_df %>% mutate(check_s = as.integer(ifelse(sensitivity >= 0, 1, 0)))
exposure_df %>% count(check_s)
exposure_df %>% filter(check_s == 1)

#save new version
write.csv(exposure_df, file.path(results_dir, "exposure_df2.csv")) 

# ==============================================================
# 08. CALCULATE DISPERSAL DISTANCE
# ==============================================================

#set input and output directories
dist_dir <- "/home/mchapman/masters_project/results/temp_dir" #input
output_dir <- "/home/mchapman/masters_project/results/distance_dir" #output
dir.create(output_dir, showWarnings = FALSE)

#filepath to subcatchment tif
subcatch_path <- "/home/mchapman/masters_project/data/raw/subcatchments.tif"

subcatch <- rast(subcatch_path)
names(subcatch) <- "subcID"

# -----------------------------
# Function to extract distances
# -----------------------------
extract_distance <- function(dist_raster_path, subcatch) {
  
  dist_r <- rast(dist_raster_path)
  names(dist_r) <- "distance"
  
  combined <- c(subcatch, dist_r)
  
  df <- as.data.frame(combined, na.rm = TRUE)
  
  df_aggr <- df %>%
    group_by(subcID) %>%
    summarise(distance = mean(distance), .groups = "drop")
  
  return(df_aggr)
}

# -----------------------------
# Get all raster files
# -----------------------------
files <- list.files(dist_dir, pattern = "\\.tif$", full.names = TRUE)

# -----------------------------
# Extract species names correctly
# aggr_present_species.tif
# -----------------------------
species_names <- unique(
  str_remove(basename(files), "^aggr_(present|ssp1|ssp3|ssp5)_") %>%
    str_remove("\\.tif$")
)

# -----------------------------
# Process each species
# -----------------------------
process_species <- function(species) {
  
  cat("Processing:", species, "\n")
  
  scenarios <- c("present", "ssp1", "ssp3", "ssp5")
  
  results <- list()
  
  for (sc in scenarios) {
    
    # Construct exact filename
    f <- file.path(dist_dir, paste0("aggr_", sc, "_", species, ".tif"))
    
    if (!file.exists(f)) {
      cat("  Missing:", sc, "\n")
      next
    }
    
    df <- extract_distance(f, subcatch)
    
    colnames(df)[2] <- paste0("dist_", sc)
    
    results[[sc]] <- df
  }
  
  # Join all scenarios
  if (length(results) == 0) return(NULL)
  
  df_final <- reduce(results, full_join, by = "subcID")
  
  df_final$species <- species
  
  # -----------------------------
  # Save per-species dataframe
  # -----------------------------
  out_file <- file.path(output_dir, paste0("dist_", species, ".csv"))
  write_csv(df_final, out_file)
  
  return(df_final)
}

# -----------------------------
# Run for all species
# -----------------------------
#process_species("Alisotrichia_alayoana")

all_results <- map(species_names, process_species)

# -----------------------------
# Test for significant differences
# -----------------------------

test_species <- function(df) {
  
  species_name <- unique(df$species)
  
  # Remove rows with NA in any scenario
  df_clean <- df %>%
    select(subcID, dist_present, dist_ssp1, dist_ssp3, dist_ssp5) %>%
    drop_na()
  
  # If too few observations, skip
  if (nrow(df_clean) < 5) {
    return(NULL)
  }
  
  # Run paired Wilcoxon tests
  tests <- list(
    ssp1 = wilcox.test(df_clean$dist_present, df_clean$dist_ssp1, paired = TRUE),
    ssp3 = wilcox.test(df_clean$dist_present, df_clean$dist_ssp3, paired = TRUE),
    ssp5 = wilcox.test(df_clean$dist_present, df_clean$dist_ssp5, paired = TRUE)
  )
  
  # Extract results
  results <- map_dfr(names(tests), function(sc) {
    tibble(
      species = species_name,
      scenario = sc,
      p_value = tests[[sc]]$p.value,
      statistic = tests[[sc]]$statistic,
      median_present = median(df_clean$dist_present),
      median_future = median(df_clean[[paste0("dist_", sc)]])
    )
  })
  
  return(results)
}

files <- list.files(output_dir, pattern = "^dist_.*\\.csv$", full.names = TRUE)

distance_df <- map_dfr(files, function(f) {
  
  # Read species dataframe
  df <- read_csv(f, show_col_types = FALSE)
  
  # Ensure species column exists
  if (!"species" %in% colnames(df)) {
    species_name <- str_remove(basename(f), "^dist_") %>%
      str_remove("\\.csv$")
    df$species <- species_name
  }
  
  # Run your existing test function
  test_species(df)
})

# Adjust p-values
distance_df <- distance_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

#check analysis
#distance_df <- read.csv(file.path(results_dir, "distance_df.csv"))
distance_df <- distance_df %>% mutate(check_d = as.integer(ifelse(p_adj < 0.05, 1, 0)))
distance_df %>% count(check_d)
distance_df %>% filter(check_d == 0)

#save df
write.csv(distance_df, file.path(results_dir, "distance_df.csv"))

# ==============================================================
# 09. PARTIAL DEPENDENCY ANALYSIS
# ==============================================================

#set directories
pdp_dir <- file.path(results_dir, "partial_dependence")

#create df to save results
pdp_df <- data.frame(
  species_name = character(),
  var1 = character(),
  var2 = character(),
  var3 = character(),
  var4 = character(),
  var5 = character(),
  stringsAsFactors = FALSE
)

#loop through species to record most important modelling variables
for (each_species in 1:nrow(working_df)) {
  
  message("Loading species ", each_species,
          " of 96")
  
  species_name <- working_df[each_species, "species"]
  
  # Select files in pdp_dir for this species
  species_files <- list.files(
    pdp_dir,
    pattern = paste0("^", species_name, "_pdp_.*\\.csv$"),
    full.names = TRUE
  )
  
  if (length(species_files) == 0) next
  
  # Order files by time created to get variable rankings
  file_info <- file.info(species_files)
  species_files <- species_files[order(file_info$mtime)]
  
  # Extract variables from filenames
  vars <- gsub(paste0(species_name, "_pdp_"), "", basename(species_files))
  vars <- gsub(".csv", "", vars)
  vars <- vars[1:min(5, length(vars))]
  vars <- c(vars, rep(NA, 5 - length(vars)))
  
  # Extract variables and fill in df
  pdp_df <- rbind(pdp_df, data.frame(
    species_name = species_name,
    var1 = vars[1],
    var2 = vars[2],
    var3 = vars[3],
    var4 = vars[4],
    var5 = vars[5],
    stringsAsFactors = FALSE
  ))
}

#convert shorthand to meaningful names
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
  "bio17" = "Precipitation of Driest Quarter",
  "c100" = "Tree and Shrub cover"
)

pdp_df <- pdp_df %>%
  mutate(across(
    starts_with("var"),
    ~ {
      code <- str_extract(.x, "bio[0-9]{2}|c100")
      ifelse(
        code %in% names(bio_lookup),
        bio_lookup[code],
        .x
      )
    }
  ))

#look at results
pdp_df %>% count(var1, sort = TRUE)

#convert format for easier visualisation
pdp_long <- pdp_df %>%
  pivot_longer(
    cols = starts_with("var"),
    names_to = "position",
    values_to = "variable"
  ) %>%
  mutate(
    position = as.numeric(gsub("var", "", position))
  )

write.csv(pdp_long, file.path(results_dir, "pdp_df.csv")) #save df

# ==============================================================
# 10. UPDATE WORKING_DF
# ==============================================================

#read in files
exposure_df <- read.csv(file.path(results_dir, "exposure_df2.csv"))
working_df <- read.csv(file.path(results_dir,  "working_df2.csv"))
distance_df <- read.csv(file.path(results_dir, "distance_df2.csv"))

#edit dfs for combination
distance_df <- distance_df %>% select(species, scenario, median_present, median_future, p_adj, check_d)
exposure_df <- exposure_df %>% mutate(species = species_name) %>% select(-species_name)

#combine 2
df <- inner_join(exposure_df, working_df, by = "species")
df <- df %>% select(family, order, class, species, scenario, variable, exposure, sensitivity, 
                    total_occurrences, mean_auc, mean_tss, cubanendemic, check, check_s)

#add 3rd
df1 <- inner_join(distance_df, df , by = c("species", "scenario"))
head(df1)

df1 %>% filter(species == "Xiphocentron_cubanum")

#filter for relevant rows
df_summary <- df1 %>%
  group_by(species, scenario) %>%
  summarise(
    exposure   = max(check, na.rm = TRUE),
    sensitivity = max(check_s, na.rm = TRUE),
    dispersal  = max(check_d, na.rm = TRUE),
    .groups = "drop"
  )

#which species are vulnerable or highly vulnerable?
venn_data <- df_summary %>%
  mutate(
    combo = case_when(
      exposure == 1 & sensitivity == 1 & dispersal == 1 ~ "all_three",
      exposure == 1 & sensitivity == 1 ~ "exp_sens",
      exposure == 1 & dispersal == 1 ~ "exp_disp",
      sensitivity == 1 & dispersal == 1 ~ "sens_disp",
      exposure == 1 ~ "exposure_only",
      sensitivity == 1 ~ "sensitivity_only",
      dispersal == 1 ~ "dispersal_only",
      TRUE ~ "none"
    ))

sp_vulnerable <- venn_data %>% filter(combo == "all_three"|combo == "exp_sens")

sp_vulnerable <- sp_vulnerable %>%
  left_join(
    working_df %>%
      select(species, order),
    by = "species"
  )

#get numbers
venn_counts <- venn_data %>%
  count(scenario, combo)

# ==============================================================
# 11. CHANGE IN HABITAT SUITABILITY
# ==============================================================

#calculate changes in habitat suitability between present and future

model_dir <- file.path(data_dir, "intermediate", "model_tables") #directory
working_df <- read.csv(file.path(results_dir, "working_df.csv")) #key df
head(working_df)

future_scenarios <- c("ssp1", "ssp3", "ssp5")

for (each_species in 1:nrow(working_df)){
  
  #each_species <- 1 #TESTING ONLY
  
  species_name <- working_df[each_species, "species"] #get species name
  
  message("Processing: ", species_name)
  
  #read in model output for present
  filename <- paste0(
    species_name, "_present_ensemble_predictions.csv"
  )
  pres_df <- read.csv(file.path(model_dir, filename))
  
  # Loop through scenarios
  for (scenario in 1:length(future_scenarios)) {
    
    #scenario <- 1 #TESTING ONLY
  
    current <- future_scenarios[scenario] #scenario name
    
    #read in model output for future
    future_model_path <- file.path(model_dir, 
                                   paste0(species_name, "_", current, "_ensemble_predictions.csv"))
    fut_df <- read.csv(future_model_path)
    
    #join dfs and calculate changes
    new_df <- pres_df %>%
      select(subcID, pres_prob = pres_prob) %>%
      inner_join(
        fut_df %>% select(subcID, fut_prob = pres_prob),
        by = "subcID"
      ) %>%
      mutate(change = fut_prob - pres_prob)
    
    #save changes
    new_filename <- paste0(species_name, "_", current, "_change.csv")
    write.csv(new_df, file.path(model_dir, new_filename), row.names = FALSE)
    
  }}

# ==============================================================
# 12. CALCULATE % RANGE CHANGE
# ==============================================================

#filter df to include only cuba subcatchments
subc_area <- read.csv(file.path(data_dir, "raw", "subcatchment_areas_sqm.csv"))
ssp1 <- read.csv(file.path(predict_dir, "ssp126_predict_table_2.csv"))

subc_area_filtered <- subc_area %>%
  semi_join(ssp1, by = c("subc_id" = "subcID")) %>%
  select(subc_id, area_sqm)

nrow(subc_area_filtered) == nrow(ssp1) #check

#save
write.csv(subc_area_filtered, file.path(data_dir, "raw", "subcatchment_areas_sqm_cuba.csv"))

#read in without repeating
subc_area <- read.csv(file.path(data_dir, "raw", "subcatchment_areas_sqm_cuba.csv"))

#create present/absence conversions

#input and output dirs
present_model_dir <- model_dir1 #input present
future_model_dir <- model_dir2 #input future
pres_abs <- file.path(data_dir, "intermediate", "pres_abs")
dir.create(pres_abs)

future_scenarios <- c("ssp1", "ssp3", "ssp5")

#create df to store range increase and decrease values
range_df <- data.frame(
  species_name = character(),
  scenario = character(),
  present_range = numeric(),
  future_range = numeric(),
  range_change = numeric()
)


for (each_species in 1:nrow(working_df)){
  
  #each_species <- 1 #TESTING ONLY
  
  species_name <- working_df[each_species, "species"] #get species name
  mean_tss <- working_df[each_species, "threshold_tss"] #select threshold
  
  message("Processing: ", species_name)
  
  #read in model output for present
  filename <- paste0(
    species_name, "_present_ensemble_predictions.csv"
  )
  pres_df <- read.csv(file.path(present_model_dir, filename))
  
  #create present and save
  reclass_pres <- pres_df %>% select(subcID, prob = pres_prob)
  reclass_pres <- reclass_pres %>% mutate(pres_abs = as.integer(ifelse(prob > mean_tss, 1, 0)))
  
  pres_filename <- paste0(species_name, "_present_pres_abs.csv")
  write.csv(reclass_pres, file.path(pres_abs, pres_filename))
  
  #calculate present_range
  present_range <- reclass_pres %>%
    left_join(subc_area,
              by = c("subcID" = "subc_id")) %>%
    filter(pres_abs == 1) %>%
    summarise(total_area = sum(area_sqm, na.rm = TRUE)) %>%
    pull(total_area)
  
  # Loop through scenarios
  for (scenario in 1:length(future_scenarios)) {
    
    #scenario <- 1 #TESTING ONLY
    
    current <- future_scenarios[scenario] #scenario name
    
    #read in model output for future
    future_model_path <- file.path(future_model_dir, 
                                   paste0(species_name, "_", current, "_ensemble_predictions.csv"))
    fut_df <- read.csv(future_model_path)
    
    #create future df and save
    reclass_fut <- fut_df %>% select(subcID, prob = pres_prob)
    reclass_fut <- reclass_fut %>% mutate(pres_abs = as.integer(ifelse(prob > mean_tss, 1, 0)))
    
    fut_filename <- paste0(species_name, "_", current, "_pres_abs.csv")
    write.csv(reclass_fut, file.path(pres_abs, fut_filename))
    
    #calculate future range
    future_range <- reclass_fut %>%
      left_join(subc_area,
                by = c("subcID" = "subc_id")) %>%
      filter(pres_abs == 1) %>%
      summarise(total_area = sum(area_sqm, na.rm = TRUE)) %>%
      pull(total_area)
    
    message("Successfully created pres/abs: ", species_name, " ", current)
    
    #calculate % range change
    range_change <- (
      (future_range - present_range) / present_range
    ) * 100
    
    #update diagram
    range_df <- rbind(
      range_df,
      data.frame(
        species_name = species_name,
        scenario = scenario,
        present_range = present_range,
        future_range = future_range,
        range_change = range_change
      )
    )
    
    message(
      "Completed: ",
      species_name,
      " ",
      scenario
    )
    
  }
  
}

write.csv(range_df, file.path(results_dir, "range_df.csv"))
