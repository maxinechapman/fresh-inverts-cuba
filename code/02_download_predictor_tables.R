# ==============================================================
# 02 Master's Project: Download Environmental Predictor Data
# Author: Maxine Chapman
# Date: 10/03/2026
# ==============================================================

# ==============================================================
# 1. SETUP
# ==============================================================

## 1.1 Load Libraries
library(hydrographr)
library(sf)
library(data.table)
library(dplyr)
library(terra)
library(tools)
library(stringr)
library(ranger)
library(kableExtra)
library(leaflet)
library(leafem)
library(htmlwidgets)
library(corrplot)
library(usdm)
library(data.table)

## 1.2 Define Directory Structure
wdir <- "/home/mchapman/masters_project" #server
setwd(wdir)

data_dir <- file.path(wdir, "data")
predict_dir <- file.path(data_dir, "raw", "predict_inputs") # create predict folder
dir.create(predict_dir, showWarnings = FALSE)

# ==============================================================
# 2. SELECTING TILES
# ==============================================================

## 2.1 Load Occurrence Data
file_raw <- file.path(wdir, "gbif_data.csv")
# file_raw <- file.path(wdir, "occurrence/gbif_data/gbif_data.csv") # alternative path
gbif <- fread(file_raw, quote = "\"", fill = Inf)
head(gbif)

## 2.2 Define Study Area (Bounding Box)
bbox <- c(-84.9749110583, 19.8554808619, -74.1780248685, 23.1886107447)

## 2.3 Identify Tiles for Download
# tile_id <- c("h08v06", "h10v06") # manual selection
tile_id <- get_tile_id(data = gbif, lon = "decimalLongitude", lat = "decimalLatitude")
tile_list <- split(tile_id, ceiling(seq_along(tile_id) / 5))

# ==============================================================
# 3. DOWNLOAD SUBCATCHMENT TILES
# ==============================================================

download_tiles(variable = c("sub_catchment"), tile_id = "h26v02", file_format = "tif",
               download_dir = wdir) # raster files

# Optional manual download
# tile_urls <- paste0(...)
# tile_files <- file.path(...)
# for(i in seq_along(tile_files)) { ... }

# ==============================================================
# 4. EXTRACT SUBCATCHMENT IDS
# ==============================================================

subc_ids <- extract_ids(data = gbif,
                        lon = "decimalLongitude",
                        lat = "decimalLatitude",
                        subc_layer = file.path(wdir, "subcatchments_gbif(1).tif")) %>%
  distinct(subcatchment_id) %>%
  pull()
head(subc_ids)
dim(subc_ids)

subc_ids_dt <- data.table(subcatchment_id = subc_ids)
fwrite(subc_ids_dt, file = file.path(data_dir, "subc_ids1.txt"))

#cuba subcatchment IDs
cuba_subc_ids <- extract_ids(subc_layer = file.path(data_dir, "raw", "subcatchments.tif"))
fwrite(cuba_subc_ids, file = file.path(data_dir, "cuba_subc_ids.txt"))

# ==============================================================
# 5. PREPARE TEMPORARY DIRECTORIES
# ==============================================================

Sys.setenv(TMPDIR = file.path(Sys.getenv("HOME"), "tmp"))
dir.create(Sys.getenv("TMPDIR"), showWarnings = FALSE, recursive = TRUE)

hydro_tmp <- file.path(wdir, "hydro_tmp")
dir.create(hydro_tmp, showWarnings = FALSE, recursive = TRUE)

# Ensure predict_dir exists and is empty
dir.create(predict_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(predict_dir, full.names = TRUE, recursive = TRUE),
       recursive = TRUE, force = TRUE)

# ==============================================================
# 6. DEFINE VARIABLES
# ==============================================================

## 6.1 Climate Variables (Observed)
climate_vars <- c(
  "bio01",  # Annual mean temperature (BIO1)
  "bio05",  # Max temperature of warmest month (BIO5)
  "bio06",  # Min temperature of coldest month (BIO6)
  "bio07",  # Temperature annual range (BIO7 = BIO5 - BIO6)
  "bio12",  # Annual precipitation (BIO12)
  "bio13",  # Precipitation of wettest month (BIO13)
  "bio14",  # Precipitation of driest month (BIO14)
  "bio15",  # Precipitation seasonality (BIO15)
  "bio17",  # Precipitation of driest quarter
  "bio18"   # Precipitation of warmest quarter (BIO18)
)

climate_observed <- paste0(climate_vars, "_1981-2010_observed")

## 6.2 Land Cover Variables
land_vars <- c(
  "c50", "c60", "c100", "c120",
  "c130", "c150", "c190", "c210"
)
land_vars_names <- paste0(land_vars, "_2020")

## 6.3 Topographical Variables
topo_vars <- c("accumulation", "gradient", "spi", "cti", "out_dist")

# ==============================================================
# 7. DOWNLOAD DATA FOR EACH TILE SECTION
# ==============================================================

tile_log <- data.frame(
  section = integer(),
  tile = character(),
  variable = character(),
  status = character(),
  stringsAsFactors = FALSE
)

for (section in seq_along(tile_list)) {
  
  tile_id <- tile_list[[section]]
  message("Processing tile group ", section, " (", length(tile_id), " tiles)")
  
  ## 7.1 Download Observed Climate Tables
  download_observed_climate_tables(
    subset = climate_observed,
    tile_ids = tile_id,
    download = TRUE,
    download_dir = predict_dir,
    file_format = "txt",
    delete_zips = TRUE,
    ignore_missing = TRUE,
    tempdir = Sys.getenv("TMPDIR"),
    quiet = FALSE
  )
  
  ## 7.2 Download Landcover Tables
  download_landcover_tables(
    base_vars = land_vars,
    years = "2020",
    tile_ids = tile_id,
    download = TRUE,
    download_dir = predict_dir,
    file_format = "txt",
    delete_zips = TRUE,
    ignore_missing = TRUE,
    tempdir = tempdir(),
    quiet = FALSE
  )
  
  ## 7.3 Download Hydrography Tables
  download_hydrography90m_tables(
    subset = topo_vars,
    tile_ids = tile_id,
    download = TRUE,
    download_dir = predict_dir,
    file_format = "txt",
    delete_zips = TRUE,
    ignore_missing = TRUE,
    tempdir = tempdir(),
    quiet = FALSE
  )
  
  ## 7.4 Create Global Prediction Table
  predict_vars <- c(clim_vars, land_vars_names, topo_vars)
  predict_table_name <- paste0("global_predict_section", section, ".csv")
  out_file_path <- file.path(wdir, predict_table_name)
  
  get_predict_table(
    variable = predict_vars,
    statistics = "mean",
    tile_id = tile_id,
    input_var_path = predict_dir,
    subcatch_id = file.path(wdir, "subc_ids.txt"),
    out_file_path = out_file_path,
    read = FALSE,
    quiet = FALSE,
    overwrite = TRUE,
    n_cores = 13
  )
  
  ## 7.5 Update Tile Log
  files <- list.files(predict_dir, pattern = "\\.txt$", recursive = TRUE, full.names = FALSE)
  for (tile in tile_id) {
    for (var in predict_vars) {
      fname <- files[grepl(var, files) & grepl(tile, files)]
      status <- ifelse(length(fname) > 0, "available", "missing")
      tile_log <- rbind(tile_log, data.frame(section = section, tile = tile, variable = var, status = status, stringsAsFactors = FALSE))
    }
  }
  
  ## 7.6 Clean Up Temporary Files
  unlink(list.files(predict_dir, full.names = TRUE, recursive = TRUE), recursive = TRUE, force = TRUE)
  unlink(list.files(hydro_tmp, full.names = TRUE, recursive = TRUE), recursive = TRUE, force = TRUE)
}

# ==============================================================
# 7. DOWNLOAD FUTURE DATA FOR CUBA
# ==============================================================

tile_id <- c("h08v06", "h10v06") #select tiles for cuba
ouput_path <- file.path(wdir, "data", "raw", "predict_tables") #output path
scenario_names <- c("ssp126", "ssp370", "ssp585")

#for loop by scenario
for (scenario in 1:length(scenario_names)) { 
  
  current <- scenario_names[scenario]
  climate_names <- paste0(
    climate_vars, "_2071-2100_mpi-esm1-2-hr_", current, "_v2_1"
  )
  table_name <- paste0(current, "_predict_table_2.csv")
  out_path <- file.path(ouput_path, table_name)
  
  # Check if file already exists
  if (file.exists(out_path)) {
    message("Skipping ", table_name, " (already exists)")
    next
  }
  
  download_projected_climate_tables(
    base_vars = c(climate_vars),
    time_periods = c("2071-2100"),
    models = c("mpi-esm1-2-hr"),
    scenarios = current,
    versions = c("v2_1"),
    subset = NULL,
    tile_ids = tile_id,
    download = TRUE,
    download_dir = predict_dir,
    file_format = "txt",
    delete_zips = TRUE,
    ignore_missing = FALSE,
    tempdir = NULL,
    quiet = FALSE
  )
  
  predict_vars <- c(climate_names, land_vars_names, topo_vars) 
  
  get_predict_table(
    variable = predict_vars,
    statistics = c("mean"),
    tile_id = tile_id,
    input_var_path = predict_dir,
    subcatch_id = file.path(data_dir, "cuba_subc_ids.txt"),
    out_file_path = out_path,
    read = FALSE,          
    quiet = FALSE,
    overwrite = TRUE,
    n_cores = 13
  )
  
}

# ==============================================================
# 8. SAVE TILE LOG
# ==============================================================
write.csv(tile_log, file.path(wdir, "tile_processing_log.csv"), row.names = FALSE)

# ==============================================================
# 9. MERGE GLOBAL PREDICT FILES
# ==============================================================
dfs <- list.files(
  path = wdir,
  pattern = "^global_predict_section[0-9]+\\.csv$",
  full.names = TRUE
)

predict_table <- rbindlist(lapply(dfs, fread), use.names = TRUE, fill = TRUE)

# ==============================================================
# 10. COLINEARITY CHECK AND FILTERING
# ==============================================================

# Drop rows with missing values (very few, at the map edges)
predict_table <- na.omit(predict_table)
colSums(is.na(predict_table))

# ------------------------------
# Select numeric predictor variables (mean values only)
# ------------------------------
numeric_vars <- predict_table %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(ends_with("_mean") | ends_with("2020"))

# ------------------------------
# Compute correlation matrix
# ------------------------------
cor_matrix <- cor(numeric_vars, use = "complete.obs")

# Visualize correlations
corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.6, diag = FALSE)

# ------------------------------
# Perform Variance Inflation Factor (VIF) analysis
# ------------------------------
vif_results <- vifstep(numeric_vars, th = 10)

# Display selected variables after removing collinear ones
vif_results@results
selected_vars <- vif_results@results$Variables
selected_vars

# ------------------------------
# Filter the prediction table to retain only non-collinear variables
# ------------------------------
predict_table <- predict_table %>%
  dplyr::select(any_of(c("subcID", selected_vars)))

## Remove Unnecessary Variables (Observed)
remove_vars_observed <- c(paste0(c("bio01","bio07","bio12"), "_1981.2010_observed_mean"), "gradient", "out_dist")
predict_table_observed <- predict_table[, !(names(predict_table) %in% remove_vars_observed)]

#Remove Unnecessary Variables (Future)

# Define scenarios and file names
scenarios <- c("ssp126", "ssp370", "ssp585")

for (sc in scenarios) {
  
  # Read file
  file_path <- file.path(predict_dir, paste0(sc, "_predict_table_2.csv"))
  df <- read.csv(file_path)
  
  # Define variables to remove (specific to scenario)
  remove_vars_future <- c(
    paste0(c("bio01","bio07","bio12"),
           "_2071.2100_mpi.esm1.2.hr_", sc, "_v2_1_mean"),
    "gradient",
    "out_dist"
  )
  
  # Remove columns
  df_clean <- df[, !(names(df) %in% remove_vars_future)]
  
  # Overwrite file
  write.csv(df_clean, file_path, row.names = FALSE)
}



# ==============================================================
# 11. SAVE PREDICTION TABLES
# ==============================================================

write.csv(predict_table_observed, "/home/mchapman/masters_project/data/present_predict_gbif.csv") #observed
write.csv(predict_table_future, file.path(wdir, paste0(scenario, "_predict_table.csv"))) #future

# ==============================================================
# 12. DELETE TEMPORARY FILES
# ==============================================================

unlink(list.files(predict_dir, full.names = TRUE, recursive = TRUE), recursive = TRUE, force = TRUE)
unlink(list.files(file.path(wdir, "calibration_areas"), full.names = TRUE, recursive = TRUE))
