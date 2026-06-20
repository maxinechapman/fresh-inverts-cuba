#Title: Master's Project modelling
#Author: Maxine Chapman
#Date: 11/03/26

# ============================================================
# 1. SETUP
# ============================================================

# ------------------------------
# 1.1 Libraries
# ------------------------------
# remotes::install_github("sjevelazco/flexsdm")
library(flexsdm)
library(dplyr)
library(terra)
library(hydrographr)
library(sf)
library(readr)
library(ranger)
library(pROC)
library(SDMtune)
library(leaflet)
library(data.table)
library(pdp)

# ------------------------------
# 1.2 Directories
# ------------------------------
# Define working directory (server or local environment)
wdir <- "/home/mchapman/masters_project" #server
# wdir <- "C:/Users/mchapman/Documents/masters_project" #local
setwd(wdir)

# Define main data and results directories
data_dir <- file.path(wdir, "data")
results_dir <- file.path(wdir, "results")
predict_dir <- file.path(data_dir, "raw", "predict_tables")
occurrence_dir <- file.path(data_dir, "raw", "occurrence_data")

# Create output directories if they do not already exist
pdp_dir <- file.path(results_dir, "partial_dependence")
dir.create(pdp_dir, showWarnings = FALSE)

model_dir <- file.path(data_dir, "intermediate","model_tables_2") #main directors
dir.create(model_dir, showWarnings = FALSE)


# ============================================================
# 2. KEY DATASETS
# ============================================================

# ------------------------------
# 2.1 Species occurrence datasets
# ------------------------------

#run code once and then just read in the df (line 88)

#read in key dataframes
gbif <- readr::read_csv(file.path(occurrence_dir, "gbif_data.csv"))
sp <- read.csv(file.path(occurrence_dir, "species_occurrences_cuba.csv")) #read in species data

#edit gbif df for compatibility
gbif$species <- gsub(" ", "_", gbif$species) #edit names for matching
gbif_subset <- gbif %>%
  select(any_of(keep_cols)) %>%
  mutate(cubanendemic = "n") # Subset GBIF to sp columns and add missing column

#get species list
species_list <- sort(unique(sp$species)) 
length(species_list) #number of species

#filter for invertebrates only (remove vertebrates)
sp_inverts <- sp %>%
  filter(!is.na(phylum) & phylum != "Chordata")
inverts_list <- sort(unique(sp_inverts$species))
length(inverts_list) #number of invertebrates

#create combined df
keep_cols <- colnames(sp_inverts)
sp_subset <- sp_inverts %>%
  select(all_of(keep_cols)) # Ensure sp has the same column order
combined_df <- bind_rows(sp_subset, gbif_subset) %>% group_by(species) # Combine datasets

#save combined dataset
write.csv(combined_df, file.path(data_dir, "global_df.csv")) 

#read in existing dataset after running code once
combined_df <- read.csv(file.path(data_dir, "raw", "occurrence_data", "global_df.csv"))
head(combined_df)


# ------------------------------
# 2.2 Model evaluation functions
# ------------------------------

#define function for mean tss to evaluate model
compute_tss <- function(obs, pred, thresholds = seq(0, 1, by = 0.01)) {
  
  tss_vals <- sapply(thresholds, function(th) {
    
    pred_bin <- ifelse(pred >= th, 1, 0)
    
    TP <- sum(pred_bin == 1 & obs == 1)
    TN <- sum(pred_bin == 0 & obs == 0)
    FP <- sum(pred_bin == 1 & obs == 0)
    FN <- sum(pred_bin == 0 & obs == 1)
    
    sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    
    sens + spec - 1
  })
  
  max(tss_vals, na.rm = TRUE)
}

#define function to calculate threshold tss to set presence/absence habitat suitability
threshold_tss <- function(obs, pred, thresholds = seq(0, 1, by = 0.01)) {
  
  tss_vals <- sapply(thresholds, function(th) {
    
    pred_bin <- ifelse(pred >= th, 1, 0)
    
    TP <- sum(pred_bin == 1 & obs == 1)
    TN <- sum(pred_bin == 0 & obs == 0)
    FP <- sum(pred_bin == 1 & obs == 0)
    FN <- sum(pred_bin == 0 & obs == 1)
    
    sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    
    sens + spec - 1
  })
  
  valid <- !is.na(tss_vals)
  
  if (!any(valid)) {
    return(list(max_tss = NA, threshold = NA))
  }
  
  max_tss <- max(tss_vals[valid])
  
  # Find first threshold giving that max value
  max_index <- which(tss_vals == max_tss)[1]
  
  list(
    threshold = thresholds[max_index]
  )
}


# ------------------------------
# 2.3 Environmental predictor datasets
# ------------------------------

#read in predict tables for different scenarios and regions
present_predict <- read.csv(file.path(predict_dir, "present_predict_table.csv"))
global_predict <- read.csv(file.path(predict_dir, "present_predict_gbif.csv"))

ssp1 <- read.csv(file.path(predict_dir, "ssp126_predict_table_2.csv"))
ssp3 <- read.csv(file.path(predict_dir, "ssp370_predict_table_2.csv"))
ssp5 <- read.csv(file.path(predict_dir, "ssp585_predict_table_2.csv"))

#create list of scenarios
scenarios <- list(
  ssp1 = ssp1,
  ssp3 = ssp3,
  ssp5 = ssp5
) 

#create list of scenario names for labelling
scenario_names <- names(scenarios)


# ------------------------------
# 2.4 Modelling progress and species lists
# ------------------------------

#create progress log to track results
progress_log <- data.frame(
  species_name = character(),
  species_index = integer(),
  total_occurrences = integer(),
  mean_auc = numeric(),
  mean_tss = numeric(),
  threshold_tss =numeric(),
  status = character(),
  stringsAsFactors = FALSE
)

#rerun only results with mean tss and auc above 0.7

# results_df <- read.csv(file.path(results_dir, "working_df.csv"))
# combined_df <- combined_df %>% filter(species %in% results_df$species) #create new versions of dfs

progress_log <- read.csv(file.path(wdir, "modelling_log.csv"))

species_list <- sort(unique(combined_df$species))
length(species_list)


# ============================================================
# 3. MODELLING LOOP
# ============================================================

# ------------------------------
# 3.1 Loop setup
# ------------------------------

for (each_species in 240:length(species_list)) {
  #for (each_species in 1:5) { # TESTING ONLY
  #each_species <- 240 #TESTING ONLY
  
  species_name <- species_list[each_species]
  #species_name
  
  message("Processing species: ", species_name)

  # ---------------------------------
  # 1. Skip if all scenario predictions exist
  # ---------------------------------
  existing_files <- list.files(
    model_dir,
    pattern = paste0("^", species_name, ".*_ensemble_predictions\\.csv$"),
    full.names = TRUE
  )

  if (length(existing_files) == length(scenarios)) {
    message("Skipping ", species_name,
            " (all scenario predictions already exist)")
    next
  }

  # ---------------------------------
  # 2. Prepare presence & pseudoabsence data
  # ---------------------------------
  sp_data <- combined_df %>% 
    dplyr::filter(species == species_name) #filter species data
  
  n_pabs <- 2 * nrow(sp_data) #set double number of pseudoabsences as occurrences
  
  occ_pts_ids <- extract_ids(
    data = sp_data,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    subc_layer = file.path(data_dir, "raw", "subcatchments_gbif(1).tif")
  ) %>% 
    distinct(subcatchment_id) %>% 
    pull() #extract subcatchment ids for this species
  
  df_pseudoabs <- global_predict %>%
    filter(!subcID %in% occ_pts_ids) %>%
    slice_sample(n = n_pabs) %>%
    mutate(PresAbs = 0) #pabs in different subcatchments to occurrences and select environmental data
  
  df_pres <- global_predict %>%
    filter(subcID %in% occ_pts_ids) %>%
    mutate(PresAbs = 1) #select environmental data for presences
  
  # ---------------------------------
  # 3. Skip species with insufficient presence data
  # ---------------------------------
  if (nrow(df_pres) < 2) {
    message("Skipping species ", species_name,
            " (insufficient presence data)")
    
    progress_log <- rbind(
      progress_log,
      data.frame(
        species_name = species_name,
        species_index = each_species,
        total_occurrences = nrow(sp_data),
        mean_auc = NA,
        mean_tss = NA,
        threshold_tss = NA,
        status = "failure",
        stringsAsFactors = FALSE
      )
    ) #record skip in progress log
    
    write.csv(progress_log, file.path(wdir, "modelling_log.csv"), row.names = FALSE)
    next
  }
  
  # ---------------------------------
  # 4. Merge presences and pseudoabsences
  # ---------------------------------
  model_fit_table <- bind_rows(df_pres, df_pseudoabs)
  model_fit_table$PresAbs <- as.factor(model_fit_table$PresAbs)
  model_fit_table <- na.omit(model_fit_table)
  
  predictors <- colnames(model_fit_table)[3:11]
  
  # ---------------------------------
  # 5. Train 10 model runs
  # ---------------------------------
  models <- list() #list for models
  metrics <- data.frame() #list for AUC and TSS
  var_imp_list <- list()  #list of variable importances
  
  for (i in 1:10) {
   #i <- 1 #TESTING ONLY
    set.seed(42 + i) #ensures reproducible but randomised results
    
    #select presences or absences
    pres_idx <- which(model_fit_table$PresAbs == 1) 
    abs_idx  <- which(model_fit_table$PresAbs == 0)
    
    #split into training and testing sets
    train_pres <- sample(pres_idx, size = floor(0.8 * length(pres_idx)))
    train_abs  <- sample(abs_idx,  size = floor(0.8 * length(abs_idx)))
    
    train_idx <- c(train_pres, train_abs)
    data_train <- model_fit_table[train_idx, ]
    data_test  <- model_fit_table[-train_idx, ]
    
    #construct model
    model <- ranger(
      x = data_train[, predictors],
      y = data_train$PresAbs,
      probability = TRUE,
      replace = TRUE,
      importance = "impurity"
    )
    
    # Sort variables, decreasing importance. Top five predictors in this example
    var_imp <- model$variable.importance
    vars <- names(sort(var_imp, decreasing = TRUE))[1:5]
    
    # Compute partial dependence for the most important predictor using the pdp package
    pd <- partial(model, pred.var = vars[1], train = data_train, prob = TRUE, which.class = 2)
    
    # Store variable importance for ensemble calculation
    var_imp_list[[i]] <- var_imp
    
    #get predictions
    test_pred <- predict(model, data_test[, predictors])$predictions[, "1"]
    
    #get metrics
    auc <- pROC::auc(data_test$PresAbs, test_pred)
    tss_value <- compute_tss(
      obs  = as.numeric(as.character(data_test$PresAbs)),
      pred = test_pred
    )
    tss_threshold <- threshold_tss(
      obs  = as.numeric(as.character(data_test$PresAbs)),
      pred = test_pred
    )
    
    #record data for this iteration
    models[[i]] <- model
    
    metrics <- rbind(
      metrics,
      data.frame(
        run = i,
        AUC = as.numeric(auc),
        TSS = as.numeric(tss_value),
        TSS_threshold = as.numeric(tss_threshold)
      )
    )
    
    message("Model complete for: ", species_name, " run ", i)
  }
  
  #get means of all runs for key metrics
  mean_auc <- mean(metrics$AUC)
  mean_tss <- mean(metrics$TSS)
  mean_threshold <- mean(metrics$TSS_threshold)
  
  # ---------------------------------
  # 5b. Ensemble PDP calculation 
  # ---------------------------------
  mean_var_imp <- Reduce("+", var_imp_list) / length(var_imp_list) #means for variables
  
  var_imp_percent <- 100 * mean_var_imp / sum(mean_var_imp) #% contributions
  
  #create and save dataframe
  var_imp_df <- data.frame(
    variable = names(var_imp_percent),
    importance_percent = var_imp_percent
  )
  
  write.csv(
    var_imp_df,
    file = file.path(pdp_dir, paste0(species_name, "_variable_importance_percent.csv")),
    row.names = FALSE
  )
  
  # ---------------------------------
  # 6. Scenario predictions
  # ---------------------------------
  for (scenario in 1:length(scenarios)) {
    predict_table <- scenarios[[scenario]] #select scenario
    current <- scenario_names[scenario] #select scenario name

    # Rename climate variables to match training data
    bio_scenario_cols <- grep("^bio", colnames(predict_table), value = TRUE)
    
    # Map to the names used in training
    bio_mapping <- predictors[grepl("^bio", predictors)]
    if (length(bio_scenario_cols) != length(bio_mapping)) {
      stop("Number of bio columns in scenario does not match training predictors")
    }
    colnames(predict_table)[match(bio_scenario_cols, colnames(predict_table))] <- bio_mapping

    # Predict 10 models
    scenario_pred_list <- vector("list", 10)
    for (i in 1:10) {
      scenario_pred_list[[i]] <- predict(
        models[[i]],
        data = predict_table[, predictors]
      )$predictions[, "1"]
    }

    # Ensemble mean
    scenario_matrix <- do.call(cbind, scenario_pred_list)
    pred_mean <- rowMeans(scenario_matrix)

    # Write CSV
    pred_df <- predict_table %>%
      select(subcID) %>%
      mutate(
        pres_prob = pred_mean,
        prob_1_int = as.integer(round(pres_prob, 2) * 100)
      )

    filename <- paste0(
      species_name, "_", current, "_ensemble_predictions.csv"
    )
    write.csv(pred_df, file.path(model_dir, filename), row.names = FALSE)

    # Create raster
    reclassTB <- pred_df %>% select(subcID, prob = prob_1_int)

    reclass_raster(
      data = reclassTB,
      rast_val = "subcID",
      new_val = "prob",
      raster_layer = file.path(data_dir, "raw", "subcatchments.tif"),
      recl_layer = file.path(model_dir, paste0(
        species_name, "_", current, "_predictions.tif"
      )),
      bigtiff = TRUE
    )

    message("Raster saved for: ", species_name, " ", current)
  }

  # ---------------------------------
  # 7. Log progress
  # ---------------------------------
  progress_log <- rbind(
    progress_log,
    data.frame(
      species_name = species_name,
      species_index = each_species,
      total_occurrences = nrow(sp_data),
      mean_auc = mean_auc,
      mean_tss = mean_tss,
      threshold_tss = mean_threshold,
      status = "success",
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(progress_log, file.path(wdir, "modelling_log.csv"), row.names = FALSE)
}

#progress_log <- read.csv(file.path(wdir, "modelling_log.csv"))

progress_log %>% count(status) #check success

# status   n
# 1 failure 101
# 2 success 353

# filter for AUC and TSS

final_results <- progress_log %>% filter(status == "success") %>% 
                                  filter(mean_auc >= 0.7) %>%
                                  filter(mean_tss >= 0.7)

failed_results <- progress_log %>% filter(mean_auc < 0.7 | is.na(mean_auc) | mean_tss < 0.7 | is.na(mean_tss))


for (each_species in 1:length(failed_results$species_name)) {
  
  #each_species <- 1 #TESTING ONLY
  
  sp <- failed_results$species_name[each_species]
  
  # Loop over each scenario name
  for (scenario in scenario_names) {
    
    current <- scenario_names[scenario] #select scenario name
    
    # Construct scenario file path (in model_dir)
    scenario_file <- file.path(
      model_dir,
      paste0(sp, "_", current, "_ensemble_predictions.csv")
    )
    
    # Delete scenario file if it exists
    if (file.exists(scenario_file)) {
      file.remove(scenario_file)
      message("Deleted: ", scenario_file)
    } else {
      message(current, " scenario file does not exist for ", sp)
    }
  }
  
  # Construct PDP file path (one per species)
  pdp_file <- file.path(
    pdp_dir,
    paste0(sp, "_variable_importance_percent.csv")
  )
  
  # Delete PDP file if it exists
  if (file.exists(pdp_file)) {
    file.remove(pdp_file)
    message("Deleted: ", pdp_file)
  } else {
    message("Variable importance file not found for ", sp)
  }
}

write.csv(final_results, file.path(results_dir, "results.csv"))

