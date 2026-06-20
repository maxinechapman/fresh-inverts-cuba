#Title: 01 Master's Project: Download GBIF data
#Author: Maxine Chapman
#Date: 10/03/2026


# Download packages
if (!require(tidyverse)) install.packages("tidyverse", dependencies = T)
if (!require(sf)) install.packages("sf", dependencies = T)
if (!require(rgbif)) install.packages("rgbif", dependencies = T)
if (!require(CoordinateCleaner)) install.packages("CoordinateCleaner", dependencies = T)
if (!require(remotes)) install.packages("remotes", dependencies = T)
remotes::install_github("glowabio/hydrographr")
library(hydrographr)
library(dplyr)
library(taxize)


# 1 Set path to import or export data
wdir <- "/home/mchapman/masters_project" #server
setwd(wdir)

data_dir <- file.path(wdir, "data")

#folder to save results
output_occur_path <- paste0(getwd(), "data/raw/occurrence_data/gbif_data")
# if(!dir.exists(output_occur_path)) dir.create(output_occur_path)

if (!dir.exists(output_occur_path)) {
  dir.create(output_occur_path, recursive = TRUE)
}


#folder with hydrography90m layers
output_hydro90m <- paste0(getwd(), "/hydrography90m")


# preliminary species list Danube
# prelim_sp_danube <- read.csv("./occurrence/JDS5_fish_list_v6.csv") 
spdata <- read.csv(file.path(data_dir, "raw", "occurrence_data", "species_occurrences_cuba.csv"), sep = ",") 
spdata <- spdata %>% filter(cubanendemic != "y")#select only non endemic species

# Search GBIF
# 2. Set up GBIF credentials
usethis::edit_r_environ()
GBIF_USER=chapmanm
GBIF_PWD=Freya123;
GBIF_EMAIL=chapmanm@student.hu-berlin.de

# 3. Get keys for each species. Keys are used to search in GBIF
taxonkey_table <- name_backbone_checklist(spdata$species) 

#### If we want to use a ploygon as a mask to download the occurrences located 
#### within the polygon:
#### 

# # 4. Polygon. Download occurrences located within this polygon
# # Import polygon
# polygon_sf <- st_read("./roi_danube.gpkg")
# 
# # Get a WKT representation. It must be counterclockwise
# 
# # Extract the coordinates of the polygon
# vertices <- st_coordinates(polygon_sf)[, 1:2]
# 
# # Function to calculate the signed area of the polygon
# calculate_signed_area <- function(vertices) {
#   n <- nrow(vertices)
#   sum <- 0
#   for (i in 1:(n-1)) {
#     sum <- sum + (vertices[i, 1] * vertices[i + 1, 2] - vertices[i + 1, 1] * vertices[i, 2])
#   }
#   sum <- sum + (vertices[n, 1] * vertices[1, 2] - vertices[1, 1] * vertices[n, 2])
#   return(sum / 2)
# }
# 
# # Calculate the signed area before transformation
# signed_area_before <- calculate_signed_area(vertices)
# 
# # Reverse the order of vertices to transform to counterclockwise
# vertices_reversed <- vertices[nrow(vertices):1, ]
# 
# # Calculate the signed area after transformation
# signed_area_after <- calculate_signed_area(vertices_reversed)
# 
# # Determine winding order based on the signed area
# winding_order_before <- ifelse(signed_area_before < 0, "clockwise", "counterclockwise")
# winding_order_after <- ifelse(signed_area_after < 0, "clockwise", "counterclockwise")
# 
# # Create new polygon with reversed vertices
# polygon_reversed_sf <- st_polygon(list(vertices_reversed))
# 
# # Print the results
# cat("Signed Area Before:", signed_area_before, "\n")
# cat("Winding Order Before:", winding_order_before, "\n")
# cat("Signed Area After:", signed_area_after, "\n")
# cat("Winding Order After:", winding_order_after, "\n")
# 
# # Print the new polygon WKT
# cat("Reversed Polygon WKT:", st_as_text(st_sfc(polygon_reversed_sf)), "\n")
# 
# # wkt with winding order counter-clockwise
# wkt <- st_as_text(st_sfc(polygon_reversed_sf))
# ####################################################################################


# 5. Download occurrences within the polygon

# start a download on GBIF servers
gbif_download <- occ_download(
  pred_and(
    # remove default geospatial issues
    pred("HAS_GEOSPATIAL_ISSUE",FALSE),
    # keep only records with coordinates
    pred("HAS_COORDINATE",TRUE),
    # remove absent records
    #pred("OCCURRENCE_STATUS","PRESENT"),
    # remove fossils and living specimens
    pred_not(pred_in("BASIS_OF_RECORD",
                     c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")))
  ),
  # only records of species list
  pred_in("taxonKey", taxonkey_table$usageKey),
  # only records from countries in the Danube basin
  #pred_within(wkt),
  # records from 1970
  # pred_gte("year", 1970),
  format = "SIMPLE_CSV")

# checks if download is finished
occ_download_wait(gbif_download)

# retrieve a download from GBIF to the computer and load the download 
# from the computer to R
raw_downloaded <- occ_download_get(gbif_download,
                                   path = output_occur_path) |>
  occ_download_import()

#6. get and save citation
gbif_citation(gbif_download[1])$download |>
  writeLines(paste0(output_occur_path,
                    "/gbif_download_citation_",Sys.Date(),".txt"))

#7. Filter records
danube_records_filtered <- raw_downloaded |>
  #rename_with(tolower)|>  # set lowercase column names to work with CoordinateCleaner
  filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) |> # below 0.01 and with missing values
  filter(coordinateUncertaintyInMeters <= 1000 | is.na(coordinateUncertaintyInMeters)) |> # below 1 km and with missing values
  filter(!coordinateUncertaintyInMeters %in% c(301,3036,999,9999)) |> # remove with known default values
  cc_cen(buffer = 1000) |> # remove country centroids within 1km 
  cc_cap(buffer = 1000) |> # remove capitals centroids within 1km
  cc_inst(buffer = 1000) |> # remove zoo and herbaria within 1km 
  distinct(decimalLongitude,decimalLatitude,speciesKey,datasetKey, .keep_all = TRUE) |> # discard location duplicates
  filter(!is.na(year)) |> # remove records without year
  filter(species != "") 


# 8 Discard duplicates, keep only one record by species at each sub-catchment

# create a column with unique rows ids. this is required by the function
# extract_ids()
# danube_records_filtered$occurrence_id <- seq_along(danube_records_filtered$gbifID)
# 
# # extract sub-catchment ids of occurrence records
# subcatchmentId <- extract_ids(data = danube_records_filtered,
#                                lon = "decimalLongitude", lat = "decimalLatitude",
#                                id = "occurrence_id",
#                               quiet = FALSE,
#                                subc_layer = paste0(output_hydro90m,
#                                                   "/sub_catchment_danube.tif"))
# #add column with sub-catchment IDs
# danube_records_filtered <- danube_records_filtered |>
#    mutate(subcatchmentId = subcatchmentId$subcatchment_id)
# 
# # delete occurrences without sub-catchment id
# danube_records_filtered2 <- danube_records_filtered |>
# slice(-which(is.na(subcatchmentId)))


# In case there are occurrences from species atlas, apply this:

#  Filter out records that have the same Latitude but different Longitude
danube_records_filtered <- danube_records_filtered |>
  group_by(decimalLatitude) %>%
  filter(n_distinct(decimalLongitude) == 1) %>%
  ungroup()

# Further filter out records that have the same Longitude but different Latitude
danube_records_filtered <- danube_records_filtered %>%
  group_by(decimalLongitude) %>%
  filter(n_distinct(decimalLatitude) == 1) %>%
  ungroup()


write.csv(danube_records_filtered, file.path(output_occur_path, "gbif_data.csv"))
