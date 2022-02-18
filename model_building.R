# This scripts include the step 4 and step 5 of the range map generate framework

# packages and functions needed ----

library(magrittr)
library(data.table)
library(raster)
library(sf)
library(dggridR)
library(dismo)
library(rJava)
library(rgdal)


lapply(list.files("functions/", full.names = TRUE),
       function(fn)source(fn)) %>% 
  invisible()


# import data needed ----
# species info (with complex tags)
species_table <- 
  fread("data/species_table.csv")

# species occurrence
occ <- 
  list.files("data/occurrence/", full.names = TRUE) %>% 
  lapply(fread) %>% 
  do.call(rbind, .)

# environment variable
# for the environment data of this analysis, please go to https://github.com/WanJyunChen/Taiwan_environmental_dataset to download the complete data set
# than unzip the GeoTIFF file, and move the tif files to folder data/environment

# environment variables used in this analysis
env_list <- readRDS("data/environment/environment_list.rds")

env <- list()
env[[1]] <- 
  list.files("data/environment") %>%
  .[. %in% env_list] %>%
  sprintf("data/environment/%s", .) %>% 
  stack

# generate grids for subsampling
# 1 km hexagonal grid
dggs_1 <- dgconstruct(spacing = 1)

# 5km hexagonal grid
dggs_5 <- dgconstruct(spacing = 5)

# occurrence data preparation ----
# add complex tags to occurrence 
occ_complex <- 
  dplyr::left_join(occ, species_table) %>% 
  dplyr::mutate(scientific_name = ifelse(complex_tags == "", 
                                  scientific_name,
                                  complex_tags))

# add grid ids to occurrence
occ_complex_grid_sf <- 
  occ_complex %>% 
  dplyr::mutate(grid_1km = dgGEO_to_SEQNUM(dggs_1, longitude, latitude)$seqnum,
                grid_5km = dgGEO_to_SEQNUM(dggs_5, longitude, latitude)$seqnum) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(3826)


occ_complex_grid_3826 <- 
  cbind(st_set_geometry(occ_complex_grid_sf, NULL),
        st_coordinates(occ_complex_grid_sf))

# MaxEnt model building (include subsampling) ----
# maxent arguments
maxent_args <- 
  c("responsecurves=false",
    "randomseed=true",
    "askoverwrite=false",
    "randomtestpoints=20",
    "replicates=5",
    "replicatetype=subsample",
    "writebackgroundpredictions=true",
    "outputgrids=false",
    "appendtoresultsfile=true")

repeat_times <- 20

# list of species
species_list <- 
  occ_complex_grid_3826$scientific_name %>% 
  unique

# Run maxent and export all the results of models
xm <- lapply(seq_along(species_list),
                function(sp){
                  
                    species_data <- 
                      occ_complex_grid_3826 %>% 
                      dplyr::filter(scientific_name == species_list[sp])
                    
                    # r for number of repeat, for this dataset we repeat 20 times for each species
                    models <- 
                      lapply(1:repeat_times, function(r){
                        
                        # step 4 subsampling --
                        # subsample of species data by 5km grids
                        sp_data <- 
                          species_data %>% 
                          # for each species or species complex, only sample one occurrence within the same 5 km grid
                          .[, .SD[sample(.N ,(min(.N, 1)))], by = grid_5km] %>% 
                          .[, list(X, Y)] %>% 
                          as.matrix
                        
                        # step 5 model building --
                        # the model will store in models/
                        
                        # create folder for model
                        path_sdm <- sprintf("%s/%s/%02d", 
                                            "models", species_list[sp], r)
                        ifelse(!dir.exists(path_sdm), 
                               dir.create(path_sdm, recursive = TRUE), 
                               FALSE)
                        
                        # run MaxEnt
                        model <- SDM2tif(predictors = env, 
                                         data = sp_data,
                                         path.sdm = path_sdm, 
                                         maxent_args = maxent_args)
                        
                        return(path_sdm)
                      })
                  }
                  )

# export binary maps by threshold suggested by experts ----
# threshold suggested by expert
threshold_table <- 
  fread("data/threshold_table.csv")

# get threshold value from maxent results
thresholds <- 
  lapply(seq_along(species_list), function(sp)
    lapply(1:repeat_times, 
           function(r)
             fread(sprintf("models/%s/%02d/maxentResults.csv", 
                           species_list[sp], r)
             )
    ) %>% 
      do.call(rbind, .) %>%
      # extract Cloglog thresholds only
      .[, colnames(.) %like% "Cloglog threshold", with = FALSE] %>%
      # average the thresholds value for each species
      .[, lapply(.SD, mean, na.rm=TRUE)] %>% 
      .[, scientific_name := species_list[sp]]
  ) %>% 
  do.call(rbind, .) %>% 
  melt(id = "scientific_name",
       measure = patterns("Cloglog threshold"),
       variable.name = "threshold",
       value.name = "value") %>% 
  setDT

# order of threshold
threshold_list <- 
  fread(sprintf("models/%s/%02d/maxentResults.csv", 
                species_list[1], 1)
  ) %>% 
  .[, colnames(.) %like% "Cloglog threshold", with = FALSE] %>% 
  names

# select threshold value that suggested by experts
threshold_use <- 
  lapply(seq_along(species_list),
         function(sp)
           thresholds[scientific_name == species_list[sp]] %>% 
           .[threshold == threshold_list[as.numeric(threshold_table$threshold_suggest[sp])]]
  ) %>%
  do.call(rbind, .)

# turn continuous model to binary model by threshold suggested by expert
binary_map <- 
  lapply(seq_along(species_list), 
            function(sp){
              
              range_map <- 
                lapply(1:repeat_times, function(r)
                  raster(sprintf("models/%s/%02d/model_predict_avg_1.tif", 
                                 species_list[sp], r))
                ) %>%
                stack %>% 
                # average the predicted value of models
                calc(fun = mean) %>% 
                # turn continuous value to bonary value
                calc(fun = function(x)ifelse(x >= threshold_use$value[sp], 1, 0)) %>% 
                rasterToPolygons()
              
              path_map <- "range_maps/"
              
              if (!dir.exists(path_map)){dir.create(path_map)}
              
              writeOGR(range_map, 
                       dsn = path_map,
                       layer = species_list[sp],
                       driver = "ESRI Shapefile",
                       delete_dsn = TRUE,
                       overwrite_layer = TRUE)
            }
            
  )
