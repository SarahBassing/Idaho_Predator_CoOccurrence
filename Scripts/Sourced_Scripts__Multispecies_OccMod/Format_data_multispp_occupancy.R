  #'  ---------------------------------------------
  #'  Format Data for Multispecies Occupancy Models 
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ---------------------------------------------
  #'  Script is sourced by 01_Run_multispp_occupancy_models.R and loads covariate 
  #'  data and detection histories for 5 predator species (black bear, bobcat,
  #'  coyote, mountain lion, and wolf) generated from camera trap data collected 
  #'  Summer 2020 and 2021 by the Northern Idaho Predator-Prey Project. Script 
  #'  formats data sets for use in multispecies occupancy models implemented with 
  #'  JAGS. 
  #'  
  #'  Required data:
  #'    1. DetectionHist_Smr20.RData and DetectionHist_Smr21.RData: each containing 
  #'       a list comprising 5 matrices. Each matrix contains the detection history of 
  #'       1 species during an 11-week sampling period (July 1 - Sept. 15, 2020/2021)
  #'       
  #'    2. Camera_Deployment_Data.csv: contains unique camera location name, GMU,
  #'       Season (summer 2020 or 2021) of deployment, height of camera from ground, 
  #'       and type of linear feature camera was facing (or a random location)
  #'       
  #'    3. Camera_Operable_Smr20.csv and Camera_Operable_Smr21.csv: each containing
  #'       the unique camera location, setup and retrieval dates, and any date
  #'       ranges when a camera was not operational
  #'       
  #'    4. Covariates_Smr20.csv and Covariates_Smr21.csv: contains site-level
  #'       covariate data including percent forest cover within 500m of each site,
  #'       elevation (m), species diversity index (H), and relative abundance 
  #'       indices (RAI) for elk, lagomorphs, moose, and white-tailed deer.
  #'  --------------------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(stringr)
  library(tidyverse)
  
  #'  Load annual detection histories and sampling effort
  load("./Data/DetectionHist_Smr20.RData")
  load("./Data/DetectionHist_Smr21.RData")
  load("./Data/SamplingEffort_Smr20.RData")
  load("./Data/SamplingEffort_Smr21.RData")
  
  #'  Read in camera deployment information and covariate data
  #'  Camera station data
  camera_stations <- read.csv("./Data/Camera_Deployment_Data.csv")
  #'  Camera operation data
  camera_op_smr20 <- read.csv("./Data/Camera_Operable_Smr20.csv")
  camera_op_smr21 <- read.csv("./Data/Camera_Operable_Smr21.csv")
  #'  Site-level covariates
  covs_smr20 <- read.csv("./Data/Covariates_Smr20.csv")
  covs_smr21 <- read.csv("./Data/Covariates_Smr21.csv")
  
  #'  Grab first detection history in each list and bind so each camera-year is represented
  dh1 <- DetectionHist_Smr20[[1]]
  dh2 <- DetectionHist_Smr21[[1]]
  newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
  colnames(dh1) <- newcols; colnames(dh2) <- newcols
  dh <- rbind(dh1, dh2) 
  
  #'  Reformat camera station data
  format_cam_station <- function(cams, season) {
    cams <- dplyr::select(cams, "NewLocationID") %>%
      mutate(GMU = str_extract(NewLocationID, "[^_]+"),
             Season = season) %>%
      left_join(camera_stations[camera_stations$Season == season,], by = c("NewLocationID", "Season")) %>%
      #'  Drop handful of duplicated rows that occur with left_join
      unique() %>%
      #'  Reduce categories to random, road, trail
      mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | 
                                   CameraFacing == "decommision" | CameraFacing == "decommission", "road", CameraFacing),
             CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing)) %>%
      arrange(NewLocationID)
    return(cams)
  }
  cams_smr20 <- format_cam_station(camera_op_smr20, season = "Smr20")
  cams_smr21 <- format_cam_station(camera_op_smr21, season = "Smr21") 
 
  #'  Scale and format site-level covariates
  format_covs <- function(cams_yr1, cams_yr2, covs_yr1, covs_yr2, rm_rows_yr1, rm_rows_yr2) {
    #'  Join camera data with extracted covariate data
    cam_covs1 <- full_join(cams_yr1, covs_yr1) %>% arrange(NewLocationID)
    cam_covs2 <- full_join(cams_yr2, covs_yr2) %>% arrange(NewLocationID)
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they are excluded
    #'  from the analysis entirely
    cam_covs1 <- cam_covs1[-rm_rows_yr1,]
    cam_covs2 <- cam_covs2[-rm_rows_yr2,]
    #'  Bind annual covariate data together
    cam_covs <- rbind(cam_covs1, cam_covs2)
    #'  Rename, format, and scale as needed
    formatted <- transmute(cam_covs,
                           NewLocationID = as.factor(NewLocationID),
                           Season = as.factor(Season),
                           GMU = as.factor(GMU),
                           CameraFacing = as.factor(CameraFacing),
                           Setup = as.factor(Setup),
                           Target = as.factor(Target),
                           Height = scale(CameraHeight_M),
                           PercForest = scale(perc_forest), 
                           Elev = scale(Elevation__10m2),
                           SppDiversity = scale(H),
                           Nelk = scale(elk_perday),  
                           Nlagomorph = scale(lagomorphs_perday),
                           Nmoose = scale(moose_perday),
                           Nwtd = scale(whitetaileddeer_perday)) 
    #'  Adjust reference category for CameraFacing factors
    order_camfacing <- c("random", "trail", "road")
    formatted <- formatted %>%
      mutate(
        CameraFacing = fct_relevel(CameraFacing, order_camfacing))
    
    return(formatted)
  }
  rm_rows_smr20 <- c(61, 79, 82, 98, 125, 157, 171, 177, 178, 181, 186, 192, 200, 214, 228, 235, 236, 259, 311, 334, 346, 361, 371, 379, 380, 385, 433, 437, 439, 458, 493)
  rm_rows_smr21 <- c(6, 106, 112, 116, 127, 145, 147, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682)
  zcovs <- format_covs(cams_yr1 = cams_smr20, cams_yr2 = cams_smr21, 
                       covs_yr1 = covs_smr20, covs_yr2 = covs_smr21, 
                       rm_rows_yr1 = rm_rows_smr20, rm_rows_yr2 = rm_rows_smr21) 
  
  #'  Save if you want
  #' write.csv(zcovs, file = "./Data/Camera_Stations.csv")
  
  #'  Correlation matrix to check for collinearity among continuous variables
  corr_matrix <- function(dat, firstcol, lastcol) {
    continuous_variables <- dat[,firstcol:lastcol]
    corr_all <- cor(continuous_variables)
    corr_all <- as.data.frame(round(corr_all, 2))
    print(corr_all)
    return(corr_all)
  }
  cov_corr_matrix <- corr_matrix(zcovs, firstcol = 7, lastcol = 14)
  
  #'  ---------------------------
  ####  Survey-level covariates  ####
  #'  ---------------------------
  #'  Bind annual sampling effort (must have matching column names)
  newcols <- c("NewLocationID", "occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11", "total_days", "total_hrs")
  colnames(SamplingEffort_Smr20) <- newcols; colnames(SamplingEffort_Smr21) <- newcols
  SamplingEffort <- rbind(SamplingEffort_Smr20, SamplingEffort_Smr21) %>%
    dplyr::select(-c(total_days, total_hrs))
  
  #'  Replace NAs in sampling effort with 0 - these sites truly were not surveyed
  #'  during those sampling occasions so survey effort = 0
  effort <- replace(SamplingEffort, is.na(SamplingEffort), 0)
  
  #'  Double check we have the same number of rows for both covariate data streams
  nrow(effort); nrow(zcovs)
  
  #'  Scale survey-level covariates
  scale_srvy_cov <- function(srvy_covs) {
    #'  Remove first column of camera names
    srvy_covs <- srvy_covs[,-1]
    
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(srvy_covs), na.rm = TRUE)
    sd <- sd(as.matrix(srvy_covs), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((srvy_covs - mu) / sd)
    scaled <- round(scaled, 3)
    scaled <- as.matrix(scaled)
    
    return(scaled)
  }
  #'  Note: rm_rows have already been removed from sampling effort data set
  zeffort <- scale_srvy_cov(srvy_covs = effort)
  
  #' #'  Create list of survey level covariates
  #' srvy_covs <- list(effort = zeffort)
  
  
  #'  -----------------------
  ####  Covariate mean & SD  ####
  #'  -----------------------
  #'  Save covariates in original format after removing specific rows (needed for plotting later on)
  unscaled_covs <- function(cams_yr1, cams_yr2, covs_yr1, covs_yr2, rm_rows_yr1, rm_rows_yr2) {
    #'  Join camera data with extracted covariate data
    cam_covs1 <- full_join(cams_yr1, covs_yr1) %>% arrange(NewLocationID)
    cam_covs2 <- full_join(cams_yr2, covs_yr2) %>% arrange(NewLocationID)
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they don't contribute
    #'  to detection data
    cam_covs1 <- cam_covs1[-rm_rows_yr1,]
    cam_covs2 <- cam_covs2[-rm_rows_yr2,]
    #'  Bind annual covariate data together
    cam_covs <- rbind(cam_covs1, cam_covs2)
    cam_covs <- transmute(cam_covs,
                          NewLocationID = as.factor(NewLocationID),
                          Season = as.factor(Season),
                          GMU = as.factor(GMU),
                          CameraFacing = as.factor(CameraFacing),
                          Setup = as.factor(Setup),
                          Target = as.factor(Target),
                          Height = CameraHeight_M,
                          PercForest = perc_forest, 
                          Elev = Elevation__10m2,
                          SppDiversity = H,
                          Nelk = elk_perday,    
                          Nlagomorph = lagomorphs_perday,
                          Nmoose = moose_perday,
                          Nwtd = whitetaileddeer_perday) %>%
      return(cam_covs)
  }
  covs <- unscaled_covs(cams_yr1 = cams_smr20, cams_yr2 = cams_smr21, 
                        covs_yr1 = covs_smr20, covs_yr2 = covs_smr21, 
                        rm_rows_yr1 = rm_rows_smr20, rm_rows_yr2 = rm_rows_smr21)
  #'  Save for later use
  save(covs, file = "./Data/covs.RData")
  