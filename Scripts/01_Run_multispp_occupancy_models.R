  #'  ---------------------------------------------
  #'  Fit Bayesian multispecies occupancy models
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ---------------------------------------------
  #'  Script sources Format_data_multispp_occupancy.R, which loads and formats data 
  #'  for 5 predator species (black bear, bobcat, coyote, mountain lion, and wolf) 
  #'  collected Summer 2020 and 2021 by the Northern Idaho Predator-Prey Project. 
  #'  Script bundles data, defines MCMC settings, and calls JAGS to fit a set of 
  #'  multispecies occupancy models to data for each predator pairing of interest. 
  #'  JAGS code for each model is sourced from a series of scripts in a folder 
  #'  named Sourced_Scripts__Multispecies_OccMod. Eight models are fit to data for
  #'  each predator dyad including:
  #'    1. Null
  #'    2. Habitat
  #'    3. Prey abundance
  #'    4. Prey diversity
  #'    5. Habitat with interaction on co-occupancy
  #'    6. Prey abundance with interaction on co-occupancy
  #'    7. Prey diversity with interaction on co-occupancy
  #'    8. Global model
  #'  Prey abundance, prey abundance with interaction, and global models have unique
  #'  scripts for the different predator dyads depending on each predator's primary
  #'  prey. e.g., script number follows 3.1, 3.2, 3.3, etc. for the different prey 
  #'  abundance models. 
  #'  The top model (selected using DIC with the 02_DIC_model_selection.R script) 
  #'  for each predator dyad was then refit with an additional interaction term 
  #'  on detection probability, resulting in 9 total models fit for each 
  #'  predator pairing. Scripts with the co-detection interaction are labeled 
  #'  with "px" following the model's number (e.g., 1_px_JAGS_null)
  #'  --------------------------------------------
  
  #'  Clean work space and load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(abind)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Call script to load and format detection and covariate data
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/Format_data_multispp_occupancy.R")
  
  #'  ----------------------------
  ####  Pair detection histories  ####
  #'  ----------------------------
  #'  Combine annual detection histories for each species
  all_detections <- function(dh1, dh2) {
    #'  Rename columns in each detection history
    newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
    colnames(dh1) <- newcols; colnames(dh2) <- newcols
    
    #'  Bind year 1 and 2 data together 
    #'  Remember that annual detection and covariate data are each arranged by  
    #'  location then stacked!
    dh <- rbind(dh1, dh2) 
    
    return(dh)
  }
  DH_bear <- all_detections(DetectionHist_Smr20[[1]], DetectionHist_Smr21[[1]])
  DH_bob <- all_detections(DetectionHist_Smr20[[2]], DetectionHist_Smr21[[2]])
  DH_coy <- all_detections(DetectionHist_Smr20[[3]], DetectionHist_Smr21[[3]])
  DH_lion <- all_detections(DetectionHist_Smr20[[4]], DetectionHist_Smr21[[4]])
  DH_wolf <- all_detections(DetectionHist_Smr20[[5]], DetectionHist_Smr21[[5]])
  
  #'  Combine species detection histories into an array (site x survey x species)
  detection_array <- function(spp1, spp2, name1, name2) {
    #'  List detection histories
    spp12_DH <- list(spp1, spp2)
    #'  Name lists based on specific species pairing
    names(spp12_DH) <- c(name1, name2)
    #'  Format list into a 3D array
    spp12_array <- abind(spp12_DH, along = 3)
    
    return(spp12_array)
  }
  wolf_bear_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_bear, name1 = "wolf", name2 = "bear")
  wolf_coy_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_coy, name1 = "wolf", name2 = "coyote")
  wolf_lion_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_lion, name1 = "wolf", name2 = "lion")
  lion_bear_DH <- detection_array(spp1 = DH_lion, spp2 = DH_bear, name1 = "lion", name2 = "bear")
  lion_bob_DH <- detection_array(spp1 = DH_lion, spp2 = DH_bob, name1 = "lion", name2 = "bobcat")
  coy_bob_DH <- detection_array(spp1 = DH_coy, spp2 = DH_bob, name1 = "coyote", name2 = "bobcat")
  
  #'  List 2-species detection arrays together for faster formatting below
  DH_array_list <- list(wolf_bear_DH, wolf_coy_DH, wolf_lion_DH, lion_bear_DH, lion_bob_DH, coy_bob_DH)
  
  #'  -------------------------
  ####  Format data for JAGS  ####
  #'  -------------------------
  #'  Define survey dimensions
  nsites <- dim(wolf_bear_DH)[1]
  nsurveys <- dim(wolf_bear_DH)[2]
  nspecies <- dim(wolf_bear_DH)[3]
  #'  Number of possible community states (species interactions): 00, 10, 01, 11
  ncat <- 2^nspecies
  #'  Number of sampling years
  nyears <- 2
  
  #####  Format detection histories  #####
  #'  --------------------------------
  #'  Function to convert 3D detection array into 2D multispecies detection history
  observation_state <- function(array_list) {
    #'  Merge species-specific observations into single 2-species observation state
    #'  (e.g., 00, 01, NANA) for each site and survey occasion
    ycat <- apply(array_list, c(1,2), paste, collapse = "")
    #'  Reclassify each observation state
    ycat[ycat == "00"] <- 1      # Neither species detected
    ycat[ycat == "10"] <- 2      # Only predator 1 (spp1) detected
    ycat[ycat == "01"] <- 3      # Only predator 2 (spp2) detected
    ycat[ycat == "11"] <- 4      # Both predator species detected
    ycat[ycat =="NANA"] <- NA    # Not sampled, no data
    #'  Make each column numeric so JAGS knows to handle it as a response variable
    ycat <- apply(ycat, 2, as.numeric)
    return(ycat)
  }
  multi_spp_DH_list <- lapply(DH_array_list, observation_state)
  names(multi_spp_DH_list) <- c("wolf_bear_DH", "wolf_coy_DH", "wolf_lion_DH", 
                                "lion_bear_DH", "lion_bob_DH", "coy_bob_DH")
  
  #####  Format covariate data  #####
  #'  ---------------------------
  #'  Format site-level covariates for detection sub-model
  final_covs <- zcovs %>%
    mutate(CameraFacing = as.factor(CameraFacing),
           Setup = ifelse(Setup == "ungulate", 0, 1),
           Height = as.numeric(Height),
           Year = ifelse(Season == "Smr20", 0, 1))
  table(final_covs[,"CameraFacing"])
  table(final_covs[,"Setup"])
  table(final_covs[, "Year"])
  
  ######  First order occupancy (psi|no second spp)  ######
  psi_cov <- matrix(NA, ncol = 10, nrow = nsites)
  psi_cov[,1] <- 1
  psi_cov[,2] <- final_covs$Setup
  psi_cov[,3] <- final_covs$Elev
  psi_cov[,4] <- final_covs$PercForest
  psi_cov[,5] <- final_covs$Year 
  psi_cov[,6] <- final_covs$SppDiversity
  psi_cov[,7] <- final_covs$Nelk
  psi_cov[,8] <- final_covs$Nmoose
  psi_cov[,9] <- final_covs$Nwtd
  psi_cov[,10] <- final_covs$Nlagomorph
  head(psi_cov); tail(psi_cov)

  ######  Second order occupancy (psi): 2-way interactions  ######
  psi_inxs_cov <- matrix(NA, ncol = 10, nrow = nsites)
  psi_inxs_cov[,1] <- 1
  psi_inxs_cov[,2] <- final_covs$Setup
  psi_inxs_cov[,3] <- final_covs$Elev
  psi_inxs_cov[,4] <- final_covs$PercForest
  psi_inxs_cov[,5] <- final_covs$Year 
  psi_inxs_cov[,6] <- final_covs$SppDiversity
  psi_inxs_cov[,7] <- final_covs$Nelk
  psi_inxs_cov[,8] <- final_covs$Nmoose
  psi_inxs_cov[,9] <- final_covs$Nwtd
  psi_inxs_cov[,10] <- final_covs$Nlagomorph
  head(psi_inxs_cov); tail(psi_inxs_cov)
  
  ######  First order detection (rho|no second spp)  ######
  rho_cov <- array(NA, dim = c(nsites, nsurveys, 3)) 
  rho_cov[,,1] <- 1
  rho_cov[,,2] <- final_covs$Setup
  rho_cov[,,3] <- zeffort
  head(rho_cov)
  
  ######  Second order detection (rho): 2-way interactions  ######
  rho_inxs_cov <- array(NA, dim = c(nsites, nsurveys, 3))
  rho_inxs_cov[,,1] <- 1
  rho_inxs_cov[,,2] <- final_covs$Setup
  rho_inxs_cov[,,3] <- zeffort
  head(rho_inxs_cov)
  
  
  #'  -----------------------------------
  ####  Data and MCMC settings for JAGS  ####
  #'  -----------------------------------
  #####  Bundle detection history and covariate data  ####
  #'  -------------------------------------------------
  #'  Function to bundle detection and covariate data
  bundle_data <- function(obs_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats, nspecies, nyears) {
    #'  list all pieces of data together
    bundled <- list(y = obs_array, psi_cov = psi_covs, psi_inxs_cov = psi_inxs,
                    rho_cov = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_1order), 
                    nsecond_order_psi = ncol(psi_2order), 
                    nfirst_order_rho = dim(rho_1order)[3], 
                    nsecond_order_rho = dim(rho_2order)[3], ncat = ncats, 
                    nspec = nspecies, nyear = nyears) 
    #'  Summarize to make sure it looks right
    str(bundled)
    return(bundled)
  }
  bundled_data_list <- lapply(multi_spp_DH_list, bundle_data, psi_covs = psi_cov, 
                               psi_inxs = psi_inxs_cov, rho_covs = rho_cov, 
                               rho_inxs = rho_inxs_cov, sites = nsites, 
                               surveys = nsurveys, psi_1order = psi_cov, 
                               psi_2order = psi_inxs_cov, rho_1order = rho_cov, 
                               rho_2order = rho_inxs_cov, ncats = ncat, 
                               nspecies = nspecies, nyears = nyears)
  #' #'  Save for later use
  #' save(bundled_data_list, file = "./Data/bundled_data_list.RData")
  
  
  #####  Initial values for model  #####
  #'  ------------------------------
  #'  Naive occupancy for each species at each site (site x spp matrix)
  initial_z <- function(bundled_dh) {
    #'  Naive occupancy for each species at each site (site x spp matrix)
    zinit <- apply(bundled_dh, c(1, 3), sum, na.rm = TRUE)
    zinit[zinit > 1] <- 1
    #'  Collapse 2-species detection state into 4 categories
    zcat <- apply(zinit, 1, paste, collapse = "")
    zcat[zcat == "00"] <- 1
    zcat[zcat == "10"] <- 2
    zcat[zcat == "01"] <- 3
    zcat[zcat == "11"] <- 4
    #'  Make z numeric again
    zcat <- as.numeric(zcat)
    
    return(zcat)
  }
  zinits <- lapply(DH_array_list, initial_z)
  names(zinits) <- c("wolf_bear_zcat", "wolf_coy_zcat", "wolf_lion_zcat", 
                     "lion_bear_zcat", "lion_bob_zcat", "coy_bob_zcat")
  
  #####  Parameters monitored  #####
  #'  --------------------------
  params <- c("betaSpp1", "betaSpp2", "alphaSpp1", "alphaSpp2", "betaSpp12", 
              "alphaSpp12", "alphaSpp21", "mean.psiSpp1", "mean.psiSpp2", 
              "mean.pSpp1", "mean.pSpp2", "z") 
  
  #####  MCMC settings  #####
  #'  -------------------
  nc <- 3
  nb <- 15000 
  nt <- 10
  na <- 1000
  #'  ni defined below for each species
  
  #'  -------------------
  ####  RUN JAGS MODELS  ####
  #'  -------------------
  #'  For each predator pairing:
  #'    1. set initial values with correct detection data
  #'    2. source and run each model in JAGS
  #'    3. visually inspect trace plots
  #'    4. save results
  #'    5. model selection (with 02_DIC_model_selection.R script)
  
  
  ####  Wolf-Bear Models  ####
  #'  ----------------------
  inits.wolf.bear <- function(){list(z = zinits[[1]])}    
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.bear.null <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null$summary)
  print(wolf.bear.null$DIC)
  which(wolf.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.null$samples)
  save(wolf.bear.null, file = "./Outputs/wolfbear_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.hab <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab$summary)
  print(wolf.bear.hab$DIC)
  which(wolf.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab$samples)
  save(wolf.bear.hab, file = "./Outputs/wolfbear_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabund <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabund$summary)
  print(wolf.bear.preyabund$DIC)
  which(wolf.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabund$samples)
  save(wolf.bear.preyabund, file = "./Outputs/wolfbear_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.preydiv <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                              "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preydiv$summary)
  print(wolf.bear.preydiv$DIC)
  which(wolf.bear.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preydiv$samples)
  save(wolf.bear.preydiv, file = "./Outputs/wolfbear_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.habx <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.habx$summary)
  print(wolf.bear.habx$DIC)
  which(wolf.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.habx$samples)
  save(wolf.bear.habx, file = "./Outputs/wolfbear_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd; psix(elk, moose, wtd); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabundx <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                         "./Outputs/06.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabundx$summary)
  print(wolf.bear.preyabundx$DIC)
  which(wolf.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabundx$samples)
  save(wolf.bear.preyabundx, file = "./Outputs/wolfbear_preyabundX.RData")
  
  #####  Prey diversity w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; psix(spp diversity); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.preydivx <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                         "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preydivx$summary)
  print(wolf.bear.preydivx$DIC)
  which(wolf.bear.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preydivx$samples)
  save(wolf.bear.preydivx, file = "./Outputs/wolfbear_preydivX.RData")
  
  #####  Global model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd, spp diversity; psix(elk, moose, wtd, spp diversity); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.1_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.global <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                             "./Outputs/08.1_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.global$summary)
  print(wolf.bear.global$DIC)
  which(wolf.bear.global$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.global$samples)
  save(wolf.bear.global, file = "./Outputs/wolfbear_global.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Prey diversity, no intx on psi
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort; px(.) 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_px_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.bear.preydiv.px <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                          "./Outputs/04_px_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort)_px(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preydiv.px$summary)
  print(wolf.bear.preydiv.px$DIC)
  which(wolf.bear.preydiv.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preydiv.px$samples)
  save(wolf.bear.preydiv.px, file = "./Outputs/wolfbear_preydiv_px.RData")
  
  
  
  #'  ----------------------
  ####  Wolf-Coyote Models  ####
  #'  ----------------------
  inits.wolf.coy <- function(){list(z = zinits[[2]])} 
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.coy.null <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.null$summary)
  print(wolf.coy.null$DIC)
  which(wolf.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.null$samples)
  save(wolf.coy.null, file = "./Outputs/wolfcoy_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                        "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary)
  print(wolf.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = "./Outputs/wolfcoy_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd, lago; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.2_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabund <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/03.2_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabund$summary)
  print(wolf.coy.preyabund$DIC)
  which(wolf.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabund$samples)
  save(wolf.coy.preyabund, file = "./Outputs/wolfcoy_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; _JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydiv <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                            "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydiv$summary)
  print(wolf.coy.preydiv$DIC)
  which(wolf.coy.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydiv$samples)
  save(wolf.coy.preydiv, file = "./Outputs/wolfcoy_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.habx <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                         "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.habx$summary)
  print(wolf.coy.habx$DIC)
  which(wolf.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.habx$samples)
  save(wolf.coy.habx, file = "./Outputs/wolfcoy_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd, lagomorph; psix(elk, moose, wtd, lagomorph); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.2_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabundx <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                               "./Outputs/06.2_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabundx$summary)
  print(wolf.coy.preyabundx$DIC)
  which(wolf.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabundx$samples)
  save(wolf.coy.preyabundx, file = "./Outputs/wolfcoy_preyabundX.RData")
  
  #####  Prey diversity w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; psix(spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydivx <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                             "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydivx$summary)
  print(wolf.coy.preydivx$DIC)
  which(wolf.coy.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydivx$samples)
  save(wolf.coy.preydivx, file = "./Outputs/wolfcoy_preydivX.RData")

  #####  Global model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd, lagomorph, spp diversity; psix(elk, moose, wtd, lagomorph, spp diversity); p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.2_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.global <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/08.2_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.global$summary)
  print(wolf.coy.global$DIC)
  which(wolf.coy.global$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.global$samples)
  save(wolf.coy.global, file = "./Outputs/wolfcoy_global.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = setup, year, elevation, forest; p = setup, effort; px(.) 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_px_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.coy.hab.px <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/02_px_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab.px$summary)
  print(wolf.coy.hab.px$DIC)
  which(wolf.coy.hab.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab.px$samples)
  save(wolf.coy.hab.px, file = "./Outputs/wolfcoy_hab_px.RData")
  
  
  #'  ---------------------
  ####  Wolf-Lion Models  ####
  #'  ---------------------
  inits.wolf.lion <- function(){list(z = zinits[[3]], mean.psiSpp1 = runif(1, 0.01, 0.15), 
                                     mean.psiSpp2 = runif(1, 0.2, 0.3),
                                     mean.pSpp1 = runif(1, 0.02, 0.12), mean.pSpp2 =  runif(1, 0.01, 0.06))}  
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.lion.null <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null$summary)
  print(wolf.lion.null$DIC)
  which(wolf.lion.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null$samples)
  save(wolf.lion.null, file = "./Outputs/wolflion_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.hab <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.hab$summary)
  print(wolf.lion.hab$DIC)
  which(wolf.lion.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.hab$samples)
  save(wolf.lion.hab, file = "./Outputs/wolflion_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabund <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                              "./Outputs/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabund$summary)
  print(wolf.lion.preyabund$DIC)
  which(wolf.lion.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabund$samples)
  save(wolf.lion.preyabund, file = "./Outputs/wolflion_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydiv <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                            "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydiv$summary)
  print(wolf.lion.preydiv$DIC)
  which(wolf.lion.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydiv$samples)
  save(wolf.lion.preydiv, file = "./Outputs/wolflion_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.habx <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.habx$summary)
  print(wolf.lion.habx$DIC)
  which(wolf.lion.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.habx$samples)
  save(wolf.lion.habx, file = "./Outputs/wolflion_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd; psix(elk, moose, wtd); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabundx <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/06.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabundx$summary)
  print(wolf.lion.preyabundx$DIC)
  which(wolf.lion.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabundx$samples)
  save(wolf.lion.preyabundx, file = "./Outputs/wolflion_preyabundX.RData")
  
  #####  Prey diversity inx model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; psix(spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydivx <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                             "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydivx$summary)
  print(wolf.lion.preydivx$DIC)
  which(wolf.lion.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydivx$samples)
  save(wolf.lion.preydivx, file = "./Outputs/wolflion_preydivX.RData")
 
  #####  Global model  #### 
  #'  psi = setup, year, elevation, forest, elk, moose, wtd, spp diversity; psix(elk, moose, wtd, spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.1_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.global <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/08.1_JAGS_global_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.global$summary)
  print(wolf.lion.global$DIC)
  which(wolf.lion.global$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.global$samples)
  save(wolf.lion.global, file = "./Outputs/wolflion_global.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_px_JAGS_null_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  wolf.lion.null.px <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                          "./Outputs/01_px_JAGS_null_psi(yr)_p(.)_px(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null.px$summary)
  print(wolf.lion.null.px$DIC)
  which(wolf.lion.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null.px$samples)
  save(wolf.lion.null.px, file = "./Outputs/wolflion_null_px.RData")
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = year; p(.); px(psi)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_px_JAGS_null_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  wolf.lion.null.px2 <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                           "./Outputs/01_px_JAGS_null_psi(yr)_p(.)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null.px2$summary)
  print(wolf.lion.null.px2$DIC)
  which(wolf.lion.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null.px2$samples)
  save(wolf.lion.null.px2, file = "./Outputs/wolflion_null_px.RData")
  
  
  
  #'  --------------------
  ####  Lion-Bear Models  ####
  #'  --------------------
  inits.lion.bear <- function(){list(z = zinits[[4]], mean.psiSpp1 = runif(1, 0.1, 0.35), 
                                     mean.psiSpp2 = runif(1, 0.6, 0.7), mean.pSpp1 = runif(1, 0.01, 0.06),
                                     mean.pSpp2 = runif(1, 0.15, 0.2))}
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bear.null <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null$summary)
  print(lion.bear.null$DIC)
  which(lion.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null$samples)
  save(lion.bear.null, file = "./Outputs/lionbear_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.hab <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.hab$summary)
  print(lion.bear.hab$DIC)
  which(lion.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.hab$samples)
  save(lion.bear.hab, file = "./Outputs/lionbear_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, elk, wtd; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.3_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabund <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                              "./Outputs/03.3_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabund$summary)
  print(lion.bear.preyabund$DIC)
  which(lion.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabund$samples)
  save(lion.bear.preyabund, file = "./Outputs/lionbear_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydiv <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydiv$summary)
  print(lion.bear.preydiv$DIC)
  which(lion.bear.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydiv$samples)
  save(lion.bear.preydiv, file = "./Outputs/lionbear_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.habx <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.habx$summary)
  print(lion.bear.habx$DIC)
  which(lion.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.habx$samples)
  save(lion.bear.habx, file = "./Outputs/lionbear_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, elk, wtd; psix(elk, wtd); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.3_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabundx <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                               "./Outputs/06.3_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabundx$summary)
  print(lion.bear.preyabundx$DIC)
  which(lion.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabundx$samples)
  save(lion.bear.preyabundx, file = "./Outputs/lionbear_preyabundX.RData")
  
  #####  Prey diversity w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; psix(spp diversity); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydivx <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                             "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydivx$summary)
  print(lion.bear.preydivx$DIC)
  which(lion.bear.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydivx$samples)
  save(lion.bear.preydivx, file = "./Outputs/lionbear_preydivX.RData")
  
  #####  Global model  ####
  #'  psi = setup, year, elevation, forest, elk, wtd, spp diversity; psix(elk, wtd, spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.3_JAGS_global_psi(global)_psix(global)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.global <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                           "./Outputs/08.3_JAGS_global_psi(global)_psix(global)_p(setup_effort)_bearlion.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.global$summary)
  print(lion.bear.global$DIC)
  which(lion.bear.global$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.global$samples)
  save(lion.bear.global, file = "./Outputs/lionbear_global.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_px_JAGS_null_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  lion.bear.null.px <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/01_px_JAGS_null_psi(yr)_p(.)_px(.).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null.px$summary)
  print(lion.bear.null.px$DIC)
  which(lion.bear.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null.px$samples)
  save(lion.bear.null.px, file = "./Outputs/lionbear_null_px.RData")
  
  
  #'  ----------------------
  ####  Lion-Bobcat Models  ####
  #'  ----------------------
  inits.lion.bob <- function(){list(z = zinits[[5]], mean.psiSpp1 = runif(1, 0.2, 0.35),
                                    mean.psiSpp2 = runif(1, 0.1, 0.3), mean.pSpp1 = runif(1, 0.01, 0.1),
                                    mean.pSpp2 = runif(1, 0.02, 0.1))}
  ni <- 75000

  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bob.null <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                       "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null$summary)
  print(lion.bob.null$DIC)
  which(lion.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null$samples)
  save(lion.bob.null, file = "./Outputs/lionbob_null.RData")
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.hab <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.hab$summary)
  print(lion.bob.hab$DIC)
  which(lion.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.hab$samples)
  save(lion.bob.hab, file = "./Outputs/lionbob_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, elk, wtd, lagomorphs; p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabund <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                              "./Outputs/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabund$summary)
  print(lion.bob.preyabund$DIC)
  which(lion.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabund$samples)
  save(lion.bob.preyabund, file = "./Outputs/lionbob_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydiv <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                            "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydiv$summary)
  print(lion.bob.preydiv$DIC)
  which(lion.bob.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydiv$samples)
  save(lion.bob.preydiv, file = "./Outputs/lionbob_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.habx <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                         "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.habx$summary)
  print(lion.bob.habx$DIC)
  which(lion.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.habx$samples)
  save(lion.bob.habx, file = "./Outputs/lionbob_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, elk, wtd, lagomorphs; psix(elk, wtd, lagomorphs); p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabundx <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                               "./Outputs/06.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabundx$summary)
  print(lion.bob.preyabundx$DIC)
  which(lion.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabundx$samples)
  save(lion.bob.preyabundx, file = "./Outputs/lionbob_abundX.RData")
  
  #####  Prey diversity w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; psix(spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydivx <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                             "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydivx$summary)
  print(lion.bob.preydivx$DIC)
  which(lion.bob.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydivx$samples)
  save(lion.bob.preydivx, file = "./Outputs/lionbob_preydivX.RData")
  
  #####  Global model  #### 
  #'  psi = setup, year, elevation, forest, elk, wtd, lagomorphs, spp diversity; psix(elk, wtd, lagomorphs, spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.4_JAGS_global_psi(global)_psix(global)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.global <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                          "./Outputs/08.4_JAGS_global_psi(global)_psix(global)_p(setup_effort)_lionbob.txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.global$summary)
  print(lion.bob.global$DIC)
  which(lion.bob.global$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.global$samples)
  save(lion.bob.global, file = "./Outputs/lionbob_global.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_px_JAGS_null_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  lion.bob.null.px <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                           "./Outputs/01_px_JAGS_null_psi(yr)_p(.)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null.px$summary)
  print(lion.bob.null.px$DIC)
  which(lion.bob.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null.px$samples)
  save(lion.bob.null.px, file = "./Outputs/lionbob_null_px.RData")
  
  
  
  #'  ------------------------
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  inits.coy.bob <- function(){list(z = zinits[[6]], mean.psiSpp1 = runif(1, 0.3, 0.4), 
                                   mean.psiSpp2 = runif(1, 0.03, 0.18),
                                   mean.pSpp1 = runif(1, 0.15, 0.2), mean.pSpp2 = runif(1, 0.03, 0.1))}
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  coy.bob.null <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/01_JAGS_null_psi(yr)_p(.).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.null$summary)
  print(coy.bob.null$DIC)
  which(coy.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.null$samples)
  save(coy.bob.null, file = "./Outputs/coybob_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, elevation, forest; p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/02_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  print(coy.bob.hab$DIC)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = "./Outputs/coybob_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, elevation, forest, wtd, lagomorphs; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.5_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabund <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                              "./Outputs/03.5_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_coybob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabund$summary)
  print(coy.bob.preyabund$DIC)
  which(coy.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabund$samples)
  save(coy.bob.preyabund, file = "./Outputs/coybob_preyabund.RData")
  
  #####  Prey diversity model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.preydiv <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/04_JAGS_preydiv_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preydiv$summary)
  print(coy.bob.preydiv$DIC)
  which(coy.bob.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preydiv$samples)
  save(coy.bob.preydiv, file = "./Outputs/coybob_preydiv.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.habx <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/05_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx$summary)
  print(coy.bob.habx$DIC)
  which(coy.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx$samples)
  save(coy.bob.habx, file = "./Outputs/coybob_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest, wtd, lagomorphs; psix(wtd, lagomorphs); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/06.5_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabundx <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                               "./Outputs/06.5_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabundx$summary)
  print(coy.bob.preyabundx$DIC)
  which(coy.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabundx$samples)
  save(coy.bob.preyabundx, file = "./Outputs/coybob_preyabundX.RData")
  
  #####  Prey diversity w/ interaction model  #### 
  #'  psi = setup, year, elevation, forest; psix(spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.preydivx <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                             "./Outputs/07_JAGS_preydivX_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preydivx$summary)
  print(coy.bob.preydivx$DIC)
  which(coy.bob.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preydivx$samples)
  save(coy.bob.preydivx, file = "./Outputs/coybob_preydivX.RData")
  
  #####  Global model  #### 
  #'  psi = setup, year, elevation, forest, wtd, lagomorphs, spp diversity; psix(wtd, lagomorphs, spp diversity); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.5_JAGS_global_psi(global)_psix(global)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.global <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                             "./Outputs/08.5_JAGS_global_psi(global)_psix(global)_p(setup_effort)_coybob.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.global$summary)
  print(coy.bob.global$DIC)
  which(coy.bob.global$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.global$samples)
  save(coy.bob.global, file = "./Outputs/coybob_global.RData")
  
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Global
  #'  psi = setup, year, elevation, forest, wtd, lagomorphs, spp diversity; psix(wtd, lagomorphs, spp diversity); p = setup, effort; px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/08.5_px_JAGS_global_psi(global)_psix(global)_p(setup_effort)_px(.)_coybob.R")
  start.time = Sys.time()
  coy.bob.global.px <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/08.5_px_JAGS_global_psi(global)_psix(global)_p(setup_effort)_px(.)_coybob.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.global.px$summary)
  print(coy.bob.global.px$DIC)
  which(coy.bob.global.px$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.global.px$samples)
  save(coy.bob.global.px, file = "./Outputs/coybob_global_px.RData")
  
  
  