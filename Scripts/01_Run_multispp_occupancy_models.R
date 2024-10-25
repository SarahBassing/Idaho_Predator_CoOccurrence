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
  #'  named Sourced_Scripts__Multispecies_OccMod. Five models are fit to data for
  #'  each predator dyad including:
  #'    1. Null
  #'    2. Habitat
  #'    3. Prey abundance
  #'    4. Habitat with interaction on co-occupancy
  #'    5. Prey abundance with interaction on co-occupancy
  #'  
  #'  Additional notes:
  #'  The prey abundance and prey abundance with interaction models have unique
  #'  scripts for the different predator dyads depending on each predator's primary
  #'  prey. e.g., script number follows 3.1, 3.2, 3.3, etc. for the different prey 
  #'  abundance models. The wolf-lion and wolf-bear dyads have the same set of primary 
  #'  prey and the same code is used for both (scripts ending in _wolfbearlion.R).
  #'  The lion-bobcat and bear-coyote dyads also have the same set of primary prey 
  #'  and the same code is used for both (scripts ending in _lionbob_bearcoy.R)
  #'  
  #'  The top model (selected using DIC with the 02_DIC_model_selection.R script) 
  #'  for each predator dyad was then refit with an additional interaction term 
  #'  on detection probability (px(.)). Goodness-of-Fit tests were conducted on the 
  #'  top model and the top model with the interaction on detection prob., resulting 
  #'  in 8 total models fit for each predator pairing. Scripts with the co-detection 
  #'  interaction are labeled w/ "px(.)" at the end of the script name. Scripts
  #'  with the Goodness-of-Fit code are labeled w/ "GoF" at the end of the script 
  #'  name. Script numbering for these added parameterizations are:
  #'    1.0 = GoF for top model
  #'    1.1 = top model with px(.)
  #'    1.2 = top model with px(.) and GoF
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
  psi_cov[,3] <- final_covs$Year
  psi_cov[,4] <- final_covs$PercForest
  psi_cov[,5] <- final_covs$Elev
  psi_cov[,6] <- final_covs$Nelk 
  psi_cov[,7] <- final_covs$Nmoose
  psi_cov[,8] <- final_covs$Nwtd
  psi_cov[,9] <- final_covs$Nlagomorph
  psi_cov[,10] <- final_covs$TRI
  head(psi_cov); tail(psi_cov)

  ######  Second order occupancy (psi): 2-way interactions  ######
  psi_inxs_cov <- matrix(NA, ncol = 10, nrow = nsites)
  psi_inxs_cov[,1] <- 1
  psi_inxs_cov[,2] <- final_covs$Setup
  psi_inxs_cov[,3] <- final_covs$Year
  psi_inxs_cov[,4] <- final_covs$PercForest
  psi_inxs_cov[,5] <- final_covs$Elev 
  psi_inxs_cov[,6] <- final_covs$Nelk
  psi_inxs_cov[,7] <- final_covs$Nmoose
  psi_inxs_cov[,8] <- final_covs$Nwtd
  psi_inxs_cov[,9] <- final_covs$Nlagomorph
  psi_inxs_cov[,10] <- final_covs$TRI
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
              "mean.pSpp1", "mean.pSpp2", "z", "chi2.obs", "chi2.sim", "ft.sims", 
              "FT.obs", "FT.sims")
  #'  FYI: monitoring chi2 & FT parameters for Goodness-of-Fit test increase 
  #'  computation time a lot (when included in the JAGS code)!
  
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
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.bear.null <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null$summary[1:15,])
  print(wolf.bear.null$DIC)
  which(wolf.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.null$samples)
  save(wolf.bear.null, file = "./Outputs/wolfbear_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.hab <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab$summary[1:30,])
  print(wolf.bear.hab$DIC)
  which(wolf.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab$samples)
  save(wolf.bear.hab, file = "./Outputs/wolfbear_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabund <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabund$summary[1:35,])
  print(wolf.bear.preyabund$DIC)
  which(wolf.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabund$samples)
  save(wolf.bear.preyabund, file = "./Outputs/wolfbear_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.habx <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.habx$summary[1:30,])
  print(wolf.bear.habx$DIC)
  which(wolf.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.habx$samples)
  save(wolf.bear.habx, file = "./Outputs/wolfbear_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd; psix(elk, moose, wtd); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabundx <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                         "./Outputs/05.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabundx$summary[1:40,])
  print(wolf.bear.preyabundx$DIC)
  which(wolf.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabundx$samples)
  save(wolf.bear.preyabundx, file = "./Outputs/wolfbear_preyabundX.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model: Habitat model
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort; px(.) 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.bear.hab.px <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                          "./Outputs/02.1.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab.px$summary[1:30,])
  print(wolf.bear.hab.px$DIC)
  which(wolf.bear.hab.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab.px$samples)
  save(wolf.bear.hab.px, file = "./Outputs/wolfbear_hab_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.0_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  wolf.bear.hab.GoF <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                            "./Outputs/02.1.0_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.bear.hab.GoF_X2 <- mean(wolf.bear.hab.GoF$sims.list$chi2.sim > wolf.bear.hab.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.bear.hab.GoF_FT <- mean(wolf.bear.hab.GoF$sims.list$FT.sims > wolf.bear.hab.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.bear.hab.GoF, file = "./Outputs/wolfbear_hab_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.2_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.)_GoF.R")
  start.time = Sys.time()
  wolf.bear.hab.px.GoF <- jags(bundled_data_list[[1]], inits = inits.wolf.bear, params,
                               "./Outputs/02.1.2_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.)_GoF.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.bear.hab.px.GoF_X2 <- mean(wolf.bear.hab.px.GoF$sims.list$chi2.sim > wolf.bear.hab.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.bear.hab.px.GoF_FT <- mean(wolf.bear.hab.px.GoF$sims.list$FT.sims > wolf.bear.hab.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.bear.hab.px.GoF, file = "./Outputs/wolfbear_hab_px_GoF.RData")
  

  #'  ----------------------
  ####  Wolf-Coyote Models  ####
  #'  ----------------------
  inits.wolf.coy <- function(){list(z = zinits[[2]])} 
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.coy.null <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.null$summary[1:15,])
  print(wolf.coy.null$DIC)
  which(wolf.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.null$samples)
  save(wolf.coy.null, file = "./Outputs/wolfcoy_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                        "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary[1:30,])
  print(wolf.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = "./Outputs/wolfcoy_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd, lago; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.2_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabund <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/03.2_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabund$summary[1:35,])
  print(wolf.coy.preyabund$DIC)
  which(wolf.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabund$samples)
  save(wolf.coy.preyabund, file = "./Outputs/wolfcoy_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.habx <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                         "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.habx$summary[1:30,])
  print(wolf.coy.habx$DIC)
  which(wolf.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.habx$samples)
  save(wolf.coy.habx, file = "./Outputs/wolfcoy_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd, lagomorph; psix(elk, moose, wtd, lagomorph); p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.2_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabundx <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                               "./Outputs/05.2_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabundx$summary[1:40,])
  print(wolf.coy.preyabundx$DIC)
  which(wolf.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabundx$samples)
  save(wolf.coy.preyabundx, file = "./Outputs/wolfcoy_preyabundX.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Habitat
  #'  psi = setup, year, elevation, forest; p = setup, effort; px(.) 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.coy.hab.px <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/02.1.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab.px$summary[1:30,])
  print(wolf.coy.hab.px$DIC)
  which(wolf.coy.hab.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab.px$samples)
  save(wolf.coy.hab.px, file = "./Outputs/wolfcoy_hab_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.0_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab.GoF <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                          "./Outputs/02.1.0_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.coy.hab.GoF_X2 <- mean(wolf.coy.hab.GoF$sims.list$chi2.sim > wolf.coy.hab.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.coy.hab.GoF_FT <- mean(wolf.coy.hab.GoF$sims.list$FT.sims > wolf.coy.hab.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.coy.hab.GoF, file = "./Outputs/wolfcoy_hab_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1.2_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.coy.hab.px.GoF <- jags(bundled_data_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/02.1.2_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.coy.hab.px.GoF_X2 <- mean(wolf.coy.hab.px.GoF$sims.list$chi2.sim > wolf.coy.hab.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.coy.hab.px.GoF_FT <- mean(wolf.coy.hab.px.GoF$sims.list$FT.sims > wolf.coy.hab.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.coy.hab.px.GoF, file = "./Outputs/wolfcoy_hab_px_GoF.RData")
  
  
  #'  ---------------------
  ####  Wolf-Lion Models  ####
  #'  ---------------------
  inits.wolf.lion <- function(){list(z = zinits[[3]], mean.psiSpp1 = runif(1, 0.01, 0.15), 
                                     mean.psiSpp2 = runif(1, 0.2, 0.3),
                                     mean.pSpp1 = runif(1, 0.02, 0.12), mean.pSpp2 =  runif(1, 0.01, 0.06))}  
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.lion.null <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null$summary[1:15,])
  print(wolf.lion.null$DIC)
  which(wolf.lion.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null$samples)
  save(wolf.lion.null, file = "./Outputs/wolflion_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.hab <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.hab$summary[1:30,])
  print(wolf.lion.hab$DIC)
  which(wolf.lion.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.hab$samples)
  save(wolf.lion.hab, file = "./Outputs/wolflion_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabund <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                              "./Outputs/03.1_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabund$summary[1:35,])
  print(wolf.lion.preyabund$DIC)
  which(wolf.lion.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabund$samples)
  save(wolf.lion.preyabund, file = "./Outputs/wolflion_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.habx <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.habx$summary[1:30,])
  print(wolf.lion.habx$DIC)
  which(wolf.lion.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.habx$samples)
  save(wolf.lion.habx, file = "./Outputs/wolflion_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, moose, wtd; psix(elk, moose, wtd); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabundx <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/05.1_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabundx$summary[1:40,])
  print(wolf.lion.preyabundx$DIC)
  which(wolf.lion.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabundx$samples)
  save(wolf.lion.preyabundx, file = "./Outputs/wolflion_preyabundX.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.1_JAGS_null_psi(yr)_px(.).R")
  start.time = Sys.time()
  wolf.lion.null.px <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                          "./Outputs/01.1.1_JAGS_null_psi(yr)_px(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null.px$summary)
  print(wolf.lion.null.px$DIC)
  which(wolf.lion.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null.px$samples)
  save(wolf.lion.null.px, file = "./Outputs/wolflion_null_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  wolf.lion.null.GoF <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                            "./Outputs/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.lion.null.GoF_X2 <- mean(wolf.lion.null.GoF$sims.list$chi2.sim > wolf.lion.null.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.lion.null.GoF_FT <- mean(wolf.lion.null.GoF$sims.list$FT.sims > wolf.lion.null.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.lion.null.GoF, file = "./Outputs/wolflion_null_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.R")
  start.time = Sys.time()
  wolf.lion.null.px.GoF <- jags(bundled_data_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (wolf.lion.null.px.GoF_X2 <- mean(wolf.lion.null.px.GoF$sims.list$chi2.sim > wolf.lion.null.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (wolf.lion.null.px.GoF_FT <- mean(wolf.lion.null.px.GoF$sims.list$FT.sims > wolf.lion.null.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(wolf.lion.null.px.GoF, file = "./Outputs/wolflion_null_px_GoF.RData")
  
  
  #'  --------------------
  ####  Lion-Bear Models  ####
  #'  --------------------
  inits.lion.bear <- function(){list(z = zinits[[4]], mean.psiSpp1 = runif(1, 0.1, 0.35), 
                                     mean.psiSpp2 = runif(1, 0.6, 0.7), mean.pSpp1 = runif(1, 0.01, 0.06),
                                     mean.pSpp2 = runif(1, 0.15, 0.2))}
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bear.null <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null$summary[1:15,])
  print(lion.bear.null$DIC)
  which(lion.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null$samples)
  save(lion.bear.null, file = "./Outputs/lionbear_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.hab <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.hab$summary[1:30,])
  print(lion.bear.hab$DIC)
  which(lion.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.hab$samples)
  save(lion.bear.hab, file = "./Outputs/lionbear_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, wtd; p = setup, effort  
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.3_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabund <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                              "./Outputs/03.3_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabund$summary[1:35,])
  print(lion.bear.preyabund$DIC)
  which(lion.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabund$samples)
  save(lion.bear.preyabund, file = "./Outputs/lionbear_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.habx <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.habx$summary[1:30,])
  print(lion.bear.habx$DIC)
  which(lion.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.habx$samples)
  save(lion.bear.habx, file = "./Outputs/lionbear_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, wtd; psix(elk, wtd); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.3_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabundx <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                               "./Outputs/05.3_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabundx$summary[1:40,])
  print(lion.bear.preyabundx$DIC)
  which(lion.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabundx$samples)
  save(lion.bear.preyabundx, file = "./Outputs/lionbear_preyabundX.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.1_JAGS_null_psi(yr)_px(.).R")
  start.time = Sys.time()
  lion.bear.null.px <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/01.1.1_JAGS_null_psi(yr)_px(.).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null.px$summary[1:15,])
  print(lion.bear.null.px$DIC)
  which(lion.bear.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null.px$samples)
  save(lion.bear.null.px, file = "./Outputs/lionbear_null_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  lion.bear.null.GoF <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                             "./Outputs/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (lion.bear.null.GoF_X2 <- mean(lion.bear.null.GoF$sims.list$chi2.sim > lion.bear.null.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (lion.bear.null.GoF_FT <- mean(lion.bear.null.GoF$sims.list$FT.sims > lion.bear.null.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(lion.bear.null.GoF, file = "./Outputs/lionbear_null_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.R")
  start.time = Sys.time()
  lion.bear.null.px.GoF <- jags(bundled_data_list[[4]], inits = inits.lion.bear, params,
                                "./Outputs/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.txt",
                                n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (lion.bear.null.px.GoF_X2 <- mean(lion.bear.null.px.GoF$sims.list$chi2.sim > lion.bear.null.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (lion.bear.null.px.GoF_FT <- mean(lion.bear.null.px.GoF$sims.list$FT.sims > lion.bear.null.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(lion.bear.null.px.GoF, file = "./Outputs/lionbear_null_px_GoF.RData")
  
  
  #'  ----------------------
  ####  Lion-Bobcat Models  ####
  #'  ----------------------
  inits.lion.bob <- function(){list(z = zinits[[5]], mean.psiSpp1 = runif(1, 0.2, 0.35),
                                    mean.psiSpp2 = runif(1, 0.1, 0.3), mean.pSpp1 = runif(1, 0.01, 0.1),
                                    mean.pSpp2 = runif(1, 0.02, 0.1))}
  ni <- 75000

  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bob.null <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                       "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null$summary[1:15,])
  print(lion.bob.null$DIC)
  which(lion.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null$samples)
  save(lion.bob.null, file = "./Outputs/lionbob_null.RData")
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.hab <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.hab$summary[1:30,])
  print(lion.bob.hab$DIC)
  which(lion.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.hab$samples)
  save(lion.bob.hab, file = "./Outputs/lionbob_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, wtd, lagomorphs; p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob_bearcoy.R")
  start.time = Sys.time()
  lion.bob.preyabund <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                              "./Outputs/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob_bearcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabund$summary[1:35,])
  print(lion.bob.preyabund$DIC)
  which(lion.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabund$samples)
  save(lion.bob.preyabund, file = "./Outputs/lionbob_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.habx <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                         "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.habx$summary[1:30,])
  print(lion.bob.habx$DIC)
  which(lion.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.habx$samples)
  save(lion.bob.habx, file = "./Outputs/lionbob_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, wtd, lagomorphs; psix(elk, wtd, lagomorphs); p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob_bearcoy.R")
  start.time = Sys.time()
  lion.bob.preyabundx <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                               "./Outputs/05.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob_bearcoy.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabundx$summary[1:40,])
  print(lion.bob.preyabundx$DIC)
  which(lion.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabundx$samples)
  save(lion.bob.preyabundx, file = "./Outputs/lionbob_abundX.RData")
  
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
  print(lion.bob.null.px$summary[1:15,])
  print(lion.bob.null.px$DIC)
  which(lion.bob.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null.px$samples)
  save(lion.bob.null.px, file = "./Outputs/lionbob_null_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  lion.bob.null.GoF <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                           "./Outputs/01.1.0_JAGS_null_psi(yr)_p(.)_GoF.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (lion.bob.null.GoF_X2 <- mean(lion.bob.null.GoF$sims.list$chi2.sim > lion.bob.null.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (lion.bob.null.GoF_FT <- mean(lion.bob.null.GoF$sims.list$FT.sims > lion.bob.null.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(lion.bob.null.GoF, file = "./Outputs/lionbob_null_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.R")
  start.time = Sys.time()
  lion.bob.null.px.GoF <- jags(bundled_data_list[[5]], inits = inits.lion.bob, params,
                            "./Outputs/01.1.2_JAGS_null_psi(yr)_px(.)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (lion.bob.null.px.GoF_X2 <- mean(lion.bob.null.px.GoF$sims.list$chi2.sim > lion.bob.null.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (lion.bob.null.px.GoF_FT <- mean(lion.bob.null.px.GoF$sims.list$FT.sims > lion.bob.null.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(lion.bob.null.px.GoF, file = "./Outputs/lionbob_null_px_GoF.RData")
  
  
  #'  ------------------------
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  inits.coy.bob <- function(){list(z = zinits[[6]], mean.psiSpp1 = runif(1, 0.3, 0.4), 
                                   mean.psiSpp2 = runif(1, 0.03, 0.18),
                                   mean.pSpp1 = runif(1, 0.15, 0.2), mean.pSpp2 = runif(1, 0.03, 0.1))}
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  coy.bob.null <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/01.1_JAGS_null_psi(yr)_p(.).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.null$summary[1:15,])
  print(coy.bob.null$DIC)
  which(coy.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.null$samples)
  save(coy.bob.null, file = "./Outputs/coybob_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort 
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary[1:30,])
  print(coy.bob.hab$DIC)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = "./Outputs/coybob_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, wtd, lagomorphs; p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/03.5_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabund <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                              "./Outputs/03.5_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_coybob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabund$summary[1:35,])
  print(coy.bob.preyabund$DIC)
  which(coy.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabund$samples)
  save(coy.bob.preyabund, file = "./Outputs/coybob_preyabund.RData")
  
  #####  Habitat w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.habx <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/04.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx$summary[1:30,])
  print(coy.bob.habx$DIC)
  which(coy.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx$samples)
  save(coy.bob.habx, file = "./Outputs/coybob_habX.RData")
  
  #####  Prey abundance w/ interaction model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, wtd, lagomorphs; psix(wtd, lagomorphs); p = setup, effort
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/05.5_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabundx <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                               "./Outputs/05.5_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabundx$summary[1:40,])
  print(coy.bob.preyabundx$DIC)
  which(coy.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabundx$samples)
  save(coy.bob.preyabundx, file = "./Outputs/coybob_preyabundX.RData")
  
  #####  Top model w/ interaction on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Habitat with an interaction
  #'  psi = setup, year, elevation, forest, ruggedness, wtd, lagomorphs; psix(.); p = setup, effort; px(.)
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x).R")
  start.time = Sys.time()
  coy.bob.habx.px <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/04.1.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx.px$summary[1:35,])
  print(coy.bob.habx.px$DIC)
  which(coy.bob.habx.px$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx.px$samples)
  save(coy.bob.habx.px, file = "./Outputs/coybob_habX_px.RData")
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  coy.bob.habx.GoF <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (coy.bob.habx.GoF_X2 <- mean(coy.bob.habx.GoF$sims.list$chi2.sim > coy.bob.habx.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (coy.bob.habx.GoF_FT <- mean(coy.bob.habx.GoF$sims.list$FT.sims > coy.bob.habx.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(coy.bob.habx.GoF, file = "./Outputs/coybob_habX_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1.2_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x)_GoF.R")
  start.time = Sys.time()
  coy.bob.habx.px.GoF <- jags(bundled_data_list[[6]], inits = inits.coy.bob, params,
                           "./Outputs/04.1.2_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x)_GoF.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (coy.bob.habx.px.GoF_X2 <- mean(coy.bob.habx.px.GoF$sims.list$chi2.sim > coy.bob.habx.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (coy.bob.habx.px.GoF_FT <- mean(coy.bob.habx.px.GoF$sims.list$FT.sims > coy.bob.habx.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(coy.bob.habx.px.GoF, file = "./Outputs/coybob_habX_px_GoF.RData")
  
  

  
  #'  ----------------------
  ####  Bear-Coyote Models  ####
  #'  ----------------------
  inits.bear.coy <- function(){list(z = zinits[[7]], mean.psiSpp1 = runif(1, 0, 1),
                                    mean.psiSpp2 = runif(1, 0, 1), mean.pSpp1 = runif(1, 0, 1),
                                    mean.pSpp2 = runif(1, 0, 1))}
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/01.1_JAGS_null_psi(yr)_p(.).R")
  start.time = Sys.time()
  bear.coy.null <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/01.1_JAGS_null_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.null$summary[1:15,])
  print(bear.coy.null$DIC)
  which(bear.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.null$samples)
  save(bear.coy.null, file = "./Outputs/bearcoy_null.RData")
  
  #####  Habitat model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  bear.coy.hab <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/02.1_JAGS_hab_psi(setup_habitat_yr)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.hab$summary[1:30,])
  print(bear.coy.hab$DIC)
  which(bear.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.hab$samples)
  save(bear.coy.hab, file = "./Outputs/bearcoy_hab.RData")
  
  #####  Prey abundance model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness, elk, wtd, lagomorphs; p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob_bearcoy.R")
  start.time = Sys.time()
  bear.coy.preyabund <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/03.4_JAGS_preyabund_psi(setup_preyabund_yr)_p(setup_effort)_lionbob_bearcoy.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.preyabund$summary[1:35,])
  print(bear.coy.preyabund$DIC)
  which(bear.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.preyabund$samples)
  save(bear.coy.preyabund, file = "./Outputs/bearcoy_preyabund.RData")
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  bear.coy.habx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.habx$summary[1:30,])
  print(bear.coy.habx$DIC)
  which(bear.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.habx$samples)
  save(bear.coy.habx, file = "./Outputs/bearcoy_habX.RData")
  
  #####  Prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, ruggedness; psix(elk, wtd, lagomorphs); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/05.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob_bearcoy.R")
  start.time = Sys.time()
  bear.coy.preyabundx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/05.4_JAGS_preyabundX_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob_bearcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.preyabundx$summary[1:40,])
  print(bear.coy.preyabundx$DIC)
  which(bear.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.preyabundx$samples)
  save(bear.coy.preyabundx, file = "./Outputs/bearcoy_preyabundX.RData")
  
  #####  Top model w/ intx on detection model  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Habitat w/ interaction
  source("./Scripts/MultiSpp_OccMod/JAGS code/04.1.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x).R")   
  start.time = Sys.time()
  bear.coy.habx.px <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/04.1.1_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_p(x).txt",    
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.habx.px$summary[1:30,])
  print(bear.coy.habx.px$DIC)
  which(bear.coy.habx.px$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.habx.px$samples)
  (bear.coy.habx.px_X2 <- mean(bear.coy.habx.px$sims.list$chi2.sim > bear.coy.habx.px$sims.list$chi2.obs)) # Bayesian p-value GOF
  (bear.coy.habx.px_FT <- mean(bear.coy.habx.px$sims.list$ft.sims > bear.coy.habx.px$sims.list$ft.obs)) # Bayesian p-value GOF
  save(bear.coy.habx.px, file = "./Outputs/bearcoy_habX_px.Rdata") 
  
  #'  -------------------------
  #####  Goodness-of-Fit tests #####
  #'  -------------------------
  ######  Top model WITHOUT co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  bear.coy.habx.GoF <- jags(bundled_data_list[[7]], inits = inits.bear.coy, params,
                           "./Outputs/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (bear.coy.habx.GoF_X2 <- mean(bear.coy.habx.GoF$sims.list$chi2.sim > bear.coy.habx.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (bear.coy.habx.GoF_FT <- mean(bear.coy.habx.GoF$sims.list$FT.sims > bear.coy.habx.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(bear.coy.habx.GoF, file = "./Outputs/bearcoy_habX_GoF.RData")
  
  ######  Top model WITH co-detection  ######
  source("./Scripts/Sourced_Scripts__Multispecies_OccMod/04.1.2_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_GoF.R")
  start.time = Sys.time()
  bear.coy.habx.px.GoF <- jags(bundled_data_list[[7]], inits = inits.bear.coy, params,
                            "./Outputs/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  (bear.coy.habx.px.GoF_X2 <- mean(bear.coy.habx.px.GoF$sims.list$chi2.sim > bear.coy.habx.px.GoF$sims.list$chi2.obs)) # Bayesian p-value GOF
  (bear.coy.habx.px.GoF_FT <- mean(bear.coy.habx.px.GoF$sims.list$FT.sims > bear.coy.habx.px.GoF$sims.list$FT.obs)) # Bayesian p-value GOF
  save(bear.coy.habx.px.GoF, file = "./Outputs/bearcoy_habX_px.GoF.RData")
  
  
  
  
  