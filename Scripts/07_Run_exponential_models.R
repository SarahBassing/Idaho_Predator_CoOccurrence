  #'  ----------------------------------------
  #'  Competition & prey effects on wait time
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Script loads time-between-detections and covaraite data, filters and formats
  #'  data for analysis, and sources JAGS code to run exponential generalized 
  #'  linear models (GLMs) to estimate effect prey availability and competition 
  #'  on wait times. JAGS code for each model is sourced from a series of scripts 
  #'  in a folder named Sourced_Scripts__Wait_Times. Nine models are fit to the 
  #'  predator - predator wait time data for each predator species including:
  #'    1. Null model
  #'    2. Species ID
  #'    3. Prey relative abundance
  #'    4. Species ID + prey relative abundance
  #'    5. Species ID x prey relative abundance
  #'    6. Prey diversity
  #'    7. Species ID + prey diversity
  #'    8. Species ID x prey diversity
  #'    9. Global model
  #'  Five additional models are fit to the prey - predator wait time data for
  #'  each predator species including:
  #'    1. Null model
  #'    2. Species ID
  #'    3. Prey relative abundance
  #'    6. Prey diversity
  #'    10. Global model
  #'  Prey abundance and global models have unique scripts for the different 
  #'  predator species depending on each predator's primary prey (e.g., script 
  #'  numbering follows 3.1, 3.2, 3.3 for the different prey abundance 
  #'  models). Species ID varied by model as well, where different predators
  #'  were included in the predator - predator wait time models and prey species
  #'  were included in the prey - predator wait time models. The global models
  #'  for the predator - predator wait time models included interactions between
  #'  species ID, prey relative abundance, and prey diversity, whereas the 
  #'  prey - predator wait time models lacked interactions.
  #'  ----------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(AICcmodavg)
  library(tidyverse)
  
  #'  Read in wait times and covariate data
  load("./Data/TBD_all_pairs.RData")
  load("./Data/covs.RData")

  #'  Join wait times & covariate data together
  tbd_w_covs_20s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr20",], covs[covs$Season == "Smr20",], by = c("NewLocationID", "GMU", "Year" = "Season")) 
  tbd_w_covs_21s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr21",], covs[covs$Season == "Smr21",], by = c("NewLocationID", "GMU", "Year" = "Season")) 
  tbd_w_covs <- rbind(tbd_w_covs_20s, tbd_w_covs_21s)
  
  #'  Filter by species
  #'  Make sure no back-to-back detections of the same species snuck in
  bear <- filter(tbd_w_covs, Focal_predator == "bear_black") %>% filter(Species_pair != "bear_black-bear_black")
  bob <- filter(tbd_w_covs, Focal_predator == "bobcat") %>% filter(Species_pair != "bobcat-bobcat")
  coy <- filter(tbd_w_covs, Focal_predator == "coyote") %>% filter(Species_pair != "coyote-coyote")
  lion <- filter(tbd_w_covs, Focal_predator == "mountain_lion") %>% filter(Species_pair != "mountain_lion-mountain_lion")
  wolf <- filter(tbd_w_covs, Focal_predator == "wolf") %>% filter(Species_pair != "wolf-wolf")
  
  #'  List data sets with extreme values removed
  pred_tbd <- list(bear, bob, coy, lion, wolf)
  
  #'  Function to identify potential outliers and remove extreme observations
  #'  Cutting off at 99% quantile for some species when 100% quantile > 20days
  tbd_summary <- function(tbd, spp, quant) {
    #'  Plot frequency of time-between-detections (should look exponential)
    hist(tbd$DaysSinceLastDet, breaks = 50, main = paste("Number of days between detections for\n", spp))
    boxplot(tbd$DaysSinceLastDet, ylab = "Days", main = paste("Number of days between detections for\n", spp))
    
    #'  Review range of TBD values
    print("Quantiles of days between sequential detections of different predators")
    print(quantile(tbd$DaysSinceLastDet))
    #'  Review 90 - 100th quantiles- where are there gaps and outliers in the distribution?
    print("90th, 95th, and 99th quantiles of days between sequentail detections of predators") 
    print(quantile(tbd$DaysSinceLastDet, c(0.9, 0.95, 0.97, 0.99, 1.0)))
    
    #'  Re-plot frequency of time-btwn-detections after removing extreme values
    short_tbd <- filter(tbd, DaysSinceLastDet <= quantile(tbd$DaysSinceLastDet, c(quant)))
    hist(short_tbd$DaysSinceLastDet, breaks = 25, main = paste("Number of days between detections for\n", spp, "up to quantile =", quant)); abline(v = 20, col = "red", lty = 2)
    boxplot(short_tbd$DaysSinceLastDet, ylab = "Days", main = paste("Number of days between detections for\n", spp, "up to quantile =", quant)); abline(h = 20, col = "red", lty = 2)
    
    #'  Summary of observations with each predator species
    print("Total TBDs with each predator species")
    print(table(short_tbd$Focal_predator))
    #'  Summary of observations in each year
    print("Total TBDs for each year")
    print(table(short_tbd$Year))
    
    #'  Actually just remove any observations over 20 days long since don't expect most 
    #'  cues from previous predator to still be detectable beyond then
    short_tbd <- filter(short_tbd, DaysSinceLastDet <= 20)
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  bear_short <- tbd_summary(bear, spp = "bear", quant = 0.99) 
  bob_short <- tbd_summary(bob, spp = "bob", quant = 0.99)
  coy_short <- tbd_summary(coy, spp = "coy", quant = 0.99)
  lion_short <- tbd_summary(lion, spp = "lion", quant = 1.0)
  wolf_short <- tbd_summary(wolf, spp = "wolf", quant = 1.0)
  
  #'  List data sets with extreme values removed
  pred_tbd_short_list <- list(bear_short, bob_short, coy_short, lion_short, wolf_short)
  
  #'  Filter to focal predator species 
  focal_species <- function(tbd_dat) {
    tbd_dat <- tbd_dat %>%
      filter(Previous_Spp == "bear_black" | Previous_Spp == "bobcat" | Previous_Spp == "coyote" | 
               Previous_Spp == "mountain_lion" | Previous_Spp == "wolf")
    return(tbd_dat)
  }
  pred_tbd_short <- lapply(pred_tbd_short_list, focal_species)
  
  #'  Filter to focal predator species and primary prey species for comparison
  prey_species <- function(tbd_dat) {
    tbd_dat <- tbd_dat %>%
      filter(Previous_Spp == "elk" | Previous_Spp == "moose" | Previous_Spp == "whitetaileddeer" | Previous_Spp == "rabbit_hare") %>% 
      mutate(Previous_Spp = ifelse(Previous_Spp == "rabbit_hare", "lagomorph", Previous_Spp)) 
    return(tbd_dat)
  }
  prey_pred_tbd_short <- lapply(pred_tbd_short_list, prey_species)
  
  #' #'  Save for permutation test
  #' save(pred_tbd_short, file = "./Data/pred_tbd_short.RData")
  #' save(prey_pred_tbd_short, file = "./Data/prey_pred_tbd_short.RData")
 
  #'  Table observations with each competitor to get a feel for sample size
  print("bear"); table(pred_tbd_short[[1]]$Previous_Spp); table(prey_pred_tbd_short[[1]]$Previous_Spp)
  print("bobcat"); table(pred_tbd_short[[2]]$Previous_Spp); table(prey_pred_tbd_short[[2]]$Previous_Spp)
  print("coyote"); table(pred_tbd_short[[3]]$Previous_Spp); table(prey_pred_tbd_short[[3]]$Previous_Spp)
  print("lion"); table(pred_tbd_short[[4]]$Previous_Spp); table(prey_pred_tbd_short[[4]]$Previous_Spp)
  print("wolf"); table(pred_tbd_short[[5]]$Previous_Spp); table(prey_pred_tbd_short[[5]]$Previous_Spp)
  #'  Limited observations for some of species combos
  #'  Using coyote as indicator variable for predator - predator models because 
  #'  has the most observations per species & because don't expect most predators 
  #'  to respond strongly to recent coyote presence
  #'  Using lagomorph as indicator variable for prey - predator models because 
  #'  don't expect most predators to respond strongly to recent lagomorph presence
  
  #'  ---------------------------------------
  ####  Set up MCMC settings and run models  ####
  #'  ---------------------------------------
  #'  Function to define and bundle data
  bundle_dat_data <- function(dat, npreyspp, species_order) {
    #'  Number of observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$NewLocationID))
    #'  Number of unique previously detected species
    nspp <- length(unique(dat$Previous_Spp))
    #'  Number of primary prey species
    npp <- npreyspp
    #'  Format covariate data
    tbd_dat <- dat %>%
      transmute(cams = as.numeric(factor(NewLocationID), levels = NewLocationID), # must be 1-n (not 0-n) for nested indexing in JAGS
                SpeciesID = as.numeric(factor(Previous_Spp), levels = species_order), # must be 1-n for nested indexing in JAGS
                GMU = as.numeric(factor(GMU), levels = c("GMU10A", "GMU6", "GMU1")),
                TBD_mins = TimeSinceLastDet,
                TBD_hrs = HoursSinceLastDet,
                TBD_days = DaysSinceLastDet,
                Elev = scale(Elev), 
                PercForest = scale(PercForest),
                SppDiversity = scale(SppDiversity),
                Nelk = scale(Nelk),
                Nmoose = scale(Nmoose),
                Nwtd = scale(Nwtd),
                Nlagomorph = scale(Nlagomorph))
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 9, nrow = ntbd)
    covs[,1] <- tbd_dat$GMU
    covs[,2] <- tbd_dat$SpeciesID
    covs[,3] <- tbd_dat$Elev
    covs[,4] <- tbd_dat$PercForest
    covs[,5] <- tbd_dat$Nelk
    covs[,6] <- tbd_dat$Nmoose
    covs[,7] <- tbd_dat$Nwtd
    covs[,8] <- tbd_dat$Nlagomorph
    covs[,9] <- tbd_dat$SppDiversity
    
    #'  Generate range of covariate values to predict across
    newElk <- seq(from = min(tbd_dat$Nelk), to = max(tbd_dat$Nelk), length.out = 100)
    newMoose <- seq(from = min(tbd_dat$Nmoose), to = max(tbd_dat$Nmoose), length.out = 100)
    newWTD <- seq(from = min(tbd_dat$Nwtd), to = max(tbd_dat$Nwtd), length.out = 100)
    newBunnies <- seq(from = min(tbd_dat$Nlagomorph), to = max(tbd_dat$Nlagomorph), length.out = 100)
    newSppDiv <- seq(from = min(tbd_dat$SppDiversity), to = max(tbd_dat$SppDiversity), length.out = 100)
    newcovs <- as.matrix(cbind(newElk, newMoose, newWTD, newBunnies, newSppDiv))
    
    #'  Number of covariates
    ncovs <- ncol(covs)
    
    #'  Time between detections
    tbd <- tbd_dat$TBD_mins
    print(summary(tbd))
    hist(tbd)
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    nspp = nspp, npp = npp, site = tbd_dat$cams, newcovs = newcovs)
    return(bundled)
    
  }
  #'  Provide specific order for SpeciesID levels - will differ for each species
  #'  Order generally goes black bear, bobcat, coyote, mountain lion, wolf
  bear_bundled <- bundle_dat_data(pred_tbd_short[[1]], npreyspp = 2, species_order = c("coyote", "bobcat", "mountain_lion", "wolf")) 
  bob_bundled <- bundle_dat_data(pred_tbd_short[[2]], npreyspp = 2, species_order = c("coyote", "bear_black", "mountain_lion", "wolf"))
  coy_bundled <- bundle_dat_data(pred_tbd_short[[3]], npreyspp = 2, species_order = c("bear_black", "bobcat", "mountain_lion", "wolf"))
  lion_bundled <- bundle_dat_data(pred_tbd_short[[4]], npreyspp = 2, species_order = c("coyote", "bear_black", "bobcat", "wolf"))
  wolf_bundled <- bundle_dat_data(pred_tbd_short[[5]], npreyspp = 3, species_order = c("coyote", "bear_black", "bobcat", "mountain_lion"))
  
  bear_bundled_nontarget <- bundle_dat_data(pred_tbd_short[[1]], npreyspp = 2, species_order = c("lagomorph", "elk", "moose", "whitetaileddeer")) 
  bob_bundled_nontarget <- bundle_dat_data(pred_tbd_short[[2]], npreyspp = 2, species_order = c("lagomorph", "elk", "moose", "whitetaileddeer"))
  coy_bundled_nontarget <- bundle_dat_data(pred_tbd_short[[3]], npreyspp = 2, species_order = c("lagomorph", "elk", "moose", "whitetaileddeer"))
  lion_bundled_nontarget <- bundle_dat_data(pred_tbd_short[[4]], npreyspp = 2, species_order = c("lagomorph", "elk", "moose", "whitetaileddeer"))
  wolf_bundled_nontarget <- bundle_dat_data(pred_tbd_short[[5]], npreyspp = 3, species_order = c("lagomorph", "elk", "moose", "whitetaileddeer"))
  
  #' #'  Save for making figures later
  #' save(bear_bundled, file = "./Data/bear_bundled.RData")
  #' save(bob_bundled, file = "./Data/bob_bundled.RData")
  #' save(coy_bundled, file = "./Data/coy_bundled.RData")
  #' save(lion_bundled, file = "./Data/lion_bundled.RData")
  #' save(wolf_bundled, file = "./Data/wolf_bundled.RData")
  #' 
  #' save(bear_bundled_nontarget, file = "./Data/bear_bundled_nontarget.RData")
  #' save(bob_bundled_nontarget, file = "./Data/bob_bundled_nontarget.RData")
  #' save(coy_bundled_nontarget, file = "./Data/coy_bundled_nontarget.RData")
  #' save(lion_bundled_nontarget, file = "./Data/lion_bundled_nontarget.RData")
  #' save(wolf_bundled_nontarget, file = "./Data/wolf_bundled_nontarget.RData")

  #'  MCMC settings
  nc <- 3
  ni <- 30000
  nb <- 5000
  nt <- 10
  na <- 1000
  
  #'  Parameters to monitor
  params <- c("alpha0", "beta.sppID", "beta.prey", "beta.div", "beta.interaction", 
              "beta.interaction.elk", "beta.interaction.wtd", "beta.interaction.moose",
              "beta.interaction.lago", "mu.tbd", "spp.tbd", "spp.tbd.elk", 
              "spp.tbd.moose", "spp.tbd.wtd", "spp.tbd.lago", "spp.tbd.div") 
  
  
  #'  Call Tyra, we need our Next Top Model!
  
  
  #'  ------------------------------
  ####  COMPETITOR - BEAR Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  bear.init <- log(aggregate(bear_bundled$y, list(bear_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bear.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.null <- jags(bear_bundled, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.null$summary)
  mcmcplot(tbd.bear.null$samples)
  save(tbd.bear.null, file = "./Outputs/tbd_bear_null.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID <- jags(bear_bundled, params, './Outputs/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID$summary)
  mcmcplot(tbd.bear.sppID$samples)
  save(tbd.bear.sppID, file = "./Outputs/tbd_bear_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.1_JAGS_tbd_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.preyabund <- jags(bear_bundled, params, './Outputs/tbd_elk_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.preyabund$summary[1:5,])
  mcmcplot(tbd.bear.preyabund$samples)
  save(tbd.bear.preyabund, file = "./Outputs/tbd_bear_preyRAI.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/04.1_JAGS_tbd_sppID_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID.preyabund <- jags(bear_bundled, params, './Outputs/tbd_sppID_elk_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID.preyabund$summary[1:15,])
  mcmcplot(tbd.bear.sppID.preyabund$samples)
  save(tbd.bear.sppID.preyabund, file = "./Outputs/tbd_bear_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bobcat [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/05.1_JAGS_tbd_sppID_X_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppIDxpreyabund <- jags(bear_bundled, params, './Outputs/tbd_sppID_X_elk_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppIDxpreyabund$summary[1:30,])
  mcmcplot(tbd.bear.sppIDxpreyabund$samples)
  save(tbd.bear.sppIDxpreyabund, file = "./Outputs/tbd_bear_sppID_X_preyRAI.RData")    
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.div <- jags(bear_bundled, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.div$summary[1:5,])
  mcmcplot(tbd.bear.div$samples)
  save(tbd.bear.div, file = "./Outputs/tbd_bear_preydiversity.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/07_JAGS_tbd_sppID_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID.div <- jags(bear_bundled, params, './Outputs/tbd_sppID_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID.div$summary[1:15,])
  mcmcplot(tbd.bear.sppID.div$samples)
  save(tbd.bear.sppID.div, file = "./Outputs/tbd_bear_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/08_JAGS_tbd_sppID_X_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppIDxdiv <- jags(bear_bundled, params, './Outputs/tbd_sppID_X_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppIDxdiv$summary[1:15,])
  mcmcplot(tbd.bear.sppIDxdiv$samples)
  save(tbd.bear.sppIDxdiv, file = "./Outputs/tbd_bear_sppID_X_preydiv.RData")

  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/09.1_JAGS_tbd_global_sppID_X_div_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.global <- jags(bear_bundled, params, './Outputs/tbd_global_sppID_X_div_elk_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.global$summary[1:21,])
  mcmcplot(tbd.bear.global$samples)
  save(tbd.bear.global, file = "./Outputs/tbd_bear_global.RData")   
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  
  #'  --------------------------------
  ####  COMPETITOR - BOBCAT Analyses  ####
  #'  --------------------------------
  #'  Setup initial values
  bob.init <- log(aggregate(bob_bundled$y, list(bob_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bob.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.null <- jags(bob_bundled, params, './Outputs/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.null$summary)
  mcmcplot(tbd.bob.null$samples)
  save(tbd.bob.null, file = "./Outputs/tbd.bob_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID <- jags(bob_bundled, params, './Outputs/02_JAGS_tbd_sppID.R', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID$summary)
  mcmcplot(tbd.bob.sppID$samples)
  save(tbd.bob.sppID, file = "./Outputs/tbd.bob_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.2_JAGS_tbd_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.preyabund <- jags(bob_bundled, params, './Outputs/tbd_wtd_lago_abundance.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.preyabund$summary[1:5,])
  mcmcplot(tbd.bob.preyabund$samples)
  save(tbd.bob.preyabund, file = "./Outputs/tbd_bob_preyRAI.RData")
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/04.2_JAGS_tbd_sppID_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID.preyabund <- jags(bob_bundled, params, './Outputs/tbd_sppID_wtd_lago_abundance.txt',
                                  inits = inits, n.chains = nc, n.iter = ni,
                                  n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID.preyabund$summary[1:15,])
  mcmcplot(tbd.bob.sppID.preyabund$samples)
  save(tbd.bob.sppID.preyabund, file = "./Outputs/tbd_bob_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/05.2_JAGS_tbd_sppID_X_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppIDxpreyabund <- jags(bob_bundled, params, './Outputs/tbd_sppID_X_wtd_lago_abundance.txt',
                                  inits = inits, n.chains = nc, n.iter = ni,
                                  n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.bob.sppIDxpreyabund$samples)
  save(tbd.bob.sppIDxpreyabund, file = "./Outputs/tbd_bob_sppID_X_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.div <- jags(bob_bundled, params, './Outputs/tbd_preydiversity.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.div$summary[1:5,])
  mcmcplot(tbd.bob.div$samples)
  save(tbd.bob.div, file = "./Outputs/tbd_bob_preydiversity.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/07_JAGS_tbd_sppID_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID.div <- jags(bob_bundled, params, './Outputs/tbd_sppID_preydiversity.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID.div$summary[1:15,])
  mcmcplot(tbd.bob.sppID.div$samples)
  save(tbd.bob.sppID.div, file = "./Outputs/tbd_bob_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/08_JAGS_tbd_sppID_X_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppIDxdiv <- jags(bob_bundled, params, './Outputs/tbd_sppID_X_preydiversity.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppIDxdiv$summary[1:18,])
  mcmcplot(tbd.bob.sppIDxdiv$samples)
  save(tbd.bob.sppIDxdiv, file = "./Outputs/tbd_bob_sppID_X_preydiv.RData")    

  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/09.2_JAGS_tbd_global_sppID_X_div_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.global <- jags(bob_bundled, params, './Outputs/tbd_global_sppID_X_div_wtd_lago_abundance.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.global$summary[1:21,])
  mcmcplot(tbd.bob.global$samples)
  save(tbd.bob.global, file = "./Outputs/tbd_bob_global.RData")   
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  
  #'  -------------------
  ####  COMPETITOR - COYOTE Analyses  ####
  #'  -------------------
  #'  Setup initial values
  coy.init <- log(aggregate(coy_bundled$y, list(coy_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = coy.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.null <- jags(coy_bundled, params, './Outputs/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.null$summary)
  mcmcplot(tbd.coy.null$samples)
  save(tbd.coy.null, file = "./Outputs/tbd.coy_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID <- jags(coy_bundled, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID$summary)
  mcmcplot(tbd.coy.sppID$samples)
  save(tbd.coy.sppID, file = "./Outputs/tbd.coy_sppID.RData") 
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.2_JAGS_tbd_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.preyabund <- jags(coy_bundled, params, './Outputs/tbd_wtd_lago_abundance.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.preyabund$summary[1:5,])
  mcmcplot(tbd.coy.preyabund$samples)
  save(tbd.coy.preyabund, file = "./Outputs/tbd_coy_preyRAI.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/04.2_JAGS_tbd_sppID_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID.preyabund <- jags(coy_bundled, params, './Outputs/tbd_sppID_wtd_lago_abundance.txt',
                                  inits = inits, n.chains = nc, n.iter = ni,
                                  n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID.preyabund$summary[1:15,])
  mcmcplot(tbd.coy.sppID.preyabund$samples)
  save(tbd.coy.sppID.preyabund, file = "./Outputs/tbd_coy_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/05.2_JAGS_tbd_sppID_X_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppIDxpreyabund <- jags(coy_bundled, params, './Outputs/tbd_sppID_X_wtd_lago_abundance.txt',
                                  inits = inits, n.chains = nc, n.iter = ni,
                                  n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.coy.sppIDxpreyabund$samples)
  save(tbd.coy.sppIDxpreyabund, file = "./Outputs/tbd_coy_sppID_X_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.div <- jags(coy_bundled, params, './Outputs/tbd_preydiversity.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.div$summary[1:5,])
  mcmcplot(tbd.coy.div$samples)
  save(tbd.coy.div, file = "./Outputs/tbd_coy_preydiversity.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/07_JAGS_tbd_sppID_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID.div <- jags(coy_bundled, params, './Outputs/tbd_sppID_preydiversity.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID.div$summary[1:15,])
  mcmcplot(tbd.coy.sppID.div$samples)
  save(tbd.coy.sppID.div, file = "./Outputs/tbd_coy_sppID_preydiv.RData") 
  
  #####  Species ID * prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/08_JAGS_tbd_sppID_X_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppIDxdiv <- jags(coy_bundled, params, './Outputs/tbd_sppID_X_preydiversity.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppIDxdiv$summary[1:18,])
  mcmcplot(tbd.coy.sppIDxdiv$samples)
  save(tbd.coy.sppIDxdiv, file = "./Outputs/tbd_coy_sppID_X_preydiv.RData")

  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/09.2_JAGS_tbd_global_sppID_X_div_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.global <- jags(coy_bundled, params, './Outputs/tbd_global_sppID_X_div_wtd_lago_abundance.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.global$summary[1:21,])
  mcmcplot(tbd.coy.global$samples)
  save(tbd.coy.global, file = "./Outputs/tbd_coy_global.RData") 
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  
  #'  ---------------------------------------
  ####  COMPETITOR - MOUNTAIN LION Analyses  ####
  #'  ---------------------------------------
  #'  Setup initial values
  lion.init <- log(aggregate(lion_bundled$y, list(lion_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = lion.init)} 
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.null <- jags(lion_bundled, params, './Outputs/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.null$summary)
  mcmcplot(tbd.lion.null$samples)
  save(tbd.lion.null, file = "./Outputs/tbd.lion_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID <- jags(lion_bundled, params, './Outputs/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID$summary)
  mcmcplot(tbd.lion.sppID$samples)
  save(tbd.lion.sppID, file = "./Outputs/tbd.lion_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.1_JAGS_tbd_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.preyabund <- jags(lion_bundled, params, './Outputs/tbd_elk_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.preyabund$summary[1:5,])
  mcmcplot(tbd.lion.preyabund$samples)
  save(tbd.lion.preyabund, file = "./Outputs/tbd_lion_preyRAI.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/04.1_JAGS_tbd_sppID_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID.preyabund <- jags(lion_bundled, params, './Outputs/tbd_sppID_elk_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID.preyabund$summary[1:15,])
  mcmcplot(tbd.lion.sppID.preyabund$samples)
  save(tbd.lion.sppID.preyabund, file = "./Outputs/tbd_lion_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Species ID * prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/05.1_JAGS_tbd_sppID_X_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppIDxpreyabund <- jags(lion_bundled, params, './Outputs/tbd_sppID_X_elk_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.lion.sppIDxpreyabund$samples)
  save(tbd.lion.sppIDxpreyabund, file = "./Outputs/tbd_lion_sppID_X_preyRAI.RData")     
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.div <- jags(lion_bundled, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.div$summary[1:5,])
  mcmcplot(tbd.lion.div$samples)
  save(tbd.lion.div, file = "./Outputs/tbd_lion_preydiversity.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/07_JAGS_tbd_sppID_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID.div <- jags(lion_bundled, params, './Outputs/tbd_sppID_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID.div$summary[1:15,])
  mcmcplot(tbd.lion.sppID.div$samples)
  save(tbd.lion.sppID.div, file = "./Outputs/tbd_lion_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/08_JAGS_tbd_sppID_X_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppIDxdiv <- jags(lion_bundled, params, './Outputs/tbd_sppID_X_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppIDxdiv$summary[1:18,])
  mcmcplot(tbd.lion.sppIDxdiv$samples)
  save(tbd.lion.sppIDxdiv, file = "./Outputs/tbd_lion_sppID_X_preydiv.RData")    
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/09.1_JAGS_tbd_global_sppID_X_div_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.global <- jags(lion_bundled, params, './Outputs/tbd_global_sppID_X_div_elk_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.global$summary[1:21,])
  mcmcplot(tbd.lion.global$samples)
  save(tbd.lion.global, file = "./Outputs/tbd_lion_global.RData")      
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  
  #'  ------------------------------
  ####  COMPETITOR - WOLF Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  wolf.init <- log(aggregate(wolf_bundled$y, list(wolf_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = wolf.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.null <- jags(wolf_bundled, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.null$summary)
  mcmcplot(tbd.wolf.null$samples)
  save(tbd.wolf.null, file = "./Outputs/tbd_wolf_null.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID <- jags(wolf_bundled, params, './Outputs/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID$summary)
  mcmcplot(tbd.wolf.sppID$samples)
  save(tbd.wolf.sppID, file = "./Outputs/tbd_wolf_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.3_JAGS_tbd_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.preyabund <- jags(wolf_bundled, params, './Outputs/tbd_elk_moose_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.preyabund$summary[1:5,])
  mcmcplot(tbd.wolf.preyabund$samples)
  save(tbd.wolf.preyabund, file = "./Outputs/tbd_wolf_preyRAI.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/04.3_JAGS_tbd_sppID_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID.preyabund <- jags(wolf_bundled, params, './Outputs/tbd_sppID_elk_moose_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID.preyabund$summary[1:25,])
  mcmcplot(tbd.wolf.sppID.preyabund$samples)
  save(tbd.wolf.sppID.preyabund, file = "./Outputs/tbd_wolf_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/05.3_JAGS_tbd_sppID_X_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppIDxpreyabund <- jags(wolf_bundled, params, './Outputs/tbd_sppID_X_elk_moose_wtd_abundance.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.wolf.sppIDxpreyabund$samples)
  save(tbd.wolf.sppIDxpreyabund, file = "./Outputs/tbd_wolf_sppID_X_preyRAI.RData")      
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.div <- jags(wolf_bundled, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.div$summary[1:5,])
  mcmcplot(tbd.wolf.div$samples)
  save(tbd.wolf.div, file = "./Outputs/tbd_wolf_preydiversity.RData") 
 
  #####  Competitor + prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/07_JAGS_tbd_sppID_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID.div <- jags(wolf_bundled, params, './Outputs/tbd_sppID_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID.div$summary[1:15,])
  mcmcplot(tbd.wolf.sppID.div$samples)
  save(tbd.wolf.sppID.div, file = "./Outputs/tbd_wolf_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/08_JAGS_tbd_sppID_X_preydiversity.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppIDxdiv <- jags(wolf_bundled, params, './Outputs/tbd_sppID_X_preydiversity.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppIDxdiv$summary[1:15,])
  mcmcplot(tbd.wolf.sppIDxdiv$samples)
  save(tbd.wolf.sppIDxdiv, file = "./Outputs/tbd_wolf_sppID_X_preydiv.RData")

  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/09.3_JAGS_tbd_global_sppID_X_div_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.global <- jags(wolf_bundled, params, './Outputs/tbd_global_sppID_X_div_elk_moose_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.global$summary[1:30,])
  mcmcplot(tbd.wolf.global$samples)
  save(tbd.wolf.global, file = "./Outputs/tbd_wolf_global.RData")    
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  
  #'  ------------------------------
  ####  NON-TARGET - BEAR Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  bear.nt.init <- log(aggregate(bear_bundled_nontarget$y, list(bear_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bear.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.null <- jags(bear_bundled_nontarget, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.null$summary)
  mcmcplot(tbd.nt.bear.null$samples)
  save(tbd.nt.bear.null, file = "./Outputs/tbd_nontarget_bear_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.sppID <- jags(bear_bundled_nontarget, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.sppID$summary)
  mcmcplot(tbd.nt.bear.sppID$samples)
  save(tbd.nt.bear.sppID, file = "./Outputs/tbd_nontarget_bear_sppID.RData") 
  #'  Keep in mind SppID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.1_JAGS_tbd_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.preyabund <- jags(bear_bundled_nontarget, params, './Outputs/tbd_elk_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.preyabund$summary[1:5,])
  mcmcplot(tbd.nt.bear.preyabund$samples)
  save(tbd.nt.bear.preyabund, file = "./Outputs/tbd_nontarget_bear_preyRAI.RData")
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.div <- jags(bear_bundled_nontarget, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.div$summary[1:5,])
  mcmcplot(tbd.nt.bear.div$samples)
  save(tbd.nt.bear.div, file = "./Outputs/tbd_nontarget_bear_preydiversity.RData") 
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/10.1_JAGS_tbd_global_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.global <- jags(bear_bundled_nontarget, params, './Outputs/tbd_global_elk_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.global$summary[1:21,])
  mcmcplot(tbd.nt.bear.global$samples)
  save(tbd.nt.bear.global, file = "./Outputs/tbd_nontarget_bear_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  
  #'  --------------------------------
  ####  NON-TARGET - BOBCAT Analyses  ####
  #'  --------------------------------
  #'  Setup initial values
  bob.nt.init <- log(aggregate(bob_bundled_nontarget$y, list(bob_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bob.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.null <- jags(bob_bundled_nontarget, params, './Outputs/tbd_null.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.null$summary)
  mcmcplot(tbd.nt.bob.null$samples)
  save(tbd.nt.bob.null, file = "./Outputs/tbd_nontarget_bob_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.sppID <- jags(bob_bundled_nontarget, params, './Outputs/tbd_sppID.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                        n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.sppID$summary)
  mcmcplot(tbd.nt.bob.sppID$samples)
  save(tbd.nt.bob.sppID, file = "./Outputs/tbd_nontarget_bob_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.2_JAGS_tbd_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.preyabund <- jags(bob_bundled_nontarget, params, './Outputs/tbd_wtd_lago_abundance.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.preyabund$summary[1:5,])
  mcmcplot(tbd.nt.bob.preyabund$samples)
  save(tbd.nt.bob.preyabund, file = "./Outputs/tbd_nontarget_bob_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.div <- jags(bob_bundled_nontarget, params, './Outputs/tbd_preydiversity.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.div$summary[1:5,])
  mcmcplot(tbd.nt.bob.div$samples)
  save(tbd.nt.bob.div, file = "./Outputs/tbd_nontarget_bob_preydiversity.RData") 
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/10.2_JAGS_tbd_global_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.global <- jags(bob_bundled_nontarget, params, './Outputs/tbd_global_wtd_lago_abundance.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.global$summary[1:21,])
  mcmcplot(tbd.nt.bob.global$samples)
  save(tbd.nt.bob.global, file = "./Outputs/tbd_nontarget_bob_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  
  #'  --------------------------------
  ####  NON-TARGET - COYOTE Analyses  ####
  #'  --------------------------------
  #'  Setup initial values
  coy.nt.init <- log(aggregate(coy_bundled_nontarget$y, list(coy_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = coy.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.null <- jags(coy_bundled_nontarget, params, './Outputs/tbd_null.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.null$summary)
  mcmcplot(tbd.nt.coy.null$samples)
  save(tbd.nt.coy.null, file = "./Outputs/tbd_nontarget_coy_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.sppID <- jags(coy_bundled_nontarget, params, './Outputs/tbd_sppID.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                        n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.sppID$summary)
  mcmcplot(tbd.nt.coy.sppID$samples)
  save(tbd.nt.coy.sppID, file = "./Outputs/tbd_nontarget_coy_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.2_JAGS_tbd_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.preyabund <- jags(coy_bundled_nontarget, params, './Outputs/tbd_wtd_lago_abundance.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.preyabund$summary[1:5,])
  mcmcplot(tbd.nt.coy.preyabund$samples)
  save(tbd.nt.coy.preyabund, file = "./Outputs/tbd_nontarget_coy_preyRAI.RData") 
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.div <- jags(coy_bundled_nontarget, params, './Outputs/tbd_preydiversity.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.div$summary[1:5,])
  mcmcplot(tbd.nt.coy.div$samples)
  save(tbd.nt.coy.div, file = "./Outputs/tbd_nontarget_coy_preydiversity.RData") 
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/10.2_JAGS_tbd_global_wtd_lago_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.global <- jags(coy_bundled_nontarget, params, './Outputs/tbd_global_wtd_lago_abundance.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.global$summary[1:21,])
  mcmcplot(tbd.nt.coy.global$samples)
  save(tbd.nt.coy.global, file = "./Outputs/tbd_nontarget_coy_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  
  #'  ---------------------------------------
  ####  NON-TARGET - MOUNTAIN LION Analyses  ####
  #'  ---------------------------------------
  #'  Setup initial values
  lion.nt.init <- log(aggregate(lion_bundled_nontarget$y, list(lion_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = lion.nt.init)} 
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.null <- jags(lion_bundled_nontarget, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.null$summary)
  mcmcplot(tbd.nt.lion.null$samples)
  save(tbd.nt.lion.null, file = "./Outputs/tbd_nontarget_lion_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.sppID <- jags(lion_bundled_nontarget, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.sppID$summary)
  mcmcplot(tbd.nt.lion.sppID$samples)
  save(tbd.nt.lion.sppID, file = "./Outputs/tbd_nontarget_lion_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.1_JAGS_tbd_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.preyabund <- jags(lion_bundled_nontarget, params, './Outputs/tbd_elk_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.preyabund$summary[1:5,])
  mcmcplot(tbd.nt.lion.preyabund$samples)
  save(tbd.nt.lion.preyabund, file = "./Outputs/tbd_nontarget_lion_preyRAI.RData")
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.div <- jags(lion_bundled_nontarget, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.div$summary[1:5,])
  mcmcplot(tbd.nt.lion.div$samples)
  save(tbd.nt.lion.div, file = "./Outputs/tbd_nontarget_lion_preydiversity.RData") 
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/10.1_JAGS_tbd_global_elk_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.global <- jags(lion_bundled_nontarget, params, './Outputs/tbd_global_elk_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.global$summary[1:21,])
  mcmcplot(tbd.nt.lion.global$samples)
  save(tbd.nt.lion.global, file = "./Outputs/tbd_nontarget_lion_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  
  #'  ------------------------------
  ####  NON-TARGET - WOLF Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  wolf.nt.init <- log(aggregate(wolf_bundled_nontarget$y, list(wolf_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = wolf.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/01_JAGS_tbd_null.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.null <- jags(wolf_bundled_nontarget, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.null$summary)
  mcmcplot(tbd.nt.wolf.null$samples)
  save(tbd.nt.wolf.null, file = "./Outputs/tbd_nontarget_wolf_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.sppID <- jags(wolf_bundled_nontarget, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.sppID$summary)
  mcmcplot(tbd.nt.wolf.sppID$samples)
  save(tbd.nt.wolf.sppID, file = "./Outputs/tbd_nontarget_wolf_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]
  
  #####  Prey relative abundance model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/03.3_JAGS_tbd_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.preyabund <- jags(wolf_bundled_nontarget, params, './Outputs/tbd_elk_moose_wtd_abundance.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.preyabund$summary[1:5,])
  mcmcplot(tbd.nt.wolf.preyabund$samples)
  save(tbd.nt.wolf.preyabund, file = "./Outputs/tbd_nontarget_wolf_preyRAI.RData") 
  
  #####  Prey diversity model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/06_JAGS_tbd_preydiversity.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.div <- jags(wolf_bundled_nontarget, params, './Outputs/tbd_preydiversity.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.div$summary[1:5,])
  mcmcplot(tbd.nt.wolf.div$samples)
  save(tbd.nt.wolf.div, file = "./Outputs/tbd_nontarget_wolf_preydiversity.RData") 
  
  #####  Global model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/10.3_JAGS_tbd_global_elk_moose_wtd_abundance.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.global <- jags(wolf_bundled_nontarget, params, './Outputs/tbd_global_elk_moose_wtd_abundance.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.global$summary[1:30,])
  mcmcplot(tbd.nt.wolf.global$samples)
  save(tbd.nt.wolf.global, file = "./Outputs/tbd_nontarget_wolf_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], elk[2], moose[3], white-tailed deer[4]

  
  #'  Fin
  #'  Next stop, 08_DIC_exponential_model_selection.R for model selection
  
  
