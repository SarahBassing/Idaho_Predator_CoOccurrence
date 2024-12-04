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
  #'  Five additional models are fit to the prey - predator wait time data for
  #'  each predator species including:
  #'    1. Null model
  #'    2. Species ID
  #'    3. Prey relative abundance
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
  # load("./Data/covs.RData")
  covs_smr20 <- read.csv("./Data/Covariates_Smr20.csv")
  covs_smr21 <- read.csv("./Data/Covariates_Smr21.csv")

  #'  Join wait times & covariate data together
  # tbd_w_covs_20s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr20",], covs[covs$Season == "Smr20",], by = c("NewLocationID", "GMU", "Year" = "Season")) 
  # tbd_w_covs_21s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr21",], covs[covs$Season == "Smr21",], by = c("NewLocationID", "GMU", "Year" = "Season")) 
  tbd_w_covs_20s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr20",], covs_smr20, by = c("NewLocationID", "GMU")) 
  tbd_w_covs_21s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr21",], covs_smr21, by = c("NewLocationID", "GMU")) 
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
    
    #'  Actually just remove any observations 7-days or longer since co-occurrence model
    #'  considered detections at summer and weekly time scales - want something finer
    #'  than the weekly time scale for this analysis
    short_tbd <- filter(short_tbd, DaysSinceLastDet < 7)
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  bear_short <- tbd_summary(bear, spp = "bear", quant = 0.99) 
  bob_short <- tbd_summary(bob, spp = "bob", quant = 0.99)
  coy_short <- tbd_summary(coy, spp = "coy", quant = 0.99)
  lion_short <- tbd_summary(lion, spp = "lion", quant = 0.99)
  wolf_short <- tbd_summary(wolf, spp = "wolf", quant = 0.99)
  
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
  #'  has the most observations per species
  
  #'  -------------
  ####  Visualize  ####
  #'  -------------
  #'  Plot histograms of raw data
  pred_tbd_short_df <- rbind(bear_short, bob_short, coy_short, lion_short, wolf_short) %>%
    mutate(Focal_predator = ifelse(Focal_predator == "bear_black", "Black bear", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "bobcat", "Bobcat", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "coyote", "Coyote", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "mountain_lion", "Mountain lion", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "wolf", "Wolf", Focal_predator))
  #'  Function to create plot for each predator species
  tbd_hist <- function(spp, pred_color) {
    dat <- filter(pred_tbd_short_df, Focal_predator == spp)
    plot_hist <- ggplot(dat, aes(x = HoursSinceLastDet)) + 
      geom_histogram(binwidth = 10, color = "black", fill = pred_color) + 
      theme_bw() + 
      xlab("Wait time (hours)") + 
      ylab("Frequency") +
      ggtitle(spp)
    print(plot_hist)
    return(plot_hist)
  }
  #'  List species and color associated with each one
  spp <- list("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")
  pred_color <- list("#98CAE1", "#A50026", "#DD3D2D", "#FDB366", "#364B9A")
  #'  Apply function to dataset
  histograms <- mapply(tbd_hist, spp, pred_color, SIMPLIFY = FALSE)
  
  #'  Create figure for publication
  (tbd_histogram_fig <- histograms[[1]] + histograms[[2]] + theme(axis.title.y = element_blank()) + 
    histograms[[3]] + theme(axis.title.y = element_blank()) + 
    histograms[[4]] + histograms[[5]]  + theme(axis.title.y = element_blank()) + 
    plot_layout(ncol = 3) + plot_annotation(title = "Predator-specific wait times following detection of a different species",
                                            tag_levels = 'a'))
  #' #'  Save
  #' ggsave("./Outputs/Figures/TBD_histograms_rawdata.tiff", tbd_histogram_fig, 
  #'        units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
 
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
      transmute(cams = as.numeric(factor(NewLocationID), levels = NewLocationID), # must be 1-n (not 0-n) for nested indexing 
                SpeciesID = as.numeric(factor(Previous_Spp), levels = species_order), # must be 1-n for nested indexing
                GMU = as.numeric(factor(GMU), levels = c("GMU10A", "GMU6", "GMU1")),
                TBD_mins = TimeSinceLastDet,
                TBD_hrs = HoursSinceLastDet,
                TBD_days = DaysSinceLastDet,
                Elev = scale(Elevation__10m2), 
                TRI = scale(TRI),
                PercForest = scale(perc_forest),
                Nelk = scale(elk_perday),  
                Nmoose = scale(moose_perday),
                Nwtd = scale(whitetaileddeer_perday),
                Nlagomorph = scale(lagomorphs_perday))
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
    covs[,9] <- tbd_dat$TRI
    
    #'  Generate range of covariate values to predict across
    newElk <- seq(from = min(tbd_dat$Nelk), to = max(tbd_dat$Nelk), length.out = 100)
    newMoose <- seq(from = min(tbd_dat$Nmoose), to = max(tbd_dat$Nmoose), length.out = 100)
    newWTD <- seq(from = min(tbd_dat$Nwtd), to = max(tbd_dat$Nwtd), length.out = 100)
    newBunnies <- seq(from = min(tbd_dat$Nlagomorph), to = max(tbd_dat$Nlagomorph), length.out = 100)
    newcovs <- as.matrix(cbind(newElk, newMoose, newWTD, newBunnies))
    
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
  
  #'  Filter to only focal prey species for nontarget analysis
  bear_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[1]][prey_pred_tbd_short[[1]]$Previous_Spp == "elk" | prey_pred_tbd_short[[1]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("elk", "whitetaileddeer")) 
  bob_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[2]][prey_pred_tbd_short[[2]]$Previous_Spp == "lagomorph" | prey_pred_tbd_short[[2]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("lagomorph", "whitetaileddeer"))
  coy_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[3]][prey_pred_tbd_short[[3]]$Previous_Spp == "lagomorph" | prey_pred_tbd_short[[3]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("lagomorph", "whitetaileddeer"))
  lion_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[4]][prey_pred_tbd_short[[4]]$Previous_Spp == "elk" | prey_pred_tbd_short[[4]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("elk", "whitetaileddeer"))
  wolf_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[5]][prey_pred_tbd_short[[5]]$Previous_Spp != "lagomorph",], npreyspp = 3, species_order = c("elk", "moose", "whitetaileddeer"))
  
  #' #'  Save for making figures later
  #' save(bear_bundled, file = "./Data/bear_bundled.RData")
  #' save(bob_bundled, file = "./Data/bob_bundled.RData")
  #' save(coy_bundled, file = "./Data/coy_bundled.RData")
  #' save(lion_bundled, file = "./Data/lion_bundled.RData")
  #' save(wolf_bundled, file = "./Data/wolf_bundled.RData")
  #' 
  # save(bear_bundled_nontarget, file = "./Data/bear_bundled_nontarget.RData")
  # save(bob_bundled_nontarget, file = "./Data/bob_bundled_nontarget.RData")
  # save(coy_bundled_nontarget, file = "./Data/coy_bundled_nontarget.RData")
  # save(lion_bundled_nontarget, file = "./Data/lion_bundled_nontarget.RData")
  # save(wolf_bundled_nontarget, file = "./Data/wolf_bundled_nontarget.RData")

  #'  MCMC settings
  nc <- 3
  ni <- 30000
  nb <- 5000
  nt <- 10
  na <- 1000
  
  #'  Parameters to monitor
  params <- c("alpha0", "beta.sppID", "beta.prey", "beta.interaction", 
              "beta.interaction.elk", "beta.interaction.wtd", "beta.interaction.moose",
              "beta.interaction.lago", "mu.tbd", "spp.tbd", "spp.tbd.elk", 
              "spp.tbd.moose", "spp.tbd.wtd", "spp.tbd.lago", "chi2.obs", "chi2.sim") 
  
  
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
  (tbd.bear.null.pval <- mean(tbd.bear.null$sims.list$chi2.sim > tbd.bear.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bear.sppID.pval <- mean(tbd.bear.sppID$sims.list$chi2.sim > tbd.bear.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bear.preyabund.pval <- mean(tbd.bear.preyabund$sims.list$chi2.sim > tbd.bear.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bear.sppID.preyabund.pval <- mean(tbd.bear.sppID.preyabund$sims.list$chi2.sim > tbd.bear.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bear.sppIDxpreyabund.pval <- mean(tbd.bear.sppIDxpreyabund$sims.list$chi2.sim > tbd.bear.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppIDxpreyabund$samples)
  save(tbd.bear.sppIDxpreyabund, file = "./Outputs/tbd_bear_sppID_X_preyRAI.RData")    
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
  tbd.bob.null <- jags(bob_bundled, params, './Outputs/tbd_null.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.null$summary)
  (tbd.bob.null.pval <- mean(tbd.bob.null$sims.list$chi2.sim > tbd.bob.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.null$samples)
  save(tbd.bob.null, file = "./Outputs/tbd.bob_null.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID <- jags(bob_bundled, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID$summary)
  (tbd.bob.sppID.pval <- mean(tbd.bob.sppID$sims.list$chi2.sim > tbd.bob.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bob.preyabund.pval <- mean(tbd.bob.preyabund$sims.list$chi2.sim > tbd.bob.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bob.sppID.preyabund.pval <- mean(tbd.bob.sppID.preyabund$sims.list$chi2.sim > tbd.bob.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.bob.sppIDxpreyabund.pval <- mean(tbd.bob.sppIDxpreyabund$sims.list$chi2.sim > tbd.bob.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppIDxpreyabund$samples)
  save(tbd.bob.sppIDxpreyabund, file = "./Outputs/tbd_bob_sppID_X_preyRAI.RData")
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
  tbd.coy.null <- jags(coy_bundled, params, './Outputs/tbd_null.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.null$summary)
  (tbd.coy.null.pval <- mean(tbd.coy.null$sims.list$chi2.sim > tbd.coy.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.null$samples)
  save(tbd.coy.null, file = "./Outputs/tbd.coy_null.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID <- jags(coy_bundled, params, './Outputs/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID$summary)
  (tbd.coy.sppID.pval <- mean(tbd.coy.sppID$sims.list$chi2.sim > tbd.coy.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.coy.preyabund.pval <- mean(tbd.coy.preyabund$sims.list$chi2.sim > tbd.coy.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.coy.sppID.preyabund.pval <- mean(tbd.coy.sppID.preyabund$sims.list$chi2.sim > tbd.coy.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.coy.sppIDxpreyabund.pval <- mean(tbd.coy.sppIDxpreyabund$sims.list$chi2.sim > tbd.coy.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppIDxpreyabund$samples)
  save(tbd.coy.sppIDxpreyabund, file = "./Outputs/tbd_coy_sppID_X_preyRAI.RData")
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
  tbd.lion.null <- jags(lion_bundled, params, './Outputs/tbd_null.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.null$summary)
  (tbd.lion.null.pval <- mean(tbd.lion.null$sims.list$chi2.sim > tbd.lion.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.null$samples)
  save(tbd.lion.null, file = "./Outputs/tbd.lion_null.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/Sourced_Scripts__Wait_Times/02_JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID <- jags(lion_bundled, params, './Outputs/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID$summary)
  (tbd.lion.sppID.pval <- mean(tbd.lion.sppID$sims.list$chi2.sim > tbd.lion.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.lion.preyabund.pval <- mean(tbd.lion.preyabund$sims.list$chi2.sim > tbd.lion.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.lion.sppID.preyabund.pval <- mean(tbd.lion.sppID.preyabund$sims.list$chi2.sim > tbd.lion.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.lion.sppIDxpreyabund.pval <- mean(tbd.lion.sppIDxpreyabund$sims.list$chi2.sim > tbd.lion.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppIDxpreyabund$samples)
  save(tbd.lion.sppIDxpreyabund, file = "./Outputs/tbd_lion_sppID_X_preyRAI.RData")     
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
  (tbd.wolf.null.pval <- mean(tbd.wolf.null$sims.list$chi2.sim > tbd.wolf.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.wolf.sppID.pval <- mean(tbd.wolf.sppID$sims.list$chi2.sim > tbd.wolf.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.wolf.preyabund.pval <- mean(tbd.wolf.preyabund$sims.list$chi2.sim > tbd.wolf.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.wolf.sppID.preyabund.pval <- mean(tbd.wolf.sppID.preyabund$sims.list$chi2.sim > tbd.wolf.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.wolf.sppIDxpreyabund.pval <- mean(tbd.wolf.sppIDxpreyabund$sims.list$chi2.sim > tbd.wolf.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppIDxpreyabund$samples)
  save(tbd.wolf.sppIDxpreyabund, file = "./Outputs/tbd_wolf_sppID_X_preyRAI.RData")      
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
  (tbd.nt.bear.null.pval <- mean(tbd.nt.bear.null$sims.list$chi2.sim > tbd.nt.bear.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.bear.sppID.pval <- mean(tbd.nt.bear.sppID$sims.list$chi2.sim > tbd.nt.bear.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.bear.preyabund.pval <- mean(tbd.nt.bear.preyabund$sims.list$chi2.sim > tbd.nt.bear.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.preyabund$samples)
  save(tbd.nt.bear.preyabund, file = "./Outputs/tbd_nontarget_bear_preyRAI.RData")
  
  
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
  (tbd.nt.bob.null.pval <- mean(tbd.nt.bob.null$sims.list$chi2.sim > tbd.nt.bob.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.bob.sppID.pval <- mean(tbd.nt.bob.sppID$sims.list$chi2.sim > tbd.nt.bob.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.bob.preyabund.pval <- mean(tbd.nt.bob.preyabund$sims.list$chi2.sim > tbd.nt.bob.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.preyabund$samples)
  save(tbd.nt.bob.preyabund, file = "./Outputs/tbd_nontarget_bob_preyRAI.RData")
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
  (tbd.nt.coy.null.pval <- mean(tbd.nt.coy.null$sims.list$chi2.sim > tbd.nt.coy.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.coy.sppID.pval <- mean(tbd.nt.coy.sppID$sims.list$chi2.sim > tbd.nt.coy.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.coy.preyabund.pval <- mean(tbd.nt.coy.preyabund$sims.list$chi2.sim > tbd.nt.coy.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.preyabund$samples)
  save(tbd.nt.coy.preyabund, file = "./Outputs/tbd_nontarget_coy_preyRAI.RData") 
  
   
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
  (tbd.nt.lion.null.pval <- mean(tbd.nt.lion.null$sims.list$chi2.sim > tbd.nt.lion.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.lion.sppID.pval <- mean(tbd.nt.lion.sppID$sims.list$chi2.sim > tbd.nt.lion.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.lion.preyabund.pval <- mean(tbd.nt.lion.preyabund$sims.list$chi2.sim > tbd.nt.lion.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.preyabund$samples)
  save(tbd.nt.lion.preyabund, file = "./Outputs/tbd_nontarget_lion_preyRAI.RData")
  
  
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
  (tbd.nt.wolf.null.pval <- mean(tbd.nt.wolf.null$sims.list$chi2.sim > tbd.nt.wolf.null$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.wolf.sppID.pval <- mean(tbd.nt.wolf.sppID$sims.list$chi2.sim > tbd.nt.wolf.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
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
  (tbd.nt.wolf.preyabund.pval <- mean(tbd.nt.wolf.preyabund$sims.list$chi2.sim > tbd.nt.wolf.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.preyabund$samples)
  save(tbd.nt.wolf.preyabund, file = "./Outputs/tbd_nontarget_wolf_preyRAI.RData") 
  
  
  #'  Fin
  #'  Next stop, 08_DIC_exponential_model_selection.R for model selection
  
  
