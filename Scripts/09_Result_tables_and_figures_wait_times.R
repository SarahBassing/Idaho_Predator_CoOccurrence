  #'  ---------------------------------
  #'  Figures and result tables for wait time analyses
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ---------------------------------
  #'  Script to create figures and result tables based on results from wait time
  #'  analyses that estimated effect of recent competitor/prey presence and prey
  #'  availability on the elapsed time before a predator is detected. 
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(ggplot2)
  library(stringr)
  library(tidyverse)
  library(khroma)
  library(patchwork)
  library(grid)
  library(png)
  library(RCurl)
  library(rphylopic)
  
  #'  Load top predator - predator models
  load("./Outputs/tbd_bear_preyRAI.RData")
  load("./Outputs/tbd_bob_sppID_preyRAI.RData")
  load("./Outputs/tbd_coy_sppID_X_preyRAI.RData")
  load("./Outputs/tbd_lion_preyRAI.RData")
  load("./Outputs/tbd_wolf_preyRAI.RData")
  
  #'  Load competitor Species ID models
  load("./Outputs/tbd_bear_sppID.RData")
  load("./Outputs/tbd_bob_sppID.RData")
  load("./Outputs/tbd_coy_sppID.RData")
  load("./Outputs/tbd_lion_sppID.RData")
  load("./Outputs/tbd_wolf_sppID.RData")
  
  #'  Load top prey - predator models
  load("./Outputs/tbd_nontarget_bear_preyRAI.RData")
  load("./Outputs/tbd_nontarget_bob_preyRAI.RData")
  load("./Outputs/tbd_nontarget_coy_preyRAI.RData")
  load("./Outputs/tbd_nontarget_lion_preyRAI.RData")
  load("./Outputs/tbd_nontarget_wolf_preyRAI.RData")
  
  #'  Load prey - predator Species ID models
  load("./Outputs/tbd_nontarget_bear_sppID.RData")
  load("./Outputs/tbd_nontarget_bob_sppID.RData")
  load("./Outputs/tbd_nontarget_coy_sppID.RData")
  load("./Outputs/tbd_nontarget_lion_sppID.RData")
  load("./Outputs/tbd_nontarget_wolf_sppID.RData")
  
  #'  Load bundled data, including new covariate data
  #'  Generated in 07_Run_exponential_models.R script
  load("./Data/bear_bundled.RData")
  load("./Data/bob_bundled.RData")
  load("./Data/coy_bundled.RData")
  load("./Data/lion_bundled.RData")
  load("./Data/wolf_bundled.RData")
  
  load("./Data/bear_bundled_nontarget_RData")
  load("./Data/bob_bundled_nontarget_RData")
  load("./Data/coy_bundled_nontarget_RData")
  load("./Data/lion_bundled_nontarget_RData")
  load("./Data/wolf_bundled_nontarget_RData")
  
  #'  ------------------------
  ####  Species silhouettes  ####
  #'  ------------------------
  #'  Silhouettes for each species from PhyloPic in two different formats (PNG & rastergrob)
  #'  White-tailed deer silhouettes created by the talented Gabriela Palomo-Munoz and uploaded to http://phylopic.org/
  wtdurlGB1 <- "https://images.phylopic.org/images/8569838c-c725-4772-b0a3-b5eb04baaada/raster/1024x850.png?v=17cfdbaf920.png"
  wtdimgGB1 <- readPNG(getURLContent(wtdurlGB1), native = T)
  wtdgrid <- rasterGrob(wtdimgGB1, interpolate = TRUE)
  wtdurlGB2 <- "https://images.phylopic.org/images/6038e80c-398d-47b2-9a69-2b9edf436f64/raster/1023x1024.png?v=17cfdb9f8b6.png"
  wtdimgGB2 <- readPNG(getURLContent(wtdurlGB2), native = T)
  wtdgridGB2 <- rasterGrob(wtdimgGB2, interpolate = TRUE)
  bunnyurl <- "https://images.phylopic.org/images/f69eb95b-3d0d-491d-9a7f-acddd419afed/raster/925x1024.png?v=177f427b3d8.png"
  bunnyimg <- readPNG(getURLContent(bunnyurl), native = T)
  bunnygrid <- rasterGrob(bunnyimg, interpolate = TRUE, scales::alpha("black", 0.5))
 
  
  #'  ------------------
  ####  Format results  ####
  #'  ------------------
  #'  Snag and reformat coefficents and predictions from each top model
  coefs <- function(mod_out, spp, prey1, prey2, prey3, comp1, comp2, comp3, comp4) {
    Species <- spp
    Estimate <- round(unlist(mod_out$mean), 2)
    lci <- round(unlist(mod_out$q2.5), 2)
    uci <- round(unlist(mod_out$q97.5), 2)
    CI <- paste(" ", lci, "-", uci) # need that extra space in front b/c excel thinks this is an equation otherwise
    overlap0 <- unlist(mod_out$overlap0)
    out <- as.data.frame(cbind(Species, Estimate, CI, lci, uci, overlap0))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "95% CI", "lci", "uci", "overlap0")
    renamed_out <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta.prey1", paste("Prey RAI:", prey1), Parameter),
             Parameter = ifelse(Parameter == "beta.prey2", paste("Prey RAI:", prey2), Parameter),
             Parameter = ifelse(Parameter == "beta.prey3", paste("Prey RAI:", prey3), Parameter),
             Parameter = ifelse(Parameter == "beta.sppID1", paste("Competitor:", comp1), Parameter), 
             Parameter = ifelse(Parameter == "beta.sppID2", paste("Competitor:", comp2), Parameter), 
             Parameter = ifelse(Parameter == "beta.sppID3", paste("Competitor:", comp3), Parameter), 
             Parameter = ifelse(Parameter == "beta.sppID4", paste("Competitor:", comp4), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd1", paste("Mean TBD:", comp1), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd2", paste("Mean TBD:", comp2), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd3", paste("Mean TBD:", comp3), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd4", paste("Mean TBD:", comp4), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd1", paste0("Competitor:Prey (", comp1, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd2", paste0("Competitor:Prey (", comp2, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd3", paste0("Competitor:Prey (", comp3, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd4", paste0("Competitor:Prey (", comp4, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago1", paste0("Competitor:Prey (", comp1, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago2", paste0("Competitor:Prey (", comp2, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago3", paste0("Competitor:Prey (", comp3, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago4", paste0("Competitor:Prey (", comp4, " x ", prey2, ")"), Parameter),
             Parameter = gsub("spp.tbd.", "tbd ", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter)) %>%
      filter(lci != 0) %>%
      filter(Parameter != "deviance") %>%
      filter(Parameter != "chi2.obs") %>%
      filter(Parameter != "chi2.sim") %>%
      mutate(Estimate = as.numeric(Estimate),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed_out)
  }
  tbd.bear.out <- coefs(tbd.bear.preyabund, spp = "Black bear", prey1 = "Elk", prey2 = "White-tailed deer")
  tbd.bob.out <- coefs(tbd.bob.sppID.preyabund, spp = "Bobcat", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.coy.out <- coefs(tbd.coy.sppIDxpreyabund, spp = "Coyote", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Black bear", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.lion.out <- coefs(tbd.lion.preyabund, spp = "Mountain lion", prey1 = "Elk", prey2 = "White-tailed deer", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Bobcat", comp4 = "Wolf")
  tbd.wolf.out <- coefs(tbd.wolf.preyabund, spp = "Wolf", prey1 = "Elk", prey2 = "Moose", prey3 = "White-tailed deer")
  
  tbd.bear.comp <- coefs(tbd.bear.sppID, spp = "Black bear", comp1 = "Coyote", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.bob.comp <- coefs(tbd.bob.sppID, spp = "Bobcat", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.coy.comp <- coefs(tbd.coy.sppID, spp = "Coyote", comp1 = "Black bear", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.lion.comp <- coefs(tbd.lion.sppID, spp = "Mountain lion", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Bobcat", comp4 = "Wolf")
  tbd.wolf.comp <- coefs(tbd.wolf.sppID, spp = "Wolf", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Bobcat", comp4 = "Mountain lion")
  
  tbd.bear.prey <- coefs(tbd.nt.bear.preyabund, spp = "Black bear", prey1 = "Elk", prey2 = "White-tailed deer") 
  tbd.bob.prey <- coefs(tbd.nt.bob.preyabund, spp = "Bobcat", prey1 = "White-tailed deer", prey2 = "Lagomorph") 
  tbd.coy.prey <- coefs(tbd.nt.coy.preyabund, spp = "Coyote", prey1 = "White-tailed deer", prey2 = "Lagomorph") 
  tbd.lion.prey <- coefs(tbd.nt.lion.preyabund, spp = "Mountain lion", prey1 = "Elk", prey2 = "White-tailed deer") 
  tbd.wolf.prey <- coefs(tbd.nt.wolf.preyabund, spp = "Wolf", prey1 = "Elk", prey2 = "Moose", prey3 = "White-tailed deer") 
  
  tbd.bear.nt <- coefs(tbd.nt.bear.sppID, spp = "Black bear", comp1 = "Elk", comp2 = "White-tailed deer") %>%
    mutate(Parameter = gsub("Competitor", "Non-target", Parameter))
  tbd.bob.nt <- coefs(tbd.nt.bob.sppID, spp = "Bobcat", comp1 = "Lagomorph", comp2 = "White-tailed deer") %>%
    mutate(Parameter = gsub("Competitor", "Non-target", Parameter))
  tbd.coy.nt <- coefs(tbd.nt.coy.sppID, spp = "Coyote", comp1 = "Lagomorph", comp2 = "White-tailed deer") %>%
    mutate(Parameter = gsub("Competitor", "Non-target", Parameter))
  tbd.lion.nt <- coefs(tbd.nt.lion.sppID, spp = "Mountain lion", comp1 = "Elk", comp2 = "White-tailed deer") %>%
    mutate(Parameter = gsub("Competitor", "Non-target", Parameter))
  tbd.wolf.nt <- coefs(tbd.nt.wolf.sppID, spp = "Wolf", comp1 = "Elk", comp2 = "Moose", comp3 = "White-tailed deer") %>%
    mutate(Parameter = gsub("Competitor", "Non-target", Parameter))
  
  #'  Pull out just coefficient estimates
  bear.coefs <- tbd.bear.out[1:3,]
  bob.coefs <- tbd.bob.out[1:6,]
  coy.coefs <- tbd.coy.out[1:12,]
  lion.coefs <- tbd.lion.out[1:3,]
  wolf.coefs <- tbd.wolf.out[1:4,]
  
  #'  Pull out mean TBD estimates
  bear.mean.tbd <- tbd.bear.out[4,1:6]
  bob.mean.tbd <- tbd.bob.out[7:11,1:6]
  coy.mean.tbd <- tbd.coy.out[13:17,1:6]
  lion.mean.tbd <- tbd.lion.out[4,1:6]
  wolf.mean.tbd <- tbd.wolf.out[5,1:6]
  
  #'  Pull out mean TBD per competitor species (note: using output from different
  #'  models depending on whether sppID was in the top model or not)
  bear.comp.tbd <- tbd.bear.comp[5:9,1:6]
  bob.comp.tbd <- tbd.bob.out[7:11,1:6]
  coy.comp.tbd <- tbd.coy.out[13:17,1:6] 
  lion.comp.tbd <- tbd.lion.comp[5:9,1:6] 
  wolf.comp.tbd <- tbd.wolf.comp[5:9,1:6]
  
  #'  Pull out predicted TBD values
  bear.tbd.predictions <- tbd.bear.out[5:204,1:6]
  bob.tbd.predictions <- tbd.bob.out[12:811,1:6]
  coy.tbd.predictions <- tbd.coy.out[18:817,1:6]
  lion.tbd.predictions <- tbd.lion.out[5:204,1:6]
  wolf.tbd.predictions <- tbd.wolf.out[6:305,1:6]
  
  #'  Pull out mean TBD estimates for prey - predator from top model
  bear.nt.mean.tbd <- tbd.bear.prey[4,1:6]
  bob.nt.mean.tbd <- tbd.bob.prey[6,1:6]
  coy.nt.mean.tbd <- tbd.coy.prey[6,1:6]
  lion.nt.mean.tbd <- tbd.lion.prey[4,1:6]
  wolf.nt.mean.tbd <- tbd.wolf.prey[5,1:6]
  
  #'  Pull out mean TBD estimates for each non-target species from prey - predator sppID model 
  bear.nt.prey.tbd <- tbd.bear.nt[4:5,1:6]
  bob.nt.prey.tbd <- tbd.bob.nt[4:5,1:6]
  coy.nt.prey.tbd <- tbd.coy.nt[4:5,1:6]
  lion.nt.prey.tbd <- tbd.lion.nt[4:5,1:6]
  wolf.nt.prey.tbd <- tbd.wolf.nt[5:7,1:6]
  
  
  #'  ------------------
  ####  Result tables  ####
  #'  -----------------
  tbd.coefs <- rbind(bear.coefs, bob.coefs, coy.coefs, lion.coefs, wolf.coefs)
  mean.tbd <- rbind(bear.mean.tbd, bob.mean.tbd, coy.mean.tbd, lion.mean.tbd, wolf.mean.tbd) %>%
    mutate(Estimate = round(Estimate/60, 2),
           lci = round(lci/60, 2),
           uci = round(uci/60, 2),
           `95% CI` = paste(" ", lci, "-", uci)) %>%
    relocate(`95% CI`, .after = "Estimate")
  competitor.tbd <- rbind(bear.comp.tbd, bob.comp.tbd, coy.comp.tbd, lion.comp.tbd, wolf.comp.tbd) %>%
    mutate(Estimate = round(Estimate/60, 2),
           lci = round(lci/60, 2),
           uci = round(uci/60, 2),
           `95% CI` = paste(" ", lci, "-", uci)) %>%
    relocate(`95% CI`, .after = "Estimate")
  nt.tbd <- rbind(bear.nt.mean.tbd, bob.nt.mean.tbd, coy.nt.mean.tbd, lion.nt.mean.tbd, wolf.nt.mean.tbd) %>%
    mutate(Parameter = "Mean TBD: Primary prey") %>%
    mutate(Estimate = round(Estimate/60, 2),
           lci = round(lci/60, 2),
           uci = round(uci/60, 2),
           `95% CI` = paste(" ", lci, "-", uci)) %>%
    relocate(`95% CI`, .after = "Estimate")
  prey.tbd <- rbind(bear.nt.prey.tbd, bob.nt.prey.tbd, coy.nt.prey.tbd, lion.nt.prey.tbd, wolf.nt.prey.tbd) %>%
    mutate(Estimate = round(Estimate/60, 2),
           lci = round(lci/60, 2),
           uci = round(uci/60, 2),
           `95% CI` = paste(" ", lci, "-", uci)) %>%
    relocate(`95% CI`, .after = "Estimate")
  predicted.tbd <- rbind(bear.tbd.predictions, bob.tbd.predictions, coy.tbd.predictions, lion.tbd.predictions, wolf.tbd.predictions) %>% 
    #'  Split out predictions based on categorical variable (mostly important for bobcat & coyote results)
    mutate(Estimate = round(Estimate/60, 2),
           lci = round(lci/60, 2),
           uci = round(uci/60, 2),
           `95% CI` = paste(" ", lci, "-", uci), Prey_species = str_replace(Parameter, "tbd ", ""),
           Prey_species = str_extract(Prey_species, "[aA-zZ]+"),
           Obs_nmbr = as.numeric(str_extract(Parameter, "[0-9]+")),
           #'  If competitor ID had no effect, use reference category (coyote)
           Competitor_ID = ifelse(Species == "Black bear" | Species == "Mountain lion" | Species == "Wolf", "Coyote", NA), 
           #'  If competitor ID had an effect, assign correct species to each data chunk
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr <101, "Coyote", Competitor_ID),
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >100, "Black bear", Competitor_ID),
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >200, "Mountain lion", Competitor_ID), 
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >300, "Wolf", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr <101, "Black bear", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >100, "Bobcat", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >200, "Mountain lion", Competitor_ID), 
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >300, "Wolf", Competitor_ID)) %>%
    relocate(`95% CI`, .after = "Estimate")
  
  #' #'  Save
  #' write.csv(tbd.coefs, "./Outputs/TBD_coefficient_estimates_allSpp.csv")
  #' write.csv(mean.tbd, "./Outputs/TBD_estimated_means_allSpp.csv")
  #' write.csv(competitor.tbd, "./Outputs/TBD_competitor_means_allSpp.csv")
  #' write.csv(nt.tbd, "./Outputs/TBD_nontarget_means_allSpp.csv")
  #' write.csv(prey.tbd, "./Outputs/TBD_nontarget_prey_means_allSpp.csv")
  #' write.csv(predicted.tbd, "./Outputs/TBD_predicted_TBD_allSpp.csv")
  
  
  #'  ---------------------------------
  ####  Plot TBD by covariate effects  ####
  #'  ---------------------------------
  #'  Format new covariate data based on length of predictions for each predator species
  bear.covs <- c(bear_bundled$newcovs[,1], bear_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  bob.covs <- c(bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("wtd", "lago"), each = 400)) %>%
    rename("cov" = ".")
  coy.covs <- c(coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("wtd", "lago"), each = 400)) %>%
    rename("cov" = ".")
  lion.covs <- c(lion_bundled$newcovs[,1], lion_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  wolf.covs <- c(wolf_bundled$newcovs[,1], wolf_bundled$newcovs[,2], wolf_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "moose", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  
  #'  Review which coefficients did NOT overlap 0 - only retain predictions for
  #'  responses to prey RAI based on these relationships
  print(bear.coefs) # Prey RAI: Elk
  print(bob.coefs) # Prey RAI: White-tailed deer
  print(coy.coefs) # Competitor:Prey (Wolf x White-tailed deer), (Wolf x Lagomorph)
  print(lion.coefs) # Prey RAI: Elk, White-tailed deer
  print(wolf.coefs) # Prey RAI: White-tailed deer
  
  #'  Combine predictions and covariate data for plotting
  #'  Filter to only predictions derived from statistically meaningful relationships (i.e., 95% CRI did NOT overlap 0)
  #'  Change Competitor_ID to "any" for species where previous competitor detection did not matter
  bear.predicted <- filter(predicted.tbd, Species == "Black bear") %>%
    bind_cols(bear.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species == "elk") %>%
    mutate(Competitor_ID = "Any predator")
  bob.predicted <- filter(predicted.tbd, Species == "Bobcat") %>%    
    bind_cols(bob.covs) %>%
    dplyr::select(-prey_spp) 
  coy.predicted <- filter(predicted.tbd, Species == "Coyote") %>%    
    bind_cols(coy.covs) %>%
    dplyr::select(-prey_spp) %>%
    #'  Excluding non-significant lion effect because predictions are on entirely 
    #'  different scale from others which makes plotting very misleading
    filter(Prey_species != "lago" | Competitor_ID != "Mountain lion")
  lion.predicted <- filter(predicted.tbd, Species == "Mountain lion") %>%   
    bind_cols(lion.covs) %>%
    dplyr::select(-prey_spp) %>%
    mutate(Competitor_ID = "Any predator")
  wolf.predicted <- filter(predicted.tbd, Species == "Wolf") %>%
    bind_cols(wolf.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species == "wtd") %>%
    mutate(Competitor_ID = "Any predator")
  
  predicted.tbd.covs <- rbind(bear.predicted, bob.predicted, coy.predicted, lion.predicted, wolf.predicted)

  #'  Format results tables for plotting
  #'  Mean TBD by competitor species
  tbd_by_competitorID <- competitor.tbd %>%
    filter(Parameter != "Mean TBD") %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Previous_Species = str_replace(Parameter, "Mean TBD: ", ""),
           Previous_Species = factor(Previous_Species, levels = c("Coyote", "Bobcat", "Black bear", "Mountain lion", "Wolf")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) %>%
    relocate(Previous_Species, .after = Species)
  
  #'  Mean TBD by prey species
  tbd_by_preyID <- prey.tbd %>%
    filter(Parameter != "Mean TBD") %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Previous_Species = str_replace(Parameter, "Mean TBD: ", ""),
           Previous_Species = factor(Previous_Species, levels = c("Lagomorph", "Elk", "Moose", "White-tailed deer")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) %>%
    relocate(Previous_Species, .after = Species)
  
  #'  Mean TBD by competitor species and mean TBD following any prey species
  tbd_meanPrey_and_competitorID <- competitor.tbd %>%
    filter(Parameter != "Mean TBD") %>%
    rbind(nt.tbd) %>% 
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Previous_Species = str_replace(Parameter, "Mean TBD: ", ""),
           Previous_Species = factor(Previous_Species, levels = c("Primary prey", "Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")), 
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) %>%
    relocate(Previous_Species, .after = Species)
  
  #'  Mean TBD by competitor species and by prey species
  tbd_by_PreyID_and_competitorID <- competitor.tbd %>%
    filter(Parameter != "Mean TBD") %>%
    rbind(prey.tbd) %>% 
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Previous_Species = str_replace(Parameter, "Mean TBD: ", ""),
           Previous_Species = factor(Previous_Species, levels = c("Lagomorph", "Elk", "Moose", "White-tailed deer", "Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")), 
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) %>%
    relocate(Previous_Species, .after = Species)

  #'  Predicted TBD over range of prey relative abundance values
  tbd_prey_prediction <- predicted.tbd.covs %>%
    # arrange(Species, Prey_species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Competitor_ID = factor(Competitor_ID, levels = c("Any predator", "Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Prey_RAI = Prey_species,
           Prey_RAI = factor(Prey_RAI, levels = c("Elk" = "elk", "Moose" = "moose", "White-tailed deer" = "wtd", "Lagomorphs" = "lago")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) 
 
  #' #'  Choose colorblind-friendly scheme
  #' plot_scheme(colour("sunset")(11))
  #' colour("sunset")(11)
  #' any.bear.bob.coy.lion.wolf_colors <- c("black", "#98CAE1", "#A50026", "#DD3D2D", "#FDB366", "#364B9A")
  #' bear.coy.lion.wolf_colors <- c("#98CAE1", "#DD3D2D", "#FDB366", "#364B9A")
  #' any.bear.coy_colors <- c("black", "#98CAE1", "#DD3D2D")
  #' bear.wolf <- c("#98CAE1", "#364B9A")
  
  #####  Competitor ID effect  #####
  #'  Effect of previous competitor detection on latency of site use
  competitorID_plot <- ggplot(tbd_by_competitorID, aes(x = Previous_Species, y = (Estimate), group = Species)) +
    geom_errorbar(aes(ymin = (lci), ymax = (uci), color = Previous_Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Previous_Species), size = 2.5, position = position_dodge(width = 0.4)) +
    theme_bw() +
    guides(color = guide_legend(title = "Previously detected species")) +
    facet_wrap(~Species, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(legend.position.inside = c(1, 0), legend.justification = c(1, 0)) +
    xlab("Previously detected species") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent competitor detection on wait time until site use")
  competitorID_plot
  
  competitorID_and_meanprey_plot <- ggplot(tbd_meanPrey_and_competitorID, aes(x = Previous_Species, y = (Estimate), group = Species)) +
    geom_errorbar(aes(ymin = (lci), ymax = (uci), color = Previous_Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Previous_Species), size = 2.5, position = position_dodge(width = 0.4)) +
    theme_bw() +
    scale_color_manual(values = any.bear.bob.coy.lion.wolf_colors) +
    guides(color = guide_legend(title = "Previously detected \nspecies")) +
    facet_wrap(~Species, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 14)) +
    theme(legend.position="none") +
    xlab("Previously detected species") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent species detection on wait time")
  competitorID_and_meanprey_plot
  
  competitorID_and_preyID_plot <- ggplot(tbd_by_PreyID_and_competitorID, aes(x = Previous_Species, y = (Estimate), group = Species)) +
    geom_errorbar(aes(ymin = (lci), ymax = (uci), color = Previous_Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Previous_Species), size = 2.5, position = position_dodge(width = 0.4)) +
    theme_bw() +
    guides(color = guide_legend(title = "Previously detected species")) +
    facet_wrap(~Species, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
    xlab("Previously detected species") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent species detection on wait time")
  competitorID_and_preyID_plot
  
  preyID_plot <- ggplot(tbd_by_preyID, aes(x = Previous_Species, y = (Estimate), group = Species)) +
    geom_errorbar(aes(ymin = (lci), ymax = (uci), color = Previous_Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Previous_Species), size = 2.5, position = position_dodge(width = 0.4)) +
    theme_bw() +
    guides(color = guide_legend(title = "Previously detected species")) +
    facet_wrap(~Species, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
    xlab("Previously detected species") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent prey detection on wait time until site use")
  preyID_plot
  
  
  #####  Prey relative abundance effects  ####
  tbd_wtdRAI_plot <- filter(tbd_prey_prediction, Prey_RAI == "wtd") %>%
    ggplot(aes(x = cov, y = (Estimate), group = Competitor_ID)) +
    geom_line(aes(color = Competitor_ID), lwd = 1.25) + 
    scale_color_manual(values = any.bear.bob.coy.lion.wolf_colors) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = (lci), ymax = (uci), fill = Competitor_ID), alpha = 0.3) +
    scale_fill_manual(values = any.bear.bob.coy.lion.wolf_colors) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(text = element_text(size = 12)) +
    xlab("White-tailed deer relative abundance (standardized)") +
    ylab("Mean number of hours between detections") +
    guides(color = guide_legend(title = "Previously detected\ncompetitor"),
           fill = guide_legend(title = "Previously detected\ncompetitor")) +
    theme(legend.position = "left",
          text = element_text(size = 12)) +
    facet_wrap(~Species, scale = "free_y") + 
    inset_element(p = bobimg, left = 0.30, bottom = 0.88, right = 0.45, top = 0.98) +
    inset_element(p = cougimgGB, left = 0.28, bottom = 0.30, right = 0.48, top = 0.45) +
    inset_element(p = coyimgGB, left = 0.82, bottom = 0.88, right = 0.97, top = 0.98) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank")) +
    inset_element(p = wolfimg, left = 0.80, bottom = 0.32, right = 0.97, top = 0.45) 
  tbd_wtdRAI_plot
  
  tbd_elkRAI_plot <- filter(tbd_prey_prediction, Prey_RAI == "elk") %>%
    filter(Species == "Mountain lion") %>%
    ggplot(aes(x = cov, y = (Estimate), group = Competitor_ID)) +
    geom_line(aes(color = Competitor_ID), lwd = 1.25) + 
    scale_color_manual(values = any.bear.coy.lion_colors) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = (lci), ymax = (uci), fill = Competitor_ID), alpha = 0.3) +
    scale_fill_manual(values = any.bear.coy.lion_colors) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "right") +
    theme(text = element_text(size = 12)) +
    theme(axis.title.y = element_blank()) +
    xlab("Relative abundance of elk \n(standardized)") +
    theme(legend.title=element_blank()) +
    facet_wrap(~Species, ncol = 1, scale = "free_x") +
    inset_element(p = cougimgGB, left = 0.65, bottom = 0.78, right = 0.98, top = 1) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank"))
  tbd_elkRAI_plot

  #'  Combine species-specific effects into a single plot for publication
  (tbd_lagoRAI_elkRAI_patchwork <- tbd_lagoRAI_plot + tbd_elkRAI_plot + plot_layout(guides = 'collect'))
  
 
  #'  Save
  ggsave("./Outputs/Figures/Mean_TBD_by_competitorID.tiff", competitorID_plot, 
         units = "in", width = 8, height = 6, dpi = 400, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/Mean_TBD_by_preyID.tiff", preyID_plot, 
         units = "in", width = 8, height = 6, dpi = 400, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/Mean_TBD_by_meanprey_competitorID.tiff", competitorID_and_meanprey_plot, 
         units = "in", width = 7, height = 5, dpi = 400, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/Mean_TBD_by_preyID_competitorID.tiff", competitorID_and_preyID_plot, 
         units = "in", width = 8, height = 8, dpi = 400, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/Predicted_TBD_by_wtdRAI.tiff", tbd_wtdRAI_plot,
         units = "in", width = 7, height = 5, dpi = 400, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/Predicted_TBD_by_lagoRAI_&_elkRAI.tiff", tbd_lagoRAI_elkRAI_patchwork,
         units = "in", width = 7, height = 4, dpi = 400, device = 'tiff', compression = 'lzw')
  

  