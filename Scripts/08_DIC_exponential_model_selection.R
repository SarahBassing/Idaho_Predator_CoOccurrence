  #'  ---------------------------------
  #'  Wait time model selection & result tables
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ---------------------------------
  #'  Script to identify most supported exponential model per predator species 
  #'  using DIC. Requires all model outputs from 07_Run_exponential_models.R are 
  #'  saved as .RData files. 
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(AICcmodavg)
  library(tidyverse)
  
  #'  -------------------------------
  ####  COMPETITOR - FOCAL PREDATOR  ####
  #'  -------------------------------
  #'  Model selection for models estimating the time that elapses before detecting
  #'  a different predator species.
  
  #####  Black bear models  ####
  load("./Outputs/tbd_bear_null.RData")
  load("./Outputs/tbd_bear_sppID.RData")
  load("./Outputs/tbd_bear_preyRAI.RData")
  load("./Outputs/tbd_bear_sppID_preyRAI.RData")
  load("./Outputs/tbd_bear_sppID_X_preyRAI.RData") 
  #'  List for model selection (some models excluded owing to poor convergence)
  bear_tbd_list <- list(tbd_bear_null, tbd_bear_sppID, tbd_bear_preyabund, tbd_bear_sppID.preyabund) #tbd_bear_sppIDxpreyabund
  bear_tbd_name <- c("tbd_bear_null", "tbd_bear_sppID", "tbd_bear_preyabund", "tbd_bear_sppID.preyabund") #"tbd_bear_sppIDxpreyabund"
  
  #####  Bobcat models  ####
  load("./Outputs/tbd_bob_null.RData")
  load("./Outputs/tbd_bob_sppID.RData")
  load("./Outputs/tbd_bob_preyRAI.RData")
  load("./Outputs/tbd_bob_sppID_preyRAI.RData")
  load("./Outputs/tbd_bob_sppID_X_preyRAI.RData")
  #'  List for model selection (some models excluded owing to poor convergence)
  bob_tbd_list <- list(tbd_bob_null, tbd_bob_sppID, tbd_bob_preyabund, tbd_bob_sppID.preyabund, tbd_bob_sppIDxpreyabund) 
  bob_tbd_name <- c("tbd_bob_null", "tbd_bob_sppID", "tbd_bob_preyabund", "tbd_bob_sppID.preyabund", "tbd_bob_sppIDxpreyabund") 
    
  #####  Coyote models  ####
  load("./Outputs/tbd_coy_null.RData")
  load("./Outputs/tbd_coy_sppID.RData")
  load("./Outputs/tbd_coy_preyRAI.RData")
  load("./Outputs/tbd_coy_sppID_preyRAI.RData")
  load("./Outputs/tbd_coy_sppID_X_preyRAI.RData")
  #'  List for model selection (some models excluded owing to poor convergence)
  coy_tbd_list <- list(tbd_coy_null, tbd_coy_sppID, tbd_coy_preyabund, tbd_coy_sppID.preyabund, tbd_coy_sppIDxpreyabund) 
  coy_tbd_name <- c("tbd_coy_null", "tbd_coy_sppID", "tbd_coy_preyabund", "tbd_coy_sppID.preyabund", "tbd_coy_sppIDxpreyabund") 
  
  #####  Mountain lion models  ####
  load("./Outputs/tbd_lion_null.RData")
  load("./Outputs/tbd_lion_sppID.RData")
  load("./Outputs/tbd_lion_preyRAI.RData")
  load("./Outputs/tbd_lion_sppID_preyRAI.RData")   
  load("./Outputs/tbd_lion_sppID_X_preyRAI.RData") 
  #'  List for model selection (some models excluded owing to poor convergence)
  lion_tbd_list <- list(tbd_lion_null, tbd_lion_sppID, tbd_lion_preyabund, tbd_lion_sppID.preyabund) #tbd_lion_sppIDxpreyabund
  lion_tbd_name <- c("tbd_lion_null", "tbd_lion_sppID", "tbd_lion_preyabund", "tbd_lion_sppID.preyabund") #"tbd_lion_sppIDxpreyabund"
    
  #####  Wolf models  ####
  load("./Outputs/tbd_wolf_null.RData")
  load("./Outputs/tbd_wolf_sppID.RData")
  load("./Outputs/tbd_wolf_preyRAI.RData")
  load("./Outputs/tbd_wolf_sppID_preyRAI.RData")
  load("./Outputs/tbd_wolf_sppID_X_preyRAI.RData") 
  #'  List for model selection (some models excluded owing to poor convergence)
  wolf_tbd_list <- list(tbd_wolf_null, tbd_wolf_sppID, tbd_wolf_preyabund, tbd_wolf_sppID.preyabund) #, tbd_wolf_sppIDxpreyabund
  wolf_tbd_name <- c("tbd_wolf_null", "tbd_wolf_sppID", "tbd_wolf_preyabund", "tbd_wolf_sppID.preyabund") #, "tbd_wolf_sppIDxpreyabund"
    
  #####  Model selection using DIC  ####
  (topmod_beartbd <- dictab(cand.set = bear_tbd_list, modnames = bear_tbd_name, sort = TRUE)) 
  (topmod_bobtbd <- dictab(cand.set = bob_tbd_list, modnames = bob_tbd_name, sort = TRUE)) 
  (topmod_coytbd <- dictab(cand.set = coy_tbd_list, modnames = coy_tbd_name, sort = TRUE)) 
  (topmod_liontbd <- dictab(cand.set = lion_tbd_list, modnames = lion_tbd_name, sort = TRUE)) 
  (topmod_wolftbd <- dictab(cand.set = wolf_tbd_list, modnames = wolf_tbd_name, sort = TRUE)) 
  
  #'  Best supported model per species
  (topmodels <- rbind(topmod_beartbd[2,], topmod_bobtbd[2,], topmod_coytbd[1,], topmod_liontbd[1,], topmod_wolftbd[1,])) 
  #'  Note: currently using 2nd most supported model for black bear since 0.23 deltaDIC from top model (null) and 95% CRI for elk effect just slightly overlaps 0
  #'  currently using 2nd most supported model for bobcat since w/in 0.90 deltaDIC of top model (preyabund)
  
  
  #'  ------------------------
  #### PREY - FOCAL PREDATOR  ####
  #'  ------------------------
  #'  Model selection for models estimating time between detections of prey followed by predators
  
  #####  Black bear models  ####
  load("./Outputs/tbd_nontarget.bear_null.RData")
  load("./Outputs/tbd_nontarget.bear_sppID.RData")
  load("./Outputs/tbd_nontarget.bear_preyRAI.RData")
  #'  List for model selection
  bear_nt_tbd_list <- list(tbd_nt_bear_null, tbd_nt_bear_sppID, tbd_nt_bear_preyabund) 
  bear_nt_tbd_name <- c("tbd_nt_bear_null", "tbd_nt_bear_sppID", "tbd_nt_bear_preyabund") 
  
  #####  Bobcat models  ####
  load("./Outputs/tbd_nontarget.bob_null.RData")
  load("./Outputs/tbd_nontarget.bob_sppID.RData")
  load("./Outputs/tbd_nontarget.bob_preyRAI.RData")
  #'  List for model selection
  bob_nt_tbd_list <- list(tbd_nt_bob_null, tbd_nt_bob_sppID, tbd_nt_bob_preyabund)  
  bob_nt_tbd_name <- c("tbd_nt_bob_null", "tbd_nt_bob_sppID", "tbd_nt_bob_preyabund")  
  
  #####  Coyote models  ####
  load("./Outputs/tbd_nontarget.coy_null.RData")
  load("./Outputs/tbd_nontarget.coy_sppID.RData")
  load("./Outputs/tbd_nontarget.coy_preyRAI.RData")
  #'  List for model selection
  coy_nt_tbd_list <- list(tbd_nt_coy_null, tbd_nt_coy_sppID, tbd_nt_coy_preyabund) 
  coy_nt_tbd_name <- c("tbd_nt_coy_null", "tbd_nt_coy_sppID", "tbd_nt_coy_preyabund") 
  
  #####  Mountain lion models  ####
  load("./Outputs/tbd_nontarget.lion_null.RData")
  load("./Outputs/tbd_nontarget.lion_sppID.RData")
  load("./Outputs/tbd_nontarget.lion_preyRAI.RData")
  #'  List for model selection
  lion_nt_tbd_list <- list(tbd_nt_lion_null, tbd_nt_lion_sppID, tbd_nt_lion_preyabund) 
  lion_nt_tbd_name <- c("tbd_nt_lion_null", "tbd_nt_lion_sppID", "tbd_nt_lion_preyabund") 
  
  #####  Wolf models  ####
  load("./Outputs/tbd_nontarget.wolf_null.RData")
  load("./Outputs/tbd_nontarget.wolf_sppID.RData")
  load("./Outputs/tbd_nontarget.wolf_preyRAI.RData")
  #'  List for model selection
  wolf_nt_tbd_list <- list(tbd_nt_wolf_null, tbd_nt_wolf_sppID, tbd_nt_wolf_preyabund) 
  wolf_nt_tbd_name <- c("tbd_nt_wolf_null", "tbd_nt_wolf_sppID", "tbd_nt_wolf_preyabund") 
  
  #####  Model selection using DIC  ####
  (topmod_nt_beartbd <- dictab(cand.set = bear_nt_tbd_list, modnames = bear_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_bobtbd <- dictab(cand.set = bob_nt_tbd_list, modnames = bob_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_coytbd <- dictab(cand.set = coy_nt_tbd_list, modnames = coy_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_liontbd <- dictab(cand.set = lion_nt_tbd_list, modnames = lion_nt_tbd_name, sort = TRUE)) 
  (topmod_nt_wolftbd <- dictab(cand.set = wolf_nt_tbd_list, modnames = wolf_nt_tbd_name, sort = TRUE)) 
  
  #'  Best supported model per species-pair
  (topmodels_nt <- rbind(topmod_nt_beartbd[1,], topmod_nt_bobtbd[1,], topmod_nt_coytbd[1,], topmod_nt_liontbd[1,], topmod_nt_wolftbd[1,]))
  
  
  #'  -----------------
  ####  RESULT TABLES  ####
  #'  -----------------
  #'  Full table of models ranked by DIC for all predator - predator analyses
  model_list_DIC <- rbind(topmod_beartbd, topmod_bobtbd, topmod_coytbd, topmod_liontbd, topmod_wolftbd) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    #'  Split model name based on placement of multiple periods
    #'  https://stackoverflow.com/questions/26265400/use-regex-in-r-to-retrieve-string-before-second-occurence-of-a-period
    mutate(Species = sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species = str_replace(Species, "tbd.", ""),
           Species = str_replace(Species, "bear", "Black bear"),
           Species = str_replace(Species, "lion", "Mountain lion"),
           Species = str_replace(Species, "coy", "Coyote"),
           Species = str_replace(Species, "bob", "Bobcat"),
           Species = str_replace(Species, "wolf", "Wolf"),
           Model_name = sub("(^[^.]+[.][^.]+)(.+$)", "\\2", Modnames),
           Model_name = str_replace(Model_name, ".", ""),
           Model = gsub("null", "Model 1", Model_name), 
           Model = gsub("sppIDxpreyabund", "Model 5", Model),
           Model = gsub("sppID.preyabund", "Model 4", Model),
           Model = gsub("sppID", "Model 2", Model),
           Model = gsub("preyabund", "Model 3", Model),
           Model_name = gsub("null", "Null", Model_name),
           Model_name = gsub("sppIDxpreyabund", "Competitor ID * prey abundance", Model_name),
           Model_name = gsub("sppID.preyabund", "Competitor ID + prey abundance", Model_name),
           Model_name = gsub("sppID", "Competitor ID", Model_name),
           Model_name = gsub("preyabund", "Prey abundance", Model_name)) %>%
    relocate(Species, .before = DIC) %>%
    relocate(Model, .after = Species) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_list_DIC) <- c("Species", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Full table of models ranked by DIC for all non-target prey - predator analyses
  model_nt_list_DIC <- rbind(topmod_nt_beartbd, topmod_nt_bobtbd, topmod_nt_coytbd, topmod_nt_liontbd, topmod_nt_wolftbd) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    mutate(Species = sub( "(^[^.]+[.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species = str_replace(Species, "tbd.nt.", ""),
           Species = str_replace(Species, "bear", "Black bear"),
           Species = str_replace(Species, "lion", "Mountain lion"),
           Species = str_replace(Species, "coy", "Coyote"),
           Species = str_replace(Species, "bob", "Bobcat"),
           Species = str_replace(Species, "wolf", "Wolf"),
           Species = gsub("\\..*", "", Species),
           Model_name = gsub("^.*\\.", "", Modnames),
           Model = gsub("null", "Model 1", Model_name),           
           Model = gsub("sppID", "Model 2", Model),
           Model = gsub("preyabund", "Model 3", Model),
           Model_name = gsub("null", "Null", Model_name),
           Model_name = gsub("sppID.preyabund", "Species ID + prey abundance", Model_name),
           Model_name = gsub("sppID", "Species ID", Model_name),
           Model_name = gsub("preyabund", "Prey abundance", Model_name)) %>%
    relocate(Species, .before = DIC) %>%
    relocate(Model, .after = Species) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_nt_list_DIC) <- c("Species", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/DIC_TBD_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/DIC_TBD_model_selection_results.csv")
  save(topmodels, file = "./Outputs/DIC_TBD_top_models.RData")
  save(model_list_DIC, file = "./Outputs/DIC_TBD_model_selection_results.RData")
  
  write.csv(topmodels_nt, file = "./Outputs/DIC_TBD_top_nt_models.csv")
  write.csv(model_nt_list_DIC, file = "./Outputs/DIC_TBD_nt_model_selection_results.csv")
  save(topmodels_nt, file = "./Outputs/DIC_TBD_top_nt_models.RData")
  save(model_nt_list_DIC, file = "./Outputs/DIC_TBD_nt_model_selection_results.RData")
  
  
  #'  Next, create figures to visualize results in 09_Result_tables_and_figures_wait_times.R
  
  
  
