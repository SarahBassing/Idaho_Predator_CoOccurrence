  #'  ---------------------------------
  #'  Model selection with DIC
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ---------------------------------------------
  #'  Script to identify most supported multispecies model per predator dyad using 
  #'  DIC. Requires all model outputs from 01_Run_multispp_occupancy_models.R are 
  #'  saved as .RData files. 
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(AICcmodavg)
  library(tidyverse)
  
  #'  ----------------------
  ####  Load model outputs  ####
  #'  ----------------------
  #####  Wolf-Bear  #####                     
  load("./Outputs/wolfbear_null.RData") 
  load("./Outputs/wolfbear_hab.RData") 
  load("./Outputs/wolfbear_preyabund.RData")
  load("./Outputs/wolfbear_preydiv.RData")
  load("./Outputs/wolfbear_habx.RData") 
  load("./Outputs/wolfbear_preyabundx.RData")
  load("./Outputs/wolfbear_preydivx.RData") 
  load("./Outputs/wolfbear_global.RData")  
  wolfbear_list <- list(wolf.bear.null, wolf.bear.hab, wolf.bear.preyabund, wolf.bear.preydiv, wolf.bear.habx, wolf.bear.preyabundx, wolf.bear.preydivx, wolf.bear.global) 
  wolfbear_name <- c("wolf.bear.null", "wolf.bear.hab", "wolf.bear.preyabund", "wolf.bear.preydiv", "wolf.bear.habx", "wolf.bear.preyabundx", "wolf.bear.preydivx", "wolf.bear.global") 
  
  #####  Wolf-Coyote  #####                     
  load("./Outputs/wolfcoy_null.RData")  
  load("./Outputs/wolfcoy_hab.RData") 
  load("./Outputs/wolfcoy_preyabund.RData") 
  load("./Outputs/wolfcoy_preydiv.RData") 
  load("./Outputs/wolfcoy_habx.RData") 
  load("./Outputs/wolfcoy_preyabundx.RData") 
  load("./Outputs/wolfcoy_preydivx.RData") 
  load("./Outputs/wolfcoy_global.RData")  
  wolfcoy_list <- list(wolf.coy.null, wolf.coy.hab, wolf.coy.preyabund, wolf.coy.preydiv, wolf.coy.habx, wolf.coy.preyabundx, wolf.coy.preydivx, wolf.coy.global) 
  wolfcoy_name <- c("wolf.coy.null", "wolf.coy.hab", "wolf.coy.preyabund", "wolf.coy.preydiv", "wolf.coy.habx", "wolf.coy.preyabundx", "wolf.coy.preydivx", "wolf.coy.global")
  
  #####  Wolf-Lion  #####                           
  load("./Outputs/wolflion_null.RData")  
  load("./Outputs/wolflion_hab.RData")
  load("./Outputs/wolflion_preyabund.RData")
  load("./Outputs/wolflion_preydiv.RData")
  load("./Outputs/wolflion_habx.RData")
  load("./Outputs/wolflion_preyabundx.RData")
  load("./Outputs/wolflion_preydivx.RData") 
  load("./Outputs/wolflion_global.RData") 
  wolflion_list <- list(wolf.lion.null, wolf.lion.hab, wolf.lion.preyabund, wolf.lion.preydiv, wolf.lion.habx, wolf.lion.preyabundx, wolf.lion.preydivx, wolf.lion.global) 
  wolflion_name <- c("wolf.lion.null", "wolf.lion.hab", "wolf.lion.preyabund", "wolf.lion.preydiv", "wolf.lion.habx", "wolf.lion.preyabundx", "wolf.lion.preydivx", "wolf.lion.global") 
  
  #####  Lion-Bear  #####                           
  load("./Outputs/lionbear_null.RData") 
  load("./Outputs/lionbear_hab.RData")
  load("./Outputs/lionbear_preyabund.RData")
  load("./Outputs/lionbear_preydiv.RData") 
  load("./Outputs/lionbear_habx.RData")
  load("./Outputs/lionbear_preyabundx.RData")
  load("./Outputs/lionbear_preydivx.RData")
  load("./Outputs/lionbear_global.RData") 
  lionbear_list <- list(lion.bear.null, lion.bear.hab, lion.bear.preyabund, lion.bear.preydiv, lion.bear.habx, lion.bear.preyabundx, lion.bear.preydivx, lion.bear.global)
  lionbear_name <- c("lion.bear.null", "lion.bear.hab", "lion.bear.preyabund", "lion.bear.preydiv", "lion.bear.habx", "lion.bear.preyabundx", "lion.bear.preydivx", "lion.bear.global")
  
  #####  Lion-Bobcat  #####                         
  load("./Outputs/lionbob_null.RData") 
  load("./Outputs/lionbob_hab.RData")
  load("./Outputs/lionbob_preyabund.RData") 
  load("./Outputs/lionbob_preydiv.RData") 
  load("./Outputs/lionbob_habx.RData")
  load("./Outputs/lionbob_preyabundx.RData")
  load("./Outputs/lionbob_preydivx.RData")
  load("./Outputs/lionbob_global.RData") 
  lionbob_list <- list(lion.bob.null, lion.bob.hab, lion.bob.preyabund, lion.bob.preydiv, lion.bob.habx, lion.bob.preyabundx, lion.bob.preydivx, lion.bob.global)
  lionbob_name <- c("lion.bob.null", "lion.bob.hab", "lion.bob.preyabund", "lion.bob.preydiv", "lion.bob.habx", "lion.bob.preyabundx", "lion.bob.preydivx", "lion.bob.global")
  
  #####  Coyote-Bobcat  #####                 
  load("./Outputs/coybob_null.RData") 
  load("./Outputs/coybob_hab.RData")
  load("./Outputs/coybob_preyabund.RData")
  load("./Outputs/coybob_preydiv.RData") 
  load("./Outputs/coybob_habx.RData")
  load("./Outputs/coybob_preyabundx.RData") 
  load("./Outputs/coybob_preydivx.RData")
  load("./Outputs/coybob_global.RData") 
  coybob_list <- list(coy.bob.null, coy.bob.hab, coy.bob.preyabund, coy.bob.preydiv, coy.bob.habx, coy.bob.preyabundx, coy.bob.preydivx, coy.bob.global) 
  coybob_name <- c("coy.bob.null", "coy.bob.hab", "coy.bob.preyabund", "coy.bob.preydiv", "coy.bob.habx", "coy.bob.preyabundx", "coy.bob.preydivx", "coy.bob.global") 
  
  #'  -------------------
  ####  Model selection  ####
  #'  -------------------
  #'  Create model selection table using DIC, deltaDIC, and model weights
  (topmod_wolfbear <- dictab(cand.set = wolfbear_list, modnames = wolfbear_name, sort = TRUE)) 
  (topmod_wolfcoy <- dictab(cand.set = wolfcoy_list, modnames = wolfcoy_name, sort = TRUE)) 
  (topmod_wolflion <- dictab(cand.set = wolflion_list, modnames = wolflion_name, sort = TRUE)) 
  (topmod_lionbear <- dictab(cand.set = lionbear_list, modnames = lionbear_name, sort = TRUE)) 
  (topmod_lionbob <- dictab(cand.set = lionbob_list, modnames = lionbob_name, sort = TRUE)) 
  (topmod_coybob <- dictab(cand.set = coybob_list, modnames = coybob_name, sort = TRUE)) 
  
  #'  Table of best supported model per species-pair
  (topmodels <- rbind(topmod_wolfbear[1,], topmod_wolfcoy[1,], topmod_wolflion[1,], topmod_lionbear[1,], topmod_lionbob[1,], topmod_coybob[1,]))
  
  #'  Full table of models ranked by DIC for all species-pairs
  model_list_DIC <- rbind(topmod_wolfbear, topmod_wolfcoy, topmod_coybob, topmod_wolflion, topmod_lionbear, topmod_lionbob) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    dplyr::select(c(Modnames, DIC, Delta_DIC, DICWt)) %>%
    #'  Rename model and species pairing
    #'  Split model name based on placement of multiple periods
    #'  https://stackoverflow.com/questions/26265400/use-regex-in-r-to-retrieve-string-before-second-occurence-of-a-period
    mutate(Species_pair = sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", Modnames), 
           Species_pair = str_replace(Species_pair, "\\.", " - "), 
           Species_pair = str_replace(Species_pair, "bear", "Black bear"),
           Species_pair = str_replace(Species_pair, "lion", "Mountain lion"),
           Species_pair = str_replace(Species_pair, "coy", "Coyote"),
           Species_pair = str_replace(Species_pair, "bob", "Bobcat"),
           Species_pair = str_replace(Species_pair, "wolf", "Wolf"),
           Model = gsub(".*null", "Model 1", Modnames), 
           Model = gsub(".*habx", "Model 5", Model), 
           Model = gsub(".*preyabundx", "Model 6", Model),
           Model = gsub(".*preydivx", "Model 7", Model),
           Model = gsub(".*hab", "Model 2", Model), 
           Model = gsub(".*preyabund", "Model 3", Model), 
           Model = gsub(".*preydiv", "Model 4", Model), 
           Model = gsub(".*global", "Model 8", Model),
           Model_name = gsub(".*null", "Null", Modnames),
           Model_name = gsub(".*habx", "Habitat with interaction", Model_name),
           Model_name = gsub(".*preyabundx", "Prey abundance with interaction", Model_name),
           Model_name = gsub(".*preydivx", "Prey diversity with interaction", Model_name),
           Model_name = gsub(".*hab", "Habitat, no interaction", Model_name),
           Model_name = gsub(".*preyabund", "Prey abundance, no interaction", Model_name),
           Model_name = gsub(".*preydiv", "Prey diversity, no interaction", Model_name),
           Model_name = gsub(".*global", "Global", Model_name)) %>%
    relocate(Species_pair, .before = DIC) %>%
    relocate(Model, .after = Species_pair) %>%
    relocate(Model_name, .after = Model) %>%
    dplyr::select(-Modnames)
  colnames(model_list_DIC) <- c("Predator pair", "Model", "Model description", "DIC", "Delta DIC", "DIC Weight")
  
  #'  Save
  write.csv(topmodels, file = "./Outputs/DIC_top_models.csv")
  write.csv(model_list_DIC, file = "./Outputs/DIC_model_selection_results.csv")
  save(topmodels, file = "./Outputs/DIC_top_models.RData")
  save(model_list_DIC, file = "./Outputs/DIC_model_selection_results.RData")
  
  
    
    
    
    
  