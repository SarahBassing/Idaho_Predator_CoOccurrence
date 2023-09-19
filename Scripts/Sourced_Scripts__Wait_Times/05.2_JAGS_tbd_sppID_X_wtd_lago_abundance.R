  #'  ----------------------------------------
  #'  Wait time species ID * prey RIA (wtd & lagomorph) model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate interactive effect between ID of the most recently detected species 
  #'  and relative abundance of wtd/lagomorph on the time that elapses before 
  #'  detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_sppID_X_wtd_lago_abundance.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous sppID
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between sppID and prey relative abundance
      beta.interaction.wtd[1] <- 0
      beta.interaction.lago[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction.wtd[spp] ~ dnorm(0, 0.01)
        beta.interaction.lago[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,7] + beta.prey[2]*covs[i,8] + 
                          beta.interaction.wtd[covs[i,2]]*covs[i,7] + beta.interaction.lago[covs[i,2]]*covs[i,8] 
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species ID at average prey RAI
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + 
                             beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*0)
      }
      
      #'  Mean TBD per species ID across range of wtd/lagomorph relative abundance values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,3] + beta.prey[2]*0 + 
                                     beta.interaction.wtd[spp]*newcovs[i,3] + beta.interaction.lago[spp]*0)
          spp.tbd.lago[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,4] + 
                                     beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*newcovs[i,4])
        }
      }
      
      
      #'  Mean TBD
      mu.tbd <- mean(spp.tbd[])
      
      } ")
