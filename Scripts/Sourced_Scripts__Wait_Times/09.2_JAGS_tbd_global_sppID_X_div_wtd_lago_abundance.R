  #'  ----------------------------------------
  #'  Wait time global (wtd & lagomorph) model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate interactive effect between ID of the most recently detected species 
  #'  relative abundance of wtd/lagomorph and prey species diversity index on the 
  #'  time that elapses before detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_global_sppID_X_div_wtd_lago_abundance.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous species
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between competitor and prey diversity/relative abundance
      beta.interaction[1] <- 0
      beta.interaction.wtd[1] <- 0
      beta.interaction.lago[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction[spp] ~ dnorm(0, 0.01)
        beta.interaction.wtd[spp] ~ dnorm(0, 0.01)
        beta.interaction.lago[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,7] + beta.prey[2]*covs[i,8] + beta.div*covs[i,9] +
                          beta.interaction.wtd[covs[i,2]]*covs[i,7] + beta.interaction.lago[covs[i,2]]*covs[i,8] + beta.interaction[covs[i,2]]*covs[i,9]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI and prey diversity
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.div*0 +
                             beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*0 + beta.interaction[spp]*0)
      }
      
      #'  Mean TBD per competitor across range of wtd/lagomorph relative abundance values & prey diversity
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,3] + beta.prey[2]*0 + beta.div*0 +
                                     beta.interaction.wtd[spp]*newcovs[i,3] + beta.interaction.lago[spp]*0 + beta.interaction[spp]*0)
          spp.tbd.lago[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,4] + beta.div*0 +
                                      beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*newcovs[i,4] + beta.interaction[spp]*0)
          spp.tbd.div[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.div*newcovs[i,5] +
                                     beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*0 + beta.interaction[spp]*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      mu.tbd <- mean(spp.tbd[])
      
      } ")
