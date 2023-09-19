  #'  ----------------------------------------
  #'  Wait time species ID + prey RIA (elk, moose, wtd) model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate effect of ID of the most recently detected species and relative
  #'  abundance of elk/moose/wtd on the time that elapses before detecting a focal 
  #'  predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_sppID_elk_moose_wtd_abundance.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous competitor
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,6] + beta.prey[3]*covs[i,7] 
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0)
      }
      
      #'  Mean TBD per competitor across range of elk/moose/wtd relative abundance values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.elk[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + beta.prey[3]*0)
          spp.tbd.moose[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,2] + beta.prey[3]*0)
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*newcovs[i,3])
        }
      }
      
      #'  Mean TBD
      mu.tbd <- mean(spp.tbd[])
      
      } ")
