  #'  ----------------------------------------
  #'  Wait time sppID + prey diversity effect model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate effect of ID of the most recently detected species and prey
  #'  species diversity index on the time that elapses before detecting a focal 
  #'  predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_sppID_preydiversity.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous sppID
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }

      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.div*covs[i,9]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.div*0)
      }
      
      #'  Mean TBD per species across range of prey diversity values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.div[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.div*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      mu.tbd <- mean(spp.tbd[])
      
      } ")
