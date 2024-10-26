  #'  ----------------------------------------
  #'  Wait time species ID model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate effect of ID of the most recently detected species on the time that
  #'  elapses before detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_sppID.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous species
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]]
      
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E ) where expected value comes from tbd_lambda
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species pair
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp])
      }
      
      #'  Mean TBD
      mu.tbd <- mean(spp.tbd[])
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
