  #'  ----------------------------------------
  #'  Wait time null model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_null.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0
      
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      
      }
      
      #'  Derived parameters
      #'  ------------------
      mu.tbd <- exp(alpha0)
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[])
      
      } ")
  