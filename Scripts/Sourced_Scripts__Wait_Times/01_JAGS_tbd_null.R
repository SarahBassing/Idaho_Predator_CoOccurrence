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
      }
      
      #'  Derived parameters
      #'  ------------------
      mu.tbd <- exp(alpha0)
      
      } ")
  