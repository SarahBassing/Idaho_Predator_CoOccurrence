  #'  ----------------------------------------
  #'  Wait time prey diversity model
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ----------------------------------------
  #'  Estimate effect of prey species diversity index on the time that elapses 
  #'  before detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/tbd_preydiversity.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)

      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.div*covs[i,9]
      }
      
    
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey diversity value
      mu.tbd <- exp(alpha0 + beta.div*0)
    
      #'  Mean TBD across range of prey diversity values
      for(i in 1:100){
        spp.tbd.div[i] <- exp(alpha0 + beta.div*newcovs[i,5])
      }
      
      } ")
