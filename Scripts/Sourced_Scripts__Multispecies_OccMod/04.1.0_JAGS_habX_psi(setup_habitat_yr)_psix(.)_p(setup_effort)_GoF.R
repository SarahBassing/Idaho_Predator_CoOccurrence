  #'  ------------------------------------
  #'  Habitat model with interaction
  #'  
  #'  Bassing et al. "Mammalian predator co-occurrence affected by prey and habitat 
  #'  more than competitor presence at multiple time scales"
  #'  ------------------------------------
  #'  Includes basic habitat variables (elevation and percent forest cover), year,
  #'  and camera setup on first-order psi. Setup and sampling effort on first-order rho.
  #'  Allows species co-occurrence to be non-independent. Second half of script 
  #'  generates Chi-squared and Freeman-Tukey test statistics for Goodness-of-Fit 
  #'  test for detection probability submodel.
  #'  ------------------------------------
  
  cat(file = './Outputs/04.1.0_JAGS_habX_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.txt', "
      model{
          
        #### Define Priors  ####
        #'  ================
        #'  First order occupancy intercerpts (psi) 
        betaSpp1[1] <- logit(mean.psiSpp1)           
        betaSpp2[1] <- logit(mean.psiSpp2)
        mean.psiSpp1 ~ dunif(0, 1)               
        mean.psiSpp2 ~ dunif(0, 1)
            
        #'  First order occupancy slopes (psi)
        for(fo_psi in 2:6){                         
          betaSpp1[fo_psi] ~ dnorm(0, 0.1)
          betaSpp2[fo_psi] ~ dnorm(0, 0.1)
        }
      
        #'  Second order occupancy intercept (psix)
        betaSpp12[1] ~ dnorm(0, 0.1)
        
        #'  First order detection intercepts (rho)
        alphaSpp1[1] <- logit(mean.pSpp1)           
        alphaSpp2[1] <- logit(mean.pSpp2)
        mean.pSpp1 ~ dunif(0, 1)                    
        mean.pSpp2 ~ dunif(0, 1)
         
        #'  First order detection slopes (rho)   
        for(fo_rho in 2:3){                         
          alphaSpp1[fo_rho] ~ dnorm(0, 0.1)  
          alphaSpp2[fo_rho] ~ dnorm(0, 0.1)
        }
      
        #'  Second order detection priors (rho)
        #'  Assumes no second-order interactions by setting these to 0
        alphaSpp12 <- 0
        alphaSpp21 <- 0
        
            
        ####  Define Likelihood  ####
        #'  =====================
        #'  1. Set up basic hierarchical model
      
        #'  Latent state model
        #'  ------------------
        #'  For each site, true occupancy (z) is drawn from a categorical distribution
        #'  with 4 mututally exclusive occupancy probabilities
            
        for(i in 1:nsites) {
          z[i] ~ dcat(lsv[i, (1:ncat)])
        }
          
        #'  Observation model
        #'  -----------------
        #'  For each site and survey occasion, the deteciton data are drawn from a
        #'  categorical distribution with 4 latent states (z)
          
        for(i in 1:nsites) {
          for(j in 1:nsurveys) {
            y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
          }
        }
          
        #'  2. Define arrays containing cell probabilities for categorical distributions
        #'  Probabilities for each latent state (z), held in latent state vector (lsv)
        for(i in 1:nsites) {
          lsv[i, 1] <- 1                   # Unoccupied
          lsv[i, 2] <- exp(psiSpp1[i])     # Pr(Spp1 present)
          lsv[i, 3] <- exp(psiSpp2[i])     # Pr(Spp2 present)
          lsv[i, 4] <- exp(psiSpp12[i])    # Pr(Spp1 & Spp2 present)
         
          #'  Probabilities for each detection array, held in rho detection matrix (rdm) 
          #'  where OS = observed state, TS = true state and each row sums to 1. 
          #'  Exponentiating log odds so rdm holds estimates on probability scale.
          for(j in 1:nsurveys) {
            #'  True state = unoccupied (z = 1 --> 00)
            rdm[i, j, 1, 1] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 1] <- 0 # ------------------------------------ OS = Spp1 present
            rdm[i, j, 3, 1] <- 0 # ------------------------------------ OS = Spp2 present
            rdm[i, j, 4, 1] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp1 present (z = 2 --> 10)
            rdm[i, j, 1, 2] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 2] <- exp(rhoSpp1[i, j]) # ------------------- OS = Spp1 present
            rdm[i, j, 3, 2] <- 0 # ------------------------------------ OS = Spp2 present
            rdm[i, j, 4, 2] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp2 present (z = 3 --> 01)
            rdm[i, j, 1, 3] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 3] <- 0 # ------------------------------------ OS = Spp1 present
            rdm[i, j, 3, 3] <- exp(rhoSpp2[i, j]) # ------------------- OS = Spp2 present
            rdm[i, j, 4, 3] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp1 & Spp2 present (z = 4 --> 11)
            rdm[i, j, 1, 4] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 4] <- exp(rhoSpp12[i, j]) # ------------------ OS = Spp1 present
            rdm[i, j, 3, 4] <- exp(rhoSpp21[i, j]) # ------------------ OS = Spp2 present
            rdm[i, j, 4, 4] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
          }
              
          #'  3. Define linear models for each fundamental parameter that governs the cell probs
          #'  These are my natural parameters (f1, f2, f12)!
          #'  Linear models for the occupancy parameters on the logit scale
              
          #'  ...for states Spp1, Spp2
          #'  Covariate order: Spp1 = Intercept[1] + Setup[2] + Year[3] + Forest[4] + Elevation[5] + TRI[10] 
          #'  Covariate order: Spp2 = Intercept[1] + Setup[2] + Year[3] + Forest[4] + Elevation[5] + TRI[10]  
          psiSpp1[i] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4] + betaSpp1[5]*psi_cov[i,5] + betaSpp1[6]*psi_cov[i,10]  
          psiSpp2[i] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4] + betaSpp2[5]*psi_cov[i,5] + betaSpp2[6]*psi_cov[i,10] 
          
          #'  ...for state Spp12
          #'  Covariate order: Intercept[1] 
          psiSpp12[i] <- psiSpp1[i] + psiSpp2[i] + betaSpp12[1]*psi_inxs_cov[i,1] 
          
          #'  Baseline linear predictors for detection
          #'  Covariate indices: Intercept[1] + Setup[2] + Sampling Effort[3]
          for(j in 1:nsurveys) {
            rhoSpp1[i, j] <- alphaSpp1[1]*rho_cov[i,j,1] + alphaSpp1[2]*rho_cov[i,j,2] + alphaSpp1[3]*rho_cov[i,j,3] 
            rhoSpp2[i, j] <- alphaSpp2[1]*rho_cov[i,j,1] + alphaSpp2[2]*rho_cov[i,j,2] + alphaSpp2[3]*rho_cov[i,j,3] 
      
            #'  Asymetric interactions between both species
            #'  Fixing to be same as species-sepcific detection probability
            rhoSpp12[i, j] <- rhoSpp1[i, j] 
            rhoSpp21[i, j] <- rhoSpp2[i, j] 
          }
        }
            
        #'  For Goodness-of-Fit test
        for(i in 1:nsites) {
          for(j in 1:nsurveys) {
            #'  Draw a replicate data set under fitted model
            y.sim[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]]) 
            
            #'  Expected detection probability based each underlying occurrence state
            #'  Need to normalized rdms because dcat can take values that sum to >1
            #'  rdm is a vector of unnormalized non-negative probability *weights* of length N
            #'  dcat pmf --> P(Y = y) = pi[y]/sum[pi[i, 1:N]] (pg. 50 JAGS manual)
            #'  True state = unoccupied (z = 1 --> 00)
            p[i, j, 1, 1] <- rdm[i, j, 1, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(unoccupied)
            p[i, j, 2, 1] <- rdm[i, j, 2, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp1 detected)
            p[i, j, 3, 1] <- rdm[i, j, 3, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp2 detected)
            p[i, j, 4, 1] <- rdm[i, j, 4, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp12 detected)
            #'  True state = Spp1 present (z = 2 --> 10)
            p[i, j, 1, 2] <- rdm[i, j, 1, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(unoccupied)
            p[i, j, 2, 2] <- rdm[i, j, 2, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp1 detected)
            p[i, j, 3, 2] <- rdm[i, j, 3, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp2 detected)
            p[i, j, 4, 2] <- rdm[i, j, 4, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp12 detected)
            #'  True state = Spp2 present (z = 3 --> 01)
            p[i, j, 1, 3] <- rdm[i, j, 1, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(unoccupied)
            p[i, j, 2, 3] <- rdm[i, j, 2, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp1 detected)
            p[i, j, 3, 3] <- rdm[i, j, 3, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp2 detected)
            p[i, j, 4, 3] <- rdm[i, j, 4, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp12 detected)
            #'  True state = Spp1 & Spp2 present (z = 4 --> 11)
            p[i, j, 1, 4] <- rdm[i, j, 1, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(unoccupied)
            p[i, j, 2, 4] <- rdm[i, j, 2, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp1 detected)
            p[i, j, 3, 4] <- rdm[i, j, 3, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp2 detected)
            p[i, j, 4, 4] <- rdm[i, j, 4, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp12 detected)
      
            y.hat[i,j,1:4] <- p[i, j, 1:4, z[i]]
            #'  Rank index position of p in descending order (i.e., identify index position of largest to smallest value of p) 
            #'  Then order data by rank (i.e, order index position based on ranking from highest to lowest value of p)
            y.hat.index[i, j, 1:4] <- order(-rank(p[i, j, 1:4, z[i]]))
            #'  Retain index position of the highest value of p
            y.hat.maxindex[i, j] <- y.hat.index[i, j, 1]

            #'  Calculate Chi-squared test statistic for observed & simulated data sets
            x2[i,j] <- pow((y[i,j] - y.hat.maxindex[i,j]), 2) / (y.hat.maxindex[i,j] + 0.0001)
            x2.sim[i,j] <- pow((y.sim[i,j] - y.hat.maxindex[i,j]), 2) / (y.hat.maxindex[i,j] + 0.0001)
            
            #'  Calculate Freeman-Tukey test statistic for observed & simulated data sets
            ft[i,j] <- pow((sqrt(y[i,j]) - sqrt(y.hat.maxindex[i,j])), 2) 
            ft.sim[i,j] <- pow((sqrt(y.sim[i,j]) - sqrt(y.hat.maxindex[i,j])), 2)
          }
      
        #'  Sum across surveys
        x2.obs[i] <- sum(x2[i,])
        x2.sims[i] <- sum(x2.sim[i,])
        
        ft.obs[i] <- sum(ft[i,])
        ft.sims[i] <- sum(ft.sim[i,])
  
        }
      
        #'  Sum across sites for final Chi-squared & Freeman-Tukey test statistics
        chi2.obs <- sum(x2.obs[])
        chi2.sim <- sum(x2.sims[])
        
        FT.obs <- sum(ft.obs[])
        FT.sims <- sum(ft.sims[])
      
      }
      ", fill=TRUE)
