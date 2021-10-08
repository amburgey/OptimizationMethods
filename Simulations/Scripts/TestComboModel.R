## need to add A, T, change K to Kcam and Kscr to jags statements


cat("
  model {
    
    ## Priors
    sigma ~ dunif(0,100)  ## both
    alpha1 <- 1/(2*sigma*sigma)  ## scr
    lam0 ~ dunif(0,5)  ## cam
    psi ~ dunif(0,1)  ## cam
    
    
    ## For SCR study:
    for(l in 1:L){   # 4 size categories
      #prior for intercept
      p0[l] ~ dunif(0,1)  ## for trapping, binomial; for visual survey, poisson
      alpha0[l] <- logit(p0[l])
      
      # Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
      Ngroup[l] <- ngroup[l] + n0[l]
    }
    
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # Gpts = number of points on integration grid
        for(j in 1:J){  # J = number of traps
          # Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),Kscr)
        } #J
        pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,]) #Prob of detect each size category across entire study area and time period
        pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
      } #G
      pstar[l] <- (sum(pdot[l,1:Gpts]*a))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
      
      # Zero trick for initial 1/pstar^n
      loglikterm[l] <- -ngroup[l] * log(pstar[l])
      lambda[l] <- -loglikterm[l] + 10000
      dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
    } #L
    
    # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    for(i in 1:n){  ## n = number of observed individuals
      ## For use when defining traps (scr and camera) on a grid cell
      s[i] ~ dcat(pi[1:Gpts])
      
      # SCR observations:
      for(j in 1:J){  ## J = number of traps
        y[i,j] ~ dbin(p[i,j],Kscr)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#n
    
    
    ## For Camera study:
    for(j in 1:T) {
      ## expected captures/occasion at camera t
      bigLambda[t] <- sum(lam[,t])
      n[t] ~ dpois(bigLambda[t]*Kcam[t]) ## camera active/inactive
    }
    
    # Camera trapping observations:
    for(i in 1:M){
      z[i] ~ dbern(psi)
      s[i] ~ dcat(pi[1:Gpts])  ## need to do again if have above? Or just do once here for all M?
      
      for(t in 1:T){  ## T = number of camera traps
        d2[i,t] <- pow(Gdist[s[i]],2) + pow(Gdist[s[i]],2)
        lam[i,t] <- lam0*exp(-(d2[i,t])/(2*sigma*sigma))*z[i]
      }#T
    }#M
    
    #derived proportion in each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
        
    Nscr <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N from SCR study
    Ncam <- sum(z[])  # population size from camera study
    Ntot <- Ncam + Nscr
    RD <- Ntot/A
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimTRAPCAM.txt")