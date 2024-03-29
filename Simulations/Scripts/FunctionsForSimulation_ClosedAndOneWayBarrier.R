#### FUNCTIONS TO SIMULATE DATA FOR BROWN TREESNAKE MONITORING OPTIMIZATION ####
## These functions are used in BaseSimulationFramework.R and create and calculate information needed to run alternate monitoring scenarios
## Study area from which real datasets were collected was the Closed Population (CP) on Guam

library(coda);library(runjags);library(sp);library(rgeos);library(raster);library(scales)

#### FUNCTION TO SPECIFY STUDY AREA TYPE.----

areatype <- function(totlocs){
  ## Define state-space of point process. (i.e., where animals live).
  ## Approximate dimensions of CP
  delta <- 11.874929
  Xl <- min(totlocs[,1]) - delta
  Xu <- max(totlocs[,1]) + delta
  Yl <- min(totlocs[,2]) - delta
  Yu <- max(totlocs[,2]) + delta
  
  sinfo <- list(Xl=Xl,Xu=Xu,Yl=Yl,Yu=Yu)
  
  return(sinfo)
}


#### FUNCTION TO FIND DISTANCES BETWEEN INTEGRATION GRID POINTS AND SAMPLING POINTS.----

e2dist <- function(x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


#### FUNCTION TO CREATE INDIVIDUAL ACTIVITY CENTERS.----
## True snake activity centers (AC), pull new values for every simulation
ActCent <- function(Gpts,N){
  s <- sample(1:Gpts,N,replace=TRUE) # pull value from integration grid locations
  
  return(s)
}


#### FUNCTION TO CREATE ONE-WAY BARRIER AND CALCULATE PROBABILITY OF LEAVING.----

ProbStay <- function(sigma, G, s, a){
  
  #### GENERATE HOME RANGE AREA FROM SIGMA.----
  #Scale parameter (sigma) can be interpreted as range size (for half-normal detection function, 95% of activity is within 2.45*sigma circle)
  rangeS <- circular.r(p=0.99,detectfn = "HHN")*sigma
  radius <- sqrt(rangeS/pi)
  locAC <- G[s,] #locations of ACs on integration grid
  
  ## Get information to plot square integration grid cells
  dfG <- as.data.frame(G)
  A2 <- cbind(dfG$Xlocs,dfG$Ylocs,rep(a,nrow(dfG)))
  colnames(A2) <- c("X","Y","Area")
  B <- SpatialPoints(A2)
  ww <- sqrt(B$Area)/2  ## Widths of buffers needed to produce desired areas    
  #Find all 5 points of each square
  pp <- list()
  for(i in seq_along(B)) {
    pp[[i]] <- gBuffer(B[i], width=ww[i], quadsegs=1, capStyle="SQUARE")
  }
  PP <- do.call(bind, pp)
  
  #The only snakes that will ever leave are those with their home ranges on the edge of the study area (the home range size is never large enough to cross the barrier from farther away)
  #Find individuals whose home ranges are outside the barrier based on centroid locations being along the edge of the study area
  circHalf <- as.numeric(row.names(subset(as.data.frame(locAC), Ylocs == G[1,2] | Ylocs == G[528,2] | Xlocs == G[1,1] | Xlocs == G[528,1])))
  
  #### PROBABILITY OF SNAKE LEAVING STUDY AREA.----
  #Find intersection with home ranges and boundaries of CP
  #Line coordinates (x1,y1) and (x2,y2); horizontal line is just y intercept, for vertical line is just x intercept
  #https://www.engineersedge.com/math/circular_segment_equation_and_calculator__13796.htm
  #For edge of CP, integration grid points pt1 = G[1,1:2], pt2 = G[22,1:2], pt3 = G[507,1:2], pt4 = G[528,1:2]
  ## Check if snake is in corner of study area as this will double the home range area potentially outside of the barrier
  doubleprob <- vector()
  for(i in 1:nrow(locAC)){
    doubleprob[i] <- ifelse(any(
      (locAC[i,1] == G[1,1] & locAC[i,2] == G[1,2]) |
        (locAC[i,1] == G[22,1] & locAC[i,2] == G[22,2]) |
        (locAC[i,1] == G[507,1] & locAC[i,2] == G[507,2]) |
        (locAC[i,1] == G[528,1] & locAC[i,2] == G[528,2])
      == TRUE), 2, 1)
  }

  #The home range will be the same size for all individuals in a simulation so only need to calculate this once per sigma value
  #Find a snake that has a home range on the edge of the study area on which to calculate area outside the barrier
  #Determine which side of the study area the example snake's home range is on
  circHalfSide <- ifelse(locAC[circHalf[1],1] == G[1,1], "1.1",              # check if overlapping circle matches left side of barrier
                         ifelse(locAC[circHalf[1],1] == G[528,1], "528.2",   # check if overlapping circle matches right side of barrier
                                ifelse(locAC[circHalf[1],2] == G[1,2], "1.3",        # check if overlapping circle matches bottom side of barrier
                                       ifelse(locAC[circHalf[1],2] == G[528,2], "528.4", "-9999"))))     # check if overlapping circle matches top side of barrier
  #split this identifier into a dataframe of which side is overlapping and the xmin (1), xmax (2), ymin (3), or ymax (4) coordinate needed
  if(is.na(circHalfSide) == FALSE){  # if no snakes near edge
    circHalfSide <- as.numeric(unlist(strsplit(circHalfSide, ".", fixed = TRUE)))
  
    #d is the distance between a home range centroid x or y-coordinate and the equation of a vertical (x) or horizontal study barrier (y) cutting through the home range
    #take absolute value as sometimes barrier is above or below the circle centroid
    d <- abs(locAC[circHalf[1],ifelse(circHalfSide[2] == 1 |circHalfSide[2] == 2,1,2)] - extent(pp[[circHalfSide[1]]])[circHalfSide[2]])
    #h is the distance between d (the barrier cutting through the home range) and the circular edge of the home range
    h <- radius - d
    #c is the length of the line segment within the home range
    c <- sqrt((radius-(h/2))*8*h)
    if(is.na(c) == TRUE){                 # if the home range is so small that a snake never would cross the barrier
      stay <- array(1, dim=c(N,K,nsims))  # fill the staying matrix with 1s to indicate snake always present
    } else {                              # if the home range is big enough to allow snake to cross the barrier and leave
      #theta is the interior angle of the segment within the circle
      theta <- 2*atan(c/(2*d))
      
      #The proportion of home range outside of the barrier to inside dictates the probability of the snake leaving
      #https://www.mathsisfun.com/geometry/circle-sector-segment.html
      #Area of circle outside of circle (when theta is in radians)
      areaseg <- ((theta-sin(theta))/2) * (radius)^2
      probleave <- areaseg/rangeS
      
      #Do Bernoulli trial of whether snake stays or snake leaves based on this probability
      #If snake leaves, it can never come back in
      stay <- array(NA, dim=c(N,K,nsims))
      for(z in 1:nsims){
        for(k in 1:K){                          #for each occasion sampled
          for(n in 1:N){                        #for each snake
            if(is.na(stay[n,K,z])){             #conditional on the snake still being in the study area
              if(any(n == circHalf) == TRUE){   #conditional on the snake's home range overlapping the fence
                stay[n,k,z] <- rbinom(1,1,(1-(probleave*doubleprob[n])))  #prob the snake stays given the proportion of the home range area outside the barrier
                if(stay[n,k,z] == 0){
                  stay[n,k:K,z] <- 0
                  next                          #If the snake leaves, fill in zeroes and skip for rest of occasions
                }
              } else {
                stay[n,k,z] <- 1
              }
            } else {
              next
            }
          }
        }
      }
    }
  } else {
    stay <- array(1, dim=c(N,K,nsims))  # fill the staying matrix with 1s to indicate snake always present
  }
  info <- list(stay=stay)
  
  return(info)
}


#### FUNCTION TO CREATE COMBINED SIGMA AND SIMULATE OBSERVATIONS DEPENDING ON THE METHOD SELECTED.----

createData <- function(){
  
  ## CREATE COMBINED SIGMA FROM MODEL RESULTS FROM VIS AND TRAP ANALYSIS
  
  #### VISUAL SURVEY REAL DATA RESULTS ####
  load("Real Data Analysis/Visual surveys/Results/NWFNVISALL_SCRpstar.Rdata")  # unified analysis of all datasets
  
  ## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
  mcmc.object <- as.mcmc.list(out$samples[[1]])
  out1 <- window(mcmc.object, start=1499)
  mcmc.object <- as.mcmc.list(out$samples[[2]])
  out2 <- window(mcmc.object, start=1499)
  mcmc.object <- as.mcmc.list(out$samples[[3]])
  out3 <- window(mcmc.object, start=1499)
  outALLV <- combine.mcmc(list(out1,out2,out3), thin = 1)
  niter(outALLV)
  postV <- as.matrix(outALLV)
  
  ## Summary info for model results
  #for median and CI
  post_sumV <- data.frame(
    mean = apply(postV, 2, function(x) mean(x)),
    med = apply(postV, 2, function(x) quantile(x, probs = 0.5, na.rm = T, names = F)),
    lower = apply(postV, 2, function(x) quantile(x, probs = 0.025, na.rm = T, names = F)),
    upper = apply(postV, 2, function(x) quantile(x, probs = 0.975, na.rm = T, names = F)),
    sd = apply(postV, 2, function(x) sd(x)))
  post_sumV$variable <- row.names(post_sumV)
  
  
  #### TRAPPING SURVEY REAL DATA RESULTS.----
  load("Real Data Analysis/Trapping/Results/NWFNTRAPALL_SCRpstar.Rdata")
  
  ## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
  mcmc.object <- as.mcmc.list(out$samples[[1]])
  out1 <- window(mcmc.object, start=1499)
  mcmc.object <- as.mcmc.list(out$samples[[2]])
  out2 <- window(mcmc.object, start=1499)
  mcmc.object <- as.mcmc.list(out$samples[[3]])
  out3 <- window(mcmc.object, start=1499)
  outALLT <- combine.mcmc(list(out1,out2,out3), thin = 1)
  niter(outALLT)
  postT <- as.matrix(outALLT)
  
  ## Summary info for model results
  #for median and CI
  post_sumT <- data.frame(
    mean = apply(postT, 2, function(x) mean(x)),
    med = apply(postT, 2, function(x) quantile(x, probs = 0.5, na.rm = T, names = F)),
    lower = apply(postT, 2, function(x) quantile(x, probs = 0.025, na.rm = T, names = F)),
    upper = apply(postT, 2, function(x) quantile(x, probs = 0.975, na.rm = T, names = F)),
    sd = apply(postT, 2, function(x) sd(x)))
  post_sumT$variable <- row.names(post_sumT)
  
  
  #### CREATE JOINT DISTRIBUTION FOR SIGMA.----
  d1 <- rnorm(1000000, subset(post_sumV,variable=="sigma")$mean, subset(post_sumV,variable=="sigma")$sd)
  d2 <- rnorm(1000000, subset(post_sumT,variable=="sigma")$mean, subset(post_sumT,variable=="sigma")$sd)
  alld <- cbind(d1,d2)
  sigma_mu <- mean(alld)
  sigma_sd <- sd(alld)
  
  
  if(stype == "closed"){
  
    if(type == "VIS"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR VIS.----
      p0MV <- list(sample(postV[,"p0[1]"]),sample(postV[,"p0[2]"]),sample(postV[,"p0[3]"]),sample(postV[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
     
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      
      ## VISUAL SURVEYS - Generate observations
      yTrueVIS <- array(NA,dim=c(N,J,nsims))
      p0V <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 5, dimnames = list(1:nsims,c("Sigma","p01","p02","p03","p04")))
      
      for(z in 1:nsims){
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        
        for(l in 1:length(Ngroup)){
          p0V[l] <- p0MV[[l]][sample(1:length(p0MV[[l]]), 1)]   ## pull value from posterior
        }
        pmatV <- p0V[Nsnsz]*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
        for(n in 1:N){
          yTrueVIS[n,,z] <- rpois(J,pmatV[n,]*K)  # observations of each snake at each trap based on encounter probability and effort
        }
        ## Save all values of sigma and p0V used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0V
      }
      
      ## yTrueVIS includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedV <- list()
      yarrV <- list()
      ## Bind sizes of snake to encounter histories
      yTrueVIS <- abind::abind(yTrueVIS,array(Nsnsz, replace(dim(yTrueVIS),2,1)), along=2)
      for(z in 1:nsims){
        capturedV[[z]] <- which(apply(yTrueVIS[,1:J,z],1,sum)>0)  # snakes that were observed at least once
        yarrV[[z]] <- yTrueVIS[capturedV[[z]],,z]  # subset to observed snakes
        write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,stype,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## VIS
  
  
    if(type == "TRAP"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR TRAP.----
      p0MT <- list(sample(postT[,"p0[1]"]),sample(postT[,"p0[2]"]),sample(postT[,"p0[3]"]),sample(postT[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      
      ## TRAPPING - Generate observations
      yTrueTRAP <- array(NA,dim=c(N,J,nsims))
      p0T <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 5, dimnames = list(1:nsims,c("Sigma","p01","p02","p03","p04")))
      
      for(z in 1:nsims){  
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        
        for(l in 1:length(Ngroup)){
          p0T[l] <- p0MT[[l]][sample(1:length(p0MT[[l]]), 1)]   ## pull value from posterior
        }
        pmatT <- p0T[Nsnsz]*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
        for(n in 1:N){
          yTrueTRAP[n,,z] <- rbinom(J,K,pmatT[n,])  # observations of each snake at each trap based on encounter probability and effort
        }
        ## Save all values of sigma and p0T used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0T
      }
      
      ## yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedT <- list()
      yarrT <- list()
      ## Bind sizes of snake to encounter histories
      yTrueTRAP <- abind::abind(yTrueTRAP,array(Nsnsz, replace(dim(yTrueTRAP),2,1)), along=2)
      for(z in 1:nsims){
        capturedT[[z]] <- which(apply(yTrueTRAP[,1:J,z],1,sum)>0)  # snakes that were observed at least once
        yarrT[[z]] <- yTrueTRAP[capturedT[[z]],,z]  # subset to observed snakes
        write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,stype,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## TRAP
  
  
    if(type == "VISTRAP"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR VIS AND TRAP.----
      p0MV <- list(sample(postV[,"p0[1]"]),sample(postV[,"p0[2]"]),sample(postV[,"p0[3]"]),sample(postV[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      p0MT <- list(sample(postT[,"p0[1]"]),sample(postT[,"p0[2]"]),sample(postT[,"p0[3]"]),sample(postT[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
      
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      
      ## VISUAL SURVEYS - Generate observations
      yTrueVIS <- array(NA,dim=c(N,J1,K,nsims))
      yTrueTRAP <- array(NA,dim=c(N,J2,K,nsims))
      yTrueVIS2 <- array(NA,dim=c(N,J1,nsims))
      yTrueTRAP2 <- array(NA,dim=c(N,J2,nsims))
      p0V <- vector()
      p0T <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 9, dimnames = list(1:nsims,c("Sigma","p0V1","p0V2","p0V3","p0V4","p0T1","p0T2","p0T3","p0T4")))
      
      for(z in 1:nsims){  
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        
        for(l in 1:length(Ngroup)){
          p0V[l] <- p0MV[[l]][sample(1:length(p0MV[[l]]), 1)]   ## pull value from posterior
          p0T[l] <- p0MT[[l]][sample(1:length(p0MT[[l]]), 1)]   ## pull value from posterior
        }
        
        pmatV <- p0V*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all surveying locations
        pmatT <- p0T*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all trapping locations
        
        for(n in 1:N){
          for(k in 1:K){
            yTrueTRAP[n,,k,z] <- rbinom(J2,1,pmatT[n,newX2])    # observations of each snake at each trap based on encounter probability and effort
            if(any(yTrueTRAP[n,,k,z] >= 1) == TRUE){
              yTrueVIS[n,,k,z] <- 0
              next
            } else {
              yTrueVIS[n,,k,z] <- rpois(J1,pmatV[n,newX1])  # observations of each snake at each survey location based on encounter probability and effort conditional on if snake was captured in a trap that evening
            }
          }
          for(j in 1:J2){
            ## Collapse observations across all occasions to give a count of observations per trap for each simulation
            yTrueTRAP2[n,j,z] <- sum(yTrueTRAP[n,j,,z])
          }
          for(j in 1:J1){
            yTrueVIS2[n,j,z] <- sum(yTrueVIS[n,j,,z])
          }
        }
        ## Save all values of sigma and p0V/T used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0V
        paramsim[z,6:9] <- p0T
      }
      
      # CHECK TO MAKE SURE NO OCCASION HAS A SNAKE TRAPPED AND VISUALLY SEEN
      for(z in 1:nsims){
        for(k in 1:K){
          ifelse(any((rowSums(yTrueTRAP[,,k,z]) & rowSums(yTrueVIS[,,k,z])) == TRUE) == TRUE , print("oh no"), print("we good"))
        }
      }
        
      ## yTrueVIS and yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedV <- list()
      capturedT <- list()
      yarrV <- list()
      yarrT <- list()
      ## Bind sizes of snake to encounter histories
      yTrueVIS3 <- abind::abind(yTrueVIS2,array(Nsnsz, replace(dim(yTrueVIS2),2,1)), along=2)
      yTrueTRAP3 <- abind::abind(yTrueTRAP2,array(Nsnsz, replace(dim(yTrueTRAP2),2,1)), along=2)
      ## Bind unique snake identity to encounter histories
      yTrueVIS3 <- abind::abind(yTrueVIS3,array(seq(1:N), replace(dim(yTrueVIS3),2,1)), along=2)
      yTrueTRAP3 <- abind::abind(yTrueTRAP3,array(seq(1:N), replace(dim(yTrueTRAP3),2,1)), along=2)
      for(z in 1:nsims){
        capturedV[[z]] <- which(apply(yTrueVIS3[,1:J1,z],1,sum)>0)  # snakes that were observed at least once
        capturedT[[z]] <- which(apply(yTrueTRAP3[,1:J2,z],1,sum)>0)  # snakes that were observed at least once
        yarrV[[z]] <- yTrueVIS3[capturedV[[z]],,z]  # subset to observed snakes
        yarrT[[z]] <- yTrueTRAP3[capturedT[[z]],,z]  # subset to observed snakes
        write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,stype,"VIS",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
        write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,stype,"TRAP",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## VISTRAP
  } ## closed
  
  
  if(stype == "oneway"){
    
      if(type == "VIS"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR VIS.----
      p0MV <- list(sample(postV[,"p0[1]"]),sample(postV[,"p0[2]"]),sample(postV[,"p0[3]"]),sample(postV[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
      
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      
      ## VISUAL SURVEYS - Generate observations
      yTrueVIS <- array(NA,dim=c(N,J,K,nsims))
      yTrueVIS2 <- array(NA,dim=c(N,J,nsims))
      p0V <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 5, dimnames = list(1:nsims,c("Sigma","p01","p02","p03","p04")))
      
      for(z in 1:nsims){
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        siminfo <- ProbStay(sigma, G, s, a)
          
        for(l in 1:length(Ngroup)){
          p0V[l] <- p0MV[[l]][sample(1:length(p0MV[[l]]), 1)]   # pull value from posterior
        }
        
        pmatV <- p0V[Nsnsz]*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
        
        for(n in 1:N){  ## for each snake
          for(k in 1:K){  ## for each occasion sampled
            yTrueVIS[n,,k,z] <- rpois(J,pmatV[n,])*siminfo$stay[n,k,z]  # observations of each snake at each survey location based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
          }
          for(j in 1:J){
            ## Collapse observations across all occasions to give a count of observations per survey location for each simulation
            yTrueVIS2[n,j,z] <- sum(yTrueVIS[n,j,,z])
          }
        }
        ## Save all values of sigma and p0V used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0V
      }
      
      ## yTrueVIS includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedV <- list()
      yarrV <- list()
      ## Bind sizes of snake to encounter histories
      yTrueVIS3 <- abind::abind(yTrueVIS2,array(Nsnsz, replace(dim(yTrueVIS2),2,1)), along=2)
      for(z in 1:nsims){
        capturedV[[z]] <- which(apply(yTrueVIS3[,1:J,z],1,sum)>0)  # snakes that were observed at least once
        yarrV[[z]] <- yTrueVIS3[capturedV[[z]],,z]  # subset to observed snakes
        write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,stype,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## VIS
  
  
    if(type == "TRAP"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR TRAP.----
      p0MT <- list(sample(postT[,"p0[1]"]),sample(postT[,"p0[2]"]),sample(postT[,"p0[3]"]),sample(postT[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      
      ## TRAPPING - Generate observations
      yTrueTRAP <- array(NA,dim=c(N,J,K,nsims))
      yTrueTRAP2 <- array(NA,dim=c(N,J,nsims))
      p0T <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 5, dimnames = list(1:nsims,c("Sigma","p01","p02","p03","p04")))
      
      for(z in 1:nsims){  
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        siminfo <- ProbStay(sigma, G, s, a)
        
        for(l in 1:length(Ngroup)){
          p0T[l] <- p0MT[[l]][sample(1:length(p0MT[[l]]), 1)]   # pull value from posterior
        }
        
        pmatT <- p0T[Nsnsz]*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
        
        for(n in 1:N){
          for(k in 1:K){  ## for each occasion sampled
            yTrueTRAP[n,,k,z] <- rbinom(J,1,pmatT[n,])*siminfo$stay[n,k,z]  # observations of each snake at each trap based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
          }
          for(j in 1:J){
            ## Collapse observations across all occasions to give a count of observations per trap for each simulation
            yTrueTRAP2[n,j,z] <- sum(yTrueTRAP[n,j,,z])
          }
        }
        ## Save all values of sigma and p0T used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0T
      }
      
      ## yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedT <- list()
      yarrT <- list()
      ## Bind sizes of snake to encounter histories
      yTrueTRAP3 <- abind::abind(yTrueTRAP2,array(Nsnsz, replace(dim(yTrueTRAP2),2,1)), along=2)
      for(z in 1:nsims){
        capturedT[[z]] <- which(apply(yTrueTRAP3[,1:J,z],1,sum)>0)  # snakes that were observed at least once
        yarrT[[z]] <- yTrueTRAP3[capturedT[[z]],,z]  # subset to observed snakes
        write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,stype,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## TRAP
  
  
    if(type == "VISTRAP"){
      
      #### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR VIS AND TRAP.----
      p0MV <- list(sample(postV[,"p0[1]"]),sample(postV[,"p0[2]"]),sample(postV[,"p0[3]"]),sample(postV[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      p0MT <- list(sample(postT[,"p0[1]"]),sample(postT[,"p0[2]"]),sample(postT[,"p0[3]"]),sample(postT[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
      
      
      #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
      ## VISUAL AND TRAPPING SURVEYS - Generate observations
      yTrueVIS <- array(NA,dim=c(N,J1,K,nsims))
      yTrueTRAP <- array(NA,dim=c(N,J2,K,nsims))
      yTrueVIS2 <- array(NA,dim=c(N,J1,nsims))
      yTrueTRAP2 <- array(NA,dim=c(N,J2,nsims))
      p0V <- vector()
      p0T <- vector()
      paramsim <- matrix(NA, nrow = nsims, ncol = 9, dimnames = list(1:nsims,c("Sigma","p0V1","p0V2","p0V3","p0V4","p0T1","p0T2","p0T3","p0T4")))
      
      for(z in 1:nsims){ 
        
        #### GENERATE ACTIVTY CENTERS
        s <- ActCent(Gpts,N)
        
        #### PULL SIGMA VALUE FIRST IN ORDER TO CALCULATE SUBSEQUENT STEPS.----
        sigma <- rnorm(1,sigma_mu,sigma_sd)   ## pull value from combined sigma distribution, save to file for reference at end of simulation runs
        alpha1 <- 1/(2*sigma*sigma)
        siminfo <- ProbStay(sigma, G, s, a)
        
        for(l in 1:length(Ngroup)){
          p0V[l] <- p0MV[[l]][sample(1:length(p0MV[[l]]), 1)]   ## pull value from posterior
          p0T[l] <- p0MT[[l]][sample(1:length(p0MT[[l]]), 1)]   ## pull value from posterior
        }
        
        pmatV <- p0V*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
        pmatT <- p0T*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
        
        for(n in 1:N){
          for(k in 1:K){
            yTrueTRAP[n,,k,z] <- rbinom(J2,1,pmatT[n,newX2])*siminfo$stay[n,k,z]    # observations of each snake at each trap based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
            if(any(yTrueTRAP[n,,k,z] >= 1) == TRUE){
              yTrueVIS[n,,k,z] <- 0
              next
            } else {
              yTrueVIS[n,,k,z] <- rpois(J1,pmatV[n,newX1])*siminfo$stay[n,k,z]  # observations of each snake at each survey location based on encounter probability and effort conditional on if snake was captured in a trap that evening multiplied by Bernoulli trial of if they leave the study area
            }
          }
          for(j in 1:J2){
            ## Collapse observations across all occasions to give a count of observations per trap for each simulation
            yTrueTRAP2[n,j,z] <- sum(yTrueTRAP[n,j,,z])
          }
          for(j in 1:J1){
            yTrueVIS2[n,j,z] <- sum(yTrueVIS[n,j,,z])
          }
        }
        ## Save all values of sigma and p0V/T used for simulations for reference
        paramsim[z,1] <- sigma
        paramsim[z,2:5] <- p0V
        paramsim[z,6:9] <- p0T
      }
      
      # CHECK TO MAKE SURE NO OCCASION HAS A SNAKE TRAPPED AND VISUALLY SEEN
      for(z in 1:nsims){
        for(k in 1:K){
          ifelse(any((rowSums(yTrueTRAP[,,k,z]) & rowSums(yTrueVIS[,,k,z])) == TRUE) == TRUE , print("oh no"), print("we good"))
        }
      }
      
      
      ## yTrueVIS and yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
      capturedV <- list()
      capturedT <- list()
      yarrV <- list()
      yarrT <- list()
      ## Bind sizes of snake to encounter histories
      yTrueVIS3 <- abind::abind(yTrueVIS2,array(Nsnsz, replace(dim(yTrueVIS2),2,1)), along=2)
      yTrueTRAP3 <- abind::abind(yTrueTRAP2,array(Nsnsz, replace(dim(yTrueTRAP2),2,1)), along=2)
      ## Bind unique snake identity to encounter histories
      yTrueVIS3 <- abind::abind(yTrueVIS3,array(seq(1:N), replace(dim(yTrueVIS3),2,1)), along=2)
      yTrueTRAP3 <- abind::abind(yTrueTRAP3,array(seq(1:N), replace(dim(yTrueTRAP3),2,1)), along=2)
      
      for(z in 1:nsims){
        capturedV[[z]] <- which(apply(yTrueVIS3[,1:J1,z],1,sum)>0)  # snakes that were observed at least once
        capturedT[[z]] <- which(apply(yTrueTRAP3[,1:J2,z],1,sum)>0)  # snakes that were observed at least once
        yarrV[[z]] <- yTrueVIS3[capturedV[[z]],,z]  # subset to observed snakes
        yarrT[[z]] <- yTrueTRAP3[capturedT[[z]],,z]  # subset to observed snakes
        write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,stype,"VIS",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
        write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,stype,"TRAP",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
      }
      write.csv(paramsim, file =paste("Simulations/simDat/TrueParameter",type,stype,N,dens,K,stde,z,".csv",sep=""))
    }  ## VISTRAP
  }  ## one-way
}

