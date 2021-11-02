#### FUNCTIONS TO SIMULATE DATA ####

#### FUNCTION TO SPECIFY STUDY AREA TYPE.----

areatype <- function(totlocs, stype, cellsize){
  if(stype == c("closed")){
    ## Define state-space of point process. (i.e., where animals live).
    ## Option 1. Slight buffer in area between transects and fence
    delta <- 11.874929
    Xl <- min(totlocs[,1]) - delta
    Xu <- max(totlocs[,1]) + delta
    Yl <- min(totlocs[,2]) - delta
    Yu <- max(totlocs[,2]) + delta
  }
  
  if(stype == c("open")){
    ## Define state-space of point process. (i.e., where animals live).
    ## Option 2. Delta is the buffer of space around the open study area
    delta <- 29
    Xl <- min(totlocs[,1]) - delta
    Xu <- max(totlocs[,1]) + delta
    Yl <- min(totlocs[,2]) - delta
    Yu <- max(totlocs[,2]) + delta
  }
  
  sinfo <- list(Xl=Xl,Xu=Xu,Yl=Yl,Yu=Yu)
  
  return(sinfo)
}


#### FUNCTION TO FIND DISTANCES BETWEEN INTEGRATION GRID POINTS AND SAMPLING POINTS.----

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


#### FUNCTION TO CREATE POSTERIOR SAMPLE FROM REAL DATA ANALYSIS AND SIMULATE OBSERVATIONS DEPENDING ON THE METHOD SELECTED.----

createData <- function(type, nsims, Ngroup, Nsnsz, stat){
  
  
  ## Read in all different model results and combine into a single posterior for each parameter
  
  if(type == "VIS"){
    
    #### VISUAL SURVEY REAL DATA RESULTS ####
    ## Only showing NWFN so far
    mdlsVIS <- c("Visual surveys/Results/NWFNVIS2_SCRpstarvisCATsizeCATdpois10GRID.Rdata",
                 "Visual surveys/Results/NWFNVISHL1_SCRpstarvisCATsizeCATdpois10GRID.Rdata",
                 "Visual surveys/Results/NWFNVISHL2_SCRpstarvisCATsizeCATdpois.Rdata",
                 "Visual surveys/Results/NWFNVISPREBT2_SCRpstarvisCATsizeCATdpoisLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTBT2_SCRpstarvisCATsizeCATdpoisLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB1_SCRpstarvisCATsizeCATdpois2sizesLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB2_SCRpstarvisCATsizeCATdpois.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB3_SCRpstarvisCATsizeCATdpois3sizesLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISTRAPVIS_SCRpstarvisCATsizeCATdpoisLONGER.Rdata")
    
    ## Parameters (p0, sigma) that influence detection of snakes will be pulled from real data posterior samples
    
    ## Initialize empty vectors for creating posterior samples
    temp <- matrix()
    temp2 <- matrix()
    sigma <- matrix()
    p01 <- matrix()
    p02 <- matrix()
    p03 <- matrix()
    p04 <- matrix()
    
    ## Currently this is contingent on the order of the models above staying the same
    for(m in 1:length(mdlsVIS)){
      load(mdlsVIS[m])
      temp <- out$sims.list$sigma
      sigma <- append(sigma, temp)
      temp2 <- out$sims.list$p0
      ## Some models couldn't estimate all 4 size categories so removed those sizes - need to add posteriors to vectors depending on what was estimated
      if(dim(temp2)[2] == 4){
        p01 <- append(p01,temp2[,1])
        p02 <- append(p02,temp2[,2])
        p03 <- append(p03,temp2[,3])
        p04 <- append(p04,temp2[,4])
      }
      if(dim(temp2)[2] == 3 & m == 8){
        p02 <- append(p02,temp2[,1])
        p03 <- append(p03,temp2[,2])
        p04 <- append(p04,temp2[,3])
      }
      if(dim(temp2)[2] == 2 & m == 6){
        p01 <- append(p01,temp2[,1])
        p02 <- append(p02,temp2[,2])
      }
    }
    
    sigmaMV <- sigma[-1]  # remove initial NA from empty vector
    p0MV <- list(sample(p01[-1]),sample(p02[-1]),sample(p03[-1]),sample(p04[-1]))  # remove initial NA from empty vector and randomize so encounter probabilities aren't from the same model iteration
    
    ## Quick plot to compare encounter probabilities of different snake categories
    c1 <- rgb(0,255,223,max = 255, alpha = 80, names = "lt.blue")
    c2 <- rgb(0,167,255,max = 255, alpha = 80, names = "lt.blue2")
    c3 <- rgb(66,129,255,max = 255, alpha = 80, names = "md.blue")
    c4 <- rgb(7,69,137,max = 255, alpha = 80, names = "dk.blue")
    
    hist(p0MV[[1]], col = c1, breaks = 25, xlim = c(0,0.021), ylim = c(0,30000))
    hist(p0MV[[2]], col = c2, add = TRUE, breaks = 25)
    hist(p0MV[[3]], col = c3, add = TRUE, breaks = 25)
    hist(p0MV[[4]], col = c3, add = TRUE, breaks = 25)
    
    
    #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
    
    ## VISUAL SURVEYS - Generate observations
    yTrueVIS <- array(NA,dim=c(N,J,nsims))
    p0V <- vector()
    
    for(z in 1:nsims){  
      sigmaV <- sigmaMV[sample(0:length(sigmaMV), 1)]   ## pull value from posterior
      alpha1V <- 1/(2*sigmaV*sigmaV)
      for(l in 1:length(Ngroup)){
        p0V[l] <- p0MV[[l]][sample(0:length(p0MV[[l]]), 1)]   ## pull value from posterior
      }
      pmatV <- p0V[Nsnsz]*exp(-alpha1V*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
      for(n in 1:N){
        yTrueVIS[n,,z] <- rpois(J,pmatV[n,]*K)  # observations of each snake at each trap based on encounter probability and effort
      }
    }
    
    ## yTrueVIS includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
    capturedV <- list()
    yarrV <- list()
    ## Bind sizes of snake to encounter histories
    yTrueVIS <- abind::abind(yTrueVIS,array(Nsnsz, replace(dim(yTrueVIS),2,1)), along=2)
    for(i in 1:nsims){
      capturedV[[i]] <- which(apply(yTrueVIS[,1:J,i],1,sum)>0)  # snakes that were observed at least once
      yarrV[[i]] <- yTrueVIS[capturedV[[i]],,i]  # subset to observed snakes
      write.csv(yarrV[[i]], file = paste("Simulations/simDat/",type,N,dens,K,stde,i,".csv",sep=""))  # write to file to keep observations
    }
  }  ## VIS
  
  
  if(type == "TRAP"){
    
    #### TRAPPING REAL DATA RESULTS ####
    ## Only showing NWFN so far
    mdlsTRAP <- c("Trapping/Results/NWFNTRAP1_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNTRAP2LINVIS_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNTRAP3_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNTRAP4LCM_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNVISTRAPTRAP_SCRpstartrapCATsizeCAT3500.Rdata",
                 "Trapping/Results/NWFNPOSTBT2TRAP_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNPOSTKBTRAP1_SCRpstartrapCATsizeCAT.Rdata",
                 "Trapping/Results/NWFNPOSTKBTRAP2_SCRpstartrapCATsizeCAT.Rdata",  ## missing size 2 (no animals caught)
                 "Trapping/Results/NWFNPREBT1TRAP_SCRpstartrapCATsizeCAT.Rdata")
    
    ## Parameters (p0, sigma) that influence detection of snakes will be pulled from real data posterior samples
    
    ## Initialize empty vectors for creating posterior samples
    temp <- matrix()
    temp2 <- matrix()
    sigma <- matrix()
    p01 <- matrix()
    p02 <- matrix()
    p03 <- matrix()
    p04 <- matrix()
    
    ## Currently this is contingent on the order of the models above staying the same
    for(m in 1:length(mdlsTRAP)){
      load(mdlsTRAP[m])
      temp <- out$sims.list$sigma
      sigma <- append(sigma, temp)
      temp2 <- out$sims.list$p0
      ## Some models could estimate all 4 size categories so removed those sizes - need to add posteriors to vectors depending on what was estimated
      if(dim(temp2)[2] == 4){
        p01 <- append(p01,temp2[,1])
        p02 <- append(p02,temp2[,2])
        p03 <- append(p03,temp2[,3])
        p04 <- append(p04,temp2[,4])
      }
      if(dim(temp2)[2] == 3 & m == 8){
        p01 <- append(p01,temp2[,1])
        p03 <- append(p03,temp2[,2])
        p04 <- append(p04,temp2[,3])
      }
    }
    
    sigmaMT <- sigma[-1]  # remove initial NA from empty vector
    p0MT <- list(sample(p01[-1]),sample(p02[-1]),sample(p03[-1]),sample(p04[-1]))  # remove initial NA from empty vector and randomize so encounter probabilities aren't from the same model iteration
    
    ## Quick plot to compare encounter probabilities of different snake categories
    c1 <- rgb(0,255,223,max = 255, alpha = 80, names = "lt.blue")
    c2 <- rgb(0,167,255,max = 255, alpha = 80, names = "lt.blue2")
    c3 <- rgb(66,129,255,max = 255, alpha = 80, names = "md.blue")
    c4 <- rgb(7,69,137,max = 255, alpha = 80, names = "dk.blue")
    
    hist(p0MT[[1]], col = c1, breaks = 25, xlim = c(0,0.021), ylim = c(0,30000))
    hist(p0MT[[2]], col = c2, add = TRUE, breaks = 25)
    hist(p0MT[[3]], col = c3, add = TRUE, breaks = 25)
    hist(p0MT[[4]], col = c3, add = TRUE, breaks = 25)
    
    #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
    
    ## TRAPPING - Generate observations
    yTrueTRAP <- array(NA,dim=c(N,J,nsims))
    p0T <- vector()
    
    for(z in 1:nsims){  
      sigmaT <- sigmaMT[sample(0:length(sigmaMT), 1)]   ## pull value from posterior
      alpha1T <- 1/(2*sigmaT*sigmaT)
      for(l in 1:length(Ngroup)){
        p0T[l] <- p0MT[[l]][sample(0:length(p0MT[[l]]), 1)]   ## pull value from posterior
      }
      pmatT <- p0T[Nsnsz]*exp(-alpha1T*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
      for(n in 1:N){
        yTrueTRAP[n,,z] <- rbinom(J,K,pmatT[n,])  # observations of each snake at each trap based on encounter probability and effort
      }
    }
    
    ## yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
    capturedT <- list()
    yarrT <- list()
    ## Bind sizes of snake to encounter histories
    yTrueTRAP <- abind::abind(yTrueTRAP,array(Nsnsz, replace(dim(yTrueTRAP),2,1)), along=2)
    for(i in 1:nsims){
      capturedT[[i]] <- which(apply(yTrueTRAP[,1:J,i],1,sum)>0)  # snakes that were observed at least once
      yarrT[[i]] <- yTrueTRAP[capturedT[[i]],,i]  # subset to observed snakes
      write.csv(yarrT[[i]], file = paste("Simulations/simDat/",type,N,dens,K,stde,i,".csv",sep=""))  # write to file to keep observations
    }
  }  ## TRAP
  
  
  if(type == "VISTRAP"){
    
    #### VISUAL SURVEY REAL DATA RESULTS ####
    mdlsVIS <- c("Visual surveys/Results/NWFNVIS2_SCRpstarvisCATsizeCATdpois10GRID.Rdata",
                 "Visual surveys/Results/NWFNVISHL1_SCRpstarvisCATsizeCATdpois10GRID.Rdata",
                 "Visual surveys/Results/NWFNVISHL2_SCRpstarvisCATsizeCATdpois.Rdata",
                 "Visual surveys/Results/NWFNVISPREBT2_SCRpstarvisCATsizeCATdpoisLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTBT2_SCRpstarvisCATsizeCATdpoisLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB1_SCRpstarvisCATsizeCATdpois2sizesLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB2_SCRpstarvisCATsizeCATdpois.Rdata",
                 "Visual surveys/Results/NWFNVISPOSTKB3_SCRpstarvisCATsizeCATdpois3sizesLONGER.Rdata",
                 "Visual surveys/Results/NWFNVISTRAPVIS_SCRpstarvisCATsizeCATdpoisLONGER.Rdata")
    
    ## Parameters (p0, sigma) that influence detection of snakes will be pulled from real data posterior samples
    
    ## Initialize empty vectors for creating posterior samples
    tempV <- matrix()
    tempV2 <- matrix()
    sigmaV <- matrix()
    pV01 <- matrix()
    pV02 <- matrix()
    pV03 <- matrix()
    pV04 <- matrix()
    
    ## Currently this is contingent on the order of the models above staying the same
    for(m in 1:length(mdlsVIS)){
      load(mdlsVIS[m])
      tempV <- out$sims.list$sigma
      sigmaV <- append(sigmaV, tempV)
      tempV2 <- out$sims.list$p0
      ## Some models couldn't estimate all 4 size categories so removed those sizes - need to add posteriors to vectors depending on what was estimated
      if(dim(tempV2)[2] == 4){
        pV01 <- append(pV01,tempV2[,1])
        pV02 <- append(pV02,tempV2[,2])
        pV03 <- append(pV03,tempV2[,3])
        pV04 <- append(pV04,tempV2[,4])
      }
      if(dim(tempV2)[2] == 3 & m == 8){
        pV02 <- append(pV02,tempV2[,1])
        pV03 <- append(pV03,tempV2[,2])
        pV04 <- append(pV04,tempV2[,3])
      }
      if(dim(tempV2)[2] == 2 & m == 6){
        pV01 <- append(pV01,tempV2[,1])
        pV02 <- append(pV02,tempV2[,2])
      }
    }
    
    sigmaMV <- sigmaV[-1]  # remove initial NA from empty vector
    p0MV <- list(sample(pV01[-1]),sample(pV02[-1]),sample(pV03[-1]),sample(pV04[-1]))  # remove initial NA from empty vector and randomize so encounter probabilities aren't from the same model iteration
  
    #### TRAPPING REAL DATA RESULTS ####
    ## Only showing NWFN so far
    mdlsTRAP <- c("Trapping/Results/NWFNTRAP1_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNTRAP2LINVIS_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNTRAP3_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNTRAP4LCM_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNVISTRAPTRAP_SCRpstartrapCATsizeCAT3500.Rdata",
                  "Trapping/Results/NWFNPOSTBT2TRAP_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNPOSTKBTRAP1_SCRpstartrapCATsizeCAT.Rdata",
                  "Trapping/Results/NWFNPOSTKBTRAP2_SCRpstartrapCATsizeCAT.Rdata",  ## missing size 2 (no animals caught)
                  "Trapping/Results/NWFNPREBT1TRAP_SCRpstartrapCATsizeCAT.Rdata")
    
    ## Parameters (p0, sigma) that influence detection of snakes will be pulled from real data posterior samples
    
    ## Initialize empty vectors for creating posterior samples
    tempT <- matrix()
    tempT2 <- matrix()
    sigmaT <- matrix()
    pT01 <- matrix()
    pT02 <- matrix()
    pT03 <- matrix()
    pT04 <- matrix()
    
    ## Currently this is contingent on the order of the models above staying the same
    for(m in 1:length(mdlsTRAP)){
      load(mdlsTRAP[m])
      tempT <- out$sims.list$sigma
      sigmaT <- append(sigmaT, tempT)
      tempT2 <- out$sims.list$p0
      ## Some models could estimate all 4 size categories so removed those sizes - need to add posteriors to vectors depending on what was estimated
      if(dim(tempT2)[2] == 4){
        pT01 <- append(pT01,tempT2[,1])
        pT02 <- append(pT02,tempT2[,2])
        pT03 <- append(pT03,tempT2[,3])
        pT04 <- append(pT04,tempT2[,4])
      }
      if(dim(tempT2)[2] == 3 & m == 8){
        pT01 <- append(pT01,tempT2[,1])
        pT03 <- append(pT03,tempT2[,2])
        pT04 <- append(pT04,tempT2[,3])
      }
    }
    
    sigmaMT <- sigmaT[-1]  # remove initial NA from empty vector
    p0MT <- list(sample(pT01[-1]),sample(pT02[-1]),sample(pT03[-1]),sample(pT04[-1]))  # remove initial NA from empty vector and randomize so encounter probabilities aren't from the same model iteration
    
    ## Quick plot to compare encounter probabilities of different snake categories
    ## Quick plot to compare encounter probabilities of different snake categories
    cV1 <- rgb(0,255,223,max = 255, alpha = 80, names = "lt.blue")
    cV2 <- rgb(0,167,255,max = 255, alpha = 80, names = "lt.blue2")
    cV3 <- rgb(66,129,255,max = 255, alpha = 80, names = "md.blue")
    cV4 <- rgb(7,69,137,max = 255, alpha = 80, names = "dk.blue")
    cT1 <- rgb(238,170,189,max = 255, alpha = 80, names = "lt.pink")
    cT2 <- rgb(239,146,173,max = 255, alpha = 80, names = "lt.pink2")
    cT3 <- rgb(244,106,146,max = 255, alpha = 80, names = "md.pink")
    cT4 <- rgb(241,61,113,max = 255, alpha = 80, names = "dk.pink")
    
    hist(p0MV[[1]], col = cV1, breaks = 25, xlim = c(0,0.021), ylim = c(0,30000))
    hist(p0MV[[2]], col = cV2, add = TRUE, breaks = 25)
    hist(p0MV[[3]], col = cV3, add = TRUE, breaks = 25)
    hist(p0MV[[4]], col = cV3, add = TRUE, breaks = 25)
    hist(p0MT[[1]], col = cT1, add = TRUE, breaks = 25)
    hist(p0MT[[2]], col = cT2, add = TRUE, breaks = 25)
    hist(p0MT[[3]], col = cT3, add = TRUE, breaks = 25)
    hist(p0MT[[4]], col = cT3, add = TRUE, breaks = 25)
    
    #### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
    
    ## VISUAL SURVEYS - Generate observations
    yTrueVIS <- array(NA,dim=c(N,J1,nsims), dimnames = list(1:N,samp1,1:nsims))
    yTrueTRAP <- array(NA,dim=c(N,J2,nsims), dimnames = list(1:N,samp2,1:nsims))
    yTrueVT <- array(NA,dim=c(N,J,nsims))
    p0V <- vector()
    p0T <- vector()
    
    for(z in 1:nsims){  
      sigmaV <- sigmaMV[sample(0:length(sigmaMV), 1)]   ## pull value from posterior
      sigmaT <- sigmaMT[sample(0:length(sigmaMT), 1)]   ## pull value from posterior
      alpha1V <- 1/(2*sigmaV*sigmaV)
      alpha1T <- 1/(2*sigmaT*sigmaT)
      
      for(l in 1:length(Ngroup)){
        p0V[l] <- p0MV[[l]][sample(0:length(p0MV[[l]]), 1)]   ## pull value from posterior
        p0T[l] <- p0MT[[l]][sample(0:length(p0MT[[l]]), 1)]   ## pull value from posterior
      }
      
      pmatV <- p0V*exp(-alpha1V*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
      pmatT <- p0T*exp(-alpha1T*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
      
      locID <- sort(c(samp1,samp2))
      '%!in%' <- function(x,y)!('%in%'(x,y))
      
      for(n in 1:N){
        for(j in 1:J){
          
          if(locID[j] %in% samp1){
            yTrueVT[n,j,z] <- rpois(1,pmatV[n,j]*K)
          }
          if(locID[j] %!in% samp2){
            yTrueVT[n,j,z] <- rbinom(1,K,pmatT[n,j])
          }
          
          # yTrueVIS[n,,z] <- rpois(J1,pmatV[n,]*K)  # observations of each snake at each trap based on encounter probability and effort
          # yTrueTRAP[n,,z] <- rbinom(J2,K,pmatT[n,])  # observations of each snake at each trap based on encounter probability and effort
          # yTrueAll <- abind::abind(yTrueVIS, yTrueTRAP, along = 2)  # combine two data types together into one dataframe
          # yTrueAll <- yTrueAll[order(yTrueAll[,as.character(c(14:39,53:78,92:117,131:156,170:195,209:234,248:273,287:312,326:351)),]),,] # order by transect location
          # yTrueAll <- sort(yTrueAll, as.character(c(14:39,53:78,92:117,131:156,170:195,209:234,248:273,287:312,326:351)))
        }
      }
    }
    
    ## yTrueVIS includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
    capturedV <- list()
    capturedT <- list()
    yarrV <- list()
    yarrT <- list()
    ## Bind sizes of snake to encounter histories
    yTrueVIS <- abind::abind(yTrueVIS,array(Nsnsz, replace(dim(yTrueVIS),2,1)), along=2)
    yTrueTRAP <- abind::abind(yTrueTRAP,array(Nsnsz, replace(dim(yTrueTRAP),2,1)), along=2)
    
    for(i in 1:nsims){
      capturedV[[i]] <- which(apply(yTrueVIS[,1:J1,i],1,sum)>0)  # snakes that were observed at least once
      yarrV[[i]] <- yTrueVIS[capturedV[[i]],,i]  # subset to observed snakes
      write.csv(yarrV[[i]], file = paste("Simulations/simDat/",type,"VIS",N,dens,K,stde,i,".csv",sep=""))  # write to file to keep observations
      capturedT[[i]] <- which(apply(yTrueTRAP[,1:J2,i],1,sum)>0)  # snakes that were observed at least once
      yarrT[[i]] <- yTrueTRAP[capturedT[[i]],,i]  # subset to observed snakes
      write.csv(yarrT[[i]], file = paste("Simulations/simDat/",type,"TRAP",N,dens,K,stde,i,".csv",sep=""))  # write to file to keep observations
    }   ## VISTRAP
  }
}