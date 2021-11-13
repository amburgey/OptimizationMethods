## Simulate capture/observation histories of snakes in different sampling scenarios
## Save data to simDat folder
## Analyze each dataset using King et al. Semicomplete Likelihood (2016) and save results in Results folder

rm(list=ls())

## Required libraries
library(jagsUI);library(secr)
## Functions for simulating data
source("Simulations/Scripts/FunctionsForSimulation.R")



#### SCENARIO DETAILS (USER SPECIFIED).----

## Question 1. Is your study area (open or closed [fenced])?
# stype <- c("open")
stype <- c("closed")

## Question 2. What type of sampling will you do?
# type <- c("VIS")
# nmeth <- 1
# type <- c("TRAP")
# nmeth <- 1
type <- c("VISTRAP")
nmeth <- 2
# type <- c("MIX")

## Question 3. How many transects will you survey?
## Full [351 transects]
# stde <- c("full")
# samp <- c(1:351)
## Half [14 transects, every other]
# stde <- c("half")
# samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
## Third [9 transects, every third]
# stde <- c("third")
# samp <- c(27:39,66:78,105:117,144:156,183:195,222:234,261:273,300:312,339:351)
## Half of two methods (e.g., VISTRAP)
# stde <- c("halfhalf")
# samp1 <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
# samp2 <- c(14:26,40:52,66:78,92:104,118:130,144:156,170:182,196:208,222:234,248:260,274:286,300:312,326:338)
## Third of two methods (e.g., VISTRAP)
stde <- c("thirdthird")
samp1 <- c(14:26,53:65,92:104,131:143,170:182,209:221,248:260,287:299,326:338)
samp2 <- c(27:39,66:78,105:117,144:156,183:195,222:234,261:273,300:312,339:351)
## CAMERA STUFF TBD

## Question 4. How many nights of sampling will you do? 
## (full [60], half [30], quarter [14])
K <- 60

## Question 5. How many snakes are there in the population?
## True number of snakes (normal [120] or low [60] density)
N <- 120

## Question 6. How many snakes are there per size category and are there more small or large snakes?
## 4 groups; <850, >=850 to <950, >=950 to <1150, >=1150
## High density - many small scenario
dens <- c("small")
Nsnsz <- sample(c(rep(1,times=50),rep(2,times=35),rep(3,times=25),rep(4,times=10)))
## High density - many large scenario
# dens <- c("large")
# Nsnsz <- sample(c(rep(1,times=10),rep(2,times=25),rep(3,times=35),rep(4,times=50)))
## Low density - many small scenario
# dens <- c("small")
# Nsnsz <- sample(c(rep(1,times=25),rep(2,times=17),rep(3,times=13),rep(4,times=5)))
## Low density - many large scenario
# dens <- c("large")
# Nsnsz <- sample(c(rep(1,times=5),rep(2,times=13),rep(3,times=17),rep(4,times=25)))
Ngroup <- as.vector(table(Nsnsz))


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area depending on if it's either a closed or open area based on "stype"
sdeets <- areatype(totlocs = totlocs, stype = stype)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (if surveying less than the full design), subset locations to only those monitored
if(type == c("VIS") | type == c("TRAP")){
  X <- X[samp,]
  ## Number of sampling points
  J <- nrow(X)
}
if(type == c("VISTRAP")){
  X1 <- X[samp1,]  ## VIS
  X2 <- X[samp2,]  ## TRAP
  X <- X[sort(c(samp1,samp2)),]
  ## Get numeric ID for grid cell subsetting for integration grid below
  newX <- cbind(X,seq(1:nrow(X)))
  newX1 <- newX[rownames(X1),][,3]
  newX2 <- newX[rownames(X2),][,3]
  ## Number of sampling points
  J1 <- nrow(X1)  ## VIS
  J2 <- nrow(X2)  ## TRAP
  J <- sum(J1 + J2)
}


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points using the general constraints of Xl, Xu, Yl, Yu
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
A <- Gpts * a                             #area of study area
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
## Create status by point matrix to indicate what kind of method implemented at which point if using multiple
if(nmeth != 1){
  GdistV <- Gdist[,newX1]
  GdistT <- Gdist[,newX2]
  x1 <- as.data.frame(samp1); colnames(x1) <- c("Loc")
  x1$Type <- 1
  x2 <- as.data.frame(samp2); colnames(x2) <- c("Loc")
  x2$Type <- 2
  stat <- rbind(x1,x2)
  stat <- stat[order(stat$Loc),][,2]
  statV <- stat
  statV[statV == 2] <- 0
  statT <- stat
  statT[statT == 1] <- 0
  statT[statT == 2] <- 1
}


#### VISUALIZE SETUP.----
plot(G, pch=16, cex=.5, col="grey")       #plot integration grid
points(X, pch=16, col="red")              #add locations of survey points
## If combined methods, see which points have which method
if(type == c("VISTRAP")){
  points(X1, col="blue", lwd=2)  ## VIS
  points(X2, col="green", lwd=2)  ## TRAP
}


#### INDIVIDUAL-SPECIFIC INFO.----
## Set seed so we get the same values each time for true information about the population
set.seed(062420)
## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE)


#### CREATE OVERALL POSTERIOR FROM WHICH TO SAMPLE AND SIMULATE OBSERVATIONS.----

nsims <- 1 #1000
## Create and save datasets matching the previously specified scenarios
set.seed(07192021)
createData(type=type,nsims=nsims,Ngroup=Ngroup,Nsnsz=Nsnsz,stat=stat)


#### READ IN DATA AND ANALYZE.----

for(i in 1:nsims){
  
  if(type != c("VISTRAP")){
    ysnsz <- read.csv(paste("Simulations/simDat/",type,N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove individual column
    y <- as.matrix(ysnsz[,-ncol(ysnsz)]) ## observations
    nind <- nrow(y)  ## number of observed individuals
    ## Categories by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
    snsz <- ysnsz[,ncol(ysnsz)]
    L <- length(unique(snsz))
    ngroup <- as.vector(table(snsz))
  }
  
  if(type == c("VISTRAP")){
    yVsnsz <- read.csv(paste("Simulations/simDat/",type,"VIS",N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove numeric column
    yTsnsz <- read.csv(paste("Simulations/simDat/",type,"TRAP",N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove numeric column
    yV <- as.matrix(yVsnsz[,-((ncol(yVsnsz)-1):ncol(yVsnsz))]) ## observations
    yT <- as.matrix(yTsnsz[,-((ncol(yTsnsz)-1):ncol(yTsnsz))]) ## observations
    nindV <- nrow(yV)  ## number of individuals visually detected
    nindT <- nrow(yT)  ## number of individuals captured in traps
    ## Categories by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
    snszV <- yVsnsz[,((ncol(yVsnsz)-1):ncol(yVsnsz))]
    snszT <- yTsnsz[,((ncol(yTsnsz)-1):ncol(yTsnsz))]
    ## All snakes ever observed and their sizes
    snsz <- merge(snszV, snszT, all = TRUE)
    snsz <- snsz[order(snsz$V119),][,1]
    snszV <- yVsnsz[,(ncol(yVsnsz)-1)]  ## snake size of visually detected
    snszT <- yTsnsz[,(ncol(yTsnsz)-1)]  ## snake size of trapped individuals
    L <- max(c(unique(snszV),unique(snszT)))
    ngroup <- as.vector(table(snsz))
    ngroupV <- as.vector(table(snszV))
    ngroupT <- as.vector(table(snszT))
  }
  
  ## Initial values for activity centers, take first location where snake found
  if(type != c("VISTRAP")){
    vsst <- list()
    for(i in 1:nrow(y)){
      vsst[i] <- apply(y,1,function(x) which(x>=1))[[i]][1]
      vsst <- unlist(vsst)
    }
  }
  
  if(type == c("VISTRAP")){
    vsstV <- list()
    vsstT <- list()
    for(i in 1:nrow(yV)){
      vsstV[i] <- apply(yV,1,function(x) which(x>=1))[[i]][1]
      vsstV <- unlist(vsstV)
    }
    for(i in 1:nrow(yT)){
      vsstT[i] <- apply(yT,1,function(x) which(x>=1))[[i]][1]
      vsstT <- unlist(vsstT)
    }
  }
  
  
  ########################################################
  ##Jags model for a King et al 2016 semicomplete likelihood
  
  ### VIS survey model using Poisson observation process ###
  
  cat("
  model {
    
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    
    for(l in 1:L){   # 4 size categories
      #prior for intercept
      p0[l] ~ dunif(0,5)
      alpha0[l] <- logit(p0[l])
      
      # Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
      Ngroup[l] <- ngroup[l] + n0[l]
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # Gpts = number of points on integration grid
        for(j in 1:J){  # J = number of traps
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K)
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
      ## For use when defining traps on a grid cell
      s[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of observed individuals:
      for(j in 1:J){  ## J = number of traps
        y[i,j] ~ dpois(p[i,j]*K)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
    
    #derived proportion in each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimVIS.txt")
  
  ### TRAP model using Binomial observation process ###
  
  cat("
  model {
    
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    
    for(l in 1:L){   # 4 size categories
      #prior for intercept
      p0[l] ~ dunif(0,1)
      alpha0[l] <- logit(p0[l])
      
      # Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
      Ngroup[l] <- ngroup[l] + n0[l]
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # Gpts = number of points on integration grid
        for(j in 1:J){  # J = number of traps
          # Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K)
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
      # For use when defining traps on a grid cell
      s[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of observed individuals:
      for(j in 1:J){  ## J = number of traps
        y[i,j] ~ dbin(p[i,j],K)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
    
    # derived proportion in each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimTRAP.txt")
  
  ### TRAP and VIS COMBO model ###
  
  cat("
  model {
    
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    
    for(l in 1:L){      # 4 size categories
      #prior for intercept
      p0V[l] ~ dunif(0,5)  # VIS
      p0T[l] ~ dunif(0,5)  # TRAP
        
      # Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
      Ngroup[l] <- ngroup[l] + n0[l]
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # Gpts = number of points on integration grid
        for(j in 1:J){  # J = number of traps (currently equal number per method)
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allKV[l,g,j] <- pow((1 - p0V[l]*exp(-alpha1*GdistV[g,j]*GdistV[g,j])),K)  # prob missed by visual searches
          miss_allKT[l,g,j] <- pow((1 - p0T[l]*exp(-alpha1*GdistT[g,j]*GdistT[g,j])),K)  # prob missed by trapping
        } #J
        pdot.temp[l,g] <- 1 - prod(miss_allKV[l,g,]*miss_allKT[l,g,]) #Prob of detect each size category across entire study area and time period
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
    
    for(i in 1:nV){  ## nV = number of individuals observed via visual surveys
      ## For use when defining traps on a grid cell
      sV[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of observed individuals from visual surveys:
      for(j in 1:J){  ## J = number of visual surveys
        yV[i,j] ~ dpois(pV[i,j]*K)
        pV[i,j] <- p0V[sizeV[i]]*exp(-alpha1*GdistV[sV[i],j]*GdistV[sV[i],j])
      }#JV
    }
    
    for(i in 1:nT){  ## nV = number of individuals observed via visual surveys
      ## For use when defining traps on a grid cell
      sT[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of individuals from traps:
      for(j in 1:J){  ## J = number of traps
        yT[i,j] ~ dpois(pT[i,j]*K)
        pT[i,j] <- p0T[sizeT[i]]*exp(-alpha1*GdistT[sT[i],j]*GdistT[sT[i],j])
      }#JT
    }#I
    
    #derived proportion in each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimVISTRAP.txt")
  
  #######################################################
  
  # MCMC settings
  nc <- 5; nAdapt=1000; nb <- 1; ni <- 10000+nb; nt <- 1 
  
  ## When only a single method used
  # Data and constants
  # jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup)
  
  ## When two methods used
  # Data and constants
  jags.data <- list (yV=yV, yT=yT, Gpts=Gpts, GdistV=GdistV, GdistT=GdistT, J=J1, A=A, K=K, a=a, nV=nindV, nT=nindT, dummy=rep(0,L), b=rep(1,Gpts), sizeV=snszV, sizeT=snszT, L=L, ngroup=ngroup)#, statT=statT, statV=statV)
  
  ## When only a single method used
  # Initial values (same as real data analysis)
  # inits <- function(){
  #   list (sigma=runif(1,50,60), n0=(ngroup+100), s=vsst, p0=runif(L,.001,.002))
  # }
  
  ## When two methods used
  # Initial values (same as real data analysis)
  inits <- function(){
    list (sigma=runif(1,50,60), n0=(ngroup+100), sV=vsstV, sT=vsstT, p0V=runif(L,.001,.002), p0T=runif(L,.001,.002))
  }
  
  ## When only a single method used
  # parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")
  
  ## When two methods used
  parameters <- c("p0V","p0T","sigma","pstar","alpha1","N","n0","Ngroup","piGroup")
  
  out <- jags(paste("Simulations/Models/SCRpstarCATsizeCAT_Sim",type,".txt",sep=""), data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")
  
  save(out, file=paste("Simulations/Results/RESULTS_",type,N,dens,K,stde,i,".Rdata",sep=""))

}
