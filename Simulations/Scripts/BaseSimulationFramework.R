## Simulate capture/observation histories of snakes in different sampling scenarios
## Save data to simDat folder
## Analyze each dataset using King et al. Semicomplete Likelihood (2016) and save results in Results folder

rm(list=ls())

## Required libraries
library(jagsUI);library(secr)
## Functions for simulating data
source("Simulations/Scripts/FunctionsForSimulation.R")


#### SCENARIO DETAILS (USER SPECIFIED).----

## Grid type (open or closed [fenced] area)
stype <- c("open")
# stype <- c("closed")

## Type of sampling design
# type <- c("VIS")
type <- c("TRAP")
# type <- c("MIX")

## Study design (full [351 transects], half [18 transects, every other], or mix)
stde <- c("full")
samp <- c(1:351)
# stde <- c("half")
# samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
# stde <- c("mixed")
## TBD mixed

## Nights of sampling (full [60], half [30], quarter [14])
K <- 60

## True number of snakes (normal [120] or low [60] density)
N <- 120

## Number of snakes per size category (4 groups; <850, >=850 to <950, >=950 to <1150, >=1150)
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
## Area
A <- sdeets$A


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (if surveying less than the full design), subset locations to only those monitored
X <- X[samp,]
## Number of sampling points
J <- nrow(X)


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")      #plot integration grid
points(X, pch=16, col="red")             #add locations of survey points


#### INDIVIDUAL-SPECIFIC INFO.----
## Set seed so we get the same values each time for true information about the population
set.seed(062420)
## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE)


#### CREATE OVERALL POSTERIOR FROM WHICH TO SAMPLE AND SIMULATE OBSERVATIONS.----

nsims <- 1 #1000
## Create and save datasets matching the previously specified scenarios
set.seed(07192021)
createData(type=type,nsims=nsims,Ngroup=Ngroup,Nsnsz=Nsnsz)


#### READ IN DATA AND ANALYZE.----

for(i in 1:nsims){
  ysnsz <- read.csv(paste("Simulations/simDat/",type,N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove individual column
  y <- as.matrix(ysnsz[,-ncol(ysnsz)]) ## observations
  nind <- nrow(y)  ## number of observed individuals
  ## Categories by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
  snsz <- ysnsz[,ncol(ysnsz)]
  L <- length(unique(snsz))
  ngroup <- as.vector(table(snsz))
  
  ## Initial values for activity centers, take first location where snake found
  vsst <- list()
  for(i in 1:nrow(y)){
    vsst[i] <- apply(y,1,function(x) which(x>=1))[[i]][1]
    vsst <- unlist(vsst)
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
      ## For use when defining traps on a grid cell
      s[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of observed individuals:
      for(j in 1:J){  ## J = number of traps
        y[i,j] ~ dbin(p[i,j],K)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
    
    #derived proportion in each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimTRAP.txt")
  
  ### COMBO model ###
  
  #######################################################
  
  # MCMC settings
  nc <- 5; nAdapt=200; nb <- 100; ni <- 25000+nb; nt <- 1 
  
  # Data and constants
  jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup)
  
  # Initial values (same as real data analysis)
  inits <- function(){
    list (sigma=runif(1,50,60), n0=(ngroup+100), s=vsst, p0=runif(L,.001,.002))
  }
  
  parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")
  
  out <- jags(paste("Simulations/Models/SCRpstarCATsizeCAT_Sim",type,".txt",sep=""), data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")
  
  save(out, file=paste("Simulations/Results/RESULTS_",type,N,dens,K,stde,i,".Rdata",sep=""))

}
