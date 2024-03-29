## Simulate observation histories of snakes from different monitoring methods in different sampling designs
## Comment out/in various components to specify design
## Analyze each dataset using King et al. Semicomplete Likelihood (2016) and save results in Results folder

rm(list=ls())

## Required libraries
library(jagsUI);library(secr)
## Functions for simulating data
source("Simulations/Scripts/FunctionsForSimulation_ClosedAndOneWayBarrier.R")



#### SCENARIO DETAILS (USER SPECIFIED, COMMENT IN OR OUT DEPENDING ON DESIGN).----

## Question 1. Is your study area (permeable [one-way movement out of area] or closed [fenced])?
stype <- c("oneway")
# stype <- c("closed")

## Question 2. What type of sampling will you do?
type <- c("VIS")
nmeth <- 1
# type <- c("TRAP")
# nmeth <- 1
# type <- c("VISTRAP")
# nmeth <- 2

## Question 3. How many transects will you survey?
## Full of ONE method [27 transects, all]
stde <- c("full")
samp <- c(1:351)
## Half of ONE method [14 transects, every other]
# stde <- c("half")
# samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
## Half of TWO methods (e.g., VISTRAP)
# stde <- c("halfhalf")
# samp1 <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
# samp2 <- c(14:26,40:52,66:78,92:104,118:130,144:156,170:182,196:208,222:234,248:260,274:286,300:312,326:338)
## Third of TWO methods (e.g., VISTRAP)
# stde <- c("thirdthird")
# samp1 <- c(14:26,53:65,92:104,131:143,170:182,209:221,248:260,287:299,326:338)
# samp2 <- c(27:39,66:78,105:117,144:156,183:195,222:234,261:273,300:312,339:351)

## Question 4. How many nights of sampling will you do? 
## (half [30], quarter [14])
K <- 14

## Question 5. How many snakes are there in the population?
## True number of snakes (normal [120] or low [60] density)
N <- 120

## Question 6. How many snakes are there per size category and are there more small or large snakes?
## 4 groups; <850, >=850 to <950, >=950 to <1150, >=1150

## High density (120 snakes) - many small scenario
dens <- c("small")
Nsnsz <- sample(c(rep(1,times=50),rep(2,times=35),rep(3,times=25),rep(4,times=10)))
## High density (120 snakes) - many large scenario
# dens <- c("large")
# Nsnsz <- sample(c(rep(1,times=10),rep(2,times=25),rep(3,times=35),rep(4,times=50)))
## Low density (60 snakes) - many small scenario
# dens <- c("small")
# Nsnsz <- sample(c(rep(1,times=25),rep(2,times=17),rep(3,times=13),rep(4,times=5)))
## Low density (60 snakes) - many large scenario
# dens <- c("large")
# Nsnsz <- sample(c(rep(1,times=5),rep(2,times=13),rep(3,times=17),rep(4,times=25)))

## Total snakes in each size category
Ngroup <- as.vector(table(Nsnsz))


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area based on rough dimensions of CP
sdeets <- areatype(totlocs = totlocs)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (i.e., if surveying less than the full design), subset locations to only those monitored
## Depending on number of monitoring methods, do the following
if(type == c("VIS") | type == c("TRAP")){
  X <- X[samp,]
  ## Number of sampling points
  J <- nrow(X)
}
if(type == c("VISTRAP")){
  X1 <- X[samp1,]  ## VIS
  X2 <- X[samp2,]  ## TRAP
  X <- X[sort(c(samp1,samp2)),]
  ## Get numeric ID of grid cell for subsetting the integration grid below
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
## Create separate distance matrices if using multiple methods
if(nmeth != 1){
  GdistV <- Gdist[,newX1]
  GdistT <- Gdist[,newX2]
}


#### VISUALIZE SETUP.----
plot(G, pch=16, cex=.5, col="grey")       #plot integration grid
points(X, pch=16, col="red")              #add locations of survey points
## If combined methods, see which points have which method
if(type == c("VISTRAP")){
  points(X1, col="blue", lwd=2)  ## VIS
  points(X2, col="green", lwd=2)  ## TRAP
}


#### CREATE COMBINED SIGMA AND SIMULATE OBSERVATIONS.----

nsims <- 100
## Create and save datasets matching the previously specified scenarios
set.seed(07192021)  # R Version 4.1.0 run on PC
createData() # warning will occur if snake never close enough to cross barrier in oneway simulations - ignore and see FunctionsForSimulation file for more information on how this is handled
## if you get "There were 50 or more warnings (use warnings() to see the first 50)" then try specifying size using code for 120 snakes vs. 60 snakes - see Question 6 above

#### READ IN DATA AND ANALYZE.----

for(i in 1:100){
  
  if(type != c("VISTRAP")){
    ysnsz <- read.csv(paste("Simulations/simDat/",type,stype,N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove individual column
    y <- as.matrix(ysnsz[,-ncol(ysnsz)]) ## observations
    nind <- nrow(y)  ## number of observed individuals
    ## Categories by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
    snsz <- ysnsz[,ncol(ysnsz)]
    L <- length(unique(snsz))
    ngroup <- as.vector(table(snsz))
    if(L == 4){  ## check if all size classes present in simulated data, will need to manually add zeros for those missing
      L
    } else {
      stop("not full length") 
    }
    ## If missing size category, check which and add as a zero
    # L <- 4
    # ngroup <- c(0,ngroup)  ##e.g., if size class 1 missing
  }
  
  if(type == c("VISTRAP")){
    yVsnsz <- read.csv(paste("Simulations/simDat/",type,stype,"VIS",N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove numeric column
    yTsnsz <- read.csv(paste("Simulations/simDat/",type,stype,"TRAP",N,dens,K,stde,i,".csv",sep=""))[,-1]  ## remove numeric column
    yV <- as.matrix(yVsnsz[,-((ncol(yVsnsz)-1):ncol(yVsnsz))]) ## observations
    yT <- as.matrix(yTsnsz[,-((ncol(yTsnsz)-1):ncol(yTsnsz))]) ## observations
    nindV <- nrow(yV)  ## number of individuals visually detected
    nindT <- nrow(yT)  ## number of individuals captured in traps
    nind <- length(unique(c(yVsnsz[,ncol(yVsnsz)],yTsnsz[,ncol(yTsnsz)])))  ## number of individuals detected overall
    ## Categories by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
    snszV <- yVsnsz[,((ncol(yVsnsz)-1):ncol(yVsnsz))];colnames(snszV) <- c("Size","ID")
    snszT <- yTsnsz[,((ncol(yTsnsz)-1):ncol(yTsnsz))];colnames(snszT) <- c("Size","ID")
    ## All snakes ever observed and their sizes
    snsz <- merge(snszV, snszT, all = TRUE)
    snsz <- snsz[order(snsz$ID),][,1]
    L <- length(unique(snsz))
    ngroup <- as.vector(table(snsz))
    if(L == 4){  ## check if all size classes present in simulated data, will need to manually add zeros for those missing
      L
    } else {
      stop("not full length") 
    }
    ## If missing size category, check which and add as a zero
    # L <- 4
    # ngroup <- c(0,ngroup)  ##e.g., if size class 1 missing
  }
  
  ## Initial values for activity centers, take first location where snake found
  if(type != c("VISTRAP")){
    vsst <- list()
    for(w in 1:nrow(y)){
      vsst[w] <- apply(y,1,function(x) which(x>=1))[[w]][1]
      vsst <- unlist(vsst)
    }
  }
  
  if(type == c("VISTRAP")){
    vsstV <- list()
    vsstT <- list()
    for(w in 1:nindV){
      vsstV[w] <- apply(yVsnsz,1,function(x) which(x>=1))[[w]][1]  ## locations from VIS
    }
    for(w in 1:nindT){
      vsstT[w] <- apply(yTsnsz,1,function(x) which(x>=1))[[w]][1]  ## locations from TRAP
    }
    vsstV <- unlist(vsstV)
    vsstT <- unlist(vsstT)
    v <- cbind(yVsnsz[,ncol(yVsnsz)],vsstV)
    t <- cbind(yTsnsz[,ncol(yTsnsz)],vsstT)
    vt <- merge(v,t, all = TRUE)   ## get locations for all snakes across both VIS and TRAP
    vsst <- apply(vt[,2:3], 1, max, na.rm = TRUE)  ## take only one location for each snake (max grid cell)
    indV <- cbind(vt[,1:2],seq(1:nrow(vt)))
    indV <- na.omit(indV)[,ncol(indV)]   ## get snake identities for ones found in VIS
    indT <- cbind(vt[,c(1,3)],seq(1:nrow(vt)))
    indT <- na.omit(indT)[,ncol(indT)]   ## get snake identities for ones found in TRAP
  }
  
  
  ########################################################
  ##Jags model for a King et al 2016 semicomplete likelihood
  
  ### VIS survey model using Poisson observation process - same as TRAP model ###
  
  cat("
  model {
    
    #Prior for spatial decay
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    #Prior for each grid cell (set b[1:Gpts] = rep(1,Gpts) below so as to be a uniform prior across all cells) 
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    for(l in 1:L){   # 4 size categories
      #Prior for intercept
      p0[l] ~ dunif(0,5)
      
      #Prior for N
      #Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures to observe by project and size category
      Ngroup[l] <- ngroup[l] + n0[l]       # number of animals observed + failures to observe by size category
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # number of points on integration grid
        for(j in 1:J){  # number of survey points
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that location
          miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K)
        } #J
        #Probability of detecting each size category across all points in the study area and time period
        pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,])
        pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  # pdot.temp is very close to zero and will lock model up with out this
      } #G
      #Probability of detecting a size category at least once in S (a=area of each integration grid, given as data)
      pstar[l] <- (sum(pdot[l,1:Gpts]*a))/A
      
      #Zero trick for initial 1/pstar^n, described in King et al. 2016
      loglikterm[l] <- -ngroup[l] * log(pstar[l])
      lambda[l] <- -loglikterm[l] + 10000
      dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
    } #L
    
    for(i in 1:n){  # number of observed individuals for each project
      #Define activity centers on a grid
      s[i] ~ dcat(pi[1:Gpts])
      
      #Model for capture histories of observed individuals:
      for(j in 1:J){
        y[i,j] ~ dpois(p[i,j]*K)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
    
    #Derived proportion of each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimVIS.txt")
  
  ### TRAP model using Poisson observation process - same as VIS model  ###
  
  cat("
  model {
    
    #Prior for spatial decay
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    #Prior for each grid cell (set b[1:Gpts] = rep(1,Gpts) below so as to be a uniform prior across all cells)
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    for(l in 1:L){   # 4 size categories
      #Prior for intercept
      p0[l] ~ dunif(0,1)
      
      #Prior for N
      #Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures to observe by project and size category
      Ngroup[l] <- ngroup[l] + n0[l]       # number of animals observed + failures to observe by size category
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # number of points on integration grid
        for(j in 1:J){  # number of trapping points
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that location
          miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K)
        } #J
        #Probability of detecting each size category across all points in the study area and time period
        pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,])
        pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  # pdot.temp is very close to zero and will lock model up with out this
      } #G
      #Probability of detecting a size category at least once in S (a=area of each integration grid, given as data)
      pstar[l] <- (sum(pdot[l,1:Gpts]*a))/A
      
      #Zero trick for initial 1/pstar^n, described in King et al. 2016
      loglikterm[l] <- -ngroup[l] * log(pstar[l])
      lambda[l] <- -loglikterm[l] + 10000
      dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
    } #L
    
    for(i in 1:n){  # number of observed individuals for each project
      #Define activity centers on a grid
      s[i] ~ dcat(pi[1:Gpts])
      
      #Model for capture histories of observed individuals:
      for(j in 1:J){
        y[i,j] ~ dpois(p[i,j]*K)
        p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
    
    #Derived proportion of each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimTRAP.txt")
  
  ### TRAP and VIS COMBO model ###
  
  cat("
  model {
    
    #Prior for spatial decay
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    #Prior for each grid cell (set b[1:Gpts] = rep(1,Gpts) below so as to be a uniform prior across all cells)
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    for(l in 1:L){      # 4 size categories
      #Prior for intercept
      p0V[l] ~ dunif(0,5)  # VIS
      p0T[l] ~ dunif(0,5)  # TRAP
        
      #Prior for N
      #Posterior conditional distribution for N-n (and hence N):
      n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures to observe by project and size category
      Ngroup[l] <- ngroup[l] + n0[l]       # number of animals observed + failures to observe by size category
    }
    
    N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
    
    for(l in 1:L){  # size category
      for(g in 1:Gpts){  # number of points on integration grid
        for(j in 1:J1){  # number of survey points
          #Probability of an individual of size i being missed at grid cell g and point j multiplied by total effort (K) at that location
          miss_allKV[l,g,j] <- pow((1 - p0V[l]*exp(-alpha1*GdistV[g,j]*GdistV[g,j])),K)  # prob missed by visual searches
        } #J1
        for(j in 1:J2){  # number of trapping points
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allKT[l,g,j] <- pow((1 - p0T[l]*exp(-alpha1*GdistT[g,j]*GdistT[g,j])),K)  # prob missed by trapping
        } #J2
        #Probability of detecting each size category across all points in the study area and time period
        pdot.temp[l,g] <- 1 - prod(miss_allKV[l,g,],miss_allKT[l,g,])
        pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  # pdot.temp is very close to zero and will lock model up with out this
      } #G
      #Probability of detecting a size category at least once in S (a=area of each integration grid, given as data)
      pstar[l] <- (sum(pdot[l,1:Gpts]*a))/A
      
      #Zero trick for initial 1/pstar^n, described in King et al. 2016
      loglikterm[l] <- -ngroup[l] * log(pstar[l])
      lambda[l] <- -loglikterm[l] + 10000
      dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
    } #L
    
    for(i in 1:n){  # number of observed individuals for each project
      #Define activity centers on a grid
      s[i] ~ dcat(pi[1:Gpts])
    }
    
    for(i in 1:nV){  # number of individuals observed via visual surveys
      #Model for capture histories of observed individuals from visual surveys:
      for(j in 1:J1){  # number of survey points
        yV[i,j] ~ dpois(pV[i,j]*K)
        pV[i,j] <- p0V[size[indV[i]]]*exp(-alpha1*GdistV[s[indV[i]],j]*GdistV[s[indV[i]],j])
      }#JV
    }#I
    
    for(i in 1:nT){  # number of individuals observed via trapping surveys
      #Model for capture histories of observed individuals from traps:
      for(j in 1:J2){  # number of trapping points
        yT[i,j] ~ dpois(pT[i,j]*K)
        pT[i,j] <- p0T[size[indT[i]]]*exp(-alpha1*GdistT[s[indT[i]],j]*GdistT[s[indT[i]],j])
      }#JT
    }#I
    
    #Derived proportion of each size class
    for(l in 1:L){
      piGroup[l] <- Ngroup[l]/N
    }
  }
  ",file = "Simulations/Models/SCRpstarCATsizeCAT_SimVISTRAP.txt")
  
  #######################################################
  
  ## MCMC settings
  nc <- 5; nAdapt=1000; nb <- 10; ni <- 5000+nb; nt <- 1 
  
  ## When only a single method used
  if(nmeth == 1){
    ## Data and constants
    jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup)
    ## Initial values
    inits <- function(){
      list (sigma=runif(1,50,60), n0=(ngroup+100), s=vsst, p0=runif(L,.001,.002))
    }
    ## Parameters to save
    parameters <- c("p0","sigma","pstar","alpha1","N","n0","Ngroup","piGroup")
  }
  
  ## When two methods used
  if(nmeth == 2){
    ## Data and constants
    jags.data <- list (yV=yV, yT=yT, Gpts=Gpts, GdistV=GdistV, GdistT=GdistT, J1=J1, J2=J2, A=A, K=K, a=a, n=nind, nV=nindV, nT=nindT, indV=indV, indT=indT, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup)
    ## Initial values
    inits <- function(){
      list (sigma=runif(1,50,60), n0=(ngroup+100), s=vsst, p0V=runif(L,.001,.002), p0T=runif(L,.001,.002))
    }
    ## Parameters to save
    parameters <- c("p0V","p0T","sigma","pstar","alpha1","N","n0","Ngroup","piGroup")
  }
  
  ## JAGS model run
  out <- jags(paste("Simulations/Models/SCRpstarCATsizeCAT_Sim",type,".txt",sep=""), data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")
  
  save(out, file=paste("Simulations/Results/RESULTS_",stype,type,N,dens,K,stde,i,".Rdata",sep=""))

}
