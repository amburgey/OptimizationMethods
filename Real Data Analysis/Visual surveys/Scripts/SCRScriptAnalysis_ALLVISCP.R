#### Unified analysis of all closed populaton (CP, aka NWFN) datasets ####
#### Visual survey datasets selected from overall database, formatted individually based on project-specific details, and formatted together for overall analysis
#### Analysis based on semi-complete likelihood spatial capture-recapture framework described by King et al. (2016)

rm(list=ls())

## Adjust default memory limit so model won't fail if it surpasses it
memory.limit(5000000000)

library(secr); library(jagsUI); library(abind)

source("Real Data Analysis/Visual surveys/Scripts/Select&PrepVisualData.R")   ## Creation of subcap and subsurv (cleaned up visual surveys from the overall combined database)
source("Real Data Analysis/Scripts/OverlayCPGrid.R")                          ## Script to define spatial information of CP


#### FORMAT DATA FOR ANALYSIS.----

## List of projects to be analyzed
projects <- c("Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VIS2.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISHL1.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISHL2.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISPREBT2.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISPOSTBT2.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISPOSTKB1.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISPOSTKB2.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISPOSTKB3.R",
        "Real Data Analysis/Visual surveys/Scripts/SCRScriptPrep_VISVISTRAP.R")

nproj <- length(projects)

## Create matrices and arrays for storing data from projects
Kall <- matrix(NA, nrow = 351, ncol = nproj)     ## maximum number of transect points is 351 so the biggest number of rows needed
nindall <- matrix(NA, nrow = nproj, ncol = 1)
yall <- array(c(NA,NA,NA), dim=c(108,351,nproj)) ## maximum number of individuals is 108 so the biggest number of rows needed
snszall <- matrix(NA, nrow = 108, ncol = nproj)
ngroupall <- matrix(NA, nrow = 4, ncol = nproj)  ## maximum number of snake size categories so the biggest number of rows needed
Lall <- 4   # maximum total snake size categories

## Loop across all the projects for this unified analysis, ignore warning about Z-dimension being discarded
for(t in 1:nproj){
  ## Read in prep code for project
  ## The first data prep script includes certain functions and information (e.g., X, A) used for all projects
  source(projects[[t]])
  ## Take prepared data and add to matrices and arrays in order to create combined datasets
  ## Survey info
  Kall[,t] <- K           # number of surveys per transect location  
  ## Capture info (will differ in number of individuals but will always be 351 columns)
  colnames(y) <- seq(1:ncol(y))
  for(r in 1:nrow(y)){
    for(c in 1:ncol(y)){
      yall[r,c,t] <- y[r,c]  # observations
    }
  }
  nindall[t,1] <- nind       # total individuals
  ## Snake size info
  for(r in 1:length(snsz)){
    snszall[r,t] <- snsz[r]  # snake sizes
  }
  for(r in 1:length(ngroup)){
    ngroupall[r,t] <- ngroup[r]  # number snakes per size category
  }
}

## Convert to vectors where needed for JAGS model
nindall <- as.vector(unlist(nindall))
snszall[is.na(snszall)] <- 0   # make NA entries (where no snakes captured) 0 instead 


#### PREPARE SPATIAL DETAILS FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS.----

## Function to find Euclidean distance between all survey points and all integration grid points
e2dist <- function (x, y) {
        i <- sort(rep(1:nrow(y), nrow(x)))
        dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
        matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Establish integration grid
Ggrid <- cellsize                         # spacing of grid, establish in OverlayCPGrid
G <- CPspecs$intgrd[,2:3]                 # location of integration grid points
Gpts <- dim(G)[1]                         # number of integration points
a <- CPspecs$area                         # area of each integration grid
Gdist <- e2dist(G, X)                     # distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##JAGS model for a King et al 2016 semicomplete likelihood

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
      
    for(t in 1:nproj){  # 9 projects
      #Prior for N
      #Posterior conditional distribution for N-n (and hence N):
      n0[l,t] ~ dnegbin(pstar[l,t],ngroup[l,t])  # number of failures to observe by project and size category
      Ngroup[l,t] <- ngroup[l,t] + n0[l,t]       # number of animals observed + failures to observe by size category
    } #T
  } #L

  for(t in 1:nproj){
    N[t] <- sum(Ngroup[1:L,t])  # successful observations plus failures to observe of each size = total N
  }

  for(t in 1:nproj){
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # number of points on integration grid
        for(j in 1:J){  # number of survey points
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that location
          miss_allK[l,g,j,t] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j,t])
        } #J
        #Probability of detecting each size category across all points in the study area and time period
        pdot.temp[l,g,t] <- 1 - prod(miss_allK[l,g,,t]) 
        pdot[l,g,t] <- max(pdot.temp[l,g,t], 1.0E-10)  # pdot.temp is very close to zero and will lock model up with out this
      } #G
      #Probability of detecting a size category at least once in S (a=area of each integration grid, given as data)
      pstar[l,t] <- (sum(pdot[l,1:Gpts,t]*a[1:Gpts]))/A
      
      #Zero trick for initial 1/pstar^n, described in King et al. 2016
      loglikterm[l,t] <- -ngroup[l,t] * log(pstar[l,t])
      lambda[l,t] <- -loglikterm[l,t] + 10000
      dummy[l,t] ~ dpois(lambda[l,t]) # dummy = 0; entered as data
    } #L

    for(i in 1:nFound[t]){  # maximum number of observed individuals for each project
      #Define activity centers on a grid
      s[i,t] ~ dcat(pi[1:Gpts])

      #Model for capture histories of observed individuals:
      for(j in 1:J){
        y[i,j,t] ~ dpois(p[i,j,t]*K[j,t])
        p[i,j,t] <- p0[size[i,t]]*exp(-alpha1*Gdist[s[i,t],j]*Gdist[s[i,t],j])
      }#J
    }#I

    #Derived proportion of each size class
    for(l in 1:L){
      piGroup[l,t] <- Ngroup[l,t]/N[t]
    }#L
  }#T
}

",file = "Real Data Analysis/Visual surveys/Models/SCRpstar_CPALL.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=1000; nb <- 1000; ni <- 10000+nb; nt <- 1  ## extend burn-in manually in simulation based on traceplots (FunctionsForSimulation_ClosedAndOneWayBarrier)

## Data and constants
jags.data <- list (y=yall, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=Kall, nFound=nindall, a=a, dummy=matrix(0,nrow=Lall,ncol=nproj), b=rep(1,Gpts), size=snszall, L=Lall, ngroup=ngroupall, nproj=nproj) 

## Initial values
inits <- function(){
  list (sigma=runif(1,30,40), n0=(ngroupall+10), p0=runif(L,.002,.003)) 
}

## Parameters to save
parameters <- c("p0","sigma","pstar","alpha1","N","n0","Ngroup","piGroup")

## JAGS model run
out <- jags("Real Data Analysis/Visual surveys/Models/SCRpstar_CPALL.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## use "factories" to keep JAGS from locking up with large categorical distribution, may speed things up a little

save(out, file="Real Data Analysis/Visual surveys/Results/NWFNVISALL_SCRpstar.Rdata")
