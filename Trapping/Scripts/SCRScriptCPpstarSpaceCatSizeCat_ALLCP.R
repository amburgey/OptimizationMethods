#### Unified analysis of all CP (NWFN) datasets ####
#### TRAPPING DATA ANALYSIS ####

rm(list=ls())

library(secr); library(reshape2); library(jagsUI); library(abind)

source("Select&PrepTrapData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Trapping/DataPrep/OverlayCPGrid.R")

projects <- c("Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_TRAP1.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_TRAP2LINVIS.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_TRAP3.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_TRAP4LCM.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISTRAPTRAP.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_POSTBT2TRAP.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_POSTKBTRAP1.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_POSTKBTRAP2.R",
        "Trapping/Scripts/SCRScriptCPpstarSpaceCatSizeCat_PREBT1TRAP.R")

nproj <- length(projects)

noccall <- matrix(NA, nrow = nproj, ncol = 1)
Kall <- matrix(NA, nrow = 169, ncol = nproj)
nindall <- matrix(NA, nrow = nproj, ncol = 1)
yall <- array(c(NA,NA,NA), dim=c(98,169,nproj)) ## maximize number of individuals is 82 so the biggest number of rows ever needed
snszall <- matrix(NA, nrow = 98, ncol = nproj)
# Lall <- as.data.frame(matrix(NA, nrow = nproj, ncol = 1))
ngroupall <- matrix(NA, nrow = 4, ncol = nproj)
vsstall <- matrix(NA, nrow = 98, ncol = nproj)

for(t in 1:nproj){
  ## Read in prep code for project
  ## The first data prep script includes loading functions and subsetting to just NWFN/CP
  source(projects[[t]])
  ## Take prepared data and add in order to create combined datasets
  ## Survey info
  noccall[t,1] <- nocc
  Kall[,t] <- K
  ## Capture info (will differ in number of individuals but will always be 351 columns)
  colnames(y) <- seq(1:ncol(y))
  for(r in 1:nrow(y)){
    for(c in 1:ncol(y)){
      yall[r,c,t] <- y[r,c]
    }
  }
  # yall <- simplify2array(lapply(list(test,test2), as.matrix))
  nindall[t,1] <- nind
  ## Snake size info
  for(r in 1:length(snsz)){
    snszall[r,t] <- snsz[r]
  }
  Lall <- L
  for(r in 1:length(ngroup)){
    ngroupall[r,t] <- ngroup[r]
  }
  ## Inits
  for(r in 1:length(vsst)){
    vsstall[r,t] <- vsst[r]
  }
}

## Convert to vectors
noccall <- as.vector(unlist(noccall)) ## don't need for model but use for error checking
nindall <- as.vector(unlist(nindall))
vsstall[is.na(vsstall)] <- 0
snszall[is.na(snszall)] <- 0
ngroupall[is.na(ngroupall)] <- 0


#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

## Function to find Euclidean distance between all survey points and all integration grid points
e2dist <- function (x, y) {
        i <- sort(rep(1:nrow(y), nrow(x)))
        dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
        matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Integration grid
Ggrid <- cellsize                         #spacing
G <- CPspecs$intgrd[,2:3]
Gpts <- dim(G)[1]                         #number of integration points
a <- CPspecs$area                         #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##Jags model for a King et al 2016 semicomplete likelihood

cat("
model {

  #prior for spatial decay
  sigma ~ dunif(0,100)
  alpha1 <- 1/(2*sigma*sigma)
  # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
  
  for(l in 1:L){   # 4 size categories
    #prior for intercept
    p0[l] ~ dunif(0,5)
    alpha0[l] <- logit(p0[l])
      
    for(t in 1:nproj){
      # Prior for N
      # Posterior conditional distribution for N-n (and hence N):
      n0[l,t] ~ dnegbin(pstar[l,t],ngroup[l,t])  # number of failures by project and size category
      Ngroup[l,t] <- ngroup[l,t] + n0[l,t]
    } #T
  } #L

  for(t in 1:nproj){
    N[t] <- sum(Ngroup[1:L,t])  # successful observations plus failures to observe of each size = total N
  }
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)

  for(t in 1:nproj){
    for(l in 1:L){  # size category
      for(g in 1:Gpts){ # Gpts = number of points on integration grid
        for(j in 1:J){  # J = number of traps
          #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
          miss_allK[l,g,j,t] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j,t])
        } #J
        pdot.temp[l,g,t] <- 1 - prod(miss_allK[l,g,,t]) #Prob of detect each size category across entire study area and time period
        pdot[l,g,t] <- max(pdot.temp[l,g,t], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
      } #G
      pstar[l,t] <- (sum(pdot[l,1:Gpts,t]*a[1:Gpts]))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)

      # Zero trick for initial 1/pstar^n
      loglikterm[l,t] <- -ngroup[l,t] * log(pstar[l,t])
      lambda[l,t] <- -loglikterm[l,t] + 10000
      dummy[l,t] ~ dpois(lambda[l,t]) # dummy = 0; entered as data
    } #L

    for(i in 1:nFound[t]){  ## nmax = maximum number of observed individuals across all projects
      ## For use when defining traps on a grid cell
      s[i,t] ~ dcat(pi[1:Gpts])

      # Model for capture histories of observed individuals:
      for(j in 1:J){  ## J = number of traps
        y[i,j,t] ~ dpois(p[i,j,t]*K[j,t])
        p[i,j,t] <- p0[size[i,t]]*exp(-alpha1*Gdist[s[i,t],j]*Gdist[s[i,t],j])
      }#J
    }#I

    #derived proportion in each size class
    for(l in 1:L){
      piGroup[l,t] <- Ngroup[l,t]/N[t]
    }#L
  }#T
}

",file = "Trapping/Models/SCRpstarCATsizeCAT_CPALL.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 10+nb; nt <- 1

## Data and constants
jags.data <- list (y=yall, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=Kall, nFound=nindall, a=a, dummy=matrix(0,nrow=Lall,ncol=nproj), b=rep(1,Gpts), size=snszall, L=Lall, ngroup=ngroupall, nproj=nproj) # ## semicomplete likelihood, #nFound=nFound, nmax=max(nindall), 

inits <- function(){
  list (sigma=runif(1,30,40), n0=(ngroupall+10), p0=runif(L,.002,.003)) #s=vsstall, 
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")

out <- jags("Trapping/Models/SCRpstarCATsizeCAT_CPALL.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little

save(out, file="Trapping/Results/NWFNVISALL_SCRpstarvisCATsizeCATdpois10GRIDnovsstALL.Rdata")

