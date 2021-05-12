##### HMU (Habitat Management Unit) is a 55-ha semi-closed (fenced to entry but not exit of snakes) study area

rm(list=ls())


source("Select&PrepVisualData.R")  ## Creation of subcap and subsurv
source("Visual surveys/DataPrep/OverlayHMUGrid.R")
source("Visual surveys/DataPrep/DataPrepHMU.R")

library(secr); library(reshape2); library(jagsUI)

## Capture data (subcap) and effort/survey data (subsurv)
HMUcaps <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR"))
HMUsurv <- subset(subsurv, (SITE == "HMUI" | SITE == "HMUR"))

## Subset to specific NWFN project (options = MWFM VIS 2, NWFN VIS HL 1, NWFN VIS HL 2, PRE NT2 VIS, POST BT2 VIS, POST KB VIS 1, POST KB VIS 2, POST KB VIS 3 EXTRA, POST KB VIS 3, NWFN VISPACE, NWFN SCENT VIS TRAIL)
HMUcaps <- subset(HMUcaps, PROJECTCODE == "EDGE EFFECT VIS")
HMUsurv <- subset(HMUsurv, PROJECTCODE == "EDGE EFFECT VIS")

##### SPECIFY DIMENSIONS AND GRID OF HMU #####
cellsize <- 10  ## dimensions of integration grid cell
HMUspecs <- overlayHMU(HMUcaps, cellsize)  ## ignore warnings, all about projections
## Area (55 ha/550,000 m2): 
A <- 550000

##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
X <- as.matrix(HMUspecs$tran[,-1])[,2:3]
J <- nrow(X)

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected or the size of snakes involved
capPROJ <- subSnk(SITEcaps=HMUspecs$snks)
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=c("05","06"))  ## specify month range
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=HMUsurv, time=c("05","06"))
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

#### FORMAT DATA FOR TRADITIONAL SCR ANALYSIS ####
# dat <- prepSCR(SCRcaps, SCReff, grid = HMUspecs$tran)  ## if error here do one below
dat <- prepSCRman(SCRcaps, SCReff, grid = HMUspecs$tran) ## MANUALLY CHANGE DISTANCES

## Observations, already in order of 1-351 CellID locations
y <- dat$y
colnames(y) <- NULL

## Uniquely marked individuals
nind <- nrow(y)

## Get sizes of individuals
snsz <- getSize(capPROJ, SCRcaps, subcap)[,2]  ## if all snakes have a measurement during that project
# snsz <- getSizeman(capPROJ, SCRcaps, subcap, time=c("2013-04-01","2013-07-31"))[,2] ## if some snake sizes are missing than expand window of time
## Categorize by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
snsz <- ifelse(snsz < 850, 1,
               ifelse(snsz >= 850 & snsz < 950, 2,
                      ifelse(snsz >= 950 & snsz < 1150, 3,
                             ifelse(snsz >= 1150, 4, -9999))))
if(max(snsz) == -9999) stop('snake size incorrect')
L <- length(unique(snsz))
ngroup <- as.vector(table(snsz))

## Active/not active for when transects run, already in order of GridID locations
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)

## Inits for activity centers, take mean grid cell location where each snake was found
locs <- HMUspecs$tran
colnames(locs)[2] <- c("CellID")
sst <- round(unlist(lapply(apply(dat$y,1,function(x) which(x==1)),function(x) max(as.numeric(names(x))))))
vlocs <- locs[locs$CellID %in% sst,]
sst <- as.data.frame(sst); colnames(sst) <- c("CellID")
vsst <- unlist(merge(sst,vlocs, by=c("CellID"))[,1])

#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Integration grid
Ggrid <- cellsize                                #spacing (check sensitivity to spacing)
G <- HMUspecs$intgrd[,2:3]
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##Jags model for a King et al 2016 semicomplete likelihood

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

  for(l in 1:4){  # size category
    for(g in 1:Gpts){ # Gpts = number of points on integration grid
      for(j in 1:J){  # J = number of traps
        #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
        one_minus_detprob[l,g,j] <- 1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])*K[j] #Gdist given as data
      } #J
      pdot.temp[l,g] <- 1 - prod(one_minus_detprob[l,g,]) #Prob of failure to detect each size category across entire study area and time period
      pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
    } #G
    pstar[l] <- (sum(pdot[l,1:Gpts])*a)/A   #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
  
    # Zero trick for initial 1/pstar^n
    loglikterm[l] <- -ngroup[l] * log(pstar[l])
    lambda[l] <- -loglikterm[l] + 1000
    dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
  } #L

  # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
 
  for(i in 1:n){  ## n = number of observed individuals
    ## For use when defining traps on a grid cell
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## J = number of traps
      y[i,j] ~ dbin(p[i,j],K[j])
      p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
    }#J
  }#I
  
  #derived proportion in each size class
  for(l in 1:L){
    piGroup[l] <- Ngroup[l]/N
  }
}
",file = "Visual surveys/Models/SCRpstarCATsizeCAT_HMU.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 10000+nb; nt <- 1
# nc <- 3; nAdapt=20; nb <- 10; ni <- 100+nb; nt <- 1

## Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=rep(0,4), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup) # ## semicomplete likelihood
#locs=X, 

inits <- function(){
  list (sigma=runif(1,45,50), n0=ngroup, s=vsst, p0=runif(L,.002,.003))
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")

out <- jags("Visual surveys/Models/SCRpstarCATsizeCAT_HMU.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little

save(out, file="Visual surveys/Results/HMUEDGE_SCRpstarvisCATsizeCAT.Rdata")  ## M = 150 (XXXXhrs)

### Ran for two days, convergence looks pretty good, cell = 5
