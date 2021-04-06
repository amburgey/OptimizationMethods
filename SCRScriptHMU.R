##### HMU (Habitat Management Unit) is a 55-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepVisualData.R")  ## Creation of subcap and subsurv
source("OverlayHMUGrid.R")
source("DataPrepHMU.R")

library(secr); library(reshape2); library(jagsUI)

## Capture data (subcap) and effort/survey data (subsurv)
HMUcaps <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR"))
HMUsurv <- subset(subsurv, (SITE == "HMUI" | SITE == "HMUR"))

## Subset to specific NWFN project (options = MWFM VIS 2, NWFN VIS HL 1, NWFN VIS HL 2, PRE NT2 VIS, POST BT2 VIS, POST KB VIS 1, POST KB VIS 2, POST KB VIS 3 EXTRA, POST KB VIS 3, NWFN VISPACE, NWFN SCENT VIS TRAIL)
HMUcaps <- subset(HMUcaps, PROJECTCODE == "EDGE EFFECT VIS")[,c("EFFORTID","PITTAG","Date","TRANSECT","LOCATION","CAPLAT","CAPLON")]
HMUsurv <- subset(HMUsurv, PROJECTCODE == "EDGE EFFECT VIS")

##### SPECIFY DIMENSIONS AND GRID OF HMU #####
HMUspecs <- overlayHMU(HMUcaps)  ## ignore warnings, all about projections
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
dat <- prepSCR(SCRcaps, SCReff, grid = HMUspecs$tran)  ## if error here do one below
dat <- prepSCRman(SCRcaps, SCReff, grid = HMUspecs$tran) ## MANUALLY CHANGE DISTANCES

## Observations
y <- dat$y
colnames(y) <- NULL

## Uniquely marked individuals
nind <- nrow(y)

## Active/not active for when transects run (one less day surveyed for TL)
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)


## Inits for activity centers, take mean grid cell location where each snake was found
locs <- HMUspecs$tran
colnames(locs)[2] <- c("CellID")
## Can't take mean location of cells here because not all cells surveyed in HMU, just pick one
sst <- unlist(lapply(apply(dat$y,1,function(x) which(x==1)),function(x) min(as.numeric(names(x)))))
vlocs <- locs[locs$CellID %in% sst,]
sst <- as.data.frame(sst); colnames(sst) <- c("CellID")
vsst <- unlist(merge(sst,vlocs, by=c("CellID"))[,1])

#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#Integration grid
Ggrid <- 10                               #spacing (verify sensitivity to spacing)
G <- HMUspecs$intgrd
# Xlocs <- seq(Yl,Yu,Ggrid)          
# G <- cbind(sort(rep(Xlocs,length(Xlocs))),rep(Xlocs,length(Xlocs))) #integration grid locations
Gpts <- dim(G)[1]                            #number of integration points
a <- Ggrid^2                                 #area of each integration grid
Gdist <- e2dist(G, X)                      #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##Jags model for a King et al 2016 semicomplete likelihood

cat("
model {

  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  sigma ~ dunif(0,100)
  alpha1 <- 1/(2*sigma*sigma)
  
  # Posterior conditional distribution for N-n (and hence N):
  n0 ~ dnegbin(pstar,n)  # number of failures
  N <- n + n0  # successful observations plus failures to observe = total N
  
  #Probability of capture for integration grid points
  #pdot = probability of being detected at least once (given location)

## Remove k loop (for decreased run time?)
  for(g in 1:Gpts){ # Gpts = number of points on integration grid
    for(j in 1:J){  # J = number of traps
      one_minus_detprob[g,j] <- 1 - p0*exp(-alpha1*Gdist[g,j]*Gdist[g,j])*K[j]
    } #J
  pdot.temp[g] <- 1 - prod(one_minus_detprob[g,])
  pdot[g] <- max(pdot.temp[g], 1.0E-10)
} #G
  
  pstar <- (sum(pdot[1:Gpts])*a)/A   #prob of detecting an individual at least once in S
  
  ### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY ### 
  # Zero trick for initial 1/pstar^n
  loglikterm <- -n * log(pstar)
  lambda <- -loglikterm + 1000
  dummy ~ dpois(lambda) # dummy = 0; entered as data
  ### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY ###
  
  # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior)   
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
  
  for(i in 1:n){  ## n = number of observed individuals
    ## For use when defining traps on a grid cell
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## J = number of traps
      y[i,j] ~ dbin(p[i,j],K[j])
      d[i,j] <- Gdist[s[i],j]
      p[i,j] <- p0*exp(-alpha1*d[i,j]*d[i,j])
    }#J
  }#n
}
",file = "SCRpstarCAT_HMU.txt")

#######################################################

## MCMC settings
# nc <- 3; nAdapt=1000; nb <- 1; ni <- 10000+nb; nt <- 1
nc <- 3; nAdapt=5; nb <- 10; ni <- 20+nb; nt <- 1

## Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=0, b=rep(1,Gpts), act=t(act)) # ## semicomplete likelihood

inits <- function(){
  list (sigma=runif(1,45,50), n0=nind, s=vsst, p0=runif(1,.002,.003)) #ran at 0.002 and 0.003 before
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N")

out <- jags("SCRpstarCAT_HMU.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use to keep JAGS from locking up with large categorical distribution, will speed things up a little

save(out, file="Results/HMUEDGE_SCRpstarvistestCAT2months10.Rdata")  ## M = 150 (XXXXhrs)

## ran for 7 hours and then hit an invalid parent value error (error in node n0)
