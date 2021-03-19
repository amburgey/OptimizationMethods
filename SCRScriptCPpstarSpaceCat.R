##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepVisualData.R")   ## Creation of subcap and subsurv (cleaned up)
source("DataPrepCP.R")              ## Functions to reshape survey and capture data

library(secr); library(reshape2); library(jagsUI)

## Subset capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project
CPcaps <- subset(CPcaps, PROJECTCODE == "NWFN VIS 2")
CPsurv <- subset(CPsurv, PROJECTCODE == "NWFN VIS 2")

##### SPECIFY DIMENSIONS OF CP #####
## Make study area grid to ensure correct size
## 27 transects (VIS and TRAP) with 13 points each, 8m from VIS to TRAP and 16m between points
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
J <- nrow(locs)

## Define state-space of point process. (i.e., where animals live).
## Don't need to estimate state-space since we know it (5 ha/50000 m2 enclosed pop) but do need this to help make integration grid below
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Area of CP
A <- (Xu-Xl)*(Yu-Yl)

##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
X <- as.matrix(locs)

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, T = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("V"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=c("02","04"))  ## this is using 3 months (Feb - April)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=c("02","04"))
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

#### FORMAT DATA FOR TRADITIONAL SCR ANALYSIS ####
dat <- prepSCR(SCRcaps, SCReff)

## Observations, already in order of 1-351 CellID locations
y <- dat$y
colnames(y) <- 1:ncol(dat$y)

## Uniquely marked individuals
nind <- nrow(y)

## Active/not active for when transects run, already in order of 1-351 CellID locations
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)

## Inits for activity centers, take mean grid cell location where each snake was found
locs <- as.data.frame(locs)
locs$CellID <- c(1:dim(locs)[1])
sst <- round(unlist(lapply(apply(y,1,function(x) which(x==1)),function(x) mean(x))))
vlocs <- locs[locs$CellID %in% sst,]
sst <- as.data.frame(sst); colnames(sst) <- c("CellID")
vsst <- merge(sst,vlocs, by=c("CellID"))[,2:3]

#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Integration grid
Ggrid <- 5                                #spacing (check sensitivity to spacing)
Xlocs <- rep(seq(Xl, Xu, Ggrid), times = 47)
Ylocs <- rep(seq(Yl, Yu, Ggrid), each = 44)
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")
points(Xl,Yu, pch=21, col="blue")         #check CP dimensions match
points(Xl,Yl, pch=21, col="blue")
points(Xu,Yu, pch=21, col="blue")
points(Xu,Yl, pch=21, col="blue")


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

  ## Removed k loop
  for(g in 1:Gpts){ # Gpts = number of points on integration grid
    for(j in 1:J){  # J = number of traps
      #Probability of being missed at grid cell g and trap j multiplied by total effort (K) at that trap
      one_minus_detprob[g,j] <- 1 - p0*exp(-alpha1*Gdist[g,j]*Gdist[g,j])*K[j] #Gdist given as data
    } #J
    pdot.temp[g] <- 1 - prod(one_minus_detprob[g,]) #Prob of failure to detect across entire study area and time period
    pdot[g] <- max(pdot.temp[g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
  } #G
  
  pstar <- (sum(pdot[1:Gpts])*a)/A   #prob of detecting an individual at least once in S (a=area of each integration grid, given as data)
  
  ##### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY ##### 
  # Zero trick for initial 1/pstar^n
  loglikterm <- -n * log(pstar)
  lambda <- -loglikterm + 1000
  dummy ~ dpois(lambda) # dummy = 0; entered as data
  ##### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY #####

  # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
 
  for(i in 1:n){  ## n = number of observed individuals
    ## For use when defining traps on a grid cell
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## J = number of traps
      y[i,j] ~ dbin(p[i,j],K[j])
      d[i,j] <- Gdist[s[i],j]  ## Doesn't Gdist already do what the line below does? As we're using categorical grid cells so s[i] are going to be one of G[i,]?
      # d[i,j] <- pow(pow(s[i,1]-locs[j,1],2) + pow(s[i,2]-locs[j,2],2),0.5)  ### traditional SCR
      p[i,j] <- p0*exp(-alpha1*d[i,j]*d[i,j])
    }#J
  }#n
}
",file = "SCRpstarCAT_CP.txt")

#######################################################

## MCMC settings
# nc <- 3; nAdapt=1000; nb <- 1; ni <- 10000+nb; nt <- 1
nc <- 3; nAdapt=50; nb <- 10; ni <- 1000+nb; nt <- 1

## Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=0, b=rep(1,Gpts), act=t(act)) # ## semicomplete likelihood
#locs=X, 

inits <- function(){
  list (sigma=runif(1,45,50), n0=nind, s=vsst, p0=runif(1,.002,.003)) #ran at 0.002 and 0.003 before
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N")

out <- jags("SCRpstarCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little

save(out, file="Results/NWFNVIS2_SCRpstarvistestCAT.Rdata")  ## M = 150 (XXXXhrs)


