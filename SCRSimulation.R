#### SIMULATION USING SCR FRAMEWORK ####
## Detect snakes based on different detection process for each method
## Use sigma info from "normal" snake density to start (23 snakes/ha; Rodda et al. 1999, Tyrrell et al. 2009, Christy et al. 2010)
## sigma ~ dgamma(274.69, 7.27) based on 2015 SCR analysis
## sigma mean (CI; SD) = 37.8 (33.6, 42.4; 2.3)

rm(list=ls())
library(secr)
library(nimble)
library(coda)
library(devtools)
set.seed(2018)

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")


## Define study area grid (random example currently)
locs <- as.matrix(secr::make.grid(nx = 10, ny = 10, spacex = 8, spacey = 8))
ntraps <- nrow(locs)

## Which parts of grid have traps
set.seed(922020)
# a=sample(100, 15)   ## remember, have to change this if changing dimensions of trapping grid above
a <- c(1:100)
X=locs[a,]
J <- nrow(X)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
delta <- 20  ## will need to play with this
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

## Establish parameters for detection
## Based on Amburgey et al. (in review) SCR (not the USCR; lam0 = 0.0036, sigma = 37.88), Gardner SSP analysis (lam0 = 0.008/0.007, sigma = ~ 38), and VIS lure/no lure analysis (lam0 = 0.0028/0.0034, sigma = 37.62)
lam0 <- 0.005
sigma <- 37.8 #, 20, 30

## Number of snakes based on probable density in Guam (Rodda, Christy, etc.)
## 23 snakes/ha -> 0.0023 snakes/m2 (see notes at beginning)
N <- round(A*0.0023)

## Number of nights trapping
K <- 20

# nsim <- 1
# 
# 
# for(iter in 1:nsim){
#   print(iter)

  ## Simulate activity centers
  sx<-runif(N,Xl,Xu)
  sy<-runif(N,Yl,Yu)
  S<-cbind(sx,sy) 
  
  ## Take a moment to plot everything to make sure it makes sense
  library(sp)
  border <- as.data.frame(matrix(c(Xl,Xl,Xu,Xu,Yl,Yu,Yl,Yu), nrow = 4, ncol = 2))  ## study area
  test <- as.data.frame(locs)  ## trapping grid
  test2 <- as.data.frame(S)  ## activity centers
  test3 <- as.data.frame(X)  ## active traps
  coordinates(border) =~ V1 + V2
  coordinates(test) =~ x+y
  coordinates(test2) =~ sx+sy
  coordinates(test3) =~ x+y
  plot(border)
  plot(test, add=TRUE, pch=21, col="blue", cex=0.5)
  plot(test2, add=TRUE, pch=22, col="red", cex=0.7)
  plot(test3, add=TRUE, pch=1, cex=2)
    
  ## Distance between each individual AC and each each trap
  ## Function to calculate distance between two sets of (x,y) locations
  e2dist <- function(A, B)  {
    xdif <- outer(A[, 1], B[, 1], "-")
    ydif <- outer(A[, 2], B[, 2], "-")
    sqrt(xdif^2 + ydif^2)
  }
  
  ####### Modify this section for each monitoring method to generate observations #######
  D <- e2dist(S,locs)
  muy <- lam0*exp(-(D*D)/(2*sigma^2))
  
  ## Simulate observations of individuals by trap
  Y<-matrix(NA,nrow=N,ncol=ntraps)
  for(i in 1:nrow(Y)){
    Y[i,]<-rpois(ntraps,K*muy[i,])
  }
  
  ## NOTE NOTE NOTE NOTE
  ## Y is a matrix of encounter frequencies of EACH individual in EACH trap
  ## As simulated here it includes the "all 0" observations.  We want
  ## to delete those to mimic real data.
  Yscr=Y[,a]
  totalcaps<-apply(Yscr,1,sum)
  Yscr<-Yscr[totalcaps>0,]
  
  # datnam<-paste('/SimDat/Dat_', iter, '.R', sep='')
  # dput(Yscr, datnam)
  
  ## Data augmentation for when N is unknown
  M <- round(N+(N*4))  ## seems like a large amount of augmentation is needed based on traceplots
  y=matrix(0, M,dim(X)[1])
  y[1:dim(Yscr)[1],]<-Yscr
  
  ## Initial values for s
  set.seed(932020)
  sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
  for(i in 1:N){
    sst[i,1] <- mean( X[y[i,]>0,1] )
    sst[i,2] <- mean( X[y[i,]>0,2] )
  }
  
  ## NIMBLE model is nearly identical to BUGS
  code <- nimbleCode({
    lam0~dunif(0,5)
    sigma~dunif(0,100) # informative prior = dgamma(274.69,7.27) based on 2015 SCR study, really bad results when using this
    psi~dunif(0,1)
    
    for(i in 1:M){
      z[i] ~ dbern(psi)
      s[i,1] ~ dunif(Xl,Xu)
      s[i,2] ~ dunif(Yl,Yu)
      
      for(j in 1:J){
        d2[i,j] <- pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2)
        p[i,j] <- z[i]*lam0*exp(-(d2[i,j])/(2*sigma*sigma))
        
        y[i,j] ~ dpois(p[i,j]*K)
      }#j
    }#i
    N <- sum(z[1:M])
    D <- N/A
  })#code
  
  
  # MCMC settings
  nc <- 3; nAdapt=5000; nb <- 20000; ni <- 80000+nb; nt <- 1
  
  # Separate data and constants (constants appear only on right-hand side of formulas)
  nim.data <- list (y=y)
  constants <- list (X=X, K=K, M=M, J=J, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A)
  
  # Initial values (same as BUGS)
  inits <- function(){
    list (z=c(rep(1, N), rep(0,M-N)), psi=runif(1), sigma=runif(1,1,50), lam0=runif(1,0.002,0.009), s=sst)
  } # lam0=runif(1,0.5,1.5)

  # Parameters (same as BUGS)
  parameters <- c("sigma","lam0","N","D")
  
  
  ## Nimble steps
  start.time <- Sys.time()  
  Rmodel <- nimbleModel(code=code, constants=constants, data=nim.data)
  conf <- configureMCMC(Rmodel,monitors=parameters,control = list(adaptInterval = nAdapt))
  Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
  Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  samplesList <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                         setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  end.time <- Sys.time()
  SCR0time<-end.time - start.time
  
  #summarize
  summaryList<-summary(samplesList)
  outSummary<-cbind(summaryList$statistics[,c("Mean","SD")],summaryList$quantiles[,c("2.5%","50%", "97.5%")],gelman.diag(samplesList,multivariate=FALSE)$psrf[,1],effectiveSize(samplesList))
  colnames(outSummary)[6:7]<-c("Rhat","n.eff")
  round(outSummary,8)
  
  #plot results
  plot(samplesList[,"sigma"])
  plot(samplesList[,"N"])
  
# }
