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

# Sys.setenv(PATH = paste("C:/rtools40/mingw64/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")


## Define study area grid (random example currently)
locs <- as.matrix(secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8))

set.seed(922020)
a=sample(351, 120, replace=TRUE)   ## remember, have to change this if changing dimensions of trapping grid above
## Which parts of grid have traps
Xt=locs[a[1:60],]
Jt <- nrow(Xt)
ntraps <- nrow(Xt)

## Which parts of grid have visual surveys
Xv=locs[a[61:120],]
Jv <- nrow(Xv)
nvis <- nrow(Xv)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
delta <- 20  ## will need to play with this
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

## Establish parameters for encounter probability and movement (equivalent right now as these data sets mixed methods)
## Trapping
## Based on Amburgey et al. (in review) SCR (not the USCR; lam0 = 0.0036, sigma = 37.88), Gardner SSP analysis (lam0 = 0.008/0.007, sigma = ~ 38), and VIS lure/no lure analysis (lam0 = 0.0028/0.0034, sigma = 37.62)
tlam0 <- 0.005
tsigma <- 37.8 #, 20, 30

## Visual searching
vlam0 <- 0.005 #0.005, 0.008
vsigma <- 37.8 #37.8 #, 20, 30

## Number of snakes based on probable density in Guam (Rodda, Christy, etc.)
## 23 snakes/ha -> 0.0023 snakes/m2 (see notes at beginning)
N <- round(A*0.0023)

## Number of nights trapping (currently have same number trapping and visual searching)
K <- 20

nsim <- 5


for(iter in 1:nsim){
  print(iter)

  ## Simulate snake activity centers
  sx<-runif(N,Xl,Xu)
  sy<-runif(N,Yl,Yu)
  S<-cbind(sx,sy) 
  
  ## Take a moment to plot everything to make sure it makes sense
  library(sp)
  border <- as.data.frame(matrix(c(Xl,Xl,Xu,Xu,Yl,Yu,Yl,Yu), nrow = 4, ncol = 2))  ## study area
  test <- as.data.frame(locs)  ## trapping grid
  test2 <- as.data.frame(S)  ## activity centers
  test3 <- as.data.frame(Xt)  ## active traps
  test4 <- as.data.frame(Xv)  ## active traps
  coordinates(border) =~ V1 + V2
  coordinates(test) =~ x+y
  coordinates(test2) =~ sx+sy
  coordinates(test3) =~ x+y
  coordinates(test4) =~ x+y
  plot(border)
  plot(test, add=TRUE, pch=21, col="black", cex=0.4)
  plot(test2, add=TRUE, pch=22, col="blue", cex=0.7)
  plot(test3, add=TRUE, pch=1, col="red", cex=2)
  plot(test4, add=TRUE, pch=1, col="orange", cex=2)
    
  ## Distance between each individual AC and each each trap
  ## Function to calculate distance between two sets of (x,y) locations
  e2dist <- function(A, B)  {
    xdif <- outer(A[, 1], B[, 1], "-")
    ydif <- outer(A[, 2], B[, 2], "-")
    sqrt(xdif^2 + ydif^2)
  }
  
  ####### Modify this section for each monitoring method to generate observations #######
  ## Trapping
  tD <- e2dist(S,Xt)
  tmuy <- tlam0*exp(-(tD*tD)/(2*tsigma^2))
  
  ## Visual searching
  vD <- e2dist(S,Xv)
  vmuy <- vlam0*exp(-(vD*vD)/(2*vsigma^2))
  
  ## Simulate observations of individuals by trap
  Ytrap<-matrix(NA,nrow=N,ncol=ntraps)
  for(i in 1:nrow(Ytrap)){
    Ytrap[i,]<-rpois(ntraps,K*tmuy[i,])
  }
  
  ## Simulate observations of individuals by visual search
  Yvis<-matrix(NA,nrow=N,ncol=nvis)
  for(i in 1:nrow(Yvis)){
    Yvis[i,]<-rpois(nvis,K*vmuy[i,])
  }
  
  ## NOTE NOTE NOTE NOTE
  ## Y is a matrix of encounter frequencies of EACH individual in EACH monitoring method
  ## As simulated here it includes the "all 0" observations.  We want
  ## to delete those to mimic real data.
  ## Trapping
  totalcaps<-apply(Ytrap,1,sum)
  Ytrap<-Ytrap[totalcaps>0,]
  
  ## Visual searching
  totalvis<-apply(Yvis,1,sum)
  Yvis<-Yvis[totalvis>0,]
  
  datnam<-paste('/Users/Staci Amburgey/Documents/Optimization/SimDat/Dat_', iter, '.R', sep='')
  dput(rbind(Ytrap, Yvis), datnam)
  
  ## Data augmentation for when N is unknown
  M <- round(N+(N*4))  ## seems like a large amount of augmentation is needed based on traceplots
  ## Trapping
  ytrap=matrix(0, M,dim(Xt)[1])
  ytrap[1:dim(Ytrap)[1],]<-Ytrap
  
  ## Visual searching
  yvis=matrix(0, M,dim(Xv)[1])
  yvis[1:dim(Yvis)[1],]<-Yvis
  
  ## Initial values for s
  set.seed(932020)
  sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
  
  ## NIMBLE model is nearly identical to BUGS
  code <- nimbleCode({
    
    lam1~dunif(0,5)  ## trapping
    lam2~dunif(0,5)  ## searching
    sigma1~dgamma(274.69,7.27) #dunif(0,100) ## trapping
    sigma2~dunif(0,100) #dunif(0,100), dgamma(274.69,7.27)  ## searching
    psi~dunif(0,1)
    
    for(i in 1:M){
      z[i] ~ dbern(psi)
      s[i,1] ~ dunif(Xl,Xu)
      s[i,2] ~ dunif(Yl,Yu)
      
      for(j in 1:Jt){  ## traps
        d2tr[i,j] <- pow(s[i,1]-Xt[j,1],2) + pow(s[i,2]-Xt[j,2],2)
        ptr[i,j] <- z[i]*lam1*exp(-(d2tr[i,j])/(2*sigma1*sigma1))
        ytrap[i,j] ~ dpois(ptr[i,j]*K)
      }#j
      
      for(v in 1:Jv){  ## visual searches
        d2v[i,v] <- pow(s[i,1]-Xv[v,1],2) + pow(s[i,2]-Xv[v,2],2)
        pv[i,v] <- z[i]*lam2*exp(-(d2v[i,v])/(2*sigma2*sigma2))
        yvis[i,v] ~ dpois(pv[i,v]*K)
      }#j
    }#i
    N <- sum(z[1:M])
    D <- N/A
  })#code
  
  
  # MCMC settings
  nc <- 3; nAdapt=3000; nb <- 10000; ni <- 70000+nb; nt <- 1
  
  # Separate data and constants (constants appear only on right-hand side of formulas)
  nim.data <- list (ytrap=ytrap, yvis=yvis)
  constants <- list (Xt=Xt, Xv=Xv, K=K, M=M, Jt=Jt, Jv=Jv, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu, A=A)
  
  # Initial values (same as BUGS)
  inits <- function(){
    list (z=c(rep(1, N), rep(0,M-N)), psi=runif(1), sigma1=runif(1,20,50), sigma2=runif(1,20,50), lam1=runif(1,0.002,0.009), lam2=runif(1,0.002,0.009), s=sst)
  }

  # Parameters (same as BUGS)
  parameters <- c("sigma1","sigma2","lam1","lam2","N","D")
  
  
  ## Nimble steps
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=code, constants=constants, data=nim.data)
  conf <- configureMCMC(Rmodel,monitors=parameters,control = list(adaptInterval = nAdapt))
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  samplesList <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                         setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  end.time <- Sys.time()
  SCR0time<-end.time - start.time
  
  # #summarize
  # summaryList<-summary(samplesList)
  # outSummary<-cbind(summaryList$statistics[,c("Mean","SD")],summaryList$quantiles[,c("2.5%","50%", "97.5%")],gelman.diag(samplesList,multivariate=FALSE)$psrf[,1],effectiveSize(samplesList))
  # colnames(outSummary)[6:7]<-c("Rhat","n.eff")
  # round(outSummary,8)
  # 
  # #plot results
  # plot(samplesList[,"sigma1"])
  # plot(samplesList[,"N"])
  
  tosave <- as.matrix(samplesList)
  
  save(tosave, file=paste("OptimSim_trapsvisINPROGinfsigmaprior",iter,".csv",sep=""))
  
}
  
