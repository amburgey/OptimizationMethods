#### Unified analysis of all CP (NWFN) datasets ####

rm(list=ls())

memory.limit(5000000000)

library(secr)
library(reshape2)
library(nimble)
library(abind)
library(coda)
library(foreach)
library(doParallel)

Sys.setenv(PATH = paste("C:/rtools40/usr/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")

source("Select&PrepVisualData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Visual surveys/DataPrep/OverlayCPGrid.R")

projects <- c("Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VIS2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL1.R")#,
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL2.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPREBT2.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTBT2.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB1.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB2.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB3.R",
        # "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISVISTRAP.R")

nproj <- length(projects)

noccall <- matrix(NA, nrow = nproj, ncol = 1)
Kall <- matrix(NA, nrow = 351, ncol = nproj)
nindall <- matrix(NA, nrow = nproj, ncol = 1)
yall <- array(c(NA,NA,NA), dim=c(108,351,nproj)) ## maximize number of individuals is 108 so the biggest number of rows ever needed
snszall <- matrix(NA, nrow = 108, ncol = nproj)
# Lall <- as.data.frame(matrix(NA, nrow = nproj, ncol = 1))
ngroupall <- matrix(NA, nrow = 4, ncol = nproj)
vsstall <- matrix(NA, nrow = 108, ncol = nproj)

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
ngroupall <- ngroupall
vsstall[is.na(vsstall)] <- 0
snszall[is.na(snszall)] <- 0
## Restrict to max number of individuals for subsetting testing
yall <- yall[1:max(nindall),,]
snszall <- snszall[1:max(nindall),]


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



## MCMC settings
# nc <- 3; nAdapt=1000; nb <- 1000; ni <- 10000+nb; nt <- 1
nc <- 3; nAdapt=10; nb <- 10; ni <- 10+nb; nt <- 1

## Data and constants
constants <- list(J=J, A=A, Gpts=Gpts, nFound=nindall, nocc=nocc, L=Lall, nproj=nproj, size=snszall, a=a, K=Kall) 

data <- list(y=yall, ngroup=ngroupall, Gdist=Gdist, dummy=matrix(0,nrow=Lall,ncol=nproj), b=rep(1,Gpts))

inits <- list (sigma=runif(1,30,40), n0=(ngroupall+10))

parameters <- c("p0","sigma","N","n0","Ngroup","alpha0","tau_p")


#### PARALLEL STUFF #####
t.start <- Sys.time() # start clock for whole process

cores=detectCores() # how many cores are available
# this line is essential on mac, not sure about pc
# does not impact performance though
used.cores = min(nc, cores) # good practice to not max out all available cores. 
# my general workflow is that in each instance of R I set up a script like this 
# that runs 3 chains in parallel (1 chain per core)
cl <- makeCluster(nc, setup_strategy = "sequential") # make cluster of cores
registerDoParallel(cl) 

foreach(i = 1:nc) %dopar% { # parallel version of for loop
  # my understanding is that each core can access your environment but not workspace 
  # i.e. you can reference data objects but you have to reload packages and functions
  # do that here
  library(nimble)


  ########################################################
  ##NIMBLE model for a King et al 2016 semicomplete likelihood
  
  NimModel <- nimbleCode({
  
    #prior for spatial decay
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    # prior for precision random effect of project
    tau_p ~ dunif(0,4)
    # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    for(l in 1:L){   # 4 size categories
      alpha0[l] ~ dnorm(0,0.1)
      
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
      # prior for random effect of project
      eta[t] ~ dnorm(0, tau_p)
      
      for(l in 1:L){  # size category
          
      # fixed intercept of size and random effect of study
      logit(p0[l,t]) <- alpha0[l] + eta[t]
      
        for(g in 1:Gpts){ # Gpts = number of points on integration grid
          for(j in 1:J){  # J = number of traps
            #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
            miss_allK[l,g,j,t] <- pow((1 - p0[l,t]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j,t])
          } #J
          pdot.temp[l,g,t] <- 1 - prod(miss_allK[l,g,1:J,t]) #Prob of detect each size category across entire study area and time period
          pdot[l,g,t] <- max(pdot.temp[l,g,t], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
        } #G
        pstar[l,t] <- (sum(pdot[l,1:Gpts,t]*a[1:Gpts]))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
  
        # Zero trick for initial 1/pstar^n
        loglikterm[l,t] <- -ngroup[l,t] * log(pstar[l,t])
        lambda[l,t] <- -loglikterm[l,t] + 10000
        dummy[l,t] ~ dpois(lambda[l,t]) # dummy = 0; entered as data
      } #L
  
      for(i in 1:nFound[t]){  ## nFound = maximum number of observed individuals in each project
        ## For use when defining traps on a grid cell
        s[i,t] ~ dcat(pi[1:Gpts])
  
        # Model for capture histories of observed individuals:
        for(j in 1:J){  ## J = number of traps
          y[i,j,t] ~ dpois(p[i,j,t]*K[j,t])
          p[i,j,t] <- p0[size[i,t],t]*exp(-alpha1*Gdist[s[i,t],j]*Gdist[s[i,t],j])
        }#J
      }#I
  
      #derived proportion in each size class
      for(l in 1:L){
        piGroup[l,t] <- Ngroup[l,t]/N[t]
      }#L
    }#T
  })
  
  
  #######################################################
  
  # ## NIMBLE functions
  # CorrectPstar <- nimbleFunction(
  #   run = function(ps = double(0)){
  #     returnType(double(1))
  #     if(ps>=1) return(1)
  #     if(ps<1) return(ps)
  #   }
  # )
  
  
  ## Compile and run in NIMBLE
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=data, inits=inits, dimensions = list(p=c(max(nindall),J,nproj), y=c(max(nindall),J,nproj), Gdist=c(Gpts,J), size=c(max(nindall),nproj), s=c(max(nindall),nproj)), check=FALSE, calculate = FALSE)
  conf <- configureMCMC(Rmodel, monitors=parameters, control = list(adaptInterval = nAdapt), thin=nt) 
  
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
  Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
  
  out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  

  saveRDS(out, paste("Visual surveys/Results/outALL-",i,".RDS", sep = "")) # save each chain with a diff name
  
}

stopCluster(cl) # important to do this to stop running things in parallel

t.end <- Sys.time() # end clock
(runTime <- t.end - t.start) # how long did whole thing take in parallel



## Diagnostics

#### LOAD AND COMBINE RESULTS ####

# for each chain, load the resulting MCMC output
for (i in 1:nc) {
  readRDS(here("Models", paste("out-", i,".RDS", sep = "")))
}

# combine them and make an mcmc object
OUT <- list(`outALL-1`, `outALL-2`, `outALL-3`) %>% as.mcmc.list()
# now you can treat it just like you ran the chains in succession
library(coda)
OUT <- out
rhat <- gelman.diag(OUT, multivariate = FALSE) 


# post <- as.matrix(out)
# 
# #for median and CI
# post_sum <- data.frame(
#   mean = apply(post, 2, function(x) mean(x)),
#   med = apply(post, 2, function(x) quantile(x, probs = 0.5, na.rm = T, names = F)),
#   lower = apply(post, 2, function(x) quantile(x, probs = 0.025, na.rm = T, names = F)),
#   upper = apply(post, 2, function(x) quantile(x, probs = 0.975, na.rm = T, names = F)))
# post_sum$variable <- row.names(post_sum)
# 
# all_pars <- colnames(post)
# 
# #traceplot
# coda::traceplot(out[,all_pars[which(grepl('alpha0', all_pars))]])
# coda::traceplot(out[,all_pars[which(grepl('tau', all_pars))]])
# coda::traceplot(out[,all_pars[which(grepl('N', all_pars))]])
# coda::traceplot(out)
# 
# #rhats
# gelman.diag(out[,c('N[1]', 'N[2]')], multivariate=F, autoburnin=F)
# gelman.diag(out, multivariate=F, autoburnin=F)
