### Toy simulation to test issues of unobservable animals for SCR abundance estimation


rm(list=ls())
library(jagsUI); library(secr); library(ggplot2)


## Function for Euclidean distance between observation locations and integration grid points
e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


## Information about true population
N <- 120
sigma <- 50
alpha1 <- 1/(2*sigma*sigma)
p0 <- 0.003


## Information about sampling
K <- 30
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
delta <- 11.874929
Xl <- min(totlocs[,1]) - delta
Xu <- max(totlocs[,1]) + delta
Yl <- min(totlocs[,2]) - delta
Yu <- max(totlocs[,2]) + delta
# Sample half the transects
samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
X <- as.matrix(totlocs)[samp,]
J <- nrow(X)


## Integration grid
Ggrid <- 10
Xlocs <- rep(seq(Xl, Xu, Ggrid), times = length(seq(Yl, Yu, Ggrid)))
Ylocs <- rep(seq(Yl, Yu, Ggrid), each = length(seq(Xl, Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
A <- Gpts * a                             #area of study area
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")       #plot integration grid
points(X, pch=16, col="red")              #add locations of survey points


## Simulate activity centers (for all)
set.seed(062420)
s <- sample(1:Gpts,N,replace=TRUE)


## Simulate observations
yTrue <- array(NA,dim=c(N,J))

pmat <- p0*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their activity centers) at all locations
for(n in 1:N){
  yTrue[n,] <- rpois(J,pmat[n,]*K)  # observations of each snake at each trap based on encounter probability and effort
}

## yTrue includes a row for every snake even if that snake was never observed. We need to check and remove these snakes to mimic real data.
captured <- which(apply(yTrue[,1:J],1,sum)>0)  # snakes that were observed at least once
yarr <- yTrue[captured,]  # subset to observed snakes (is all snakes still here)


## Simulate snakes that are unobservable (e.g., lack PIT tag)
rsnks <- sample(1:nrow(yarr),40,replace=FALSE)  ## snakes we can't observed
yobs <- yarr[-rsnks,]
nind <- nrow(yobs)
vsst <- list()
for(w in 1:nrow(yobs)){
  vsst[w] <- apply(yobs,1,function(x) which(x>=1))[[w]][1]
  vsst <- unlist(vsst)
}



########################################################
##Jags model for a King et al 2016 semicomplete likelihood

### VIS survey model using Poisson observation process ###

cat("
  model {
    
    #prior for spatial decay
    sigma ~ dunif(0,100)
    alpha1 <- 1/(2*sigma*sigma)
    #prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    
    #prior for intercept
    p0 ~ dunif(0,5)
      
    # Posterior conditional distribution for N-n (and hence N):
    n0 ~ dnegbin(pstar,n)  # number of failures by size category
    N <- n + n0
    
    #Probability of capture for integration grid points
    #pdot = probability of being detected at least once (given location)
  
    for(g in 1:Gpts){ # Gpts = number of points on integration grid
      for(j in 1:J){  # J = number of traps
        #Probability of an individual being missed at grid cell g and trap j multiplied by total effort (K) at that trap
        miss_allK[g,j] <- pow((1 - p0*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K)
      } #J
      pdot.temp[g] <- 1 - prod(miss_allK[g,]) #Prob of detect each size category across entire study area and time period
      pdot[g] <- max(pdot.temp[g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
    } #G
    pstar <- (sum(pdot[1:Gpts]*a))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
      
    # Zero trick for initial 1/pstar^n
    loglikterm <- -n * log(pstar)
    lambda <- -loglikterm + 10000
    dummy ~ dpois(lambda) # dummy = 0; entered as data
    
    for(i in 1:n){  ## n = number of observed individuals
      ## For use when defining traps on a grid cell
      s[i] ~ dcat(pi[1:Gpts])
      
      # Model for capture histories of observed individuals:
      for(j in 1:J){  ## J = number of traps
        y[i,j] ~ dpois(p[i,j]*K)
        p[i,j] <- p0*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
      }#J
    }#I
  }
",file = "Simulations/Models/SCRpstarCAT_SimUNMKtoMK.txt")


## MCMC settings
nc <- 3; nAdapt=100; nb <- 1; ni <- 3000+nb; nt <- 1 

jags.data <- list (y=yobs, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=0, b=rep(1,Gpts))

inits <- function(){
  list (sigma=runif(1,50,60), n0=nind, s=vsst, p0=runif(1,.001,.002))
}

parameters <- c("p0","sigma","pstar","alpha1","N","n0")

out <- jags("Simulations/Models/SCRpstarCAT_SimUNMKtoMK.txt", data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")

save(out, file="Simulations/Results/RESULTS_SimUNMKtoMK40.Rdata")



## Simulate snakes that are unobservable (e.g., lack PIT tag)
rsnks <- sample(1:nrow(yarr),20,replace=FALSE)  ## snakes we can't observed
yobs <- yarr[-rsnks,]
nind <- nrow(yobs)
vsst <- list()
for(w in 1:nrow(yobs)){
  vsst[w] <- apply(yobs,1,function(x) which(x>=1))[[w]][1]
  vsst <- unlist(vsst)
}

jags.data <- list (y=yobs, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=0, b=rep(1,Gpts))

inits <- function(){
  list (sigma=runif(1,50,60), n0=nind, s=vsst, p0=runif(1,.001,.002))
}

out <- jags("Simulations/Models/SCRpstarCAT_SimUNMKtoMK.txt", data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")

save(out, file="Simulations/Results/RESULTS_SimUNMKtoMK20.Rdata")



## Simulate snakes that are unobservable (e.g., lack PIT tag)
rsnks <- sample(1:nrow(yarr),60,replace=FALSE)  ## snakes we can't observed
yobs <- yarr[-rsnks,]
nind <- nrow(yobs)
vsst <- list()
for(w in 1:nrow(yobs)){
  vsst[w] <- apply(yobs,1,function(x) which(x>=1))[[w]][1]
  vsst <- unlist(vsst)
}

jags.data <- list (y=yobs, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=0, b=rep(1,Gpts))

inits <- function(){
  list (sigma=runif(1,50,60), n0=nind, s=vsst, p0=runif(1,.001,.002))
}

out <- jags("Simulations/Models/SCRpstarCAT_SimUNMKtoMK.txt", data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")

save(out, file="Simulations/Results/RESULTS_SimUNMKtoMK60.Rdata")



## Simulate snakes that are unobservable (e.g., lack PIT tag)
rsnks <- sample(1:nrow(yarr),80,replace=FALSE)  ## snakes we can't observed
yobs <- yarr[-rsnks,]
nind <- nrow(yobs)
vsst <- list()
for(w in 1:nrow(yobs)){
  vsst[w] <- apply(yobs,1,function(x) which(x>=1))[[w]][1]
  vsst <- unlist(vsst)
}

jags.data <- list (y=yobs, Gpts=Gpts, Gdist=Gdist, J=J, A=A, K=K, a=a, n=nind, dummy=0, b=rep(1,Gpts))

inits <- function(){
  list (sigma=runif(1,50,60), n0=nind, s=vsst, p0=runif(1,.001,.002))
}

out <- jags("Simulations/Models/SCRpstarCAT_SimUNMKtoMK.txt", data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")

save(out, file="Simulations/Results/RESULTS_SimUNMKtoMK80.Rdata")




#### Load results and plot to see trends
inst <- c(20,40,60,80)
vals <- data.frame(matrix(NA, nrow = 4, ncol = 4, dimnames = list(1:4, c("EstN","Q2.5","Q97.5","Type"))))

for(i in 1:4){
  load(paste("Simulations/Results/RESULTS_SimUNMKtoMK",inst[i],".Rdata", sep=""))
  vals[i,1] <- out$mean$N
  vals[i,2] <- out$q2.5$N
  vals[i,3] <- out$q97.5$N
  vals[i,4] <- paste("Removed",inst[i], sep="")
}

p1 <- ggplot(vals, aes(x = Type, y = EstN)) +
  geom_point() +
  geom_pointrange(aes(ymin = `Q2.5`, ymax = `Q97.5`)) +
  geom_hline(yintercept = 120, col = "blue") +
  annotate("text", x = 2.5, y = 118, label = "Truth = 120 snakes", col = "blue")

