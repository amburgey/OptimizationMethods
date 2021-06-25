## Simulate capture/observation histories of snakes in different sampling scenarios
## Save data to SimDat folder
## Analyze each dataset and save results in Results folder

rm(list=ls())

library(jagsUI)


#### STUDY AREA INFORMATION ----

# CP-sized, fenced area
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)

# Define state-space of point process. (i.e., where animals live).
# Option 1. Delta is the buffer of space between the end of transects and the fence
delta <- 11.874929
Xl <- min(totlocs[,1]) - delta
Xu <- max(totlocs[,1]) + delta
Yl <- min(totlocs[,2]) - delta
Yu <- max(totlocs[,2]) + delta
## Area of CP
A <- (Xu-Xl)*(Yu-Yl)

# Define state-space of point process. (i.e., where animals live).
# Option 2. Delta is the buffer of space around the open study area
# delta <- 20
# Xl <- min(totlocs[,1]) - delta
# Xu <- max(totlocs[,1]) + delta
# Yl <- min(totlocs[,2]) - delta
# Yu <- max(totlocs[,2]) + delta
# ## Area of CP
# A <- (Xu-Xl)*(Yu-Yl)


#### SURVEY INFORMATION ----

X <- as.matrix(totlocs)
# If applicable, subset locations to those monitored
# E.g., every 2nd transect
samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
X <- X[samp,]
J <- nrow(X)


#### INTEGRATION GRID ----

# Function to find distance between each integration grid point and each location of trap/visual survey/camera/etc.
e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

Ggrid <- 10                                #spacing
Xlocs <- rep(seq(Xl, Xu, Ggrid), times = length(seq(Yl, Yu, Ggrid)))
Ylocs <- rep(seq(Yl, Yu, Ggrid), each = length(seq(Xl, Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid^2                              #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")      #plot integration grid
points(X, pch=16, col="red")             #add locations of survey points


#### INDIVIDUAL-SPECIFIC INFO ----

# Set seed so we get the same values each time for true information about the population
set.seed(062420)
# True number of snakes
N <- 120
# Number of snakes per size category (4 groups; <850, >=850 to <950, >=950 to <1150, >=1150)
Nsnsz <- sample(c(1:4),N,replace=TRUE)
Ngroup <- as.vector(table(Nsnsz))
# True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE)


#### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----

# Parameters (p0, sigma) that influence detections of snakes will be pulled from real data posterior samples

yTrue <- matrix(NA,dim=c(N,J))
for(k in 1:K){
  D<- e2dist(S,locs)
  pmat<- plogis(alpha0)*exp(-alpha1*D*D)
  yTrue[,,k]<-rbinom(prod(dim(pmat)),1,pmat)
}

captured<-which(apply(yTrue,1,sum)>0)
yarr<-yTrue[captured,,]
y<- apply(yarr,c(1,2),sum)
nind <- length(captured)




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
  
  for(l in 1:L){  # size category
    for(g in 1:Gpts){ # Gpts = number of points on integration grid
      for(j in 1:J){  # J = number of traps
        #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
        miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j])
      } #J
      pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,]) #Prob of detect each size category across entire study area and time period
      pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
    } #G
    pstar[l] <- (sum(pdot[l,1:Gpts]*a[1:Gpts]))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
    
    # Zero trick for initial 1/pstar^n
    loglikterm[l] <- -ngroup[l] * log(pstar[l])
    lambda[l] <- -loglikterm[l] + 10000
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
",file = "Visual surveys/Models/SCRpstarCATsizeCAT_CP.txt")

#######################################################

# MCMC settings
nc <- 5; nAdapt=100; nb <- 10; ni <- 1000+nb; nt <- 1 

# Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup)

# Initial values (same as real data analysis)
inits <- function(){
  list (sigma=runif(1,30,40), n0=(ngroup+10), s=vsst, p0=runif(L,.002,.003))
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")

out <- jags("Simulations/Models/SCRpstarCATsizeCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE, n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE")

save(out, file="Simulations/Results/DataSim1TEST.Rdata")
