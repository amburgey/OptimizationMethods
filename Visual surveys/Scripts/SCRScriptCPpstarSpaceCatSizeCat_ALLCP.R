#### Unified analysis of all CP (NWFN) datasets ####

rm(list=ls())

library(secr); library(reshape2); library(jagsUI)

source("Select&PrepVisualData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Visual surveys/DataPrep/OverlayCPGrid.R")

projects <- c("Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VIS2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL1.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPREBT2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTBT2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB1.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB3.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISVISTRAP.R")

noccall <- as.data.frame(matrix(NA, nrow = length(projects), ncol = 1))
Kall <- as.data.frame(matrix(NA, nrow = 351, ncol = length(projects)))
nindall <- as.data.frame(matrix(NA, nrow = length(projects), ncol = 1))
yall <- as.data.frame(matrix(NA, nrow = 1, ncol = 351));colnames(yall)=seq(1:ncol(yall))
snszall <- as.data.frame(matrix(NA, nrow = 1, ncol = 1));colnames(snszall)=seq(1:ncol(snszall))
Lall <- as.data.frame(matrix(NA, nrow = length(projects), ncol = 1))
ngroupall <- as.data.frame(matrix(NA, nrow = 1, ncol = 1));colnames(ngroupall)=seq(1:ncol(ngroupall))
vsstall <- as.data.frame(matrix(NA, nrow = 1, ncol = 1));colnames(vsstall)=seq(1:ncol(vsstall))

for(p in 1:length(projects)){
        ## Read in prep code for project
        ## The first data prep script includes loading functions and subsetting to just NWFN/CP
        source(projects[[p]])
        ## Take prepared data and add in order to create combined datasets
        ## Survey info
        noccall[p,1] <- nocc
        Kall[,p] <- K
        ## Capture info (will differ in number of individuals but will always be 351 columns)
        colnames(y) <- seq(1:ncol(y))
        yall <- rbind(yall, y)
        nindall[p,1] <- nind
        ## Snake size info
        snsz <- as.data.frame(snsz); colnames(snsz) = seq(1:ncol(snsz))
        snszall <- rbind(snszall, snsz)
        Lall[p,1] <- L
        ngroup <- as.data.frame(ngroup); colnames(ngroup) = seq(1:ncol(ngroup))
        ngroupall <- rbind(ngroupall, ngroup)
        ## Inits
        vsst <- as.data.frame(vsst); colnames(vsst) = seq(1:ncol(vsst))
        vsstall <- rbind(vsstall, vsst)

}

## Remove initial row of NAs needed to create dataframe
yall <- yall[-1,]
snszall <- snszall[-1,]
ngroupall <- ngroupall[-1,]
vsstall <- vsstall[-1,]
## Convert single column dataframes to vectors


#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

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

# cat("
# model {
#   
#   sigma ~ dunif(0,100)
#   alpha1 <- 1/(2*sigma*sigma)
#   
#   for(l in 1:L){   # 4 size categories
#     #prior for intercept
#     p0[l] ~ dunif(0,5)
#     alpha0[l] <- logit(p0[l])
#     
#     # Posterior conditional distribution for N-n (and hence N):
#     n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
#     Ngroup[l] <- ngroup[l] + n0[l]
#   }
#   
#   N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
#   
#   #Probability of capture for integration grid points
#   #pdot = probability of being detected at least once (given location)
#   
#   for(l in 1:L){  # size category
#     for(g in 1:Gpts){ # Gpts = number of points on integration grid
#       for(j in 1:J){  # J = number of traps
#         #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
#         miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j])
#       } #J
#       pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,]) #Prob of detect each size category across entire study area and time period
#       pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
#     } #G
#     pstar[l] <- (sum(pdot[l,1:Gpts]*a[1:Gpts]))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
#     
#     # Zero trick for initial 1/pstar^n
#     loglikterm[l] <- -ngroup[l] * log(pstar[l])
#     lambda[l] <- -loglikterm[l] + 10000
#     dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
#   } #L
#   
#   # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
#   pi[1:Gpts] ~ ddirch(b[1:Gpts])
#   
#   for(i in 1:n){  ## n = number of observed individuals
#     ## For use when defining traps on a grid cell
#     s[i] ~ dcat(pi[1:Gpts])
#     
#     # Model for capture histories of observed individuals:
#     for(j in 1:J){  ## J = number of traps
#       y[i,j] ~ dpois(p[i,j]*K[j])
#       p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
#     }#J
#   }#I
#   
#   #derived proportion in each size class
#   for(l in 1:L){
#     piGroup[l] <- Ngroup[l]/N
#   }
# }
# ",file = "Visual surveys/Models/SCRpstarCATsizeCAT_CP.txt")
# 
# #######################################################
# 
# ## MCMC settings
# nc <- 3; nAdapt=200; nb <- 100; ni <- 2500+nb; nt <- 1
# 
# ## Data and constants
# jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup) # ## semicomplete likelihood
# # jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=0, b=rep(1,Gpts))
# 
# 
# inits <- function(){
#   list (sigma=runif(1,30,40), n0=(ngroup+10), s=vsst, p0=runif(L,.002,.003))
# }
# # inits <- function(){
# #   list (sigma=runif(1,30,40), n0=(nind+30), s=vsst, p0=runif(1,.002,.003))
# # }
# 
# parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")
# # parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0")
# 
# out <- jags("Visual surveys/Models/SCRpstarCATsizeCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
#             n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little
# # out <- jags("Archive/SCRpstarCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
#             # n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little
# 
# 
# save(out, file="Visual surveys/Results/NWFNVIS2_SCRpstarvisCATsizeCATdpois10GRID.Rdata")
# # save(out, file="Visual surveys/Results/NWFNVIS2_SCRpstarvisCATNOSIZEgrid10.Rdata")

