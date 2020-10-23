### USING CP, TRY OUT DIFFERENT PERIODS OF SAMPLING TO SEE WHEN ESTIMATES PLATEAU AND VARIANCE ASYMPTOTES ###
### THIS WILL HELP SET THE WINDOW OF TIME THAT WE CAN LIMIT SAMPLING TO SO WE DON'T NEED TO INCLUDE A CORRELATED ERROR STRUCTURE ###

rm(list=ls())

library(abind); library(nimble)

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

## VISUAL DATA ##

# Survey information
allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VISinfo.csv")[,-c(21:62)]
CPvissurv <- subset(allsurv, SITE == "NWFN" & PROJECTCODE == "NWFN VISTRAP VIS")
CPvissurv$Date2 <- as.Date(as.character(CPvissurv$Date), format = "%d-%b-%y")

# Capture information
allcap <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VIScaptures.csv")[,-c(6,15,22:25,27,29:32,34:36)]
CPviscap <- subset(allcap, SITE == "NWFN" & PROJECTCODE == "NWFN VISTRAP VIS")
CPviscap$Date2 <- as.Date(as.character(CPviscap$Date), format = "%d-%b-%y")
## Add transect ID
CPviscap$TRANID <- paste(CPviscap$TRANSECT,CPviscap$LOCATION, sep="")


#### CHECK FOR MISSING LOCATION INFO ####
remov <- subset(CPvissurv, is.na(TRANSECT) | is.na(STARTNUMBER) | TRANSECT == "" | STARTNUMBER == "")  ## all are incidental
'%!in%' <- Negate('%in%')
subsurv <- CPvissurv[CPvissurv$TRANSID %!in% as.vector(remov$TRANSID) & CPvissurv$EFFORTID %!in% as.vector(remov$EFFORTID), ]
subcap <- CPviscap[CPviscap$TRANSID %!in% as.vector(remov$TRANSID) & CPviscap$EFFORTID %!in% as.vector(remov$EFFORTID), ]

#### REMOVE EDGE TRANSECTS AS NOT ENOUGH DATA TO INCLUDE AND MODEL SEPARATELY ADEQUATELY ####
subsurv <- subsurv[!(subsurv$TRANSECT=="1PR"| subsurv$TRANSECT=="SWE" | subsurv$TRANSECT=="NEE"),]
subcap <- subcap[!(subcap$TRANSECT=="1PR"| subcap$TRANSECT=="SWE" | subcap$TRANSECT=="NEE"),]

#### CHECK FOR MISSING SURVEYS ####
###### Check that dates of captures match dates when surveys were conducted, this will flag for captures without survey dates but not the other way around (because snakes weren't caught at every survey so wouldn't match) ######
test <- ifelse(CPviscap$TRANSID %in% CPvissurv$TRANSID & CPviscap$EFFORTID %in% CPvissurv$EFFORTID, 1, 0)
for(i in 1:length(test)){
  if(test[i] == 0){
    stop('survey info does not exist for all captures')  ## originally flagged two records, fixed above so shouldn't result in error now
  }
}


#### SHAPE DATA FOR ANALYSIS ####
ord <- paste(rep(c(LETTERS[1:26],"AA"), each=13),seq(1:13),sep="")

survs <- subsurv[,c("Date2","TRANSECT")]
survs <- unique(survs)# Remove extra survey from second searcher
survs <- survs[rep(seq_len(nrow(survs)), each = 13), ]
survs$TRANID <- paste(survs$TRANSECT, rep(1:13, times = 422), sep = "")
survs$Active <- c(1)
caps <- subcap[,c(1:2,5,6,15:16,23:24)]
eff <- merge(survs,caps, by = c("TRANSECT","TRANID","Date2"), all = TRUE)
## TRANID by Date matrix of effort (when traps were active)
eff2 <- reshape2::dcast(eff, TRANID ~ Date2, fun.aggregate = length, value.var = "Active")
eff2 <- eff2[order(match(eff2$TRANID, ord)), ]
## PITTAG by Date matrix of captures
snks <- reshape2::dcast(data = eff, formula = PITTAG ~ TRANID, fun.aggregate = length, value.var = "LOCATION")[-(length(unique(subcap$PITTAG))+1),]  # adds an extra row
snks <- snks[ord]


#### SPECIFY STATE SPACE ####
# Create trapping grid of CP dimensions (5 ha, 50,000 m2)
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
X <- as.matrix(locs)

## Number of transect grid locations
J <- nrow(X)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
## Don't need to estimate state-space since we know it (5 ha enclosed pop)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
# Check area: 
A <- (Xu-Xl)*(Yu-Yl)

# Observations
y <- as.matrix(snks)
colnames(y) <- NULL
rownames(y) <- NULL
# Uniquely marked individuals
nind <- nrow(y)

# Active/not active for when transects run (one less day surveyed for TL)
act <- as.matrix(eff2[,-1])
colnames(act) <- NULL
## number of occasions a camera was active
K <- as.vector(apply(act, 1, sum))  

## Date augmentation
M <- 300
y <-rbind(y,array(0,dim=c((M-nrow(y)),ncol(y))))

## Starting values for activity centers
## set inits of AC to mean location of individuals and then to some area within stat space for augmented
set.seed(10232020)
sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
for(i in 1:nind){
  sst[i,1] <- mean( X[y[i,]>0,1] )
  sst[i,2] <- mean( X[y[i,]>0,2] )
}

## NIMBLE model is nearly identical to BUGS
code <- nimbleCode({
  lam0~dunif(0,5)
  sigma~dunif(0,100) # informative prior = dgamma(274.69,7.27) based on 2015 SCR study
  psi~dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
    
    for(j in 1:J){
      d2[i,j] <- pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2)
      p[i,j] <- z[i]*lam0*exp(-(d2[i,j])/(2*sigma*sigma))
      
      y[i,j] ~ dpois(p[i,j]*K[i])
    }#j
  }#i
  N <- sum(z[1:M])
  D <- N/A
})#code


# MCMC settings
nc <- 3; nAdapt=1000; nb <- 1000; ni <- 3000+nb; nt <- 1

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

tosave <- as.matrix(samplesList)

save(tosave, file=paste("OptimSim_15traps",iter,".csv",sep=""))

