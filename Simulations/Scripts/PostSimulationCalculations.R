### Post simulation calculations
### Calculate RMSE, CV, coverage

rm(list=ls())


library(HDInterval); library(dplyr); library(stringr)


nsims <- 100

#### Scenario.----

scen <- c(
  "closedVIS120small30full","closedVIS120small30half",
  "closedVIS120small14full","closedVIS120small14half",
  "closedVIS60small30full","closedVIS60small30half",
  "closedVIS60small14full","closedVIS60small14half",
  "closedVIS120large30full","closedVIS120large30half",
  "closedVIS120large14full","closedVIS120large14half",
#   "closedVIS60large30full","closedVIS60large30half",
#   "closedVIS60large14full","closedVIS60large14half",
  "closedTRAP120small30full","closedTRAP120small30half",
  "closedTRAP120small14full",
# "closedTRAP120small14half",
#   "closedTRAP60small30full","closedTRAP60small30half",
#   "closedTRAP60small14full","closedTRAP60small14half",
  "closedTRAP120large30full",
# "closedTRAP120large30half",
  "closedTRAP120large14full",
# "closedTRAP120large14half",
#   "closedTRAP60large30full","closedTRAP60large30half",
#   "closedTRAP60large14full","closedTRAP60large14half",
#   "closedVISTRAP120small30halfhalf","closedVISTRAP120small30thirdthird",
#   "closedVISTRAP120small14halfhalf","closedVISTRAP120small14thirdthird",
#   "closedVISTRAP60small30halfhalf","closedVISTRAP60small30thirdthird",
#   "closedVISTRAP60small14halfhalf","closedVISTRAP60small14thirdthird",
#   "closedVISTRAP120large30halfhalf","closedVISTRAP120large30thirdthird",
#   "closedVISTRAP120large14halfhalf","closedVISTRAP120large14thirdthird",
#   "closedVISTRAP60large30halfhalf","closedVISTRAP60large30thirdthird",
#   "closedVISTRAP60large14halfhalf","closedVISTRAP60large14thirdthird",
"onewayVIS120small30full","onewayVIS120small30half",
"onewayVIS120small14full",
# "onewayVIS120small14half",
"onewayVIS60small30full","onewayVIS60small30half",
"onewayVIS60small14full","onewayVIS60small14half",
"onewayVIS120large30full","onewayVIS120large30half",
"onewayVIS120large14full","onewayVIS120large14half"
#   "onewayVIS60large30full","onewayVIS60large30half",
#   "onewayVIS60large14full","onewayVIS60large14half",
# "onewayTRAP120small30full","onewayTRAP120small30half",
# "onewayTRAP120small14full","onewayTRAP120small14half",
#   "onewayTRAP60small30full","onewayTRAP60small30half",
#   "onewayTRAP60small14full","onewayTRAP60small14half",
# "onewayTRAP120large30full","onewayTRAP120large30half",
# "onewayTRAP120large14full","onewayTRAP120large14half",
#   "onewayTRAP60large30full","onewayTRAP60large30half",
#   "onewayTRAP60large14full","onewayTRAP60large14half",
#   "onewayVISTRAP120small30halfhalf","onewayVISTRAP120small30thirdthird",
#   "onewayVISTRAP120small14halfhalf","onewayVISTRAP120small14thirdthird",
#   "onewayVISTRAP60small30halfhalf","onewayVISTRAP60small30thirdthird",
#   "onewayVISTRAP60small14halfhalf","onewayVISTRAP60small14thirdthird",
#   "onewayVISTRAP120large30halfhalf","onewayVISTRAP120large30thirdthird",
#   "onewayVISTRAP120large14halfhalf","onewayVISTRAP120large14thirdthird",
#   "onewayVISTRAP60large30halfhalf","onewayVISTRAP60large30thirdthird",
#   "onewayVISTRAP60large14halfhalf","onewayVISTRAP60large14thirdthird"
)

type <- length(scen)

#### Study area size.----

A <- 52800
# for(t in 1:length(scen)){
#   if(str_detect(scen[t], "closed")){  # if study area is closed, use this area
#     A <- 52800
#   } else{                             # if study area is open, use this area
#     A <- 52800   ### UPDATE WITH FINAL OPEN AREA FROM SIMULATION
#   }
# }

#### Function to find mode for continuous data.----

mode <- function(data) {
  # Function for mode estimation of a continuous variable
  # Kernel density estimation by Ted Harding & Douglas Bates (found on RSiteSearch)	
  
  x<-data
  lim.inf=min(x)-1; lim.sup=max(x)+1
  
  # hist(x,freq=FALSE,breaks=seq(lim.inf,lim.sup,0.2))
  s<-density(x,from=lim.inf,to=lim.sup,bw=0.2)
  n<-length(s$y)
  v1<-s$y[1:(n-2)];
  v2<-s$y[2:(n-1)];
  v3<-s$y[3:n]
  ix<-1+which((v1<v2)&(v2>v3))
  
  # lines(s$x,s$y,col="red")
  # points(s$x[ix],s$y[ix],col="blue")
  
  md <- s$x[which(s$y==max(s$y))] 
  
  md
}

#### Create dataframe to populate with simulation results calculations.----

simSumm <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("Sim","Nmn","Nmod","bias", "Q2.5", "Q97.5","CV","Rhat","RMSE"))

#### Pull info from model results.----

count <- 1  ## incremental key to store data in correct dataframe position

for(t in 17){ #1:length(scen)
  for(i in 1:100){ #1:nsims
  
  load(file=paste("D:/Results/RESULTS_",scen[t],i,".Rdata",sep=""))
      
  print(t)
  print(i)
      
  simSumm[count,1] <- scen[t]                                     # Scenario
  simSumm[count,2] <- out$mean$N                                  # Mean Abundance
  simSumm[count,3] <- mode(out$sims.list$N)                       # Mode Abundance
  if(str_detect(scen[t], "120")){                                 # If scenario used 120 snakes in simulation
    simSumm[count,4] <- (1/120)*(out$mean$N - 120)                # Calculate bias
  } else{                                                         # If scenario used 60 snakes in simulation
    simSumm[count,4] <- (1/60)*(out$mean$N - 60)                  # Calculate bias
  }
  simSumm[count,5] <- hdi(out, credMass = 0.95)[,"N"][1]          # 2.5% credible interval of N
  simSumm[count,6] <- hdi(out, credMass = 0.95)[,"N"][2]          # 97.5% credible interval of N
  simSumm[count,7] <- out$sd$N/out$mean$N                         # Coefficient of variation
  simSumm[count,8] <- out$Rhat$N                                  # Rhat for N
  simSumm[count,9] <- mean(sqrt(var(out$sims.list$N) + (out$sims.list$N - 120)^2))  # RMSE
      
  count <- count + 1
      
  }
}


s#### Check if any simulation didn't converge.----
for(i in 1:nrow(simSumm)){
  try(if(simSumm$Rhat[i] > 1.1) stop("did not converge"))
}

#### Write individual results to file.----
write.csv(simSumm, file = "Simulations/Results/SimIndividualResults.csv")

############ TEMPORARY #################
## remove subset for right now that didn't converge
simSumm <- subset(simSumm, Rhat < 1.1)

#### Calculate metrics.----

## Calculate mean snake abundance by each monitoring scenario
snakenums <- aggregate(simSumm$Nmn, list(simSumm$Sim), mean)
snakenums <- snakenums[order(snakenums[,1]),]

## Calculate mode snake abundance by each monitoring scenario
snakemod <- aggregate(simSumm$Nmod, list(simSumm$Sim), mean)
snakemod <- snakemod[order(snakemod[,1]),]

## Calculate coverage (true value [120 or 60] contained within 95% credible intervals [1] or not [0])
simSumm$Cov <- NA
for(t in 1:nrow(simSumm)){
  if(str_detect(simSumm$Sim[t], "120")){                                                   # If scenario used 120 snakes in simulation
    simSumm[t,ncol(simSumm)] <- ifelse(simSumm[t,5] <= 120 & simSumm[t,6] >= 120, 1, 0)     # Calculate coverage
  } else{                                                                                  # If scenario used 60 snakes in simulation
    simSumm[t,ncol(simSumm)] <- ifelse(simSumm[t,5] <= 60 & simSumm[t,6] >= 60, 1, 0)      # Calculate coverage
  }
}
## Percent of time real value included within credible intervals (assuming 100 simulations used)
snakecov <- aggregate(simSumm$Cov, list(simSumm$Sim), sum)
snakecov <- snakecov[order(snakecov[,1]),]

## Calculate RMSE of snake abundance by each monitoring scenario
# mse <- as.data.frame(matrix(NA, nrow = length(scen), ncol = 2, dimnames = list(seq(1:length(scen)),c("Group.1","x"))))
# rmse <- as.data.frame(matrix(NA, nrow = length(scen), ncol = 2, dimnames = list(seq(1:length(scen)),c("Group.1","x"))))
# for(t in 1:nrow(simSumm)){
#   if(str_detect(simSumm$Sim[t], "120")){                                  # If scenario used 120 snakes in simulation
#     mse[t,1] <- simSumm$Sim[t]                                            # Scenario type
#     mse[t,2] <- apply((as.data.frame(snakenums[t,2]) - 120)^2 ,1,mean)    # mean squared error
#     rmse[t,1] <- simSumm$Sim[t]                                           # Scenario type
#     rmse[t,2] <- sqrt(mse[t,2])                                           # root mean squared error
#   } else{                                                                 # If scenario used 60 snakes in simulation
#     mse[t,1] <- simSumm$Sim[t]                                            # Scenario type
#     mse[t,2] <- apply((as.data.frame(snakenums[t,2]) - 60)^2 ,1,mean)     # mean squared error
#     rmse[t,1] <- simSumm$Sim[t]                                           # Scenario type
#     rmse[t,2] <- sqrt(mse[t,2])                                           # root mean squared error
#   }
# }
snakermse <- aggregate(simSumm$RMSE, list(simSumm$Sim), mean)
snakermse <- snakermse[order(snakermse[,1]),]

## Calculate % coefficient of variation
snakeCV <- aggregate(simSumm$CV, list(simSumm$Sim), mean)
snakeCV$x <- snakeCV$x*100
snakeCV <- snakeCV[order(snakeCV[,1]),]

## Calculate % bias
snakebias <- cbind(aggregate(simSumm$bias, list(simSumm$Sim), mean))
snakebias$x <- snakebias$x*100
snakebias <- snakebias[order(snakebias[,1]),]

#### Create results table and write to file.----
all <- cbind(snakenums,snakemod[,2],snakecov[,2],snakeCV[,2],snakebias[,2],snakermse[,2])
colnames(all) <- c("Sim","Nmn","Nmod","%Coverage","%CV","%Bias","RMSE")

write.csv(all,"Simulations/Results/SimResultsSummary.csv")
