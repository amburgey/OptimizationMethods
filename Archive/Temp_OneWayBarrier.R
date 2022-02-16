### Simulation of one-way barrier surrounding study area
### Snakes can leave the study area but cannot enter

library(coda);library(runjags);library(phonTools);library(EnvStats);library(bda);library(sp);library(rgeos);library(raster)

#### VISUAL SURVEY REAL DATA RESULTS.----
load("Visual surveys/Results/NWFNVISALL_SCRpstarvisCATsizeCATdpois10GRIDnovsstALL.Rdata")#RE.Rdata")  # unified analysis of all datasets

## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
mcmc.object <- as.mcmc.list(out$samples[[1]])
out1 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[2]])
out2 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[3]])
out3 <- window(mcmc.object, start=1499)
outALLV <- combine.mcmc(list(out1,out2,out3), thin = 1)
niter(outALLV)
postV <- as.matrix(outALLV)

#for median and CI
post_sumV <- data.frame(
  mean = apply(postV, 2, function(x) mean(x)),
  med = apply(postV, 2, function(x) quantile(x, probs = 0.5, na.rm = T, names = F)),
  lower = apply(postV, 2, function(x) quantile(x, probs = 0.025, na.rm = T, names = F)),
  upper = apply(postV, 2, function(x) quantile(x, probs = 0.975, na.rm = T, names = F)),
  sd = apply(postV, 2, function(x) sd(x)))
post_sumV$variable <- row.names(post_sumV)

all_parsV <- colnames(postV)

# traceplots
# coda::traceplot(outALLV[,all_parsV[which(grepl('sigma', all_parsV))]])
# coda::traceplot(outALLV[,all_parsV[which(grepl('N', all_parsV))]])
# coda::traceplot(outALLV)


#### TRAPPING SURVEY REAL DATA RESULTS.----
load("Trapping/Results/NWFNVISALL_SCRpstarvisCATsizeCATdpois10GRIDnovsstALL.Rdata")

## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
mcmc.object <- as.mcmc.list(out$samples[[1]])
out1 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[2]])
out2 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[3]])
out3 <- window(mcmc.object, start=1499)
outALLT <- combine.mcmc(list(out1,out2,out3), thin = 1)
niter(outALLT)
postT <- as.matrix(outALLT)

#for median and CI
post_sumT <- data.frame(
  mean = apply(postT, 2, function(x) mean(x)),
  med = apply(postT, 2, function(x) quantile(x, probs = 0.5, na.rm = T, names = F)),
  lower = apply(postT, 2, function(x) quantile(x, probs = 0.025, na.rm = T, names = F)),
  upper = apply(postT, 2, function(x) quantile(x, probs = 0.975, na.rm = T, names = F)),
  sd = apply(postT, 2, function(x) sd(x)))
post_sumT$variable <- row.names(post_sumT)

all_parsT <- colnames(postT)

# # traceplots
# coda::traceplot(outALLT[,all_parsT[which(grepl('sigma', all_parsT))]])
# coda::traceplot(outALLT[,all_parsT[which(grepl('N', all_parsT))]])
# coda::traceplot(outALLT)


#### CREATE JOINT DISTRIBUTION FOR SIGMA.----
d1 <- rnorm(1000000, subset(post_sumV,variable=="sigma")$mean, subset(post_sumV,variable=="sigma")$sd)
d2 <- rnorm(1000000, subset(post_sumT,variable=="sigma")$mean, subset(post_sumT,variable=="sigma")$sd)
alld <- cbind(d1,d2)
sigma_mu <- mean(alld)
sigma_sd <- sd(alld)

dfV <- as.data.frame((postV[,c("sigma")]));colnames(dfV) <- c("sigma")
dfT <- as.data.frame(postT[,c("sigma")]);colnames(dfT) <- c("sigma")
hist(rbind(dfV,dfT)$sigma)
abline(v = c(mean(alld)), col="red", lwd=2)
abline(v = c(mean(alld)-sd(alld), mean(alld)+sd(alld)), col="red", lwd=2, lty=2)

combosigma <- hist(rnorm(1000,44.64,7.35))


#### GET POSTERIOR FOR P0 FOR EACH SNAKE SIZE FOR EACH MONITORING METHOD.----
p0MV <- list(sample(postV[,"p0[1]"]),sample(postV[,"p0[2]"]),sample(postV[,"p0[3]"]),sample(postV[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration
p0MT <- list(sample(postT[,"p0[1]"]),sample(postT[,"p0[2]"]),sample(postT[,"p0[3]"]),sample(postT[,"p0[4]"])) # randomize so encounter probabilities aren't from the same model iteration


#### GENERATE HOME RANGE AREA FROM SIGMA.----
#Scale parameter (sigma) can be interpreted as range size (for half-normal detection function, 95% of activity is within 2.45*sigma circle)
rangeS <- circular.r(p=0.99,detectfn = "HHN")*44.64
radius <- sqrt(rangeS/pi)
locAC <- G[s,] #locations of ACs on integration grid

## Get information to plot square integration grid cells
dfG <- as.data.frame(G)
A2 <- cbind(dfG$Xlocs,dfG$Ylocs,rep(a,nrow(dfG)))
colnames(A2) <- c("X","Y","Area")
B <- SpatialPoints(A2)
ww <- sqrt(B$Area)/2  ## Widths of buffers needed to produce desired areas    
#Find all 5 points of each square
pp <- list()
for(i in seq_along(B)) {
  pp[i] <- gBuffer(B[i], width=ww[i], quadsegs=1, capStyle="SQUARE")  #ignore warnings, deprecated approach but still works
}
PP <- do.call(bind, pp)

#Visualize
plot(PP)
# plot(PP, xlim = c(152.4605,197.3189), ylim = c(210.1056,229.1364))  #close up view
# plot(PP, xlim = c(-4.208262,38.5987), ylim = c(-24.31353, 0.247878))
plot(B, add=TRUE, pch=21, col="red", cex=0.5)
points(X, pch=18, col="orange")               #add locations of survey points
points(locAC, pch=16, col="black")            #add locations of ACs
#Plot home range around each AC (study area is 50,000m2 so home ranges are fairly small)
for(i in 1:nrow(locAC)){
  symbols(x=locAC[i,1], y=locAC[i,2], add=TRUE, inches = FALSE, circles = radius)
}
#The only snakes that will ever leave are those with their home ranges right on the edge of the study area (the home range size is never large enough to cross the barrier from farther away)
#Find individuals whose home ranges are outside the barrier based on centroid locations being right along the edge of the study area
circHalf <- as.numeric(row.names(subset(as.data.frame(locAC), Ylocs == G[1,2] | Ylocs == G[528,2] | Xlocs == G[1,1] | Xlocs == G[528,1])))
for(i in 1:length(circHalf)){
  symbols(x=locAC[circHalf[i],1], y=locAC[circHalf[i],2], add=TRUE, bg="#aaf0d1", inches = FALSE, circles = radius)
}
#Centroids of grid cells on the edge
# abline(v = G[1,1])
# abline(v = G[528,1])
# abline(h = G[1,2])
# abline(h = G[528,2])
#Anchor points of grid cells
# plot(pp[[1]], add=TRUE, col="red")
# plot(pp[[22]], add=TRUE, col="red")
# plot(pp[[507]], add=TRUE, col="red")
# plot(pp[[528]], add=TRUE, col="red")
#XY coordinates of farthest points of study area
# points(extent(pp[[1]])[1],extent(pp[[1]])[3], col="blue", cex=3, pch=21)
# points(extent(pp[[22]])[2],extent(pp[[22]])[3], col="blue", cex=3, pch=21)
# points(extent(pp[[507]])[1],extent(pp[[507]])[4], col="blue", cex=3, pch=21)
# points(extent(pp[[528]])[2],extent(pp[[528]])[4], col="blue", cex=3, pch=21)


#### PROBABILITY OF SNAKE LEAVING STUDY AREA.----
#Find intersection with home ranges and boundaries of CP
#https://www.engineersedge.com/math/circular_segment_equation_and_calculator__13796.htm
#For edge of CP, integration grid points pt1 = G[1,1:2], pt2 = G[22,1:2], pt3 = G[507,1:2], pt4 = G[528,1:2]
#Line coordinates (x1,y1) and (x2,y2); horizontal line is just y intercept, for vertical line is just x intercept
#The home range will be the same size for all individuals in a simulation so only need to calculate this once per sigma value
#d is the distance between a home range centroid y-coordinate and the equation of a horizontal study barrier (y) cutting through the home range
d <- locAC[circHalf[1],2] - extent(pp[[1]])[3]
#h is the distance between d (the barrier cutting through the home range) and the circular edge of the home range
h <- radius - d
#c is the length of the line segment within the home range
c <- sqrt((radius-(h/2))*8*h)
#theta is the interior angle of the segment within the circle
theta <- 2*atan(c/(2*d))

#The proportion of home range outside of the barrier to inside dictates the probability of the snake leaving
#https://www.mathsisfun.com/geometry/circle-sector-segment.html
#Area of circle outside of circle (when theta is in radians)
areaseg <- ((theta-sin(theta))/2) * (radius)^2
probleave <- areaseg/rangeS

#Do Bernoulli trial of whether snake stays or snake leaves based on this probability
#If snake leaves, it can never come back in
stay <- array(NA, dim=c(N,K,nsims))
for(z in 1:nsims){
  for(k in 1:K){                          #for each occasion sampled
    for(n in 1:N){                        #for each snake
      if(is.na(stay[n,K,z])){               #conditional on the snake still being in the study area
        if(any(n == circHalf) == TRUE){   #conditional on the snake's home range overlapping the fence
          stay[n,k,z] <- rbinom(1,1,(1-probleave))
          if(stay[n,k,z] == 0){
            stay[n,k:K,z] <- 0
            next                          #If the snake leaves, fill in zeroes and skip for rest of occasions
          }
        } else {
          stay[n,k,z] <- 1
        }
      } else {
        next
      }
    }
  }
}


#### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----

## VISUAL SURVEYS - Generate observations
yTrueVIS <- array(NA,dim=c(N,J,K,nsims))
yTrueVIS2 <- array(NA,dim=c(N,J,nsims))
p0V <- vector()

for(z in 1:nsims){
  sigmaV <- rnorm(1,44.64,7.35)   ## pull value from combined sigma distribution
  alpha1V <- 1/(2*sigmaV*sigmaV)
  for(l in 1:length(Ngroup)){
    p0V[l] <- p0MV[[l]][sample(0:length(p0MV[[l]]), 1)]   ## pull value from posterior
  }
  print(p0V)
  
  pmatV <- p0V[Nsnsz]*exp(-alpha1V*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
  
  for(n in 1:N){  ## for each snake
    for(k in 1:K){  ## for each occasion sampled
      yTrueVIS[n,,k,z] <- rpois(J,pmatV[n,])*stay[n,k,z]  # observations of each snake at each survey location based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
    }
    for(j in 1:J){
      ## Collapse observations across all occasions to give a count of observations per survey location for each simulation
      yTrueVIS2[n,j,z] <- sum(yTrueVIS[n,j,,z])
    }
  }
}


## yTrueVIS includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
capturedV <- list()
yarrV <- list()
## Bind sizes of snake to encounter histories
yTrueVIS3 <- abind::abind(yTrueVIS2,array(Nsnsz, replace(dim(yTrueVIS2),2,1)), along=2)
for(z in 1:nsims){
  capturedV[[z]] <- which(apply(yTrueVIS3[,1:J,z],1,sum)>0)  # snakes that were observed at least once
  yarrV[[z]] <- yTrueVIS3[capturedV[[z]],,z]  # subset to observed snakes
  write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
}



#### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----

## TRAPPING - Generate observations
yTrueTRAP <- array(NA,dim=c(N,J,K,nsims))
yTrueTRAP2 <- array(NA,dim=c(N,J,nsims))
p0T <- vector()

for(z in 1:nsims){  
  sigmaT <- rnorm(1,44.64,7.35)   ## pull value from combined sigma distribution
  alpha1T <- 1/(2*sigmaT*sigmaT)
  for(l in 1:length(Ngroup)){
    p0T[l] <- p0MT[[l]][sample(0:length(p0MT[[l]]), 1)]   ## pull value from posterior
  }
  print(p0T)
  
  pmatT <- p0T[Nsnsz]*exp(-alpha1T*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all locations
  
  for(n in 1:N){
    for(k in 1:K){  ## for each occasion sampled
      yTrueTRAP[n,,k,z] <- rbinom(J,1,pmatT[n,])*stay[n,k,z]  # observations of each snake at each trap based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
    }
    for(j in 1:J){
      ## Collapse observations across all occasions to give a count of observations per trap for each simulation
      yTrueTRAP2[n,j,z] <- sum(yTrueTRAP[n,j,,z])
    }
  }
}

## yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
capturedT <- list()
yarrT <- list()
## Bind sizes of snake to encounter histories
yTrueTRAP3 <- abind::abind(yTrueTRAP2,array(Nsnsz, replace(dim(yTrueTRAP2),2,1)), along=2)
for(z in 1:nsims){
  capturedT[[z]] <- which(apply(yTrueTRAP3[,1:J,z],1,sum)>0)  # snakes that were observed at least once
  yarrT[[z]] <- yTrueTRAP3[capturedT[[z]],,z]  # subset to observed snakes
  write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
}



#### SIMULATE OBSERVATIONS OF SNAKES BASED ON THIS DESIGN ----
## VISUAL AND TRAPPING SURVEYS - Generate observations
yTrueVIS <- array(NA,dim=c(N,J1,K,nsims))
yTrueTRAP <- array(NA,dim=c(N,J2,K,nsims))
yTrueVIS2 <- array(NA,dim=c(N,J1,nsims))
yTrueTRAP2 <- array(NA,dim=c(N,J2,nsims))
p0V <- vector()
p0T <- vector()

for(z in 1:nsims){  
  sigma <- rnorm(1,44.64,7.35)   ## pull value from combined sigma distribution
  alpha1 <- 1/(2*sigma*sigma)
  
  for(l in 1:length(Ngroup)){
    p0V[l] <- p0MV[[l]][sample(0:length(p0MV[[l]]), 1)]   ## pull value from posterior
    p0T[l] <- p0MT[[l]][sample(0:length(p0MT[[l]]), 1)]   ## pull value from posterior
  }
  print(p0V)
  print(p0T)
  
  pmatV <- p0V*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
  pmatT <- p0T*exp(-alpha1*Gdist[s,]*Gdist[s,])  # encounter probabilities of all snakes (based on their size and activity centers) at all location
  
  for(n in 1:N){
    for(k in 1:K){
      yTrueTRAP[n,,k,z] <- rbinom(J2,1,pmatT[n,newX2])*stay[n,k,z]    # observations of each snake at each trap based on encounter probability and effort multiplied by Bernoulli trial of if they leave the study area
      if(any(yTrueTRAP[n,,k,z] >= 1) == TRUE){
        yTrueVIS[n,,k,z] <- 0
        next
        } else {
        yTrueVIS[n,,k,z] <- rpois(J1,pmatV[n,newX1])*stay[n,k,z]  # observations of each snake at each survey location based on encounter probability and effort conditional on if snake was captured in a trap that evening multiplied by Bernoulli trial of if they leave the study area
      }
    }
    for(j in 1:J2){
      ## Collapse observations across all occasions to give a count of observations per trap for each simulation
      yTrueTRAP2[n,j,z] <- sum(yTrueTRAP[n,j,,z])
    }
    for(j in 1:J1){
      yTrueVIS2[n,j,z] <- sum(yTrueVIS[n,j,,z])
    }
  }
}

# # CHECK TO MAKE SURE NO OCCASION HAS A SNAKE TRAPPED AND VISUALLY SEEN
# for(z in 1:nsims){
#   for(k in 1:K){
#     ifelse(any((rowSums(yTrueTRAP[,,k,z]) & rowSums(yTrueVIS[,,k,z])) == TRUE) == TRUE , print("oh no"), print("we good"))
#   }
# }


## yTrueVIS and yTrueTRAP includes a row for every snake even if that snake was never observed. We need to remove these snakes to mimic real data.
capturedV <- list()
capturedT <- list()
yarrV <- list()
yarrT <- list()
## Bind sizes of snake to encounter histories
yTrueVIS3 <- abind::abind(yTrueVIS2,array(Nsnsz, replace(dim(yTrueVIS2),2,1)), along=2)
yTrueTRAP3 <- abind::abind(yTrueTRAP2,array(Nsnsz, replace(dim(yTrueTRAP2),2,1)), along=2)
## Bind unique snake identity to encounter histories
yTrueVIS3 <- abind::abind(yTrueVIS3,array(seq(1:N), replace(dim(yTrueVIS3),2,1)), along=2)
yTrueTRAP3 <- abind::abind(yTrueTRAP3,array(seq(1:N), replace(dim(yTrueTRAP3),2,1)), along=2)

for(z in 1:nsims){
  capturedV[[z]] <- which(apply(yTrueVIS3[,1:J1,z],1,sum)>0)  # snakes that were observed at least once
  capturedT[[z]] <- which(apply(yTrueTRAP3[,1:J2,z],1,sum)>0)  # snakes that were observed at least once
  yarrV[[z]] <- yTrueVIS3[capturedV[[z]],,z]  # subset to observed snakes
  yarrT[[z]] <- yTrueTRAP3[capturedT[[z]],,z]  # subset to observed snakes
  write.csv(yarrV[[z]], file = paste("Simulations/simDat/",type,"VIS",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
  write.csv(yarrT[[z]], file = paste("Simulations/simDat/",type,"TRAP",N,dens,K,stde,z,".csv",sep=""))  # write to file to keep observations
}




##REJECT CODE
#Plot integration grid and activity centers
# xyG <- as.data.frame(G)
# xyAC <- as.data.frame(locAC)
# coordinates(xyG) =~ Xlocs + Ylocs
# coordinates(xyAC) =~ Xlocs + Ylocs
# plot(xyG, pch = 1)
# plot(xyAC, add=TRUE, col="red", pch=21)
# 
# #Find top left and bottom right corners of bounding box around each AC
# bb <- as.data.frame(matrix(NA, nrow = nrow(locAC), ncol = 5))
# for(i in 1:nrow(locAC)){
#   bb[i,1:2] <- c(locAC[i,1] - radius, locAC[i,2] + radius)
#   bb[i,3:4] <- c(locAC[i,1] + radius, locAC[i,2] - radius)
#   bb[i,5] <- i
# }
#Plot bounding box edges
# xyTL <- as.data.frame(bb[,1:2])
# xyBR <- as.data.frame(bb[,3:4])
# coordinates(xyTL) =~ V1 + V2
# coordinates(xyBR) =~ V3 + V4
# plot(xyTL, add=TRUE, col="blue", pch=18)
# plot(xyBR, add=TRUE, col="blue", pch=18)


# overlap <- vector()
# for(i in 1:nrow(locAC)){
#   center <- locAC[i,]
#   print(i)
#   #(x−point_x)^2+(y−point_y)^2=rangeS
#   ## Test lower boundary
#   yL = G[1,2] #Horizontal, lower boundary
#   xL1 <- locAC[i,1] + sqrt(rangeS-(yL-locAC[i,2])^2)
#   xL2 <- locAC[i,1] - sqrt(rangeS-(yL-locAC[i,2])^2)
#   if(is.na(xL1) | is.na(xL2)){
#     print("No intersection with lower")
#   } else {
#     print("Intersection with lower")
#   }
#   ## Test upper boundary
#   yU = G[528,2] #Horiztonal, upper boundary
#   xU1 <- locAC[i,1] + sqrt(rangeS-(yU-locAC[i,2])^2)
#   xU2 <- locAC[i,1] - sqrt(rangeS-(yU-locAC[i,2])^2)
#   if(is.na(xU1) | is.na(xU2)){
#     print("No intersection with upper")
#   } else {
#     print("Intersection with upper")
#   }
#   ## Test left boundary
#   xL = G[1,1] #Vertical, left boundary
#   yL1 <- center[2] + sqrt((radius)^2-(xL-center[1])^2)
#   yL2 <- center[2] - sqrt((radius)^2-(xL-center[1])^2)
#   # yL1 <- locAC[i,1] + sqrt(rangeS-(xL-locAC[i,2])^2)
#   # yL2 <- locAC[i,1] - sqrt(rangeS-(xL-locAC[i,2])^2)
#   if(is.na(yL1) | is.na(yL2)){
#     print("No intersection with left")
#   } else {
#     print("Intersection with left")
#   }
#   ## Test right boundary
#   xR = G[528,1] #Vertical, right boundary
#   yR1 <- center[2] + sqrt((radius)^2-(xR-center[1])^2)
#   yR2 <- center[2] - sqrt((radius)^2-(xR-center[1])^2)
#   # yR1 <- locAC[i,1] + sqrt(rangeS-(xR-locAC[i,2])^2)
#   # yR2 <- locAC[i,1] - sqrt(rangeS-(xR-locAC[i,2])^2)
#   if(is.na(yR1) | is.na(yR2)){
#     print("No intersection with right")
#   } else {
#     print("Intersection with right")
#   }
#   
#   if(!is.na(xL1|xL2|xU1|xU2|yL1|yL2|yR1|yR2) == TRUE){
#     overlap[i] <- i
#   }
# }

