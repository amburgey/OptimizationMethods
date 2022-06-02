### Brown treesnake simulation - Visualize different transect sampling designs (visual surveys [VIS], trapping surveys [TRAP], and combined surveys [VISTRAP])

rm(list=ls())

## Required libraries
library(jagsUI);library(secr);library(ggspatial);library(patchwork);library(ggplot2)
## Functions for simulating data
source("Simulations/Scripts/FunctionsForSimulation_ClosedAndOneWayBarrier.R")



#### SCENARIO ONE.----

## Type of sampling
type <- c("VIS")  # Same locations for TRAP so just do one to plot pattern
nmeth <- 1

## Full of ONE method [27 transects, all]
stde <- c("full")
samp <- c(1:351)

## True number of snakes (low density)
N <- 60


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area based on rough dimensions of CP
sdeets <- areatype(totlocs = totlocs)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (i.e., if surveying less than the full design), subset locations to only those monitored
X <- X[samp,]


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points using the general constraints of Xl, Xu, Yl, Yu
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points


#### SIMULATION INFO.----

## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE) #pull value from integration grid locations


#### PLOT ONE.----

## Either VIS or TRAP
BOTHfull <- ggplot() + 
  geom_point(data = as.data.frame(G), aes(x = Xlocs, y = Ylocs), cex = 0.8, color = "darkgrey") + 
  geom_point(data = as.data.frame(X), aes(x = x, y = y), cex = 3, pch = 21, fill = "darkgrey", color = "black") +
  scale_y_continuous(breaks = seq(min(G[,2])-5, max(G[,2])+5, by = 10)) +
  scale_x_continuous(breaks = seq(min(G[,1])-5, max(G[,1])+5, by = 10)) +
  theme(panel.grid = element_line(color = "lightgrey",
                                  size = 0.75,
                                  linetype = 1),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank()) +
  geom_line(aes(x=c(min(G[,1])-5,min(G[,1])-5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(max(G[,1])+5,max(G[,1])+5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(min(G[,2])-5, min(G[,2])-5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(max(G[,2])+5, max(G[,2])+5))) +
  ylab("Y locations")



#### SCENARIO TWO.----

## Type of sampling
type <- c("VIS")  # Same locations for TRAP so just do one to plot pattern
nmeth <- 1

## Half of ONE method [14 transects, every other]
stde <- c("half")
samp <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)

## True number of snakes (low density)
N <- 60


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area based on rough dimensions of CP
sdeets <- areatype(totlocs = totlocs)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
X <- X[samp,]


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points using the general constraints of Xl, Xu, Yl, Yu
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points


#### SIMULATION INFO.----

## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE) #pull value from integration grid locations


#### PLOT TWO.----
  
## Either VIS or TRAP
BOTHhalf <- ggplot() + 
  geom_point(data = as.data.frame(G), aes(x = Xlocs, y = Ylocs), cex = 0.8, color = "darkgrey") + 
  geom_point(data = as.data.frame(X), aes(x = x, y = y), cex = 3, pch = 21, fill = "darkgrey", color = "black") +
  scale_y_continuous(breaks = seq(min(G[,2])-5, max(G[,2])+5, by = 10)) +
  scale_x_continuous(breaks = seq(min(G[,1])-5, max(G[,1])+5, by = 10)) +
  theme(panel.grid = element_line(color = "lightgrey",
                                  size = 0.75,
                                  linetype = 1),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_line(aes(x=c(min(G[,1])-5,min(G[,1])-5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(max(G[,1])+5,max(G[,1])+5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(min(G[,2])-5, min(G[,2])-5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(max(G[,2])+5, max(G[,2])+5)))



#### SCENARIO THREE.----

## Type of sampling
type <- c("VISTRAP")
nmeth <- 2

## Half of TWO methods (e.g., VISTRAP)
stde <- c("halfhalf")
samp1 <- c(1:13,27:39,53:65,79:91,105:117,131:143,157:169,183:195,209:221,235:247,261:273,287:299,313:325,339:351)
samp2 <- c(14:26,40:52,66:78,92:104,118:130,144:156,170:182,196:208,222:234,248:260,274:286,300:312,326:338)

## True number of snakes (low density)
N <- 60


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area based on rough dimensions of CP
sdeets <- areatype(totlocs = totlocs)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (i.e., if surveying less than the full design), subset locations to only those monitored
X1 <- X[samp1,]  ## VIS
X2 <- X[samp2,]  ## TRAP
X <- X[sort(c(samp1,samp2)),]


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points using the general constraints of Xl, Xu, Yl, Yu
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points


#### SIMULATION INFO.----

## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE) #pull value from integration grid locations


#### PLOT THREE.----

## VISTRAP
VISTRAPfull <- ggplot() + 
  geom_point(data = as.data.frame(G), aes(x = Xlocs, y = Ylocs), cex = 0.8, color = "darkgrey") + 
  geom_point(data = as.data.frame(X1), aes(x = x, y = y), cex = 3, pch = 21, fill = "#de993e", color = "black") +
  geom_point(data = as.data.frame(X2), aes(x = x, y = y), cex = 3, pch = 21, fill = "#84ad32", color = "black") +
  scale_y_continuous(breaks = seq(min(G[,2])-5, max(G[,2])+5, by = 10)) +
  scale_x_continuous(breaks = seq(min(G[,1])-5, max(G[,1])+5, by = 10)) +
  theme(panel.grid = element_line(color = "lightgrey",
                                  size = 0.75,
                                  linetype = 1),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_line(aes(x=c(min(G[,1])-5,min(G[,1])-5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(max(G[,1])+5,max(G[,1])+5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(min(G[,2])-5, min(G[,2])-5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(max(G[,2])+5, max(G[,2])+5))) +
  ylab("Y locations") + xlab("X locations")



#### SCENARIO FOUR.----

## Type of sampling
type <- c("VISTRAP")
nmeth <- 2

## Third of TWO methods (e.g., VISTRAP)
stde <- c("thirdthird")
samp1 <- c(14:26,53:65,92:104,131:143,170:182,209:221,248:260,287:299,326:338)
samp2 <- c(27:39,66:78,105:117,144:156,183:195,222:234,261:273,300:312,339:351)

## True number of snakes (low density)
N <- 60


#### STUDY AREA INFORMATION.----

## Study area based on the Closed Population (CP; 5-ha) with transects spaced every 8-m from each other and points on the transects spaced every 16-m
## Create location points and get coordinates
totlocs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
## Get dimensions of the study area based on rough dimensions of CP
sdeets <- areatype(totlocs = totlocs)


#### SURVEY INFORMATION.----

## Create matrix of sampling location options
X <- as.matrix(totlocs)
## If applicable (i.e., if surveying less than the full design), subset locations to only those monitored
X1 <- X[samp1,]  ## VIS
X2 <- X[samp2,]  ## TRAP
X <- X[sort(c(samp1,samp2)),]


#### INTEGRATION GRID.----

## Spacing of grid cells
Ggrid <- 10
## Find XY locations of all integration grid cell points using the general constraints of Xl, Xu, Yl, Yu
Xlocs <- rep(seq(sdeets$Xl, sdeets$Xu, Ggrid), times = length(seq(sdeets$Yl, sdeets$Yu, Ggrid)))
Ylocs <- rep(seq(sdeets$Yl, sdeets$Yu, Ggrid), each = length(seq(sdeets$Xl, sdeets$Xu, Ggrid)))
G <- cbind(Xlocs, Ylocs)
Gpts <- dim(G)[1]                         #number of integration points


#### SIMULATION INFO.----

## True snake activity centers (AC)
s <- sample(1:Gpts,N,replace=TRUE) #pull value from integration grid locations


#### PLOT FOUR.----

## VISTRAP
VISTRAPhalf <- ggplot() + 
  geom_point(data = as.data.frame(G), aes(x = Xlocs, y = Ylocs), cex = 0.8, color = "darkgrey") + 
  geom_point(data = as.data.frame(X1), aes(x = x, y = y), cex = 3, pch = 21, fill = "#de993e", color = "black") +
  geom_point(data = as.data.frame(X2), aes(x = x, y = y), cex = 3, pch = 21, fill = "#84ad32", color = "black") +
  scale_y_continuous(breaks = seq(min(G[,2])-5, max(G[,2])+5, by = 10)) +
  scale_x_continuous(breaks = seq(min(G[,1])-5, max(G[,1])+5, by = 10)) +
  theme(panel.grid = element_line(color = "lightgrey",
                                  size = 0.75,
                                  linetype = 1),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank()) +
  geom_line(aes(x=c(min(G[,1])-5,min(G[,1])-5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(max(G[,1])+5,max(G[,1])+5),y=c(min(G[,2])-5, max(G[,2])+5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(min(G[,2])-5, min(G[,2])-5))) +
  geom_line(aes(x=c(min(G[,1])-5,max(G[,1])+5),y=c(max(G[,2])+5, max(G[,2])+5))) +
  xlab("X locations")


#### Create combined panel figure of different simulation sampling designs

png(file="Simulations/Figures/SimDesignScenarios.png",width=12,height=9,units="in",res=600)
((BOTHfull + BOTHhalf) / (VISTRAPfull + VISTRAPhalf))
dev.off()



#### PLOT LEGEND TO COMBINE WITH ABOVE FIGURE.----

X2 <- as.data.frame(X)
X2$Type <- c(rep("VIS or TRAP",78),rep("VIS",78),rep("TRAP",78))

Legend <- ggplot(data = X2, aes(x = x, y = y, fill = Type)) +
  geom_point(cex = 3, pch = 21, color = "black") +
  scale_fill_manual(values = c("darkgrey","#de993e","#84ad32"), breaks = c("VIS or TRAP","VIS","TRAP")) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 15), legend.position = "top")

png(file="Simulations/Figures/SimDesignScenariosLEGEND.png",width=12,height=9,units="in",res=600)
Legend
dev.off()

