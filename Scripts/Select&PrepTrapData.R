#### Nov 2, 2020; S. Amburgey ####
## Read in data (survey info and captures) from numerous trapping projects
## Select what projects to use
## Format for SCR analysis

library(lubridate); library(reshape2); library(dplyr); library(sp); library(ggplot2); library(rgdal); library(raster); library(sf); library(tidyverse); library(plyr); library(plotKML)

# source("RenamingGrid.R") # don't need as no traps set on edges or perimeter



###### SURVEY DATA ######
allsurv <- read.csv("TRAPinfo.csv")

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("ASOI","NWFN","HMUI","HMUR","HMU","RTSI"))
## remove for right now: "NWFO"

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("TOX DROP TRAP","POST KB TRAP 1","POST KB TRAP 2", "NWFN TRAP 1","NWFN TRAP 2 LINVIS","NWFN TRAP 3","NWFN TRAP 4 LCM","NWFN VISTRAP TRAP","PRE BT1 TRAP","POST BT2 TRAP","GNWR SSP MOUSE", "SWIFT TRAP ARRAY"))
## removed for right now: "NWFO TRAP 1","NWFO TRAP 2","NWFO TRAP 3","NWFO TRAP 4","GNWR SSP BIRD","HMU TOX DROP 2 TRAP","KB TRAP"

## Specify date and separate year
subsurv$Date <- dmy(subsurv$Date)
subsurv$Year <- year(subsurv$Date)

#### CLEAN SURVEY DATA ####

subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "51100" & subsurv$EFFORTID == "9638"),] # Survey marked as occurring on RIM (transect) with location blank


###### CAPTURE DATA ######
allcap <- read.csv("TRAPcaptures.csv")

## Study Areas deemed suitable
subcap <- subset(allcap, SITE %in% c("ASOI","NWFN","HMUI","HMUR","HMU","RTSI"))
#remove for right now: "NWFO"

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subcap <- subset(subcap, PROJECTCODE %in% c("TOX DROP TRAP","POST KB TRAP 1","POST KB TRAP 2", "NWFN TRAP 1","NWFN TRAP 2 LINVIS","NWFN TRAP 3","NWFN TRAP 4 LCM","NWFN VISTRAP TRAP","PRE BT1 TRAP","POST BT2 TRAP","GNWR SSP MOUSE", "SWIFT TRAP ARRAY"))
#remove for right now: "NWFO TRAP 1","NWFO TRAP 2","NWFO TRAP 3","NWFO TRAP 4","GNWR SSP BIRD","HMU TOX DROP 2 TRAP","KB TRAP"

## Specify date and separate year
subcap$Date <- dmy(subcap$Date)
subcap$Year <- year(subcap$Date)

#### CLEAN CAPTURE DATA ####
subcap <- subcap[!(subcap$SITE == "HMUR" & subcap$TRANSID == "76886" & subcap$EFFORTID == "15039"),] # Snake survey missing, but the rest of the captures for this animal are part of HMU KNOWN FATE (including its death)
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "57216" & subcap$EFFORTID == "12309"),] # Snake found incidentally dead, died during telemtry project likely due to force feeding
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "29784" & subcap$EFFORTID == "5310"),]  # removed from unknown transect
# subcap <- subcap[!(subcap$SITE == "NWFO" & subcap$TRANSID == "29783" & subcap$EFFORTID == "5151"),]  # removed from unknown transect
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "55997" & subcap$EFFORTID == "12223"),] # Snake found incidentally
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "56296" & subcap$EFFORTID == "12295"),] # Snake left in box, incidental death
subcap[subcap$EFFORTID == "11794" & subcap$TRANSID == "53624" & subcap$LOCATION == "F05","LOCATION"] <- "F5" # capture listed at location F05 when others listed as F5
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "51100" & subcap$EFFORTID == "9638"),] # Survey marked as occurring on RIM (transect) with location blank
subcap[subcap$EFFORTID == "5618" & subcap$TRANSID == "5136" & subcap$TRANSECT == "UNKN" & subcap$LOCATION == "",c("TRANSECT","LOCATION")] <- c("X","10") # capture listed as unknown transect and blank location, but comments say it looks like a poorly written X10 (and surveys did happen on X that evening)
subcap[subcap$EFFORTID == "7793" & subcap$TRANSID == "37059" & subcap$TRANSECT == "" & subcap$LOCATION == "T4",c("TRANSECT","LOCATION")] <- c("T","4") # capture listed as unknown transect and blank location, but comments say it looks like a poorly written X10 (and surveys did happen on X that evening)
subcap <- subcap[!(subcap$EFFORTID == "9968" & subcap$TRANSID == "48687" & subcap$TRANSECT == "" & subcap$LOCATION == ""),] # Snake found incidentally while walking the transect (didn't remove survey as this was the trapping transect being checked)

## Transect + location = Point
subcap$Point <- paste(subcap$TRANSECT, subcap$LOCATION, sep = "")

#### CHECK FOR MISSING SURVEYS ####
###### Check that dates of captures match dates when surveys were conducted, this will flag for captures without survey dates but not the other way around (because snakes weren't caught at every survey so wouldn't match) ######
test <- ifelse(subcap$TRANSID %in% subsurv$TRANSID & subcap$EFFORTID %in% subsurv$EFFORTID, 1, 0)
for(i in 1:length(test)){
  if(test[i] == 0){
    stop('survey info does not exist for all captures')  ## originally flagged two records, fixed above so shouldn't result in error now
  }
}

## further investigate to fix any errors
# which(test == 0)


##### CORRECT CAPTURE LOCATION INFORMATION #####

#### SPECIFY STUDY AREA FOR EACH (HMU, NCR, and NWFN) and each project ####
ToCheck <- unique(subsurv[,c("SITE","PROJECTCODE")])
ToCheck$checked <- NA


#### ASOI SWIFT TRAP ARRAY ####

asoi <- read.csv("Data/TrapRingWeb_coords.csv")
# asoi2 <- asoi[,4:5]
# coordinates(asoi2)=~Lon+Lat
# proj4string(asoi2) <- CRS("+proj=longlat +datum=WGS84")
# 
# asoiUTM <- spTransform(asoi2, CRS("+proj=utm +zone=55 +ellps=WGS84"))
asoi2 <- asoi[,6:7]
coordinates(asoi2)=~UTM.X+UTM.Y
proj4string(asoi2) <- CRS("+proj=utm +zone=55 +ellps=WGS84 +datum=WGS84")

asoisw <- subset(subcap, SITE == "ASOI" & PROJECTCODE == "SWIFT TRAP ARRAY")[,c("TRANSECT","LOCATION")]

acaps <- merge(asoisw, asoi, by=c("LOCATION"))
coordinates(acaps) =~ UTM.X + UTM.Y
proj4string(acaps) <- CRS("+proj=utm +zone=55 +ellps=WGS84 +datum=WGS84")

asoi2 <- as.data.frame(asoi2)
acaps <- as.data.frame(acaps)

p1 <- ggplot() + geom_point(data=asoi2, aes(x=UTM.X, y=UTM.Y), col="black", cex=1) +
  geom_count(data=acaps, aes(x=UTM.X, y=UTM.Y), col="purple", alpha=0.5) +
  scale_size_area(c(0,10))



ToCheck[ToCheck$SITE == "ASOI" & ToCheck$PROJECTCODE == "SWIFT TRAP ARRAY" & is.na(ToCheck$checked),"checked"] <- 1


#### RTSI GNWR SSP BIRD/MOUSE ####

# bird <- readGPX("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_SSP/SSP Bird Trap Locations.gpx") #AB-EF, 1-16
# mouse <- readGPX("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_SSP/SSP Mouse Trap Locations.gpx") #A-F, 1-18 transects
# fence <- readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_SSP/Boundary Fence.kml")

# rtsi <- rbind(as.data.frame(mouse$waypoints[1:3]), as.data.frame(bird$waypoints[1:3]))
# rtsi$Type <- as.factor(c(rep(c("Mouse"),nrow(mouse$waypoints[1:2])), rep(c("Bird"),nrow(bird$waypoints[1:2]))))
# colnames(rtsi) <- c("Lon","Lat","LOCATION","Type")
# 
# # rtsicaps <- subset(subcap, SITE == "RTSI" & (PROJECTCODE == "GNWR SSP BIRD" | PROJECTCODE == "GNWR SSP MOUSE"))[,c("TRANSECT","LOCATION")]
# rtsicaps <- subset(subcap, SITE == "RTSI" & (PROJECTCODE == "GNWR SSP MOUSE"))[,c("TRANSECT","LOCATION")]
# 
# rcaps <- merge(rtsicaps, rtsi, by=c("LOCATION"))
# coordinates(rcaps) =~ Lon + Lat
# proj4string(rcaps) <- CRS("+proj=longlat +datum=WGS84")
# rcapsUTM <- spTransform(rcaps, CRS("+proj=utm +zone=55 +ellps=WGS84"))
# 
# coordinates(rtsi) =~ Lon + Lat
# proj4string(rtsi) <- CRS("+proj=longlat +datum=WGS84")
# rallUTM <- spTransform(rtsi, CRS("+proj=utm +zone=55 +ellps=WGS84"))
# 
# rcapsUTM <- as.data.frame(rcapsUTM)
# rallUTM <- as.data.frame(rallUTM)
# 
# p1 <- ggplot() + geom_point(data=rallUTM, aes(x=Lat, y=Lon), col="black", cex=1) + ## labeled Lat Lon still but actually UTM
#   geom_count(data=rcapsUTM, aes(x=Lat, y=Lon), col="purple", alpha=0.5) +
#   scale_size_area(c(0,10))

# ToCheck[ToCheck$SITE == "RTSI" & (ToCheck$PROJECTCODE == "GNWR SSP BIRD" | ToCheck$PROJECTCODE == "GNWR SSP MOUSE") & is.na(ToCheck$checked),"checked"] <- 1
ToCheck[ToCheck$SITE == "RTSI" & (ToCheck$PROJECTCODE == "GNWR SSP MOUSE") & is.na(ToCheck$checked),"checked"] <- 1


#### NWFN KB TRAP, NWFN TRAP 1, TRAP 2 LINVIS, TRAP 3, VISTRAP TRAP, POST BT2 TRAP, POST KB TRAP 1, POST KB TRAP 2, PRE BT1 TRAP ####
## Check on specifics of KB and BT methods

nwfncaps <- subset(subcap, SITE == "NWFN")[,c("TRANSECT","LOCATION")]

## Run CPlocs in specifystatespace.R
# source("SpecifyStateSpace.R")  ## FIGOURE OUT HOW TO CROSSWALK ALL FILES

# nwfncaps$Location <- paste(nwfncaps$TRANSECT,nwfncaps$LOCATION, sep="")
# nwfnallcaps <- merge(nwfncaps, CPlocs, by = c("Location"))
# 
# p1 <- ggplot() + geom_point(data = CPlocs, aes(x=x, y=y), pch=19, cex=1) +
#   geom_count(data=nwfnallcaps, aes(x=x, y=y), col="green", alpha=0.5)

ToCheck[ToCheck$SITE == "NWFN" & is.na(ToCheck$checked),"checked"] <- 1



#### HMUI TOX DROP TRAP ####
## KMZ file - HMUI EDGE EFFECTS VIS or EDGE EFFECTS PEK -> HMU VIS EDGE and HMU VIS INTERIOR
## KML HMU VIS EDGE/INTERIOR
# HMUItoxpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU TRAPS.kml","HMU TRAPS", require_geomType = "wkbPoint")
# HMUItoxline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU VIS EDGE.kml","HMU VIS EDGE", require_geomType = "wkbLineString")

# hmutxdrop <- subset(subcap, SITE == "HMUI" & PROJECTCODE == "TOX DROP TRAP")[,c("TRANSECT","LOCATION")]
# 
# ## Convert points to UTM
# hmuTOX <- as.data.frame(spTransform(HMUItoxpts, CRS("+proj=utm +zone=55")))
# hmuTOX$utmE <- as.vector(apply(as.data.frame(hmuTOX[,3]), 2, function(x) x-mean(x)))
# hmuTOX$utmN <- as.vector(apply(as.data.frame(hmuTOX[,4]), 2, function(x) x-mean(x)))
# hmuTOX <- hmuTOX %>%
#   separate(Name, 
#            into = c("TRANSECT", "LOCATION"),
#            sep = "(?<=[A-Za-z])(?=[0-9])"
#   )
# colnames(hmuTOX) <- c("TRANSECT","LOCATION","Description","UTME","UTMN","NA","utmE","utmN")
# hmuTOX$LOCATION <- sub("^0+", "", hmuTOX$LOCATION)
# 
# hmutoxcaps <- merge(hmutxdrop, hmuTOX[,-c(3,6:8)], by = c("TRANSECT","LOCATION"))
# 
# ## Plot to check
# p1 <- ggplot() + geom_point(data=hmuTOX, aes(x=UTME,y=UTMN), pch=19, cex=1) +
#   geom_count(data=hmutoxcaps, aes(x=UTME,y=UTMN), col="red", alpha=0.5)


ToCheck[ToCheck$SITE == "HMUI" & ToCheck$PROJECTCODE == "TOX DROP TRAP" & is.na(ToCheck$checked),"checked"] <- 1



#### NWFO NWFO TRAP 1, 2, 3, 4 ####

# nwfotrap <- subset(subcap, SITE == "NWFO" & (PROJECTCODE == "NWFO TRAP 1" | PROJECTCODE == "NWFO TRAP 2" | PROJECTCODE == "NWFO TRAP 3" | PROJECTCODE == "NWFO TRAP 4"))[,c("TRANSECT","LOCATION")]
# 
# ## Figure out how to crosswalk with SpecifyStateSpace.R
# 
# CPO <- cbind(CPOlocs,rownames(CPOlocs))
# colnames(CPO) <- c("x","y","Location")
# CPO <- CPO %>%
#   separate(Location, into = c("TRANSECT", "LOCATION"), sep = "(?<=[A-Za-z])(?=[0-9])")
# 
# nwfocaps <- merge(nwfotrap, CPO, by = c("TRANSECT","LOCATION"))
# 
# p1 <- ggplot() + geom_point(data=nwfocaps, aes(x=x,y=y), pch=19, cex=1) +
#   geom_count(data=nwfocaps, aes(x=x,y=y), col="orange", alpha=0.5)
# 
# 
# ToCheck[ToCheck$SITE == "NWFO" & (ToCheck$PROJECTCODE == "NWFO TRAP 1" | ToCheck$PROJECTCODE == "NWFO TRAP 2" | ToCheck$PROJECTCODE == "NWFO TRAP 3" | ToCheck$PROJECTCODE == "NWFO TRAP 4") & is.na(ToCheck$checked),"checked"] <- 1
#### KIND OF - NOT CHECKED EXACT LOCATION YET!!!! 


#### HMU/I/R TOX DROP 2 TRAP ####

# hmutxdrop2 <- subset(subcap, (SITE == "HMU" | SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "HMU TOX DROP 2 TRAP")[,c("TRANSECT","LOCATION")]

## TRANSECT: 1, blank, HMUI
## LOCATION: HMU59, HMU45, HMU12, HMU71, HMU91, HMU52, HMU78, blank, HMU67, HMU8, HMU94, HMU70, HMU37, HMU86, HMU69, HMU99, HMU26, HMU53, HMU73, HU31, HMU6, HMU93, HMU87, HMU79, HMU96, HMU82, HMU9, HMU60, HMU84, HMU27, HMU38, HMU48, HMU29, HMU10, HMU54, HMU2, HMU57, HMU30, HMU75, HMU19, HMU50, HMU51, HMU49, HMU18, HMU44, HMU90, HMU7, HMU23, HMU62, HMU34, HMU11


