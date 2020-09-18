#### July 29, 2020; S. Amburgey ####
## Read in data (survey info and captures) from visual surveys from numerous projects
## Select what projects to use
## Format for SCR analysis

rm(list=ls())

library(lubridate); library(reshape2); library(dplyr); library(sp); library(ggplot2); library(rgdal)


###### SURVEY DATA ######
allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VISinfo.csv")[,-c(21:62)]

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("NWFN","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("NWFN SCENT VIS TRAIL","NWFN TOXDROP VIS","NWFN VIS 1","NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2","NWFN VISPACE", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3","POST KB VIS 3 EXTRA","EDGE EFFECT VIS","LOWDENS SUPPVIS","LOWDENS VIS","TOX DROP VIS 1","TOX DROP VIS 2","TOX DROP VIS 3","HMU TOX DROP 2 VIS"))

## Rename two sites for simplification
subsurv$SITE[subsurv$SITE == "HMUI2B"] <- "HMUI2"
subsurv$SITE[subsurv$SITE == "HMUI5B"] <- "HMUI5"

## Specify date and separate year
subsurv$Date <- dmy(subsurv$Date)
subsurv$Year <- year(subsurv$Date)




###### CAPTURE DATA ######
allcap <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VIScaptures.csv")[,-c(6,15,22:25,27,29:32,34:36)]

## Study Areas deemed suitable
subcap <- subset(allcap, SITE %in% c("NWFN","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subcap <- subset(subcap, PROJECTCODE %in% c("NWFN SCENT VIS TRAIL","NWFN TOXDROP VIS","NWFN VIS 1","NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2","NWFN VISPACE", "NWFN VISTRAP VIS","PRE BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3","POST KB VIS 3 EXTRA","EDGE EFFECT VIS","LOWDENS SUPPVIS","LOWDENS VIS","TOX DROP VIS 1","TOX DROP VIS 2","TOX DROP VIS 3","HMU TOX DROP 2 VIS"))

## Specify date and separate year
subcap$Date <- dmy(subcap$Date)
subcap$Year <- year(subcap$Date)

## Transect + location = Point
subcap$Point <- paste(subcap$TRANSECT, subcap$LOCATION, sep = "")


#### REMOVE ALL HMUR TOX DROP 2 VIS?????



##### CREATE EFFORT MATRICES #####
## By person-hours
## By person-distance

## Clean data
subsurv$TRANSECT[subsurv$TRANSECT==""] <- "-9999"  ## denote surveys currently lacking transect info (NWFN, HMUI, HMUR, NCRI, NCRR, HMUI4, HMUI3)
## remove survey that looks like incidental (same beginning and end time) but no capture recorded
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "30992" & subsurv$EFFORTID == "6078" & subsurv$TRANSECT == "-9999"),]
## remove survey with no end time that looks like incidental but can't tell
subsurv <- subsurv[!(subsurv$SITE == "HMUI" & subsurv$TRANSID == "45499" & subsurv$EFFORTID == "9521" & subsurv$TRANSECT == "-9999"),]
## remove survey with no end time that looks like incidental but can't tell
subsurv <- subsurv[!(subsurv$SITE == "HMUR" & subsurv$TRANSID == "55518" & subsurv$EFFORTID == "12117" & subsurv$TRANSECT == "-9999"),]

## remove surveys that are incidental (same beginning and end time); removes many of the blank transect entries but also others with transect identified, some incidental surveys are listed as TRANSECT == INCID but get removed here as well using this criterion
remov <- subset(subsurv, BEGIN==END)  ## keep record of what was removed
surveys <- subsurv %>% anti_join(remov)  ## remove from surveys
surveys <- droplevels(surveys)

## Cast into dataframes showing elapsed time and distance by survey
efftime <- dcast(surveys, PROJECTCODE + SITE + Date + EFFORTID + SEARCHER ~ TRANSECT, mean, value.var = "ELAPSED", na.rm = TRUE)
effdist <- dcast(surveys, PROJECTCODE + SITE + Date + EFFORTID + SEARCHER ~ TRANSECT, mean, value.var = "DISTANCE", na.rm = TRUE)

efftime <- efftime[with(efftime, order(SITE, PROJECTCODE, Date)),]
effdist <- effdist[with(effdist, order(SITE, PROJECTCODE, Date)),]




##### REMOVE INCIDENTAL SURVEYS FROM CAPTURE DATA #####
## Remove animals from incidental surveys above and ones where TRAPTYPE == "I" for incidental
## some incidental surveys are listed as TRANSECT == INCI
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "30992" & subcap$EFFORTID == "6078"),]  ## no capture but just making sure
subcap <- subcap[!(subcap$SITE == "HMUI" & subcap$TRANSID == "45499" & subcap$EFFORTID == "9521"),]
subcap <- subcap[!(subcap$SITE == "HMUR" & subcap$TRANSID == "55518" & subcap$EFFORTID == "12117"),]
'%!in%' <- Negate('%in%')
subcap <- subcap[subcap$TRANSID %!in% as.vector(remov$TRANSID) & subcap$EFFORTID %!in% as.vector(remov$EFFORTID), ]

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




##### PREPARE DATA FOR ANALYSIS #####

#### SPECIFY STUDY AREA FOR EACH (HMU, NCR, and NWFN) and each project ####

## HMU EDGE EFFECT VIS
## KMZ file - HMUI EDGE EFFECTS VIS or EDGE EFFECTS PEK -> HMU VIS EDGE and HMU VIS INTERIOR
hmuEdgeT <- "/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMUI EDGE EFFECT VIS.kml"
readOGR(hmuEdgeT,"HMU VIS INTERIOR")   ### Problem reading
hmuEdge <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "EDGE EFFECT VIS")[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]
## several odd points in coordinates; coordinates are fairly imprecise at times
## HE04 location 1 (13.59550, 144.8643)
## HE08 location 3 (13.59433, 144.86)
## HE01 location 5 (13.59913, 144.8625)
## Solutions:
## replace HE04 coordinates with other HE04 location 1 coordinates (13.59690, 144.8668)
## HE08 3 can just usecoordinates for 2 as precision is low between coordinates within the same transect e.g. HE08 1 (13.59613, 144.8622) & 2 (13.59567, 144.8625) & 2.5 is (13.59613, 144.8622)
## replace HE01 location 5 with coordinates for 4.3 (13.60048, 144.8648)
hmuEdge[hmuEdge$TRANSECT == "HE04" & hmuEdge$CAPLAT == 13.5955,"CAPLAT"] <- 13.5969
hmuEdge[hmuEdge$TRANSECT == "HE04" & hmuEdge$CAPLON == 144.86432,"CAPLON"] <- 144.86678   ## won't show fifth decimal place, had to go look
hmuEdge[hmuEdge$TRANSECT == "HE08" & hmuEdge$CAPLAT == 13.59433,"CAPLAT"] <- 13.59567
hmuEdge[hmuEdge$TRANSECT == "HE08" & hmuEdge$CAPLON == 144.86,"CAPLON"] <- 144.8625
hmuEdge[hmuEdge$TRANSECT == "HE01" & hmuEdge$CAPLAT == 13.59913,"CAPLAT"] <- 13.60048
hmuEdge[hmuEdge$TRANSECT == "HE01" & hmuEdge$CAPLON == 144.86247,"CAPLON"] <-  144.86485

## Convert to UTM for easier handling in SCR
coordinates(hmuEdge)=~CAPLON+CAPLAT
proj4string(hmuEdge) <- CRS("+proj=longlat +datum=WGS84")
## Convert points to UTM
hmuEdge <- as.data.frame(spTransform(hmuEdge, CRS("+proj=utm +zone=55")))
hmuEdge$utmE <- as.vector(apply(as.data.frame(hmuEdge[,3]), 2, function(x) x-mean(x)))
hmuEdge$utmN <- as.vector(apply(as.data.frame(hmuEdge[,4]), 2, function(x) x-mean(x)))
colnames(hmuEdge) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")

## Plot to check
p1 <- ggplot(hmuEdge, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
  scale_fill_manual(values = c("#A3E4D7","#16A085","#58D68D","#82E0AA","#1E5C50","#19987F","#06FECC","#A6EADD","#E0EFEC","#D35400","#F5B041","#AA6030","#F5A06A","#F2D5C3","#E8C080","#996F2C","#DD8A05","#933C02"))

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
hmuEdelta <- 5  ## will need to play with this
XlhmuE<-min(hmuEdge[,3]) - hmuEdelta
XuhmuE<-max(hmuEdge[,3]) + hmuEdelta
YlhmuE<-min(hmuEdge[,4]) - hmuEdelta
YuhmuE<-max(hmuEdge[,4]) + hmuEdelta
AhmuE <- (XuhmuE-XlhmuE)*(YuhmuE-YlhmuE)



## HMU LOWDENS/SUPP VIS
## KMZ file - HMUI EDGE EFFECTS VIS or EDGE EFFECTS PEK -> HMU VIS EDGE and HMU VIS INTERIOR
hmuLD <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & (PROJECTCODE == "LOWDENS VIS" | PROJECTCODE == "LOWDENS SUPPVIS"))[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]

## Several odd points in coordinates
## Missing latitude for HMUI, HM TRANSECT and STARTNUMBER of 85 on 2015-06-15; the longitude that is there looks off as well
## Point on H6 looks incorrect but actually the fence does go at an angle around a utility building
## Solutions:
## Find 85m on Google Earth and replace lat/lon with that value
# library(geosphere)
# midPoint(c(144.862,13.59793),c(144.8625,13.59822))
hmuLD[hmuLD$TRANSECT == "HM" & is.na(hmuLD$CAPLAT),"CAPLAT"] <- 13.59808
hmuLD[hmuLD$TRANSECT == "HM" & hmuLD$CAPLON == 144.51724,"CAPLON"] <- 144.8622

## Convert to UTM for easier handling in SCR
coordinates(hmuLD)=~CAPLON+CAPLAT
proj4string(hmuLD) <- CRS("+proj=longlat +datum=WGS84")
## Convert points to UTM
hmuLD <- as.data.frame(spTransform(hmuLD, CRS("+proj=utm +zone=55")))
hmuLD$utmE <- as.vector(apply(as.data.frame(hmuLD[,3]), 2, function(x) x-mean(x)))
hmuLD$utmN <- as.vector(apply(as.data.frame(hmuLD[,4]), 2, function(x) x-mean(x)))
colnames(hmuLD) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")

## Plot to check
p1 <- ggplot(hmuLD, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
  scale_fill_manual(values = c("#A3E4D7","#16A085","#58D68D","#82E0AA","#1E5C50","#19987F","#06FECC","#A6EADD","#E0EFEC","#D35400","#F5B041","#AA6030","#F5A06A","#F2D5C3","#E8C080","#996F2C","#DD8A05","#933C02"))





## For troubleshooting
# x <- subset(hmuLD, TRANSECT == "H1") #& LOCATION == "85")
# coordinates(x) =~ CAPLON + CAPLAT
# proj4string(x) <- CRS("+proj=longlat +datum=WGS84")
# plot(hmuLD, pch=21, cex=0.5)
# plot(x, add=TRUE, pch=21, cex=1, col="red")


## NWFN
## Define study area grid (random example currently)
locs <- as.matrix(secr::make.grid(nx = 10, ny = 10, spacex = 8, spacey = 8))

## Which parts of grid have traps
set.seed(922020)
a=sample(100, 15)
X=locs[a,]
J <- nrow(X)
# Yscr <- locs
# ntraps <- nrow(locs)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
delta <- 5  ## will need to play with this
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)





###### COMBINE PROJECTS THAT CAN BE TREATED AS ONE DATASET (SAME TIME, AREA, AND METHODS) > check time still
subsurv$NEWPROJ <- ifelse(subsurv, (SITE == "HMUI" | SITE == "HMUR") & (PROJECTCODE == "LOWDENS SUPPVIS" | PROJECTCODE == "LOWDENS VIS"), 1,
                          ifelse(subsurv, (SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "EDGE EFFECT VIS", 2, 
                          ifelse(subsurv, (SITE == "HMUI" | SITE == "HMUR") & (PROJECTCODE == "TOX DROP VIS 1" | PROJECTCODE == "TOX DROP VIS 2" | PROJECTCODE == "TOX DROP VIS 3"), 3, 
                          ifelse(subsurv, (SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "HMU TOX DROP 2 VIS", 4, 
                          ifelse(subsurv, (SITE == "NCRI" | SITE == "NCRR") & PROJECTCODE == "EDGE EFFECT VIS", 5, 
                          ifelse(subsurv, SITE == "NWFN" & PROJECTCODE == "NWFN SCENT VIS TRAIL", 6, 
                          ifelse(subsurv, SITE == "NWFN" & PROJECTCODE == "NWFN TOXDROP VIS", 7,   ## Check this group
                          ifelse(subsurv, SITE == "NWFN" & (PROJECTCODE == "NWFN VIS 1" | PROJECTCODE == "NWFN VIS 2"), 8,
                          ifelse(subsurv, SITE == "NWFN" & (PROJECTCODE == "NWFN VIS HL 1" | PROJECTCODE == "NWFN VIS HL 2"), 9, 
                          ifelse(subsurv, SITE == "NWFN" & PROJECTCODE == "NWFN VISPACE", 10,
                          ifelse(subsurv, SITE == "NWFN" & PROJECTCODE == "NWFN VISTRAP VIS", 11,
                          ifelse(subsurv, SITE == "NWFN" & PROJECTCODE == "PRE BT2 VIS", 12,
                          ifelse(subsurv, SITE == "NWFN" & (PROJECTCODE == "POST BT2 VIS" | PROJECTCODE == "POST KB VIS 1" | PROJECTCODE == "POST KB VIS 2" | PROJECTCODE == "POST KB VIS 3" | PROJECTCODE == "POST KB VIS 3 EXTRA"), 13, 9999))))))))))))) ## check this group



