## The purpose of this code is to: 
## 1) read in data (survey info and captures) from trapping surveys from numerous projects
## 2) select what projects to use
## 3) check capture and survey information and fix errors

library(lubridate); library(dplyr)



###### READ IN SURVEY DATA.----
allsurv <- read.csv("Real Data Analysis/Trapping/Data/TRAPinfo.csv")

## Check that only study area represented is NWFN (aka, CP = Closed Population)
subsurv <- subset(allsurv, SITE %in% c("NWFN"))

## Check using only projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("POST KB TRAP 1","POST KB TRAP 2", "NWFN TRAP 1","NWFN TRAP 2 LINVIS","NWFN TRAP 3","NWFN TRAP 4 LCM","NWFN VISTRAP TRAP","PRE BT1 TRAP","POST BT2 TRAP"))

#### CLEAN SURVEY DATA.----
# Survey marked as occurring on RIM (transect) with location blank
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "51100" & subsurv$EFFORTID == "9638"),]

## Specify date and separate year
subsurv$Date <- dmy(subsurv$Date)
subsurv$Year <- year(subsurv$Date)



###### READ IN CAPTURE DATA.----
allcap <- read.csv("Real Data Analysis/Trapping/Data/TRAPcaptures.csv")

## Check that only study area represented is NWFN (aka, CP = Closed Population)
subcap <- subset(allcap, SITE %in% c("NWFN"))

## Check using only projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subcap <- subset(subcap, PROJECTCODE %in% c("POST KB TRAP 1","POST KB TRAP 2", "NWFN TRAP 1","NWFN TRAP 2 LINVIS","NWFN TRAP 3","NWFN TRAP 4 LCM","NWFN VISTRAP TRAP","PRE BT1 TRAP","POST BT2 TRAP"))

## Specify date and separate year
subcap$Date <- dmy(subcap$Date)
subcap$Year <- year(subcap$Date)

#### CLEAN CAPTURE DATA.----
## Snake found incidentally dead, died during telemetry project likely due to feeding
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "57216" & subcap$EFFORTID == "12309"),]
## Removed from unknown transect
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "29784" & subcap$EFFORTID == "5310"),]  
## Snake found incidentally
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "55997" & subcap$EFFORTID == "12223"),]
## Snake incidental death
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "56296" & subcap$EFFORTID == "12295"),]
## Survey marked as occurring on RIM (transect) with location blank
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "51100" & subcap$EFFORTID == "9638"),]
## Capture listed on unknown transect and blank location, but comments say it looks like a poorly written X10 (and surveys did happen on X that evening)
subcap[subcap$EFFORTID == "5618" & subcap$TRANSID == "5136" & subcap$TRANSECT == "UNKN" & is.na(subcap$LOCATION),c("TRANSECT","LOCATION")] <- c("X","10") 

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




##### CORRECT CAPTURE LOCATION INFORMATION.----

## List of each project to check
ToCheck <- unique(subsurv[,c("SITE","PROJECTCODE")])
ToCheck$checked <- NA


## All NWFN projects

nwfncaps <- subset(subcap, SITE == "NWFN")[,c("TRANSECT","LOCATION")]

## Check that all captures have locations
if(nrow(subset(nwfncaps, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "")) > 0){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & is.na(ToCheck$checked),"checked"] <- 1


