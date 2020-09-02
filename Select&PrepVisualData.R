#### July 29, 2020; S. Amburgey ####
## Read in data (survey info and captures) from visual surveys from numerous projects
## Select what projects to use
## Format for SCR analysis

rm(list=ls())

library(lubridate)


###### SURVEY DATA ######
allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VISinfo.csv")[,-c(21:62)]

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("NWFN","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("NWFN SCENT VIS TRAIL","NWFN TOXDROP VIS","NWFN VIS 1","NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2","NWFN VISPACE", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3","POST KB VIS 3 EXTRA","EDGE EFFECT VIS","LOWDENS SUPPVIS","LOWDENS VIS","TOX DROP VIS 1","TOX DROP VIS 2","TOX DROP VIS 3","HMU TOX DROP 2 VIS"))

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


##### THINGS STILL TO-DO! #####
## Remove locations that aren't part of grid in CP (but that might be specific to project)



###### Check that dates of captures match dates when surveys were conducted ######
checkDims <- function(subsurv, subcap){
  x <- unique(subsurv[,c("Date","EFFORTID","PROJECTCODE")])
  y <- unique(subcap[,c("Date","EFFORTID")])
  
  new <- merge(x, y, by=c("EFFORTID","Date"), all = TRUE)
  
  if(nrow(new) != nrow(unique(subsurv[,c("Date","EFFORTID","PROJECTCODE")]))) stop('mismatch in effort and capture dimensions')
  
}

checkDims(subsurv, subcap)


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



