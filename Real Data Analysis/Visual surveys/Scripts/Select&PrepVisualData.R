#### July 29, 2020; Code written by Staci Amburgey ####
## Read in data (survey info and captures) from visual surveys from numerous projects
## Select what projects to use

library(lubridate); library(dplyr)



###### READ IN SURVEY DATA.----
allsurv <- read.csv("Real Data Analysis/Visual surveys/Data/VISinfo.csv")

## Check that only study area represented is NWFN (aka, CP = Closed Population)
subsurv <- subset(allsurv, SITE %in% c("NWFN"))

## Check using only projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3"))

#### CLEAN SURVEY DATA.----
## Remove random NWFN survey that seems incorrect for NWFN VIS 2
subsurv <- subsurv[!(subsurv$TRANSID=="29700" & subsurv$EFFORTID=="4160"),]
## Remove random NWFN record where TRANSECT is missing (is 6 in capture file) and LOCATION is 440
subsurv <- subsurv[!(subsurv$TRANSID=="38317" & subsurv$EFFORTID=="8101"),]
## Rename typo
subsurv[subsurv$SITE == "NWFN" & subsurv$TRANSECT == "OPR", "TRANSECT"] <- "0PR"
## Rename typo
subsurv[subsurv$SITE == "NWFN" & subsurv$TRANSECT == "NEW", "TRANSECT"] <- "NWE"
## Denote surveys currently lacking transect info
subsurv$TRANSECT[subsurv$TRANSECT==""] <- "-9999"
## Remove survey that looks like incidental (same beginning and end time) but no capture recorded
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "30992" & subsurv$EFFORTID == "6078" & subsurv$TRANSECT == "-9999"),]
## Remove survey where no transect info is available (-9999)
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "35242" & subsurv$EFFORTID == "7289" & subsurv$TRANSECT == "-9999"),]
## Decide whether to include surveys on NWFN perimeter fence (except 0PR as that seems like an error or training survey) and edge transects (NEE, SWE). Incidentals (INCID, NSP, NSV, PR) are removed by time check lower down. Random transects (e.g., 0) are also removed.
subsurv <- subsurv[!(subsurv$SITE=="NWFN" & (subsurv$TRANSECT=="NWE" | subsurv$TRANSECT=="0" | subsurv$TRANSECT=="0PR")),]
## Rename 1PR and 2PR to just PR > 1 means a perimeter search at the beginning of the evening, 2 means at the end
subsurv[subsurv$TRANSECT == "2PR" | subsurv$TRANSECT == "1PR", "TRANSECT"] <- "PR"

## Specify date and separate year
subsurv$Date <- dmy(subsurv$Date)
subsurv$Year <- year(subsurv$Date)

## remove surveys that are incidental (same beginning and end time OR NA begin and end and elapsed of 0 or NA); removes many of the blank transect entries but also others with transect identified, some incidental surveys are listed as TRANSECT == INCID but get removed here as well using this criterion
remov <- subset(subsurv, BEGIN==END | (is.na(BEGIN) & is.na(END)))  ## keep record of what was removed
subsurv <- subsurv %>% anti_join(remov)  ## remove from surveys
subsurv <- droplevels(subsurv)
rownames(subsurv) <- 1:nrow(subsurv)




###### READ IN CAPTURE DATA.----
allcap <- read.csv("Real Data Analysis/Visual surveys/Data/VIScaptures.csv")

## Check that only study area represented is NWFN (aka, CP = Closed Population)
subcap <- subset(allcap, SITE %in% c("NWFN"))

## Check using only projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subcap <- subset(subcap, PROJECTCODE %in% c("NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3"))

#### CLEAN CAPTURE DATA.----
## NWFN VIS 2 record listed as N13.5, making N13
subcap[subcap$SITE == "NWFN" & subcap$CAPID == "4012","LOCATION"] <- "13"
## NWFN record matches a survey of SWE and another recorded listed as AA7 (but likely should be SWE)
subcap[subcap$EFFORTID == "14901" & subcap$TRANSID == "76172" & subcap$TRANSECT == "UNKN","TRANSECT"] <- "SWE"
subcap[subcap$EFFORTID == "14901" & subcap$TRANSID == "76172" & subcap$TRANSECT == "AA","TRANSECT"] <- "SWE"
## Remove NWFN capture with missing (blank) transect info (listed as -9999 from subsurv)
subcap <- subcap[!(subcap$TRANSID=="35242" & subcap$EFFORTID=="7289" & subcap$TRANSECT==""),]
## All NWFN others with transect UNKN are actually 2PR (removed from surveys but in case we're using need to rename)
check <- cbind(c("4034","4721","4342","4412","4474","4765","4566","4647","4648","4659","4296","4405","4385","4497","4090","4153","4422","4523","4567","4634","4158","4058","4670","4770","4123","4143","4189","4361","4387","4003","4462","4748","4763","4646","4086","4138","4026","4008","4546","4634","4562","242"), c("381","4238","2110","2515","2839","4488","3365","3830","3837","3889","1851","2476","2357","2967","697","1041","2573","3103","3374","3745","1069","513","3956","4520","873","982","1252","2224","2368","216","2768","4392","4478","3820","671","955","338","242","3244","3745","3341","4008"))
# test <- subset(subcap, EFFORTID %in% check[,1] & TRANSID %in% check[,2])
subcap[subcap$EFFORTID %in% check & subcap$TRANSID %in% check & subcap$TRANSECT == "UNKN","TRANSECT"] <- "2PR"
## Keep animals detected on NWFN perimeter fence (except 0PR, shouldn't exist but is an error in survey data) and edge transects (NEE, SWE). Incidentals (INCID, NSP, NSV) are removed by time check lower down. Random transects (e.g., 0) are also removed.
subcap <- subcap[!(subcap$SITE=="NWFN" & (subcap$TRANSECT=="0"| subcap$TRANSECT=="6" | subcap$TRANSECT=="0PR")),]
## Rename 1PR and 2PR to just PR > 1 means a perimeter search at the beginning of the evening, 2 means at the end
subcap[subcap$TRANSECT == "2PR" | subcap$TRANSECT == "1PR", "TRANSECT"] <- "PR"

## Specify date and separate year
subcap$Date <- dmy(subcap$Date)
subcap$Year <- year(subcap$Date)

##### REMOVE INCIDENTAL SURVEYS FROM CAPTURE DATA #####
## Remove animals from incidental surveys above and ones where TRAPTYPE == "I" for incidental
## some incidental surveys are listed as TRANSECT == INCI
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSID == "30992" & subcap$EFFORTID == "6078"),]  ## no capture but just making sure
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

subcap <- droplevels(subcap)
row.names(subcap) <- 1:nrow(subcap)

## Transect + location = Point
subcap$Point <- paste(subcap$TRANSECT, subcap$LOCATION, sep = "")




##### CORRECT CAPTURE LOCATION INFORMATION.----

## List of each project to check
ToCheck <- unique(subsurv[,c("SITE","PROJECTCODE")])
ToCheck$checked <- NA

## NWFN HL 1 and 2
nwfnVISHL12 <- subset(subcap, SITE == "NWFN" & (PROJECTCODE == "NWFN VIS HL 1" | PROJECTCODE == "NWFN VIS HL 2"))[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVISHL12, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "")) > 0){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & (ToCheck$PROJECTCODE == "NWFN VIS HL 1" | ToCheck$PROJECTCODE == "NWFN VIS HL 2") & is.na(ToCheck$checked),"checked"] <- 1


## NWFN 2
## One capture has a blank LOCATION (TRANSECT F, EFFORTID == 4089 & TRANSID == 688) but no details to figure out the location, remove capture (but not survey)
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSECT == "F" & subcap$LOCATION == ""),]
## One capture has a blank LOCATION (TRANSECT U, EFFORTID == 4235 & TRANSID == 1506) but no details to figure out the location, remove capture (but not survey)
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSECT == "U" & subcap$LOCATION == ""),]
subcap <- subset(subcap, !is.na(TRANSECT) & !is.na(SITE)) ## some reason the above rows leave empty (NA) rows

subcap <- subset(subcap, !is.na(TRANSECT) & !is.na(SITE)) ## removing the data above leaves the rows empty (NA)

nwfnVIS2 <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN VIS 2")[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVIS2, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN VIS 2" & is.na(ToCheck$checked),"checked"] <- 1


## PRE BT2 VIS
nwfnPREKB <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "PRE BT2 VIS")[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPREKB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "PRE BT2 VIS" & is.na(ToCheck$checked),"checked"] <- 1


## POST BT2 VIS
nwfnPOSTKB <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "POST BT2 VIS")[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPOSTKB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "POST BT2 VIS" & is.na(ToCheck$checked),"checked"] <- 1


## POST KB VIS 1-3
nwfnPOST13KB <- subset(subcap, SITE == "NWFN" & (PROJECTCODE == "POST KB VIS 1" | PROJECTCODE == "POST KB VIS 2" | PROJECTCODE == "POST KB VIS 3"))[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPOST13KB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & (ToCheck$PROJECTCODE == "POST KB VIS 1" | ToCheck$PROJECTCODE == "POST KB VIS 2" | ToCheck$PROJECTCODE == "POST KB VIS 3") & is.na(ToCheck$checked),"checked"] <- 1


## NWFN VISTRAP VIS
nwfnVISTRAP <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN VISTRAP VIS")[,c("TRANSECT","LOCATION","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVISTRAP, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN VISTRAP VIS" & is.na(ToCheck$checked),"checked"] <- 1

