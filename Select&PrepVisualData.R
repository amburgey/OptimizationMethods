#### July 29, 2020; S. Amburgey ####
## Read in data (survey info and captures) from visual surveys from numerous projects
## Select what projects to use
## Format for SCR analysis

library(lubridate); library(reshape2); library(dplyr); library(sp); library(ggplot2); library(rgdal); library(raster); library(sf); library(tidyverse); library(plyr); library(plotKML)

source("RenamingGrid.R")


###### SURVEY DATA ######
allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VISinfo.csv")[,-c(21:62)]

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("NWFN","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("NWFN SCENT VIS TRAIL","NWFN VIS 1","NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2","NWFN VISPACE", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3","POST KB VIS 3 EXTRA","EDGE EFFECT VIS","LOWDENS SUPPVIS","LOWDENS VIS","TOX DROP VIS 1","TOX DROP VIS 2","TOX DROP VIS 3","HMU TOX DROP 2 VIS"))
## removed "NWFN TOXDROP VIS" as need permission specifically from Melia to use

#### CLEAN SURVEY DATA ####
## Rename two sites to correct
subsurv$SITE[subsurv$SITE == "HMUI2B"] <- "HMUI2"
subsurv$SITE[subsurv$SITE == "HMUI5B"] <- "HMUI5"
## Remove random NWFN survey that seems incorrect and mistmatched for NWFN VIS 1
subsurv <- subsurv[!(subsurv$TRANSID=="13788" & subsurv$EFFORTID=="45"),]
## Remove random NWFN survey that seems incorrect for NWFN VIS 2
subsurv <- subsurv[!(subsurv$TRANSID=="29700" & subsurv$EFFORTID=="4160"),]
## Remove random HMUI and NWFN TOXDROP VIS record
subsurv <- subsurv[!(subsurv$TRANSID=="80081" & subsurv$EFFORTID=="15761"),]
## Remove random NWFN record where TRANSECT is missing (is 6 in capture file) and LOCATION is 440
subsurv <- subsurv[!(subsurv$TRANSID=="38317" & subsurv$EFFORTID=="8101"),]
## Remove survey done in HMU that has no start point defined and no capture can be found for the animal that was caught
subsurv <- subsurv[!(subsurv$TRANSID=="51871" & subsurv$EFFORTID=="9832"),]
## Remove one survey that is shorter than all others and doesn't record end time plus elapsed time seems suspect > doesn't match walking pace for all other surveys in this project
subsurv <- subsurv[!(subsurv$SITE == "HMUI" & subsurv$TRANSID=="45703" & subsurv$EFFORTID=="9572"),]
## Change one record that has SITE specified as NWFN when it should be HMUI
subsurv[subsurv$SITE == "NWFN" & subsurv$PROJECTCODE == "HMU TOX DROP 2 VIS", "SITE"] <- "HMUI"
## Rename typo
subsurv[subsurv$SITE == "NWFN" & subsurv$TRANSECT == "OPR", "TRANSECT"] <- "0PR"
## Rename typo
subsurv[subsurv$SITE == "NWFN" & subsurv$TRANSECT == "NEW", "TRANSECT"] <- "NWE"
## Subset to not include NWFN VIS 1 as it's missing survey info
subsurv <- subset(subsurv, PROJECTCODE != "NWFN VIS 1")
## Denote surveys currently lacking transect info (NWFN, HMUI, HMUR, NCRI, NCRR, HMUI4, HMUI3)
subsurv$TRANSECT[subsurv$TRANSECT==""] <- "-9999"  
## Remove survey that looks like incidental (same beginning and end time) but no capture recorded
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "30992" & subsurv$EFFORTID == "6078" & subsurv$TRANSECT == "-9999"),]
## Remove survey with no end time that looks like incidental but can't tell
subsurv <- subsurv[!(subsurv$SITE == "HMUI" & subsurv$TRANSID == "45499" & subsurv$EFFORTID == "9521" & subsurv$TRANSECT == "-9999"),]
## Remove survey with no end time that looks like incidental but can't tell
subsurv <- subsurv[!(subsurv$SITE == "HMUR" & subsurv$TRANSID == "55518" & subsurv$EFFORTID == "12117" & subsurv$TRANSECT == "-9999"),]
## Remove survey where no transect info is available (-9999)
subsurv <- subsurv[!(subsurv$SITE == "NWFN" & subsurv$TRANSID == "35242" & subsurv$EFFORTID == "7289" & subsurv$TRANSECT == "-9999"),]
## Include surveys on NWFN perimeter fence (except 0PR as that seems like an error or training survey) and edge transects (NEE, SWE). Incidentals (INCID, NSP, NSV, PR) are removed by time check lower down. Random transects (e.g., 0) are also removed.
subsurv <- subsurv[!(subsurv$SITE=="NWFN" & (subsurv$TRANSECT=="NWE" | subsurv$TRANSECT=="0" | subsurv$TRANSECT=="0PR")),]
# subsurv$TRANSECT=="1PR" | subsurv$TRANSECT=="PR" | subsurv$TRANSECT=="2PR" | subsurv$TRANSECT=="SWE" | subsurv$TRANSECT=="NEE" |
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

## RENAME NWFN PERIMETER SITES TO MATCH GRID ESTABLISHED
## Not a perfect match to CP but the dimensions and general guidelines are followed when splitting sections into grid cells
## Split into multiple ifelse statements as was exceeding stack memory
subsurv$STARTNUMBER <- as.numeric(as.character(subsurv$STARTNUMBER))  ## ignore warning, blank entries becomes NA
subsurv$STARTNUMBER <- gridrenam(data=subsurv, type=c("STARTNUMBER")) 




###### CAPTURE DATA ######
allcap <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VIScaptures.csv")[,-c(6,15,22:25,27,29:32,34:36)]

## Study Areas deemed suitable
subcap <- subset(allcap, SITE %in% c("NWFN","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subcap <- subset(subcap, PROJECTCODE %in% c("NWFN SCENT VIS TRAIL","NWFN VIS 1","NWFN VIS 2","NWFN VIS HL 1","NWFN VIS HL 2","NWFN VISPACE", "NWFN VISTRAP VIS","PRE BT2 VIS","POST BT2 VIS","POST KB VIS 1","POST KB VIS 2", "POST KB VIS 3","POST KB VIS 3 EXTRA","EDGE EFFECT VIS","LOWDENS SUPPVIS","LOWDENS VIS","TOX DROP VIS 1","TOX DROP VIS 2","TOX DROP VIS 3","HMU TOX DROP 2 VIS"))
## Remove "NWFN TOXDROP VIS" as need specific permission from Melia to use

#### CLEAN CAPTURE DATA ####
## Subset to not include NWFN VIS 1 as it's missing survey info
subcap <- subset(subcap, PROJECTCODE != "NWFN VIS 1")
## Remove random NWFN survey that seems incorrect and mistmatched for NWFN VIS 1
subcap <- subcap[!(subcap$TRANSID=="13788" & subcap$EFFORTID==45),]
## Mis-named NCRR records
## NCRR record matches a survey of RE09
subcap[subcap$SITE == "NCRR" & subcap$TRANSECT == "NCRR","TRANSECT"] <- "RE09"
## HE07 records matches a survey of RE07
subcap[subcap$SITE == "NCRR" & subcap$TRANSECT == "HE07","TRANSECT"] <- "RE07"
## NE08 records matches a survey of RE08
subcap[subcap$SITE == "NCRR" & subcap$TRANSECT == "NE08","TRANSECT"] <- "RE08"
## Mis-named NWFN records
## NWFN VIS 2 record listed as N13.5, making N13
subcap[subcap$SITE == "NWFN" & subcap$CAPID == "4012","LOCATION"] <- "13"
## NWFN record matches a survey of SWE and another recorded listed as AA7 (but likely should be SWE)
subcap[subcap$EFFORTID == "14901" & subcap$TRANSID == "76172" & subcap$TRANSECT == "UNKN","TRANSECT"] <- "SWE"
subcap[subcap$EFFORTID == "14901" & subcap$TRANSID == "76172" & subcap$TRANSECT == "AA","TRANSECT"] <- "SWE"
## Remove random NWFN survey with TRANSECT SE, no matching survey record
subcap <- subcap[!(subcap$TRANSID=="16503" & subcap$EFFORTID=="2760"),]
## Remove NWFN capture with missing (blank) transect info (listed as -9999 from subsurv)
subcap <- subcap[!(subcap$TRANSID=="35242" & subcap$EFFORTID=="7289" & subcap$TRANSECT==""),]
## Remove NWFN captures with missing transect (blank) and location (NA) info. According to survey info, should be transects Y and S but no idea where on it
subcap <- subcap[!(subcap$TRANSID=="42931" & subcap$EFFORTID=="9079" & subcap$TRANSECT==""),]
subcap <- subcap[!(subcap$TRANSID=="42932" & subcap$EFFORTID=="9079" & subcap$TRANSECT==""),]
## Correct typos with capture type (NEW vs. OLD)
## Blank is OLD and ODL should be OLD
subcap[subcap$EFFORTID == "13985" & subcap$TRANSID == "68400" & subcap$NEW == "","NEW"] <- "OLD"
subcap[subcap$EFFORTID == "13041" & subcap$TRANSID == "59308" & subcap$NEW == "ODL","NEW"] <- "OLD"
## All NWFN others with transect UNKN are actually 2PR (removed from surveys but in case we're using need to rename)
check <- cbind(c("4034","4721","4342","4412","4474","4765","4566","4647","4648","4659","4296","4405","4385","4497","4090","4153","4422","4523","4567","4634","4158","4058","4670","4770","4123","4143","4189","4361","4387","4003","4462","4748","4763","4646","4086","4138","4026","4008","4546","4634","4562","242"), c("381","4238","2110","2515","2839","4488","3365","3830","3837","3889","1851","2476","2357","2967","697","1041","2573","3103","3374","3745","1069","513","3956","4520","873","982","1252","2224","2368","216","2768","4392","4478","3820","671","955","338","242","3244","3745","3341","4008"))
# test <- subset(subcap, EFFORTID %in% check[,1] & TRANSID %in% check[,2])
subcap[subcap$EFFORTID %in% check & subcap$TRANSID %in% check & subcap$TRANSECT == "UNKN","TRANSECT"] <- "2PR"
## Keep animals detected on NWFN perimeter fence (except 0PR, shouldn't exist but is an error in survey data) and edge transects (NEE, SWE). Incidentals (INCID, NSP, NSV) are removed by time check lower down. Random transects (e.g., 0) are also removed.
subcap <- subcap[!(subcap$SITE=="NWFN" & (subcap$TRANSECT=="0"| subcap$TRANSECT=="6" | subcap$TRANSECT=="0PR")),]
# subcap$TRANSECT=="1PR" | subcap$TRANSECT=="PR" | subcap$TRANSECT=="2PR" | subcap$TRANSECT=="SWE" | subcap$TRANSECT=="NEE" | subcap$TRANSECT=="NWE" | 
## Rename 1PR and 2PR to just PR > 1 means a perimeter search at the beginning of the evening, 2 means at the end
subcap[subcap$TRANSECT == "2PR" | subcap$TRANSECT == "1PR", "TRANSECT"] <- "PR"

## Specify date and separate year
subcap$Date <- dmy(subcap$Date)
subcap$Year <- year(subcap$Date)

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

subcap <- droplevels(subcap)
row.names(subcap) <- 1:nrow(subcap)

## RENAME NWFN PERIMETER LOCATIONS TO MATCH GRID ESTABLISHED
## Not a perfect match to CP but the dimensions and general guidelines are followed when splitting sections into grid cells
## Split into multiple ifelse statements as was exceeding stack memory
subcap$LOCATION <- as.numeric(as.character(subcap$LOCATION))  ## ignore warning, blank entries becomes NA
subcap$LOCATION <- gridrenam(data=subcap, type= c("LOCATION"))

## Transect + location = Point
subcap$Point <- paste(subcap$TRANSECT, subcap$LOCATION, sep = "")



##### CREATE EFFORT MATRICES #####
## By person-hours
## By person-distance
## Cast into dataframes showing elapsed time and distance by survey
efftime <- dcast(subsurv, PROJECTCODE + SITE + Date + EFFORTID + SEARCHER ~ TRANSECT, mean, value.var = "ELAPSED", na.rm = TRUE)
effdist <- dcast(subsurv, PROJECTCODE + SITE + Date + EFFORTID + SEARCHER ~ TRANSECT, mean, value.var = "DISTANCE", na.rm = TRUE)

efftime <- efftime[with(efftime, order(SITE, PROJECTCODE, Date)),] ## -9999 will be dropped when HMU TOX DROP 2 VIS and POST BT2 VIS fixed
effdist <- effdist[with(effdist, order(SITE, PROJECTCODE, Date)),]




##### CORRECT CAPTURE LOCATION INFORMATION #####

#### SPECIFY STUDY AREA FOR EACH (HMU, NCR, and NWFN) and each project ####
ToCheck <- unique(subsurv[,c("SITE","PROJECTCODE")])
ToCheck$checked <- NA

## HMU EDGE EFFECT VIS
## KMZ file - HMUI EDGE EFFECTS VIS or EDGE EFFECTS PEK -> HMU VIS EDGE and HMU VIS INTERIOR
## KML HMU VIS EDGE/INTERIOR
# HMUREEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU VIS EDGE.kml","HMU VIS EDGE", require_geomType = "wkbLineString")
# HMUREEpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU VIS EDGE.kml","HMU VIS EDGE", require_geomType = "wkbPoint")
# HMUIEEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU VIS INTERIOR.kml","HMUI EDGE EFFECT VIS", require_geomType = "wkbLineString")
# HMUIEEpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMU VIS INTERIOR.kml","HMUI EDGE EFFECT VIS", require_geomType = "wkbPoint")

## several odd points in coordinates; coordinates are fairly imprecise at times
## HE04 location 1 (13.59550, 144.8643)
## HE08 location 3 (13.59433, 144.86)
## HE01 location 5 (13.59913, 144.8625)
## Solutions:
## replace HE04 coordinates with other HE04 location 1 coordinates (13.59690, 144.8668)
## HE08 3 can just usecoordinates for 2 as precision is low between coordinates within the same transect e.g. HE08 1 (13.59613, 144.8622) & 2 (13.59567, 144.8625) & 2.5 is (13.59613, 144.8622)
## replace HE01 location 5 with coordinates for 4.3 (13.60048, 144.8648)
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE04" & subcap$CAPLAT == 13.5955,"CAPLAT"] <- 13.5969
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE04" & subcap$CAPLON == 144.86432,"CAPLON"] <- 144.86678   ## won't show fifth decimal place, had to go look
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE08" & subcap$CAPLAT == 13.59433,"CAPLAT"] <- 13.59567
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE08" & subcap$CAPLON == 144.86,"CAPLON"] <- 144.8625
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE01" & subcap$CAPLAT == 13.59913,"CAPLAT"] <- 13.60048
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HE01" & subcap$CAPLON == 144.86247,"CAPLON"] <-  144.86485

hmuEdge <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "EDGE EFFECT VIS")[,c("PITTAG","Date","TRANSECT","LOCATION","CAPLAT","CAPLON")]

## Convert to UTM for easier handling in SCR
# coordinates(hmuEdge)=~CAPLON+CAPLAT
# proj4string(hmuEdge) <- CRS("+proj=longlat +datum=WGS84")
# ## Convert points to UTM
# hmuEdge <- as.data.frame(spTransform(hmuEdge, CRS("+proj=utm +zone=55")))
# hmuEdge$utmE <- as.vector(apply(as.data.frame(hmuEdge[,3]), 2, function(x) x-mean(x)))
# hmuEdge$utmN <- as.vector(apply(as.data.frame(hmuEdge[,4]), 2, function(x) x-mean(x)))
# colnames(hmuEdge) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# HMUREEbe <- as.data.frame(spTransform(HMUREEpts, CRS("+proj=utm +zone=55")))
# HMUREEbe <- cbind(c(rep(c("HE01be"),2),rep(c("HE02be"),2),rep(c("HE03be"),2),rep(c("HE04be"),2),rep(c("HE05be"),2),rep(c("HE06be"),2),rep(c("HE07be"),2),rep(c("HE08be"),2),rep(c("HE09be"),2)),c(1),HMUREEbe[,3:4],c(NA),c(NA))
# colnames(HMUREEbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# HMUIEEbe <- as.data.frame(spTransform(HMUIEEpts, CRS("+proj=utm +zone=55")))
# HMUIEEbe <- cbind(c("HP01b","HP02b","HP03b","HP04b","HP05b","HP06b","HP07b","HP08b","HP09b"),c(1),HMUIEEbe[10:18,3:4],c(NA),c(NA))
# colnames(HMUIEEbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# hmuEdgeAll <- rbind(hmuEdge,HMUREEbe,HMUIEEbe)

# ## Plot to check
# p1 <- ggplot(hmuEdgeAll, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) +
#   scale_fill_manual(values = c("#A3E4D7","black","#16A085","black","#58D68D","black","#82E0AA","black","#1E5C50","black","#19987F","black","#06FECC","black","#A6EADD","black","#E0EFEC","black","#D35400","black","#F5B041","black","#AA6030","black","#F5A06A","black","#F2D5C3","black","#E8C080","black","#996F2C","black","#DD8A05","black","#933C02","black"))



ToCheck[(ToCheck$SITE == "HMUI" | ToCheck$SITE == "HMUR") & ToCheck$PROJECTCODE == "EDGE EFFECT VIS" & is.na(ToCheck$checked),"checked"] <- 1


## HMU LOWDENS/SUPP VIS
## KMZ file - HMUI EDGE EFFECTS VIS or EDGE EFFECTS PEK -> HMU VIS EDGE and HMU VIS INTERIOR
## KML HMUI/HMUR LOW DENS
# HMURLDline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMUR_transects LOW DENS.kml","HMUR", require_geomType = "wkbLineString")
# HMURLDpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMUR_transects LOW DENS.kml","HMUR", require_geomType = "wkbPoint")
# HMUILDline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMUI_transects LOW DENS.kml","HMUI_transects", require_geomType = "wkbLineString")
# HMUILDpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/HMUI_transects LOW DENS.kml","HMUI_transects", require_geomType = "wkbPoint")

## Several odd points in coordinates
## Missing latitude for HMUI, HM TRANSECT and STARTNUMBER of 85 on 2015-06-15; the longitude that is there looks off as well
## Point on H6 looks incorrect but actually the fence does go at an angle around a utility building
## Solutions:
## Find 85m on Google Earth and replace lat/lon with that value
# library(geosphere)
# midPoint(c(144.862,13.59793),c(144.8625,13.59822))
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HM" & is.na(subcap$CAPLAT),"CAPLAT"] <- 13.59808
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HM" & subcap$CAPLON == 144.51724,"CAPLON"] <- 144.8622

hmuLD <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & (PROJECTCODE == "LOWDENS VIS" | PROJECTCODE == "LOWDENS SUPPVIS"))[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]

## Convert to UTM for easier handling in SCR
# coordinates(hmuLD)=~CAPLON+CAPLAT
# proj4string(hmuLD) <- CRS("+proj=longlat +datum=WGS84")
# ## Convert points to UTM
# hmuLD <- as.data.frame(spTransform(hmuLD, CRS("+proj=utm +zone=55")))
# hmuLD$utmE <- as.vector(apply(as.data.frame(hmuLD[,3]), 2, function(x) x-mean(x)))
# hmuLD$utmN <- as.vector(apply(as.data.frame(hmuLD[,4]), 2, function(x) x-mean(x)))
# colnames(hmuLD) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# HMURLDbe <- as.data.frame(spTransform(HMURLDpts, CRS("+proj=utm +zone=55")))
# HMURLDbe <- cbind(c(rep(c("H1be"),2),rep(c("H2be"),2),rep(c("H3be"),2),rep(c("H4be"),2),rep(c("H5be"),2),rep(c("H6be"),2),rep(c("H7be"),2),rep(c("H8be"),2)),c(1),HMURLDbe[,3:4],c(NA),c(NA))
# colnames(HMURLDbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# HMUILDbe <- as.data.frame(spTransform(HMUILDpts, CRS("+proj=utm +zone=55")))
# HMUILDbe <- cbind(c(rep(c("HKbe"),2),rep(c("HLbe"),2),rep(c("HMbe"),2),rep(c("HNbe"),2),rep(c("HObe"),2),rep(c("HPbe"),2),rep(c("HQbe"),2),rep(c("HRbe"),2),rep(c("HSbe"),2),rep(c("HTbe"),2)),c(1),HMUILDbe[,3:4],c(NA),c(NA))
# colnames(HMUILDbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# hmuLDAll <- rbind(hmuLD,HMURLDbe,HMUILDbe)
# 
# ## Plot to check
# p1 <- ggplot(hmuLDAll, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
#   scale_fill_manual(values = c("#A3E4D7","black","#16A085","black","#58D68D","black","#82E0AA","black","#1E5C50","black","#19987F","black","#06FECC","black","#A6EADD","black","#E0EFEC","black","#D35400","black","#F5B041","black","#AA6030","black","#F5A06A","black","#F2D5C3","black","#E8C080","black","#996F2C","black","#DD8A05","black","#933C02","black"))


ToCheck[(ToCheck$SITE == "HMUI" | ToCheck$SITE == "HMUR") & (ToCheck$PROJECTCODE == "LOWDENS VIS" | ToCheck$PROJECTCODE == "LOWDENS SUPPVIS") & is.na(ToCheck$checked),"checked"] <- 1




## HMUI/R TOX DROP VIS 1 & 2
## Uses same transects as EDGE (HE01-HE09, HP01-HP09)
## Use kml file from above; HMUREEline, HMUREEpts, HMUIEEline, HMUIEEpts
## Several odd points in coordinates
## The end points of HP03 are outside HMU (location 25; 13.59730, 144.8588) and nearly to other edge (72.5; 13.60032, 144.8638)
## Used Google Earth to map other points on HP03 and then measure out distance for these two points
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HP03" & subcap$CAPLAT == 13.59730,"CAPLAT"] <- 13.598772
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HP03" & subcap$CAPLON == 144.85876,"CAPLON"] <- 144.861086
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HP03" & subcap$CAPLAT == 13.60032,"CAPLAT"] <- 13.599015
subcap[(subcap$SITE == "HMUI" | subcap$SITE == "HMUR") & subcap$TRANSECT == "HP03" & subcap$CAPLON == 144.86384,"CAPLON"] <- 144.861619

hmuTD12 <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & (PROJECTCODE == "TOX DROP VIS 1" | PROJECTCODE == "TOX DROP VIS 2"))[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]

## Convert to UTM for easier handling in SCR
# coordinates(hmuTD12)=~CAPLON+CAPLAT
# proj4string(hmuTD12) <- CRS("+proj=longlat +datum=WGS84")
# ## Convert points to UTM
# hmuTD12 <- as.data.frame(spTransform(hmuTD12, CRS("+proj=utm +zone=55")))
# hmuTD12$utmE <- as.vector(apply(as.data.frame(hmuTD12[,3]), 2, function(x) x-mean(x)))
# hmuTD12$utmN <- as.vector(apply(as.data.frame(hmuTD12[,4]), 2, function(x) x-mean(x)))
# colnames(hmuTD12) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# hmuTD12All <- rbind(hmuTD12,HMUREEbe,HMUIEEbe)
# 
# ## Plot to check
# p1 <- ggplot(hmuTD12All, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
#   scale_fill_manual(values = c("#A3E4D7","black","#16A085","black","#58D68D","black","#82E0AA","black","#1E5C50","black","#19987F","black","#06FECC","black","#A6EADD","black","#E0EFEC","black","#D35400","black","#F5B041","black","#AA6030","black","#F5A06A","black","#F2D5C3","black","#E8C080","black","#996F2C","black","#DD8A05","black","#933C02","black"))


ToCheck[(ToCheck$SITE == "HMUI" | ToCheck$SITE == "HMUR") & (ToCheck$PROJECTCODE == "TOX DROP VIS 1" | ToCheck$PROJECTCODE == "TOX DROP VIS 2") & is.na(ToCheck$checked),"checked"] <- 1




## HMUI/R TOX DROP VIS 3
## Uses same transects as LOWDENS (H1-9, HK-HT)
## Use kml file from above; HMURLDline, HMURLDpts, HMUILDline, HMUILDpts
hmuTD3 <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR") & PROJECTCODE == "TOX DROP VIS 3")[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]
## No coordinate issue (remember, point on H6 looks incorrect but actually the fence does go at an angle around a utility building)
## Convert to UTM for easier handling in SCR
# coordinates(hmuTD3)=~CAPLON+CAPLAT
# proj4string(hmuTD3) <- CRS("+proj=longlat +datum=WGS84")
# ## Convert points to UTM
# hmuTD3 <- as.data.frame(spTransform(hmuTD3, CRS("+proj=utm +zone=55")))
# hmuTD3$utmE <- as.vector(apply(as.data.frame(hmuTD3[,3]), 2, function(x) x-mean(x)))
# hmuTD3$utmN <- as.vector(apply(as.data.frame(hmuTD3[,4]), 2, function(x) x-mean(x)))
# colnames(hmuTD3) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# hmuTD3All <- rbind(hmuTD3,HMURLDbe,HMUILDbe)
# 
# ## Plot to check
# p1 <- ggplot(hmuTD3All, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
#   scale_fill_manual(values = c("#A3E4D7","black","#16A085","black","#58D68D","black","#82E0AA","black","#1E5C50","black","#19987F","black","#06FECC","black","#A6EADD","black","#E0EFEC","black","#D35400","black","#F5B041","black","#AA6030","black","#F5A06A","black","#F2D5C3","black","#E8C080","black","#996F2C","black","#DD8A05","black","#933C02","black"))

ToCheck[(ToCheck$SITE == "HMUI" | ToCheck$SITE == "HMUR") & ToCheck$PROJECTCODE == "TOX DROP VIS 3" & is.na(ToCheck$checked),"checked"] <- 1




##### NEED TO DECIDE IF ADDING OR NOT #####
## HMUI/R 1-5 HMU TOX DROP 2 VIS
## KMZ file - missing
## KML file - missing
# hmu15TD2 <- subset(subcap, (SITE == "HMUI" | SITE == "HMUR" |  SITE == "HMUI1" | SITE == "HMUI2" | SITE == "HMUI3" | SITE == "HMUI4" | SITE == "HMUI5") & PROJECTCODE == "HMU TOX DROP 2 VIS")[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]


# ## For troubleshooting
# x <- subset(ncrEE, TRANSECT == "RE01" )# & LOCATION == "2") #  "RE01"
# coordinates(x) =~ CAPLON + CAPLAT
# proj4string(x) <- CRS("+proj=longlat +datum=WGS84")
# coordinates(ncrEE) =~ CAPLON + CAPLAT
# proj4string(ncrEE) <- CRS("+proj=longlat +datum=WGS84")
# plot(ncrEE, pch=21, cex=0.5)
# plot(x, add=TRUE, pch=21, cex=1, col="red")




## NCRI/R EDGE EFFECT VIS
## Uses kmz and kml NCRI EDGE EFFECTS VIS
# NCRREEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRR EDGE EFFECT VIS.kml","NCRR EDGE EFFECT VIS", require_geomType = "wkbLineString")
# NCRREEpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRR EDGE EFFECT VIS.kml","NCRR EDGE EFFECT VIS", require_geomType = "wkbPoint")
# NCRIEEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRI EDGE EFFECT VIS.kml","NCRI EDGE EFFECT VIS", require_geomType = "wkbLineString")
# NCRIEEpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRI EDGE EFFECT VIS.kml","NCRI EDGE EFFECT VIS", require_geomType = "wkbPoint")

## Several odd points in coordinates
## First, where is the line and pts info for RP10-14?
## RP14 and Location 121
## RP14 and Location 170
## HE07, NCRR, and NE08 transect seems incorrect or unclear > fixed above so reflected in effort matrix
## RP13 has three problem locations, 37, 4, 11, 69
## RP07 and Location 20
## RP11 and Locations 93, 116, 125
## RP10 and Locations 7, 40, 40, 46, 175
## RP09 and Locations 133, 170, 196
## RE08 and Location 1 (13.59905, 144.8601 from 2013-08-28) looks to be on different transect but survey info matches for RE08
## RE04 and Location 2
## RE01 and Location 2 (13.59042, 144.8616)
## Solutions:
## NP14 Location 170 shifted back based on Location 162 and 176 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP14" & subcap$CAPLAT == 13.58819,"CAPLAT"] <- 13.589716
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP14" & subcap$CAPLON == 144.85945,"CAPLON"] <- 144.861897
## NP14 Location 121 shifted back based on Location 115 and 130 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP14" & subcap$CAPLAT == 13.58789,"CAPLAT"] <-  13.589413
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP14" & subcap$CAPLON == 144.85965,"CAPLON"] <- 144.862106
## RP13 Location 37 shifted back based on Location 30 and 48 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLAT == 13.58795,"CAPLAT"] <- 13.589342
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLON == 144.86130,"CAPLON"] <- 144.863682
## RP13 Location 4 shifted back based on Location 3 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLAT == 13.58914,"CAPLAT"] <-  13.589057
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLON == 144.86286,"CAPLON"] <- 144.863894
## RP13 Location 11 shifted back based on Location 3 and 22 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLAT == 13.58769,"CAPLAT"] <- 13.589113
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLON == 144.86147,"CAPLON"] <- 144.863848
## RP13 Location 69 shifted back based on Location 53 and 77 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLAT == 13.59017,"CAPLAT"] <- 13.589567
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP13" & subcap$CAPLON == 144.86281,"CAPLON"] <- 144.863549
## RP07 Location 20 shifted back based on Location 30 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP07" & subcap$CAPLAT == 13.58617,"CAPLAT"] <- 13.587806
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP07" & subcap$CAPLON == 144.85852,"CAPLON"] <- 144.860846
## RP11 Location 93 shifted back based on Location 95 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLAT == 13.58889,"CAPLAT"] <- 13.587496
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLON == 144.86183,"CAPLON"] <- 144.859392
## RP11 Location 116 shifted back based on Location 100 and 140 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLAT == 13.58750,"CAPLAT"] <- 13.587655
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLON == 144.85891,"CAPLON"] <- 144.859309
## RP11 Location 125 shifted back based on Location 100 and 140 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLAT == 13.58637,"CAPLAT"] <- 13.587729
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP11" & subcap$CAPLON == 144.85681,"CAPLON"] <- 144.859274
## RP10 Location 7 shifted back based on Location 15 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLAT == 13.58559,"CAPLAT"] <- 13.586943
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLON == 144.85750,"CAPLON"] <- 144.859945
## RP10 both Location 40s shifted back based on Location 40 coordinates that are correct (on transect)
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLAT == 13.58725,"CAPLAT"] <- 13.58724
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLON == 144.85976,"CAPLON"] <- 144.85980
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLAT == 13.58866,"CAPLAT"] <- 13.58724
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLON == 144.86225,"CAPLON"] <- 144.85980
## RP10 Location 175 shifted back based on Location 175 coordinates that are correct (on transect)
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLAT == 13.58696,"CAPLAT"] <- 13.58828
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLON == 144.85684,"CAPLON"] <- 144.85928
## RP10 Location 46 shifted back based on Location 40 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLAT == 13.58221,"CAPLAT"] <- 13.587250
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP10" & subcap$CAPLON == 144.85974,"CAPLON"] <- 144.859796
## RP09 Location 133 shifted back based on Location 135 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$CAPLAT == 13.58954,"CAPLAT"] <- 13.588005
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$CAPLON == 144.86213,"CAPLON"] <- 144.859845
## RP09 Location 170 shifted back based on Location 154 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$CAPLAT == 13.58702,"CAPLAT"] <- 13.588344
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$CAPLON == 144.85719,"CAPLON"] <- 144.859558
## RP09 Location 196 shifted back based on Location 154 and 170 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$LOCATION == "196" & subcap$CAPLAT == 13.58727,"CAPLAT"] <-  13.588539  ## another with same latitude
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RP09" & subcap$CAPLON == 144.85709,"CAPLON"] <- 144.859423
## RE08 and Location 1 shifted back based on Location of other 1s on this transect
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE08" & subcap$CAPLAT == 13.59905,"CAPLAT"] <-  13.59993
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE08" & subcap$CAPLON == 144.86009,"CAPLON"] <- 144.85976
## RE04 and Location 2 shifted back based on Location 2.1 coordinates plotted in Google Earth
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE04" & subcap$CAPLAT == 13.59417,"CAPLAT"] <- 13.595095
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE04" & subcap$CAPLON == 144.85954,"CAPLON"] <- 144.862101
## RE01 and Location 2 shifted back based on other Location 2 coordinates
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE01" & subcap$CAPLAT == 13.59042,"CAPLAT"] <- 13.59155
subcap[(subcap$SITE == "NCRI" | subcap$SITE == "NCRR") & subcap$TRANSECT == "RE01" & subcap$CAPLON == 144.86162,"CAPLON"] <- 144.86421

ncrEE <- subset(subcap, (SITE == "NCRI" | SITE == "NCRR") & PROJECTCODE == "EDGE EFFECT VIS")[,c("TRANSECT","LOCATION","CAPLAT","CAPLON")]

## Convert to UTM for easier handling in SCR
# coordinates(ncrEE)=~CAPLON+CAPLAT
# proj4string(ncrEE) <- CRS("+proj=longlat +datum=WGS84")
# ## Convert points to UTM
# ncrEE <- as.data.frame(spTransform(ncrEE, CRS("+proj=utm +zone=55")))
# ncrEE$utmE <- as.vector(apply(as.data.frame(ncrEE[,3]), 2, function(x) x-mean(x)))
# ncrEE$utmN <- as.vector(apply(as.data.frame(ncrEE[,4]), 2, function(x) x-mean(x)))
# colnames(ncrEE) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# NCRREEbe <- as.data.frame(spTransform(NCRREEpts, CRS("+proj=utm +zone=55")))
# NCRREEbe <- cbind(c(rep(c("RE01be"),2),rep(c("RE02be"),2),rep(c("RE03be"),2),rep(c("RE04be"),2),rep(c("RE05be"),2),rep(c("RE06be"),2),rep(c("RE07be"),2),rep(c("RE08be"),2),rep(c("RE09be"),2)),c(1),NCRREEbe[,3:4],c(NA),c(NA))
# colnames(NCRREEbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# NCRIEEbe <- as.data.frame(spTransform(NCRIEEpts, CRS("+proj=utm +zone=55")))[10:18,]
# NCRIEEbe <- cbind(c(rep(c("RP01be"),2),rep(c("RP02be"),2),rep(c("RP03be"),2),rep(c("RP04be"),2),rep(c("RP05be"),2),rep(c("RP06be"),2),rep(c("RP07be"),2),rep(c("RP08be"),2),rep(c("RP09be"),2)),c(1),NCRIEEbe[,3:4],c(NA),c(NA))
# colnames(NCRIEEbe) <- c("TRANSECT","LOCATION","UTME","UTMN","utmE","utmN")
# 
# ncrEEAll <- rbind(ncrEE,NCRREEbe,NCRIEEbe)
# 
# ## Plot to check
# p1 <- ggplot(ncrEEAll, aes(x=UTME, y=UTMN,fill=TRANSECT)) + geom_point(pch=21) + 
#   scale_fill_manual(values = c("#A3E4D7","black","#16A085","black","#58D68D","black","#82E0AA","black","#1E5C50","black","#19987F","black","#06FECC","black","#A6EADD","black","#E0EFEC","black","#D35400","black","#F5B041","black","#AA6030","black","#F5A06A","black","#F2D5C3","black","#E8C080","black","#996F2C","black","#DD8A05","black","#933C02","black","blue","black","red","black","orange","black","green",'black',"purple","black","yellow","black"))

ToCheck[(ToCheck$SITE == "NCRI" | ToCheck$SITE == "NCRR") & ToCheck$PROJECTCODE == "EDGE EFFECT VIS" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN NWFN HL 1 and 2
## Only KMZ or KML files show only the TRAP stations but this area uses a standard grid of 13 x 27 cells
# NWFNtrpts <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/CP TRAP STATIONS.kml","Waypoints", require_geomType = "wkbPoint")

nwfnVISHL12 <- subset(subcap, SITE == "NWFN" & (PROJECTCODE == "NWFN VIS HL 1" | PROJECTCODE == "NWFN VIS HL 2"))[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVISHL12, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "")) > 0){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & (ToCheck$PROJECTCODE == "NWFN VIS HL 1" | ToCheck$PROJECTCODE == "NWFN VIS HL 2") & is.na(ToCheck$checked),"checked"] <- 1




## NWFN VIS 1 and 2; survey is missing from all VIS 1 so removed for now
## Same KMZ file showing only the TRAP stations
## One capture has a blank LOCATION (TRANSECT F, EFFORTID == 4089 & TRANSID == 688) but no details to figure out the location, remove capture (but not survey)
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSECT == "F" & subcap$LOCATION == ""),]
## One capture has a blank LOCATION (TRANSECT U, EFFORTID == 4235 & TRANSID == 1506) but no details to figure out the location, remove capture (but not survey)
subcap <- subcap[!(subcap$SITE == "NWFN" & subcap$TRANSECT == "U" & subcap$LOCATION == ""),]
subcap <- subset(subcap, !is.na(TRANSECT) & !is.na(SITE)) ## some reason the above rows leave empty (NA) rows

nwfnVIS12 <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN VIS 2")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVIS12, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN VIS 2" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN PRE BT2 VIS
## Same KMZ file showing only the TRAP stations
nwfnPREKB <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "PRE BT2 VIS")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPREKB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "PRE BT2 VIS" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN POST BT2 VIS
## Same KMZ file showing only the TRAP stations
nwfnPOSTKB <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "POST BT2 VIS")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPOSTKB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "POST BT2 VIS" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN POST KB VIS 1-3 and 3 EXTRA
## Same KMZ file showing only the TRAP stations
nwfnPOST13KB <- subset(subcap, SITE == "NWFN" & (PROJECTCODE == "POST KB VIS 1" | PROJECTCODE == "POST KB VIS 2" | PROJECTCODE == "POST KB VIS 3" | PROJECTCODE == "POST KB VIS 3 EXTRA"))[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnPOST13KB, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & (ToCheck$PROJECTCODE == "POST KB VIS 1" | ToCheck$PROJECTCODE == "POST KB VIS 2" | ToCheck$PROJECTCODE == "POST KB VIS 3" | ToCheck$PROJECTCODE == "POST KB VIS 3 EXTRA") & is.na(ToCheck$checked),"checked"] <- 1




## NWFN VISTRAP VIS
## Same KMZ file showing only the TRAP stations
nwfnVISTRAP <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN VISTRAP VIS")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVISTRAP, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN VISTRAP VIS" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN VISPACE
## Same KMZ file showing only the TRAP stations
nwfnVISPACE <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN VISPACE")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnVISPACE, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN VISPACE" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN SCENT VIS TRAIL
## Same KMZ file showing only the TRAP stations
nwfnSCENT <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN SCENT VIS TRAIL")[,c("TRANSECT","LOCATION","COMMENT","Point")]

## Check that all captures have locations
if(nrow(subset(nwfnSCENT, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
  stop('info missing for grid capture location or transect')
}

ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN SCENT VIS TRAIL" & is.na(ToCheck$checked),"checked"] <- 1




## NWFN TOXDROP VIS
## Same KMZ file showing only the TRAP stations
# nwfnTOX <- subset(subcap, SITE == "NWFN" & PROJECTCODE == "NWFN TOXDROP VIS")[,c("TRANSECT","LOCATION","COMMENT","Point")]
# 
# ## Check that all captures have locations
# if(nrow(subset(nwfnTOX, is.na(TRANSECT) | is.na(LOCATION) | TRANSECT == "" | LOCATION == "") > 0)){
#   stop('info missing for grid capture location or transect')
# }
# 
# ToCheck[ToCheck$SITE == "NWFN" & ToCheck$PROJECTCODE == "NWFN TOXDROP VIS" & is.na(ToCheck$checked),"checked"] <- 1


### SUM INFO BY EACH STUDY ###
vissub <- subset(subsurv, SITE == "NWFN" & PROJECTCODE == "POST BT2 VIS" & TARGET == "BI")
ggplot(vissub, aes(x=Date, y=BI)) + geom_point()
sumdat <- as.data.frame(table(vissub$Date, vissub$TARGET))
ggplot(sumdat, aes(x=Var1,y=Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Date") + ylab("Count Across All Transects")


#### SUMMARIZE INFO (project length, captures, recaptures) ABOUT VARIOUS STUDIES ####

## Days each project was conducted and date range
subsurv <- subsurv[subsurv$PROJECTCODE != "HMU TOX DROP 2 VIS",]

times <- aggregate(data=subsurv, Date ~ PROJECTCODE, function(x) length(unique(x)))
times$min <- aggregate(data=subsurv, Date ~ PROJECTCODE, function(x) min(x))[,2]
times$max <- aggregate(data=subsurv, Date ~ PROJECTCODE, function(x) max(x))[,2]

tplot <- ggplot(times, aes(x=Date, y=PROJECTCODE)) + geom_bar(stat = "identity") + xlim(0,425) +
  xlab("Days Surveyed") +
  annotate(geom="text", x=(times$Date[1]+100), y=times$PROJECTCODE[1], label=paste(times$min[1]," to ",times$max[1],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[2]+100), y=times$PROJECTCODE[2], label=paste(times$min[2]," to ",times$max[2],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[3]+100), y=times$PROJECTCODE[3], label=paste(times$min[3]," to ",times$max[3],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[4]+100), y=times$PROJECTCODE[4], label=paste(times$min[4]," to ",times$max[4],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[5]+100), y=times$PROJECTCODE[5], label=paste(times$min[5]," to ",times$max[5],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[6]+100), y=times$PROJECTCODE[6], label=paste(times$min[6]," to ",times$max[6],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[7]+100), y=times$PROJECTCODE[7], label=paste(times$min[7]," to ",times$max[7],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[8]+100), y=times$PROJECTCODE[8], label=paste(times$min[8]," to ",times$max[8],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[9]+100), y=times$PROJECTCODE[9], label=paste(times$min[9]," to ",times$max[9],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[10]+100), y=times$PROJECTCODE[10], label=paste(times$min[10]," to ",times$max[10],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[11]+100), y=times$PROJECTCODE[11], label=paste(times$min[11]," to ",times$max[11],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[12]+100), y=times$PROJECTCODE[12], label=paste(times$min[12]," to ",times$max[12],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[13]+100), y=times$PROJECTCODE[13], label=paste(times$min[13]," to ",times$max[13],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[14]+100), y=times$PROJECTCODE[14], label=paste(times$min[14]," to ",times$max[14],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[15]+100), y=times$PROJECTCODE[15], label=paste(times$min[15]," to ",times$max[15],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[16]+100), y=times$PROJECTCODE[16], label=paste(times$min[16]," to ",times$max[16],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[17]+100), y=times$PROJECTCODE[17], label=paste(times$min[17]," to ",times$max[17],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[18]+100), y=times$PROJECTCODE[18], label=paste(times$min[18]," to ",times$max[18],sep=""), color="red", size=3.5) +
  annotate(geom="text", x=(times$Date[19]+100), y=times$PROJECTCODE[19], label=paste(times$min[19]," to ",times$max[19],sep=""), color="red", size=3.5)


## The count of previously marked and new snakes per project
caps <- subcap
caps[caps$NEW == "NEW","NEW"] <- 1
caps[caps$NEW == "OLD","NEW"] <- 0
caps2 <- aggregate(data=caps, as.numeric(NEW) ~ PROJECTCODE, function(x) sum(x))
caps2$TYPE <- c("NEW")
caps <- subcap
caps[caps$NEW == "NEW","NEW"] <- 0
caps[caps$NEW == "OLD","NEW"] <- 1
caps3 <- aggregate(data=caps, as.numeric(NEW) ~ PROJECTCODE, function(x) sum(x))
caps3$TYPE <- c("OLD")
caps <- rbind(caps2,caps3)
colnames(caps) <- c("PROJECTCODE","COUNT","TYPE")

cplot <- ggplot(caps, aes(x=COUNT, y=PROJECTCODE, fill=factor(TYPE))) + geom_bar(stat = "identity", position="dodge")


## The number of unique individuals per project
snks <- aggregate(data=subcap, PITTAG ~ PROJECTCODE, function(x) length(unique(x)))

## The projects under which each snake was caught
proj <- aggregate(data=subcap, PROJECTCODE ~ PITTAG, function(x) length(unique(x)))  ## snakes caught in 1.6 projects on av.
proj2 <- aggregate(data=subcap, PROJECTCODE ~ PITTAG, function(x) unique(x))  ## list of projects each caught in
projsn <- proj2 %>% unnest(PROJECTCODE) %>% group_by(PITTAG) %>% mutate(col=seq_along(PROJECTCODE)) %>% spread(key=unique(PROJECTCODE), value=PROJECTCODE)

## Number of times each individual was caught per project
csnk <- count(subcap, c("PITTAG","PROJECTCODE"))  ## snakes caught about 3 times per project (though some projects will be grouped so increase captures)

splot <- ggplot(subset(subcap, (SITE == "HMUR" | SITE == "HMUI") & PROJECTCODE == "EDGE EFFECT VIS"), aes(x=PITTAG)) +
  geom_histogram(stat="count") + theme(axis.text.x = element_text(angle=90))

sallplot <- ggplot(subcap, aes(x=PITTAG)) +
  geom_histogram(stat="count") + theme(axis.text.x = element_text(angle=90))

## Number of times each individual was caught over all projects
asnk <- count(subcap, c("PITTAG"))  ## snakes caught about 5 times


