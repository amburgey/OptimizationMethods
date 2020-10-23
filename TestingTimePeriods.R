### USING CP, TRY OUT DIFFERENT PERIODS OF SAMPLING TO SEE WHEN ESTIMATES PLATEAU AND VARIANCE ASYMPTOTES ###
### THIS WILL HELP SET THE WINDOW OF TIME THAT WE CAN LIMIT SAMPLING TO SO WE DON'T NEED TO INCLUDE A CORRELATED ERROR STRUCTURE ###

rm(list=ls())

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
snks <- reshape2::dcast(data = eff, formula = PITTAG ~ Date2, fun.aggregate = length, value.var = "LOCATION")[-(length(unique(subcap$PITTAG))+1),]  # adds an extra row


#### SPECIFY STATE SPACE ####
# Create trapping grid of CP dimensions (5 ha, 50,000 m2)
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
locs$TRANID <- ord
pts <- as.matrix(locs)

## Number of transect grid locations
J <- nrow(pts)

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

