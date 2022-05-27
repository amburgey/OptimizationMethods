##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

## Subset down to captures needed for project:
subSnk <- function(SITEcaps, type, info){

  ##### READ IN AND FORMAT CAPTURES #####
  ## Remove edge (SWE, SE, NWE, NEE, SW, PR) and one-sided (letter-0 and letter-14) captures
  SITEcaps <- droplevels(subset(SITEcaps, TRANSECT %in% c(LETTERS,"AA")))
  SITEcaps <- droplevels(subset(SITEcaps, LOCATION %in% c(1:13)))
  ## Set levels for correct ordering
  SITEcaps$Point <- as.factor(SITEcaps$Point)
  SITEcaps$Point <- factor(SITEcaps$Point, levels = paste(rep(c(LETTERS,"AA"),each=13),rep(1:13,times=27),sep=""))
  ## Order by Point (composite of TRANSECT and LOCATION)
  SITEcaps <- SITEcaps[order(SITEcaps$Point),]
  ## Separate year and month out for organizing
  SITEcaps$YEAR <- format(as.Date(SITEcaps$Date, format="%m/%d/%Y"),"%Y")
  SITEcaps$MONTH <- format(as.Date(SITEcaps$Date, format="%m/%d/%Y"),"%m")
    
  SITEcaps2 <- matrix(NA)
  
  ## Make sure only using visual survey data (should be determined by PROJECTCODE but double check)
  if(type == c("TRAPTYPE")){
    if(info == c("V")){
      ## Subset captures to those from traps only, no visual surveys
      SITEcaps2 <- subset(SITEcaps, TRAPTYPE == "V")
    }
  }
  
  return(SITEcaps2)
  
}

## Specify time and subset snake captures to that timeframe
subYr <- function(SITEcaps, time){
  
  capsyr <- subset(SITEcaps, MONTH >= time[1] & MONTH <= time[2])
  capsyr <- capsyr[,c("EFFORTID","Date","PITTAG","Point")]
  
  return(capsyr)
}


## Determine the sampling effort (when transects were surveyed)
effSnk <- function(eff, time){
  ##### READ IN TRAPPING EFFORT #####
  ## Remove edge (SWE, SE, NWE, NEE, SW, PR) and one-sided (letter-0 and letter-14) captures
  eff <- droplevels(subset(eff, TRANSECT %in% c(LETTERS,"AA")))
  ## Set levels for correct ordering
  eff$TRANSECT <- as.factor(eff$TRANSECT)
  eff$TRANSECT <- factor(eff$TRANSECT, levels = c(LETTERS,"AA"))
  ## Separate month out for organizing
  eff$MONTH <- format(as.Date(eff$Date, format="%m/%d/%Y"),"%m")
  ## Subset to timespan
  effyr <- subset(eff, MONTH >= time[1] & MONTH <= time[2])
  effyr <- effyr[,c("EFFORTID","Date","TRANSECT","DISTANCE")]
  effyr <- effyr[order(effyr$TRANSECT, effyr$Date),]
  
  return(effyr)
}


##### CHECK FOR DUPLICATE SURVEYS OR SNAKES THAT SHOULD NOT HAVE BEEN CAPTURED TWICE IN THE SAME NIGHT.----

checkSnks <- function(SCRcaps){
  ## Check if there are erroneous duplicates
  if(length(unique(duplicated(SCRcaps[,-2]))) >= 2) stop('fix duplicates')
  
  return(SCRcaps)
}


##### COMPARE DATES BETWEEN CAPTURES AND EFFORT #####
## Make sure that the captures of snakes match the dates when surveys occurred
checkDims <- function(SCRseff, SCRcaps){
  x <- unique(SCReff[,c("Date","EFFORTID")])
  y <- unique(SCRcaps[,c("Date","EFFORTID")])
  
  new <- merge(x, y, by=c("EFFORTID","Date"), all = TRUE)
  
  if(nrow(new) != nrow(unique(SCReff[,c("Date","EFFORTID")]))) stop('mismatch in effort and capture dimensions')

}


##### SETUP DATA FOR SCR ANALYSIS FORMAT #####
prepSCR <- function(SCRcaps, SCReff, grid){
  
  ## set up observations as count of individuals (rows) and traps (columns) over study
  y <- dcast(data=SCRcaps, formula=PITTAG ~ GridID, length, fill=0, value.var = "GridID")
  ## Make sure snakes ordered by PITTAG in order to match to size
  y <- y[order(y$PITTAG),]
  y <- y[,-1]
  ## Find names of missing columns
  Missing <- setdiff(as.character(unique(grid$GridID)),colnames(y))
  
  if(length(Missing)!= 0){
    ## Add them, filled with '0's
    y[Missing] <- 0
    ## Put columns in desired order
    y <- y[as.character(grid$GridID)]  
    ## Prep for model
    y <- as.matrix(y)
  }
  
  ## Set up effort matrix (Grid cell by Date and indicate if active or not)
  ## Check if effort should be scaled due to different survey lengths
  if(length(unique(SCReff$DISTANCE)) != 1) stop('mismatch in effort and capture dimensions')
  
  ## Subset to which transects were done on which dates
  act <- SCReff[,c("Date","TRANSECT")]
  ## Indicate these were the active times
  act$Active <- c(1)
  ## Create dataframe of all transect-location combinations
  allact <- melt(data.frame(rep(rev(unique(act$TRANSECT)), each=13), unique(fullX$Point))
                 , id.vars = c("rep.rev.unique.act.TRANSECT....each...13.", "unique.fullX.Point."))
  colnames(allact) <- c("TRANSECT","Point")
  ## Add GridID to Points
  allact <- merge(allact, grid[,1:2], by=c("Point"))
  ## Expand dataframe to be all points and not just at the broad transect level
  allact <- merge(act, allact, by = c("TRANSECT"))
  ## Reshape to be all points by dates and 1=surveyed, 0=not surveyed
  act2 <- reshape2::dcast(allact, GridID ~ Date, fun.aggregate = sum, value.var = "Active")
  ## Two surveyors at each survey so change to 1 and factor in two people later in cost model
  act2 <- cbind(act2[,1], act2[,2:ncol(act2)] %>% mutate_if(is.numeric, ~1 * (. > 0))); colnames(act2)[1] <- c("GridID")
  ## Prep for model
  
  both <- list(y=y,act=act2)
  
  return(both)
}


prepSCRman <- function(SCRcaps, SCReff, grid){
  
  ## set up observations as count of individuals (rows) and traps (columns) over study
  y <- dcast(data=SCRcaps, formula=PITTAG ~ GridID, length, fill=0, value.var = "GridID")
  ## Make sure snakes ordered by PITTAG in order to match to size
  y <- y[order(y$PITTAG),]
  y <- y[,-1]
  ## Find names of missing columns
  Missing <- setdiff(as.character(unique(grid$GridID)),colnames(y))
  
  if(length(Missing)!= 0){
    ## Add them, filled with '0's
    y[Missing] <- 0
    ## Put columns in desired order
    y <- y[as.character(grid$GridID)]  
    ## Prep for model
    y <- as.matrix(y)
  }
  
  ## Set up effort matrix (Grid cell by Date and indicate if active or not)
  ## Check if effort should be scaled due to different survey lengths
  if(length(unique(SCReff$DISTANCE)) == 1) stop('no mismatch in effort and capture dimensions')
  
  ## Subset to which transects were done on which dates
  act <- SCReff[,c("Date","TRANSECT")]
  ## Indicate these were the active times
  act$Active <- c(1)
  ## Create dataframe of all transect-location combinations
  allact <- melt(data.frame(rep(rev(unique(act$TRANSECT)), each=13), unique(fullX$Point))
                 , id.vars = c("rep.rev.unique.act.TRANSECT....each...13.", "unique.fullX.Point."))
  colnames(allact) <- c("TRANSECT","Point")
  ## Add GridID to Points
  allact <- merge(allact, grid[,1:2], by=c("Point"))
  ## Expand dataframe to be all points and not just at the broad transect level
  allact <- merge(act, allact, by = c("TRANSECT"))
  ## Reshape to be all points by dates and 1=surveyed, 0=not surveyed
  act2 <- reshape2::dcast(allact, GridID ~ Date, fun.aggregate = sum, value.var = "Active")
  ## Two surveyors at each survey so change to 1 and factor in two people later in cost model
  act2 <- cbind(act2[,1], act2[,2:ncol(act2)] %>% mutate_if(is.numeric, ~1 * (. > 0))); colnames(act2)[1] <- c("GridID")
  
  
  both <- list(y=y,act=act2)
  
  return(both)
}


getSize <- function(capPROJ, SCRcaps, subcap){
  
  ## Find body sizes from current dataset
  bod <- capPROJ[,c("PITTAG","SVL")]
  ## Eliminate ones not retained for analysis above
  bod <- bod[(bod$PITTAG %in% SCRcaps$PITTAG),]
  ## Find mean body size in case multiple measures taken during this study
  bod <- aggregate(bod[,c("SVL")],list(bod$PITTAG),mean,na.rm = TRUE)
  bod <- bod[order(bod$Group.1),]
  if(dim(subset(bod, is.na(x)))[1] != 0) stop("Not all snakes have size")
  
  return(bod)
  
}

getSizeman <- function(capPROJ, SCRcaps, subcap, time){
  
  ## Find body sizes from current dataset
  bod <- capPROJ[,c("PITTAG","SVL")]
  ## Eliminate ones not retained for analysis above
  bod <- bod[(bod$PITTAG %in% SCRcaps$PITTAG),]
  ## Find mean body size in case multiple measures taken during this study
  bod <- aggregate(bod[,c("SVL")],list(bod$PITTAG),mean,na.rm = TRUE)
  
  ## If some snakes are missing body size so see if there is another record from same area and ballpark time frame
  Missing <- subset(bod, is.na(x))
  tf <- subset(subcap, SITE == "NWFN" & Date >= as.Date(time[1]) & Date <= as.Date(time[2]))
  tf <- tf[tf$PITTAG %in% Missing$Group.1,]
  tfbod <- aggregate(tf[,c("SVL")],list(tf$PITTAG),mean,na.rm = TRUE)
  ## When manual solutions are needed
  tfbod[tfbod$Group.1 == "45296F4049" & is.na(tfbod$x), "x"] <- 961.5  ## HL 2
  ## Insert found values
  for(i in 1:nrow(tfbod)){
    bod[bod$Group.1 == tfbod$Group.1[i] & is.na(bod[2]), "x"] <- tfbod[i,2]
  }
  bod <- bod[order(bod$Group.1),]
  ## If error triggered again then...other creative solutions to follow
  if(dim(subset(bod, is.na(x)))[1] != 0) stop("Not all snakes have size")
  
  return(bod)
}
