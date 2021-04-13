##### NCR (Naval Computer Telecommunications Area Master Station) is an open (no fence) study area

## NCR was was only surveyed for this Edge Effects project

## Subset down to captures needed for project:
subSnk <- function(SITEcaps){

  ##### READ IN AND FORMAT CAPTURES #####
  ## Order by GridID
  SITEcaps <- SITEcaps[order(SITEcaps$GridID),]
  ## Separate year and month out for organizing
  SITEcaps$YEAR <- format(as.Date(SITEcaps$Date, format="%m/%d/%Y"),"%Y")
  SITEcaps$MONTH <- format(as.Date(SITEcaps$Date, format="%m/%d/%Y"),"%m")
  
  return(SITEcaps)
  
}


## Specify time and subset snake captures to that timeframe
subYr <- function(SITEcaps, time){
  
  capsyr <- subset(SITEcaps, MONTH >= time[1] & MONTH <= time[2])
  capsyr <- capsyr[,c("EFFORTID","Date","PITTAG","TRANSECT","GridID")]
  capsyr <- capsyr[order(capsyr$TRANSECT, capsyr$Date),]
  
  return(capsyr)
}


## Determine the sampling effort (when transects were surveyed)
effSnk <- function(eff, time){
  ##### READ IN TRAPPING EFFORT #####
  ## Separate month out for organizing
  eff$MONTH <- format(as.Date(eff$Date, format="%m/%d/%Y"),"%m")
  ## Subset to timespan
  effyr <- subset(eff, MONTH >= time[1] & MONTH <= time[2])
  effyr <- effyr[,c("EFFORTID","Date","BI","TRANSECT","DISTANCE")]
  effyr <- effyr[order(effyr$TRANSECT, effyr$Date),]
  
  return(effyr)
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
  y <- y[,-1]
  
  ## Find names of missing columns
  Missing <- as.character(setdiff(grid$GridID,colnames(y)))
  
  if(length(Missing)!= 0){
    ## Add them, filled with '0's
    y[Missing] <- 0
    ## Put columns in desired order
    y <- y[as.character(grid$GridID)]  
    ## Prep for model
    y <- as.matrix(y)
  }
  
  ## Transform any multiple captures at a single location to just 1
  y <- ifelse(y>=1,1,y)
  
  ## Set up effort matrix (Grid cell by Date and indicate if active or not)
  ## Check if effort should be scaled due to different survey lengths
  if(length(unique(SCReff$DISTANCE)) != 1) stop('mismatch in effort and capture dimensions')
  
  ## Subset to which transects were done on which dates
  act <- SCReff[,c("Date","TRANSECT")]
  ## Indicate these were the active times
  act$Active <- c(1)
  ## Add missing transect - date combinations
  act <- act %>%
    complete(Date, nesting(TRANSECT), fill = list(Active = 0))
  ## Name grid to match act
  colnames(grid) <- c("TRANSECT","GridID","x","y")
  ## Expand dataframe to be all grid cells and not just at the broad transect level
  allact <- merge(grid, act, by = c("TRANSECT"), all=TRUE)
  allact <- allact[order(allact$TRANSECT,allact$Date,allact$GridID),]
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
  y <- y[,-1]
  
  ## Find names of missing columns
  Missing <- as.character(setdiff(grid$GridID,colnames(y)))
  
  if(length(Missing)!= 0){
    ## Add them, filled with '0's
    y[Missing] <- 0
    ## Put columns in desired order
    y <- y[as.character(grid$GridID)]  
    ## Prep for model
    y <- as.matrix(y)
  }
  
  ## Transform any multiple captures at a single location to just 1
  y <- ifelse(y>=1,1,y)
  
  ## Set up effort matrix (Grid cell by Date and indicate if active or not)
  ## Check if effort should be scaled due to different survey lengths
  if(length(unique(SCReff$DISTANCE)) == 1) stop('no mismatch in effort and capture dimensions')
  
  ## Subset to which transects were done on which dates
  act <- SCReff[,c("Date","TRANSECT")]
  ## Indicate these were the active times
  act$Active <- c(1)
  ## Add missing transect - date combinations
  act <- act %>%
    complete(Date, nesting(TRANSECT), fill = list(Active = 0))
  ## Name grid to match act
  colnames(grid) <- c("TRANSECT","GridID","x","y")
  ## Expand dataframe to be all grid cells and not just at the broad transect level
  allact <- merge(grid, act, by = c("TRANSECT"), all=TRUE)
  allact <- allact[order(allact$TRANSECT,allact$Date,allact$GridID),]
  ## Reshape to be all points by dates and 1=surveyed, 0=not surveyed
  act2 <- reshape2::dcast(allact, GridID ~ Date, fun.aggregate = sum, value.var = "Active")
  ## Two surveyors at each survey so change to 1 and factor in two people later in cost model
  act2 <- cbind(act2[,1], act2[,2:ncol(act2)] %>% mutate_if(is.numeric, ~1 * (. > 0))); colnames(act2)[1] <- c("GridID")
  
  ##### DO THIS STEP MANUALLY, HAVE TO SET WHICH POINTS BASED ON STARTINGNUMBER AND DISTANCE TRAVELED
  act2[act2$GridID == "2849",7] <- 0
  act2[act2$GridID == "2868",7] <- 0
  act2[act2$GridID == "2869",7] <- 0
  act2[act2$GridID == "2870",7] <- 0 
  act2[act2$GridID == "2871",7] <- 0
  act2[act2$GridID == "2872",7] <- 0
  act2[act2$GridID == "2873",7] <- 0
  act2[act2$GridID == "2874",7] <- 0
  act2[act2$GridID == "2875",7] <- 0
  act2[act2$GridID == "2876",7] <- 0
  
  both <- list(y=y,act=act2)
  
  return(both)
}