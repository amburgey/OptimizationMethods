### Import observations
### Overlay 16m x 8m grid over area
### Assign animal observations to grid cell
### Determine where survey started and ended
### Record which grid cells were surveyed per evening

library(maptools)

## Snake dataframe (PITTAG by Date)
# hmuEsnks <- as.data.frame(caps %>% 
#   filter(SITE == c("HMUI", "HMUR") & PROJECTCODE == c("EDGE EFFECT VIS")) %>%
#   dplyr::select(PITTAG,Date,TRANSECT, LOCATION, seen) %>% 
#   pivot_wider(names_from = Date, values_from = seen, values_fill = 0))  ## add traps

#### HMU ####

PrepDat <- function(hmuEdge, subsurv, hmuSpace){

  HMUtran <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/HMU_EElocations.csv")
  # coordinates(HMUtran) =~ UTME + UTMN
  # proj4string(HMUtran) <- CRS("+proj=utm +zone=55")
  # raster::extract(hmuSpace, HMUtran)  ## check that all points lie somewhere on raster within HMU (shifted two points that were slighlty outside)
  ## Create transects lines between start and end points
  hmutran <- list()
  hmutran[[1]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[1,4],HMUtran[2,4]),c(HMUtran[1,5],HMUtran[2,5]))), ID = "HE01")))
  hmutran[[2]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[3,4],HMUtran[4,4]),c(HMUtran[3,5],HMUtran[4,5]))), ID = "HE02")))
  hmutran[[3]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[5,4],HMUtran[6,4]),c(HMUtran[5,5],HMUtran[6,5]))), ID = "HE03")))
  hmutran[[4]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[7,4],HMUtran[8,4]),c(HMUtran[7,5],HMUtran[8,5]))), ID = "HE04")))
  hmutran[[5]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[9,4],HMUtran[10,4]),c(HMUtran[9,5],HMUtran[10,5]))), ID = "HE05")))
  hmutran[[6]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[11,4],HMUtran[12,4]),c(HMUtran[11,5],HMUtran[12,5]))), ID = "HE06")))
  hmutran[[7]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[13,4],HMUtran[14,4]),c(HMUtran[13,5],HMUtran[14,5]))), ID = "HE07")))
  hmutran[[8]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[15,4],HMUtran[16,4]),c(HMUtran[15,5],HMUtran[16,5]))), ID = "HE08")))
  hmutran[[9]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[17,4],HMUtran[18,4]),c(HMUtran[17,5],HMUtran[18,5]))), ID = "HE09")))
  hmutran[[10]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[19,4],HMUtran[20,4]),c(HMUtran[19,5],HMUtran[20,5]))), ID = "HP01")))
  hmutran[[11]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[21,4],HMUtran[22,4]),c(HMUtran[21,5],HMUtran[22,5]))), ID = "HP02")))
  hmutran[[12]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[23,4],HMUtran[24,4]),c(HMUtran[23,5],HMUtran[24,5]))), ID = "HP03")))
  hmutran[[13]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[25,4],HMUtran[26,4]),c(HMUtran[25,5],HMUtran[26,5]))), ID = "HP04")))
  hmutran[[14]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[27,4],HMUtran[28,4]),c(HMUtran[27,5],HMUtran[28,5]))), ID = "HP05")))
  hmutran[[15]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[29,4],HMUtran[30,4]),c(HMUtran[29,5],HMUtran[30,5]))), ID = "HP06")))
  hmutran[[16]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[31,4],HMUtran[32,4]),c(HMUtran[31,5],HMUtran[32,5]))), ID = "HP07")))
  hmutran[[17]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[33,4],HMUtran[34,4]),c(HMUtran[33,5],HMUtran[34,5]))), ID = "HP08")))
  hmutran[[18]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[35,4],HMUtran[36,4]),c(HMUtran[35,5],HMUtran[36,5]))), ID = "HP09")))
  
  
  ## Import capture locations and convert to UTM
  hmuEdge <- hmuEdge[order(hmuEdge$TRANSECT, hmuEdge$LOCATION),]
  hmucaps <- hmuEdge
  coordinates(hmucaps)=~CAPLON+CAPLAT
  proj4string(hmucaps) <- CRS("+proj=longlat +datum=WGS84")
  ## Convert points to UTM
  hmucaps <- spTransform(hmucaps, CRS("+proj=utm +zone=55"))
  
  ## snap these capture locations to nearest transect line to correct for GPS inaccuracy
  names <- c("HE01","HE02","HE03","HE04","HE05","HE06","HE07","HE08","HE09","HP01","HP02","HP03","HP04","HP05","HP06","HP07","HP08","HP09")
  
  snaplist <- list()
  for (i in 1:length(hmutran)){
    snaplist[[i]] <- as.data.frame(snapPointsToLines(subset(hmucaps, TRANSECT == names[i]), hmutran[[i]]))
  }
  
  allsnap <- do.call(rbind,snaplist)
  hmucaps <- cbind(as.data.frame(hmucaps),allsnap)
  coordinates(hmucaps)=~X+Y  ## use snapped utms
  proj4string(hmucaps) <- CRS("+proj=utm +zone=55")
  
  ## Plot onto raster and figure out cell IDs for all captures
  nr <- dim(HMUhabmat)[1]
  nc <- dim(HMUhabmat)[2]
  hmuEdge$cellIDs <- cellFromXY(hmuSpace, hmucaps)
  hmuEdge$row <- NA
  hmuEdge$col <- NA
  for(i in 1:nrow(hmuEdge)){
    hmuEdge[i,(ncol(hmuEdge)-1)] <- floor(hmuEdge[i,ncol(hmuEdge)-2]/nc) # get row
    hmuEdge[i,ncol(hmuEdge)] <- hmuEdge[i,ncol(hmuEdge)-2] - (floor(hmuEdge[i,ncol(hmuEdge)-2]/nc)*nc) # get col
  }
  
  ## Where did each survey start and end
  hmuEsurv <- tibble(subsurv) %>% 
    filter(SITE == c("HMUI","HMUR") & PROJECTCODE == "EDGE EFFECT VIS") %>%
    dplyr::select(Date, TRANSECT, STARTNUMBER, DISTANCE)
  
  
  ## Create and fill list for all HMU cells that each transect passess through
  cellIDs <- vector("list", length(unique(HMUtran[,2])))
  names(cellIDs) <- unique(HMUtran[,2])
  
  for(i in 1:18){
    cellIDs[[i]] <- unlist(cellFromLine(hmuSpace, hmutran[[i]]))
  }
  
  ## Convert HMU cell IDs to row and column
  hmuvis <- as.data.frame(matrix(NA, ncol = 3))
  
  key <- 0
  
  for(i in 1:length(cellIDs)){
    for(j in 1:length(cellIDs[[i]])){
      hmuvis[key+j,1] <- names(cellIDs)[[i]]
      hmuvis[key+j,2] <- floor(cellIDs[[i]][j]/nc) # get row
      hmuvis[key+j,3] <- cellIDs[[i]][j] - (floor(cellIDs[[i]][j]/nc)*nc) # get col
    } 
    key <- nrow(hmuvis)
  }
  
  hmuvis <- hmuvis %>%
    # dplyr::rename(TRANSECT = V1, row = V2, col = V3) %>%
    group_by(V1) %>%
    dplyr::mutate(row_number = dplyr::row_number())
  colnames(hmuvis) <- c("TRANSECT","row","col","Site")
  hmuvis$LOC <- paste(hmuvis$TRANSECT,hmuvis$Site,sep=".")
  locs <- cbind(hmuvis[,5],hmuvis[,2:3])
  
  ## Find Transect-Site location of captured snakes
  snks <- merge(hmuvis, hmuEdge[,c(1:3,8:9)], by=c("TRANSECT","row","col"))
  snks$LOC <- paste(snks$TRANSECT,snks$Site,sep=".")
  snks$Cap <- 1
  allsnks <- reshape::cast(snks, PITTAG ~ LOC, value = "Cap", fun.aggregate = "length")
  miss <- anti_join(hmuvis,snks,by=c("TRANSECT","Site","LOC","row","col"))  ## grid cells that never had a capture
  allsnks[miss$LOC] <- 0                    # Add missing grid cells, filled with '0's
  tags <- allsnks[,1]
  allsnks <- allsnks[hmuvis$LOC]
  allsnks <- add_column(allsnks, PITTAG = tags, .after = 0)
  
  ## Create a dataframe to indicate when each transect (and cells within each transect) were "active" (aka, surveyed)
  hmuEsurv$Act <- 1
  act <- reshape::cast(hmuEsurv, TRANSECT ~ Date, value = "Act", fun.aggregate = "length")
  act[,2:70] <- as.numeric(act[,2:70]>0)  ## double surveyor days = 2 so simplify everything that's active to 1
  # act$LOC <- paste(act$TRANSECT,1,sep=".")
  
  ## Just remove X from before column headers below so easier to see
  destroyX = function(es) {
    f = es
    for (col in c(1:ncol(f))){ #for each column in dataframe
      if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
        colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
      }
    }
    assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
  }
  
  tempdat <- as.data.frame(matrix(NA, ncol = 71, nrow = nrow(hmuvis), dimnames = list(1:nrow(hmuvis), c("TRANSECT","LOC",colnames(act[,2:70])))))
  tempdat <- destroyX(tempdat)
  key <- 0
  for (i in 1:nrow(act)){
    tra <- subset(hmuvis, TRANSECT == act[i,1])
    for (j in 1:nrow(tra)){
      tempdat[j+key,1] <- tra[j,1]
      tempdat[j+key,2] <- tra[j,5]
      tempdat[j+key,3:71] <- act[i,2:70]
    }
    key <- sum(!is.na(tempdat[,1]))
  }
  
  tempdat$LOC <- factor(tempdat$LOC, levels=c(hmuvis$LOC))
  tempdat <- tempdat[order(tempdat$LOC),]
  
  ## CORRECT HP08, Dist = 0.02 (only did 20 m starting at 200 m side), TRANSID == 45711 & EFFORTID == 9576 by hand on 2013-05-13
  # test <- tempdat$`2013-05-13`[which(tempdat$`2013-05-13` & tempdat$TRANSECT == "HP08" & tempdat$LOC != )]
  
  ## Plot and check, pick a random transect to check (e.g., cellIDs[[6]])
  # xyFromCell(hmuSpace, 64452)
  # xyFromCell(hmuSpace, 67207)
  # xyFromCell(hmuSpace, 70237)
  # 
  # plot(hmuSpace)
  # plot(hmutran[[6]], add=TRUE, cex=5)
  # points(x=268823.7, y=1503811, type="p", cex=1)
  # points(x=268844.2, y=1503769, type="p", cex=1)
  # points(x=268864.8, y=1503723, type="p", cex=1)
  
  hmuEdgedat <- list(act = tempdat, locs = locs, caps = allsnks)
  
  return(hmuEdgedat)


}


