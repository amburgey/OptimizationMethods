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




PrepDat2 <- function(){
  
  NCRIEEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRI EDGE EFFECT VIS.kml","NCRI EDGE EFFECT VIS", require_geomType = "wkbLineString")
  NCRREEline <- rgdal::readOGR("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/NCRR EDGE EFFECT VIS.kml","NCRR EDGE EFFECT VIS", require_geomType = "wkbLineString")
  
}
















plot(hmuSpace)
## 1318 long, 491 wide, 55 ha but if just making rectangle should be ~62 ha
## Grid should be 8 m by 16 m

cellsize = c(16,8)
bbox <- st_sfc(st_polygon(list(rbind(c(xmin(hmuSpace)+250,ymin(hmuSpace)-80), c(xmax(hmuSpace)-250,ymin(hmuSpace)-80), c(xmax(hmuSpace)-250,ymax(hmuSpace)-90), c(xmin(hmuSpace)+250,ymin(hmuSpace)-80)))))
grd <- sf::st_make_grid(bbox, cellsize = cellsize, square = TRUE)
rotang = -27
rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
grd_rot <- (grd - st_centroid(st_union(grd))) * rot(rotang * pi / 180) +
  st_centroid(st_union(grd))
plot(grd_rot, add=TRUE)

## Add CRS to grd_rot
grd_rot <- grd_rot %>%
  st_set_crs(st_crs("+proj=utm +zone=55"))

## Get centroids of grid cells
grd_cts <- grd_rot %>%
  st_set_crs(st_crs("+proj=utm +zone=55")) %>%
  st_centroid()


HMUtran <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/HMU_EElocations.csv")
# colnames(HMUtran)[7:8] <- c("UTME","UTMN")
# coordinates(HMUtran) =~ UTME + UTMN
# proj4string(HMUtran) <- CRS("+proj=utm +zone=55")
# raster::extract(hmuSpace, HMUtran)  ## check that all points lie somewhere on raster within HMU (shifted two points that were slighlty outside)
## Create transects lines between start and end points

# Create list of simple feature geometries (linestrings)
# l_sf <- vector("list", nrow(HMUtran))
# for (i in seq_along(l_sf)){
#   l_sf[[i]] <- st_linestring(as.matrix(rbind(HMUtran[i,4:5], HMUtran[i,7:8])))
# }
# Create simple feature geometry list column
# l_sfc <- st_sfc(l_sf, crs = "+proj=utm +zone=55")

hmutran <- list()
hmutran[[1]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[1,4],HMUtran[1,7]),c(HMUtran[1,5],HMUtran[1,8]))), ID = "HE01")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[2]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[2,4],HMUtran[2,7]),c(HMUtran[2,5],HMUtran[2,8]))), ID = "HE02")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[3]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[3,4],HMUtran[3,7]),c(HMUtran[3,5],HMUtran[3,8]))), ID = "HE03")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[4]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[4,4],HMUtran[4,7]),c(HMUtran[4,5],HMUtran[4,8]))), ID = "HE04")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[5]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[5,4],HMUtran[5,7]),c(HMUtran[5,5],HMUtran[5,8]))), ID = "HE05")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[6]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[6,4],HMUtran[6,7]),c(HMUtran[6,5],HMUtran[6,8]))), ID = "HE06")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[7]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[7,4],HMUtran[7,7]),c(HMUtran[7,5],HMUtran[7,8]))), ID = "HE07")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[8]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[8,4],HMUtran[8,7]),c(HMUtran[8,5],HMUtran[8,8]))), ID = "HE08")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[9]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[9,4],HMUtran[9,7]),c(HMUtran[9,5],HMUtran[9,8]))), ID = "HE09")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[10]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[10,4],HMUtran[10,7]),c(HMUtran[10,5],HMUtran[10,8]))), ID = "HP01")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[11]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[11,4],HMUtran[11,7]),c(HMUtran[11,5],HMUtran[11,8]))), ID = "HP02")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[12]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[12,4],HMUtran[12,7]),c(HMUtran[12,5],HMUtran[12,8]))), ID = "HP03")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[13]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[13,4],HMUtran[13,7]),c(HMUtran[13,5],HMUtran[13,8]))), ID = "HP04")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[14]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[14,4],HMUtran[14,7]),c(HMUtran[14,5],HMUtran[14,8]))), ID = "HP05")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[15]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[15,4],HMUtran[15,7]),c(HMUtran[15,5],HMUtran[15,8]))), ID = "HP06")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[16]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[16,4],HMUtran[16,7]),c(HMUtran[16,5],HMUtran[16,8]))), ID = "HP07")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[17]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[17,4],HMUtran[17,7]),c(HMUtran[17,5],HMUtran[17,8]))), ID = "HP08")), proj4string = CRS("+proj=utm +zone=55"))
hmutran[[18]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[18,4],HMUtran[18,7]),c(HMUtran[18,5],HMUtran[18,8]))), ID = "HP09")), proj4string = CRS("+proj=utm +zone=55"))


## Import capture locations and convert to UTM
hmucaps <- hmuEdge[order(hmuEdge$TRANSECT, hmuEdge$LOCATION),]
coordinates(hmucaps)=~CAPLON+CAPLAT
proj4string(hmucaps) <- CRS("+proj=longlat +datum=WGS84")
## Convert points to UTM
hmucaps <- spTransform(hmucaps, CRS("+proj=utm +zone=55"))
# hmucaps <- st_as_sf(as.data.frame(hmucaps), coords = c("CAPLON","CAPLAT"))
# hmucaps <- hmucaps %>%
#   st_set_crs(st_crs("+proj=utm +zone=55"))

## snap these capture locations to nearest transect line to correct for GPS inaccuracy
names <- c("HE01","HE02","HE03","HE04","HE05","HE06","HE07","HE08","HE09","HP01","HP02","HP03","HP04","HP05","HP06","HP07","HP08","HP09")

snaplist <- list()
for (i in 1:length(hmutran)){
  snaplist[[i]] <- as.data.frame(snapPointsToLines(subset(hmucaps, TRANSECT == names[i]), hmutran[[i]]))
}

# snaplist <- list()
# for (i in 1:length(l_sf)){
#   snaplist[[i]] <- as.data.frame(st_nearest_points(subset(hmucaps$geometry, TRANSECT == names[i]), l_sf[[i]], tolerance=8))
# }

allsnap <- do.call(rbind,snaplist)
hmucaps <- cbind(as.data.frame(hmucaps),allsnap)
coordinates(hmucaps)=~X+Y  ## use snapped utms
proj4string(hmucaps) <- CRS("+proj=utm +zone=55")

## use drawExtent() to find spot on hmuSpace
plot(grd_rot, xlim=c(268614.2, 268677.9), ylim=c(1504381,1504436))
plot(grd_cts, xlim=c(268614.2, 268677.9), ylim=c(1504381,1504436), add=TRUE, pch=21, col="blue", cex=0.5)
plot(l_sf[[12]], add=TRUE, col="red", lty=2)
plot(hmucaps, add=TRUE, cex=1, pch=20)

## Create and fill list for all HMU cells that each transect passess through
## Find proportion of line that goes through the cell to scale effort (e.g., if only 50% line goes through cell, than effort in that cell should only be 50% of regular)

## Lines have to be linestrings instead of spatiallines
## Transform spatiallines to dataframe
datlin <- as.data.frame(matrix(NA, nrow = length(hmutran), ncol = 5))
for(i in 1:8){  ## the N-S oriented transects need to be done separately from the E-W oriented transects
  datlin[i,1] <- names[i]
  datlin[i,2] <- xmin(hmutran[[i]])
  datlin[i,3] <- ymax(hmutran[[i]])
  datlin[i,4] <- xmax(hmutran[[i]])
  datlin[i,5] <- ymin(hmutran[[i]])
}

for(i in 9:18){  ## the N-S oriented transects need to be done separately from the E-W oriented transects
  datlin[i,1] <- names[i]
  datlin[i,2] <- xmin(hmutran[[i]])
  datlin[i,3] <- ymin(hmutran[[i]])
  datlin[i,4] <- xmax(hmutran[[i]])
  datlin[i,5] <- ymax(hmutran[[i]])
}

## Convert data frame to linestrings
out <- pmap(datlin[-1], ~c(...) %>%
              matrix(., , ncol=2, byrow = TRUE) %>% 
              st_linestring)

test <- st_as_sfc(out, crs = "+proj=utm +zone=55")

grd_rotID <- sf::st_sf(grd_rot, 'ID' = seq(length(grd_rot)), grd_rot)

test2 <- grd_rotID[which(st_intersects(test[1], grd_rotID, sparse = FALSE)), ]



