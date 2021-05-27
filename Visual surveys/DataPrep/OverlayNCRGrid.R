## Code to specify grid locations for NCR (Naval Computer Telecommunications Area Master Station) visual surveys

library(rgeos); library(maptools); library(spatialEco)


overlayNCR <- function(NCRcaps, cellsize){
  # Section 1. Load coordinates for NCRI/NCRR Edge Effects project ----
  ## Use kml transects and convert to UTM
  edgeline <- rgdal::readOGR("Data/NCRR EDGE EFFECT VIS.kml","NCRR EDGE EFFECT VIS", require_geomType = "wkbLineString")
  edgeline <- spTransform(edgeline, CRS("+proj=utm +zone=55"))
  intline <- rgdal::readOGR("Data/NCRI EDGE EFFECT VIS.kml","NCRI EDGE EFFECT VIS", require_geomType = "wkbLineString")
  intline <- spTransform(intline, CRS("+proj=utm +zone=55"))
  ## Combine all transects into one object
  lines <- rbind(edgeline,intline)
  lines <- lines[,!(names(lines) %in% c("Description"))]
  ## Convert to dataframe in order to separate lines
  df <- as.data.frame(geom(lines))
  ## Name transects to simpler names that match capture locations
  new_IDs <- rep(paste0(c(rep("RE0", times=9),rep("RP0", times=9)), 1:9), each=2)
  df$Transect <- new_IDs
  ## Convert all to spatial lines
  lines2 <- list()
  key <- 1
  for(i in 1:length(unique(df$Transect))){
    lines2[[i]] <- SpatialLines(list(Lines(Line(cbind(c(df[key,4],df[key+1,4]),c(df[key,5],df[key+1,5]))), ID = new_IDs[i])), proj4string = CRS("+proj=utm +zone=55"))
    key <- key + 2
  }
  
  ## Some transects needed to be replaced with new locations due to construction
  missing <- read.csv("Data/NCRI_missing.csv")
  coordinates(missing) =~ Lon + Lat
  proj4string(missing) <- CRS("+proj=longlat +datum=WGS84")
  ## Convert missing points to UTM and create lines between beginning and end points
  M <- as.data.frame(spTransform(missing, CRS("+proj=utm +zone=55")))
  miss <- list()
  new_IDs <- paste0("RP", 10:14)
  key <- 1
  for(i in 1:length(unique(M$Transect))){
    miss[[i]] <- SpatialLines(list(Lines(Line(cbind(c(M[key,3],M[key+1,3]),c(M[key,4],M[key+1,4]))), ID = new_IDs[i])), proj4string = CRS("+proj=utm +zone=55"))
    key <- key + 2
  }
  ## Combine all transects into one object
  lines <- c(lines2,miss)
  # ## Create combined dataframe of all transects for later use
  # M2 <- cbind(as.data.frame(matrix(c(rep(19:23, each=2),rep(1, times=10),rep(19:23, each=2)),nrow=10, ncol = 3)),M[,3:4],M[,1]); colnames(M2) <- c("object","part","cump","x","y","Transect")
  # dfall <- rbind(df,M2)
  # write.csv(dfall,"Data/NCRCoordinates.csv")
  tran <- read.csv("Data/NCRCoordinates.csv")[,-1]
  
  
  # Section 2. Create Grid Polygon
  ## Grid should be 8 m by 16 m to match surveying grid of Closed Population
  cellsize1 = c(8,16)
  ## Overlay a grid of these dimensions across the space of the NCR (but with room for rotation) and then rotate
  ## drawExtent()
  ## What grid size is good to use? Probably the same as the state space
  bbox <- st_sfc(st_polygon(list(rbind(c(267700,1503020), c(269000,1503020), c(269000,1504763), c(267700,1503020)))))
  grd <- sf::st_make_grid(bbox, cellsize = cellsize1, square = TRUE)
  rotang = -26.5
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  grd_rot <- (grd - st_centroid(st_union(grd))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd))
  ## Specify projection again
  st_crs(grd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  grd_rot <- as_Spatial(grd_rot, cast=TRUE, IDs=paste0("ID", seq_along(grd_rot)))
  plot(grd_rot, col="red")
  
  
  # Section 3. Convert Transect Survey Locations to Grid Cells ----
  ## Make sure all transects pass through cells
  # test <- c()
  # for(i in 1:length(lines)){
  #   test[i] <- gIntersects(lines[[i]],grd_rot)
  # }
  # unique(test)  ## takes a while to do so only do when checking first time
  
  ## Find which grid cells a transect passes through
  ## Lines have to be linestrings instead of spatial lines
  ## Convert data frame to linestrings
  ltran <- pmap(tran[,-c(1:3,8)], ~c(...) %>%
                  matrix(., , ncol=2, byrow = TRUE) %>% 
                  st_linestring)
  
  ltran <- st_as_sfc(ltran, crs = "+proj=utm +zone=55")
  
  ## Convert to sf object and assign ID to each grid cell
  sfNCR <- sf::st_as_sf(grd_rot)
  grd_rotID <- sf::st_sf(sfNCR$geometry, 'ID' = seq(length(sfNCR$geometry)), sfNCR$geometry)
  
  ## Identify each grid cell that each line goes through, don't do in for loop so maintains projection
  polylin1 <- grd_rotID[which(st_intersects(ltran[1], grd_rotID, sparse = FALSE)), ]
  polylin2 <- grd_rotID[which(st_intersects(ltran[2], grd_rotID, sparse = FALSE)), ]
  polylin3 <- grd_rotID[which(st_intersects(ltran[3], grd_rotID, sparse = FALSE)), ]
  polylin4 <- grd_rotID[which(st_intersects(ltran[4], grd_rotID, sparse = FALSE)), ]
  polylin5 <- grd_rotID[which(st_intersects(ltran[5], grd_rotID, sparse = FALSE)), ]
  polylin6 <- grd_rotID[which(st_intersects(ltran[6], grd_rotID, sparse = FALSE)), ]
  polylin7 <- grd_rotID[which(st_intersects(ltran[7], grd_rotID, sparse = FALSE)), ]
  polylin8 <- grd_rotID[which(st_intersects(ltran[8], grd_rotID, sparse = FALSE)), ]
  polylin9 <- grd_rotID[which(st_intersects(ltran[9], grd_rotID, sparse = FALSE)), ]
  polylin10 <- grd_rotID[which(st_intersects(ltran[10], grd_rotID, sparse = FALSE)), ]
  polylin11 <- grd_rotID[which(st_intersects(ltran[11], grd_rotID, sparse = FALSE)), ]
  polylin12 <- grd_rotID[which(st_intersects(ltran[12], grd_rotID, sparse = FALSE)), ]
  polylin13 <- grd_rotID[which(st_intersects(ltran[13], grd_rotID, sparse = FALSE)), ]
  polylin14 <- grd_rotID[which(st_intersects(ltran[14], grd_rotID, sparse = FALSE)), ]
  polylin15 <- grd_rotID[which(st_intersects(ltran[15], grd_rotID, sparse = FALSE)), ]
  polylin16 <- grd_rotID[which(st_intersects(ltran[16], grd_rotID, sparse = FALSE)), ]
  polylin17 <- grd_rotID[which(st_intersects(ltran[17], grd_rotID, sparse = FALSE)), ]
  polylin18 <- grd_rotID[which(st_intersects(ltran[18], grd_rotID, sparse = FALSE)), ]
  polylin19 <- grd_rotID[which(st_intersects(ltran[19], grd_rotID, sparse = FALSE)), ]
  polylin20 <- grd_rotID[which(st_intersects(ltran[20], grd_rotID, sparse = FALSE)), ]
  polylin21 <- grd_rotID[which(st_intersects(ltran[21], grd_rotID, sparse = FALSE)), ]
  polylin22 <- grd_rotID[which(st_intersects(ltran[22], grd_rotID, sparse = FALSE)), ]
  polylin23 <- grd_rotID[which(st_intersects(ltran[23], grd_rotID, sparse = FALSE)), ]
  
  my.list <- list(polylin1, polylin2, polylin3, polylin4, polylin5, polylin6, polylin7, polylin8, polylin9, polylin10, polylin11, polylin12, polylin13, polylin14, polylin15, polylin16, polylin17, polylin18, polylin19, polylin20, polylin21, polylin22, polylin23)
  
  ## Find centroid of all grid cells
  eff_cts <- list()
  for(i in 1:length(my.list)){
    eff_cts[[i]] <- st_centroid(my.list[[i]])  ## Warnings about geometries, ignore
  }
  
  ## Names (in order) of transects
  tnames <- tran$Transect
  names(eff_cts) <- tnames
  
  ## Create a matrix of transect name, grid cell ID, and X-Y coordinate information
  ## Number of rows needed based on elements of the list
  entries <- nrow(eff_cts[[1]]) + nrow(eff_cts[[2]]) + nrow(eff_cts[[3]]) + nrow(eff_cts[[4]]) + nrow(eff_cts[[5]]) + nrow(eff_cts[[6]]) + nrow(eff_cts[[7]]) + nrow(eff_cts[[8]]) + nrow(eff_cts[[9]]) + nrow(eff_cts[[10]]) + nrow(eff_cts[[11]]) + nrow(eff_cts[[12]]) + nrow(eff_cts[[13]]) + nrow(eff_cts[[14]]) + nrow(eff_cts[[15]]) + nrow(eff_cts[[16]]) + nrow(eff_cts[[17]]) + nrow(eff_cts[[18]]) + nrow(eff_cts[[19]]) + nrow(eff_cts[[20]]) + nrow(eff_cts[[21]]) + nrow(eff_cts[[22]]) + nrow(eff_cts[[23]])
  
  vistran <- setNames(as.data.frame(matrix(NA, nrow = entries, ncol = 4)), c("TranID","GridID","x","y"))
  key <- 1
  for(i in 1:length(eff_cts)){
    for(j in 1:nrow(eff_cts[[i]])){
      vistran[key,1] <- names(eff_cts)[i]
      vistran[key,2] <- unlist(map(eff_cts[[i]],j))[1]  
      vistran[key,3] <- unlist(map(eff_cts[[i]],j))[2]
      vistran[key,4] <- unlist(map(eff_cts[[i]],j))[3]
      
      key <- key + 1
    }
  }
  
  vistran <- vistran[order(vistran$GridID),]
  
  
  # Section 4. Convert BTS Captures to Grid Cells ----
  ## Import capture locations (from Select&Prep) and convert to UTM
  ncrcaps <- NCRcaps[order(NCRcaps$TRANSECT, NCRcaps$LOCATION),]
  coordinates(ncrcaps)=~CAPLON+CAPLAT
  ## Specify projection
  proj4string(ncrcaps) <- CRS("+proj=longlat +datum=WGS84")
  ## Convert points to UTM
  ncrcaps <- spTransform(ncrcaps, CRS("+proj=utm +zone=55"))
  ## snap these capture locations to nearest transect line to correct for GPS inaccuracy
  names <- c("RE01","RE02","RE03","RE04","RE05","RE06","RE07","RE08","RE09","RP01","RP02","RP03","RP04","RP05","RP06","RP07","RP08","RP09","RP10","RP11","RP12","RP13","RP14")
  
  snaplist <- list()
  for (i in 1:length(lines)){
    snaplist[[i]] <- as.data.frame(snapPointsToLines(subset(ncrcaps, TRANSECT == names[i]), lines[[i]]))  ## warning about CRS, ignore
  }
  ## Convert from list of dataframes to a single dataframe with transect info
  allsnap <- do.call(rbind,snaplist)
  capncr <- cbind(as.data.frame(ncrcaps),allsnap)
  ## Convert to spatial points dataframe
  spd <- capncr
  coordinates(spd) =~ X + Y
  proj4string(spd) <- CRS("+proj=utm +zone=55")
  capncr$GridID <- NA
  for(i in 1:length(spd)){
    capncr[i,29] <- point.in.poly(spd[i,],grd_rot)@data$poly.ids  ## takes A LONG time because so large a space
  }
  
  ## Check that no grid cells were identified as having snakes captured that weren't part of transect grid cells surveyed
  ifelse(unique(unique(capncr$GridID) %in% unique(vistran$GridID)) == TRUE, "GridID matched", "error: GridID mismatch")
  
  ## Dataframe of snake capture info and the grid cell on which it occurred
  snkcap <- capncr[,c("EFFORTID","PITTAG","SVL","Date","TRANSECT","X","Y","GridID")]
  
  
  # Section 5. Create rotated integration grid ----
  ## Desired grid cell
  cellsize2 = cellsize
  ## Overlay a grid of these dimensions across the space of the HMU (but with room for rotation) and then rotate
  bbox2 <- bbox
  grd2 <- sf::st_make_grid(bbox2, cellsize = cellsize2, square = TRUE)
  rotang = rotang
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  intgrd_rot <- (grd2 - st_centroid(st_union(grd2))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd2))
  ## Specify projection again
  st_crs(intgrd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  intgrd_rot <- as_Spatial(intgrd_rot, cast=TRUE, IDs=paste0("ID", seq_along(intgrd_rot)))
  ## Find centroid of all grid cells
  intgrd_cts <- gCentroid(intgrd_rot, byid = TRUE)
  plot(intgrd_cts, col="red", pch=21, cex=0.2)
  intgrd <- geom(intgrd_cts)
  
  # Check. Plot all transects and grid----
  # par(mar = c(1,1,1,1))
  # plot(lines, col="red")
  # plot(edgeline, col="blue", add=TRUE)
  # plot(intline, col="green", add=TRUE)
  # plot(miss, col="orange", add=TRUE)
  # plot(grd_rot, col="red")
  # for(i in 1:length(ltran)){
  #   plot(ltran[[i]], add=TRUE, col="white", lwd=2)
  # }
  
  # ## Checking. Plot transects up close ----
  # plot(grd_rot, xlim=c(268849,268377.9), ylim=c(1503555,1504715))
  # plot(lines[[1]], add=TRUE, col="red")
  # plot(lines[[2]], add=TRUE, col="green")
  # plot(lines[[3]], add=TRUE, col="blue")
  # plot(lines[[4]], add=TRUE, col="orange")
  # plot(lines[[5]], add=TRUE, col="red")
  # plot(lines[[6]], add=TRUE, col="green")
  # plot(lines[[7]], add=TRUE, col="blue")
  # plot(lines[[8]], add=TRUE, col="orange")
  # plot(lines[[9]], add=TRUE, col="red")
  # t <- extent(lines[[23]])
  # plot(grd_rot, xlim=c(as.numeric(xmin(t))-0.5,as.numeric(xmax(t))+0.5), ylim=c(as.numeric(ymin(t))-5,as.numeric(ymax(t))+5))
  # plot(lines[[23]], add=TRUE, col="red")
  # plot(my.list[[23]]$sfNCR.geometry, add=TRUE, col="red")
  
  
  ## Checking. Plot snapped captures to transects ----
  # t <- extent(lines[[23]])
  # plot(grd_rot, xlim=c(as.numeric(xmin(t))-0.5,as.numeric(xmax(t))+0.5), ylim=c(as.numeric(ymin(t))-5,as.numeric(ymax(t))+5))
  # plot(my.list[[23]]$sfNCR.geometry, add=TRUE, col="red")
  # snaps1 <- snaplist[[23]]
  # coordinates(snaps1) =~ X + Y
  # plot(ncrcaps, add=TRUE, pch=21, cex=0.5) ## original
  # plot(snaps1, add=TRUE, pch=21, cex=0.5, col="yellow") ## adjusted
  
  dat <- list(tran = vistran, snks = snkcap, grd = grd, intgrd = intgrd)
  
  return(dat)
}
