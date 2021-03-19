## Code to specify grid locations for HMU visual surveys

library(rgeos); library(maptools); library(spatialEco)

overlayHMU <- function(HMUcaps){
  # Section 1 Create HMU Polygon ----
  ## Read in shapefile of HMU fence
  hmuline <- st_read("hmuline/hmu.shp")
  ## Specify projection
  hmuline <- st_transform(hmuline, CRS("+proj=utm +zone=55 +units=m +datum=WGS84"))
  ## Convert from a spatial line to a sfc polygon
  hmuline <- st_cast(hmuline, to = "POLYGON")
  ## Create spatial polygon
  ## get polygon geometry from hmuline but have to copy and paste below
  HMU <- readWKT("POLYGON ((269022.5 1503602, 269066.6 1503608, 269092.1 1503607, 269119.9 1503602, 269139.6 1503597, 269163 1503588, 269177.5 1503579, 269188.4 1503571, 269212 1503549, 269236.4 1503533, 269271.1 1503518, 269290.2 1503514, 269316.2 1503509, 269339.3 1503509, 269359.6 1503511, 269370 1503513, 269389.8 1503518, 269395.9 1503521, 269417.8 1503530, 269472.6 1503556, 269474.3 1503562, 269468.7 1503572, 269347.9 1503819, 269338.7 1503821, 269285.9 1503795, 269280.3 1503796, 269255.7 1503826, 269248.5 1503829, 269239.1 1503831, 269234.5 1503834, 269229.3 1503840, 269225.7 1503848, 269223.8 1503857, 269223.7 1503865, 269224.8 1503875, 269227.3 1503884, 269292.5 1503916, 269294.8 1503925, 269260.1 1503995, 269252.4 1503996, 269174.2 1503957, 269168.9 1503960, 269151.2 1504014, 269156.2 1504036, 269227.8 1504071, 269230 1504078, 269220.4 1504093, 269156.5 1504221, 269066.3 1504405, 269065.6 1504408, 269060.2 1504419, 268963.4 1504614, 268899.9 1504741, 268891.6 1504742, 268834.7 1504707, 268828.5 1504695, 268716.8 1504641, 268469.9 1504520, 268467.3 1504514, 268484.8 1504487, 268499.2 1504461, 268507.7 1504444, 268590.2 1504274, 268686.5 1504081, 268733.5 1503985, 268800.7 1503849, 268850.1 1503748, 268892.6 1503663, 268900.5 1503648, 268953.5 1503580, 268978.7 1503592, 268981.4 1503597, 268971.8 1503623, 268995.9 1503635, 269005.6 1503636, 269010.6 1503632, 269015.8 1503619, 269022.5 1503602))")
  ## Specify projection
  proj4string(HMU) <- CRS("+proj=utm +zone=55 +units=m +datum=WGS84")
  plot(HMU)
  
  
  # Section 2. Create Grid Polygon ----
  ## Grid should be 8 m by 16 m
  hmuSpace <- raster("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/Reclass_hmu_71.tif")
  ## Specify projection
  hmuSpace <- projectRaster(hmuSpace, crs="+proj=utm +zone=55 +units=m +datum=WGS84")
  ## Desired grid cell to match surveying grid of Closed Population
  cellsize = c(16,8)
  ## Overlay a grid of these dimensions across the space of the HMU (but with room for rotation) and then rotate
  bbox <- st_sfc(st_polygon(list(rbind(c(xmin(hmuSpace)+320,ymin(hmuSpace)-95), c(xmax(hmuSpace)-320,ymin(hmuSpace)-95), c(xmax(hmuSpace)-320,ymax(hmuSpace)-95), c(xmin(hmuSpace)+320,ymin(hmuSpace)-95)))))
  grd <- sf::st_make_grid(bbox, cellsize = cellsize, square = TRUE)
  rotang = -26.5
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  grd_rot <- (grd - st_centroid(st_union(grd))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd))
  ## Specify projection again
  st_crs(grd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  grd_rot <- as_Spatial(grd_rot, cast=TRUE, IDs=paste0("ID", seq_along(grd_rot)))
  plot(grd_rot, add=TRUE, col="red")
  
  
  # Section 3. Create HMU-clipped Survey Grid ----
  clipHMU <- gIntersection(HMU, grd_rot, byid = TRUE, drop_lower_td = TRUE)
  # par(mar = c(1,1,1,1))  ## if for some reason the margins are weird
  plot(clipHMU)
  
  
  # Section 4. Convert Transect Survey Locations to Grid Cells ----
  HMUtran <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/HMU_EElocations.csv")
  
  ## Create spatial lines from dataframe
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
  
  ## Make sure all transects pass through cells
  for(i in 1:length(hmutran)){
    test[i] <- gIntersects(hmutran[[i]],clipHMU)
  }
  unique(test)
  
  ## Find which grid cells a transect passes through
  ## Lines have to be linestrings instead of spatial lines
  ## Convert data frame to linestrings
  ltran <- pmap(HMUtran[,-c(1:3,6)], ~c(...) %>%
                matrix(., , ncol=2, byrow = TRUE) %>% 
                st_linestring)
  
  ltran <- st_as_sfc(ltran, crs = "+proj=utm +zone=55")
  
  ## Convert to sf object and assign ID to each grid cell
  sfHMU <- sf::st_as_sf(clipHMU)
  grd_rotID <- sf::st_sf(sfHMU$geometry, 'ID' = seq(length(sfHMU$geometry)), sfHMU$geometry)
  
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
  
  my.list <- list(polylin1, polylin2, polylin3, polylin4, polylin5, polylin6, polylin7, polylin8, polylin9, polylin10, polylin11, polylin12, polylin13, polylin14, polylin15, polylin16, polylin17, polylin18)
  
  ## Find centroid of all grid cells
  eff_cts <- list()
  for(i in 1:length(my.list)){
    eff_cts[[i]] <- st_centroid(my.list[[i]])  ## Warnings about geometries, ignore
  }
  
  ## Names (in order) of transects
  tnames <- c("HE01","HE02","HE03","HE04","HE05","HE06","HE07","HE08","HE09","HP01","HP02","HP03","HP04","HP05","HP06","HP07","HP08","HP09")
  
  names(eff_cts) <- tnames
  
  ## Create a matrix of transect name, grid cell ID, and X-Y coordinate information
  ## Number of rows needed based on elements of the list
  entries <- nrow(eff_cts[[1]]) + nrow(eff_cts[[2]]) + nrow(eff_cts[[3]]) + nrow(eff_cts[[4]]) + nrow(eff_cts[[5]]) + nrow(eff_cts[[6]]) + nrow(eff_cts[[7]]) + nrow(eff_cts[[8]]) + nrow(eff_cts[[9]]) + nrow(eff_cts[[10]]) + nrow(eff_cts[[11]]) + nrow(eff_cts[[12]]) + nrow(eff_cts[[13]]) + nrow(eff_cts[[14]]) + nrow(eff_cts[[15]]) + nrow(eff_cts[[16]]) + nrow(eff_cts[[17]]) + nrow(eff_cts[[18]])
  
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
  
  # Section 5. Convert BTS Captures to Grid Cells ----
  ## Import capture locations (from Select&Prep) and convert to UTM
  hmucaps <- HMUcaps[order(HMUcaps$TRANSECT, HMUcaps$LOCATION),]
  coordinates(hmucaps)=~CAPLON+CAPLAT
  ## Specify projection
  proj4string(hmucaps) <- CRS("+proj=longlat +datum=WGS84")
  ## Convert points to UTM
  hmucaps <- spTransform(hmucaps, CRS("+proj=utm +zone=55"))
  ## snap these capture locations to nearest transect line to correct for GPS inaccuracy
  names <- c("HE01","HE02","HE03","HE04","HE05","HE06","HE07","HE08","HE09","HP01","HP02","HP03","HP04","HP05","HP06","HP07","HP08","HP09")
  
  snaplist <- list()
  for (i in 1:length(hmutran)){
    snaplist[[i]] <- as.data.frame(snapPointsToLines(subset(hmucaps, TRANSECT == names[i]), hmutran[[i]]))  ## warning about CRS, ignore
  }
  ## Convert from list of dataframes to a single dataframe with transect info
  allsnap <- do.call(rbind,snaplist)
  caphmu <- cbind(as.data.frame(hmucaps),allsnap)
  ## Convert to spatial points dataframe
  spd <- caphmu
  coordinates(spd) =~ X + Y
  proj4string(spd) <- CRS("+proj=utm +zone=55")
  caphmu$GridID <- NA
  for(i in 1:length(spd)){
    caphmu[i,12] <- point.in.poly(spd[i,],clipHMU)@data$poly.ids
  }
  
  ## Check that no grid cells were identified as having snakes captured that weren't part of transect grid cells surveyed
  ifelse(unique(unique(caphmu$GridID) %in% unique(vistran$GridID)) == TRUE, "GridID matched", "error: GridID mismatch")
  
  ## Dataframe of snake capture info and the grid cell on which it occurred
  snkcap <- caphmu[,c("EFFORTID","PITTAG","Date","TRANSECT","X","Y","GridID")]
  
  



  # ## Checking. Plot Interior transects ----
  # plot(clipHMU, xlim=c(268550,268860), ylim=c(1504100,1504490))
  # plot(hmutran[[18]], add=TRUE, col="red")
  # plot(hmutran[[17]], add=TRUE, col="green")
  # plot(hmutran[[16]], add=TRUE, col="blue")
  # plot(hmutran[[15]], add=TRUE, col="orange")
  # plot(hmutran[[14]], add=TRUE, col="brown")
  # plot(hmutran[[13]], add=TRUE, col="purple")
  # plot(hmutran[[12]], add=TRUE, col="yellow")
  # plot(hmutran[[11]], add=TRUE, col="pink")
  # plot(hmutran[[10]], add=TRUE, col="lightblue")
  # plot(my.list[[18]], add=TRUE, col="red")
  # plot(my.list[[17]], add=TRUE, col="green")
  # plot(my.list[[16]], add=TRUE, col="blue")
  # plot(my.list[[15]], add=TRUE, col="orange")
  # plot(my.list[[14]], add=TRUE, col="brown")
  # plot(my.list[[13]], add=TRUE, col="purple")
  # plot(my.list[[12]], add=TRUE, col="yellow")
  # plot(my.list[[11]], add=TRUE, col="pink")
  # plot(my.list[[10]], add=TRUE, col="lightblue")
  # 
  # ## Checking. Plot Edge transects ----
  # plot(clipHMU, xlim=c(268735.1,268737.6), ylim=c(1503600,1504100))
  # plot(hmutran[[8]], add=TRUE, col="red")
  # plot(hmutran[[7]], add=TRUE, col="green")
  # plot(hmutran[[6]], add=TRUE, col="blue")
  # plot(hmutran[[5]], add=TRUE, col="orange")
  # plot(my.list[[8]], add=TRUE, col="red")
  # plot(my.list[[7]], add=TRUE, col="green")
  # plot(my.list[[6]], add=TRUE, col="blue")
  # plot(my.list[[5]], add=TRUE, col="orange")
  # ## This one seems too short but that's because edge transects were only 0.1km and this was the only one oriented longways (so fewer grid cells are part of the transect)
  # plot(clipHMU, xlim=c(268690.7,268819.6), ylim=c(1504643,1504685))
  # plot(hmutran[[9]], add=TRUE, col="red")
  # plot(my.list[[9]], add=TRUE, col="red")
  # plot(clipHMU, xlim=c(268945.6,269187.2), ylim=c(1504147,1504634))
  # plot(hmutran[[4]], add=TRUE, col="red")
  # plot(hmutran[[3]], add=TRUE, col="green")
  # plot(hmutran[[2]], add=TRUE, col="blue")
  # plot(hmutran[[1]], add=TRUE, col="orange")
  # plot(my.list[[4]], add=TRUE, col="red")
  # plot(my.list[[3]], add=TRUE, col="green")
  # plot(my.list[[2]], add=TRUE, col="blue")
  # plot(my.list[[1]], add=TRUE, col="orange")
  # 
  # ## Checking. Plot snapped captures to transects ----
  # snaps1 <- snaplist[[10]]
  # snaps2 <- snaplist[[11]]
  # snaps3 <- snaplist[[12]]
  # coordinates(snaps1) =~ X + Y
  # coordinates(snaps2) =~ X + Y
  # coordinates(snaps3) =~ X + Y
  # plot(clipHMU, xlim=c(268550,268860), ylim=c(1504100,1504490))
  # plot(hmucaps, add=TRUE, pch=21, cex=0.5) ## original
  # plot(snaps1, add=TRUE, pch=21, cex=0.5, col="red") ## adjusted
  # plot(snaps2, add=TRUE, pch=21, cex=0.5, col="red") ## adjusted
  # plot(snaps3, add=TRUE, pch=21, cex=0.5, col="red") ## adjusted
  # 
  # ## use drawExtent() to find spot on hmuSpace
  # plot(clipHMU, xlim=c(268614.2, 268677.9), ylim=c(1504381,1504436))
  # plot(hmucaps, add=TRUE, cex=1, pch=20)
  #### #####
  
  dat <- list(tran = vistran, snks = snkcap)
  
  return(dat)
}