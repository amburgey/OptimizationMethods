## Code to specify grid locations for CP surveys

library(rgeos); library(maptools); library(spatialEco); library(raster); library(sf); library(reshape2)

overlayCP <- function(CPcaps, cellsize){
  # Section 1. Load coordinates for CP Perimeter fence ----
  ## Use kml and convert to UTM
  edgeline <- rgdal::readOGR("Real Data Analysis/Data/CPUperimeter.kml",require_geomType = "wkbLineString")
  ## Specify projection
  edgeline <- spTransform(edgeline, CRS("+proj=utm +zone=55"))
  ## Create spatial polygon
  ## get polygon geometry from edgeline (geom()) but have to copy and paste below
  CP <- readWKT("POLYGON ((268881.7 1508943, 269033.4 1508778, 269198.8 1508928, 269047.5 1509092, 268881.7 1508943))")
  ## Specify projection
  proj4string(CP) <- CRS("+proj=utm +zone=55 +units=m +datum=WGS84")
  plot(CP)
  
  
  # Section 2. Create Grid Polygon ----
  ## Grid should be 8 m by 16 m
  cpSpace <- raster(edgeline)
  ## Desired grid cell to match surveying grid of Closed Population
  cellsize1 = c(16,8)
  ## Overlay a grid of these dimensions across the space of the HMU (but with room for rotation) and then rotate
  bbox <- st_sfc(st_polygon(list(rbind(c(xmin(cpSpace)+7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymax(cpSpace)-2), c(xmin(cpSpace)+7,ymax(cpSpace)-2), c(xmin(cpSpace)+7,ymin(cpSpace)-2)))))
  grd <- sf::st_make_grid(bbox, cellsize = cellsize1, square = TRUE)
  rotang = -42.5
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  grd_rot <- (grd - st_centroid(st_union(grd))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd))
  ## Specify projection again
  st_crs(grd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  grd_rot <- as_Spatial(grd_rot, cast=TRUE, IDs=paste0("ID", seq_along(grd_rot)))
  plot(grd_rot, add=TRUE, col="red")
  
  
  # Section 3. Create CPU-clipped Survey Grid ----
  clipCP <- gIntersection(CP, grd_rot, byid = TRUE, drop_lower_td = TRUE)
  # par(mar = c(1,1,1,1))  ## if for some reason the margins are weird
  plot(clipCP)
  
  
  # Section 4. Convert Transect Survey Locations to Grid Cells ----
  ## Find centroid of all grid cells
  eff_cts <- list()
  for(i in 1:length(clipCP)){
    eff_cts[[i]] <- gCentroid(clipCP[i])  ## Warnings about geometries, ignore
  }
  
  ## Identify transects for CP visual surveys (all used at various times so 1-13 and A-AA)
  vistran <- setNames(data.frame(matrix(NA, nrow=(27*13), ncol=4)), c("TranID","GridID","x","y"))
  vistran$TranID <- as.vector(paste(rep(c(LETTERS,"AA"),each=13),1:13,sep=""))
  vistran$GridID <- c(407:419,392:404,377:389,362:374,347:359,332:344,317:329,302:314,287:299,272:284,257:269,242:254,227:239,212:224,197:209,182:194,167:179,152:164,137:149,122:134,107:119,92:104,77:89,62:74,47:59,32:44,17:29)
  for(i in 1:nrow(vistran)){
    vistran[i,3] <- xmin(extent(eff_cts[[vistran[i,2]]]))
    vistran[i,4] <- ymin(extent(eff_cts[[vistran[i,2]]]))
  }
  
  vistran <- vistran[order(vistran$GridID),]
  
  
# Section 5. Create rotated integration grid ----
  ## Desired grid cell
  cellsize2 = cellsize
  ## Overlay a grid of these dimensions across the space of the HMU (but with room for rotation) and then rotate
  bbox2 <- st_sfc(st_polygon(list(rbind(c(xmin(cpSpace)+7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymax(cpSpace)-2), c(xmin(cpSpace)+7,ymax(cpSpace)-2), c(xmin(cpSpace)+7,ymin(cpSpace)-2)))))
  grd2 <- sf::st_make_grid(bbox2, cellsize = cellsize2, square = TRUE)
  rotang = -42.5
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  intgrd_rot <- (grd2 - st_centroid(st_union(grd2))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd2))
  ## Specify projection again
  st_crs(intgrd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  intgrd_rot <- as_Spatial(intgrd_rot, cast=TRUE, IDs=paste0("ID", seq_along(intgrd_rot)))
  
# Section 7. Create CP-clipped Integration Grid ----
  clipG <- gIntersection(CP, intgrd_rot, byid = TRUE, drop_lower_td = TRUE)
  # par(mar = c(1,1,1,1))  ## if for some reason the margins are weird
  plot(clipG, col="red")
  ## Determine area of each integration grid cell
  area <- area(clipG)
  ## Find centroid of all grid cells
  intgrd_cts <- gCentroid(clipG, byid = TRUE)
  plot(intgrd_cts, col="red", pch=21, cex=0.2)
  intgrd <- geom(intgrd_cts)

  ## Plotting to visualize spatial information for troubleshooting transects
  # test <- as.data.frame(X)
  # coordinates(test) =~ x+y
  # # plot(clipCP, xlim=c(269028.4,269069.5), ylim=c(1508799,1508820))
  # plot(clipCP, xlim=c(269006.2,269214.2), ylim=c(1508796,1508972))
  # plot(intgrd_cts, add=TRUE, col="red", pch=21, cex=0.5)
  # plot(test, add=TRUE, col="pink", pch=5, cex=0.5)
  # ## Problem grid points
  # points(x=269050.2, y=1508808.4, col="blue", pch=3)
  # points(x=269109.2, y=1508862.5, col="blue", pch=3)
  # points(x=269168.2, y=1508916.5, col="blue", pch=3)
  # points(x=268996.2, y=1508867.4, col="blue", pch=3)
  # points(x=269055.2, y=1508921.5, col="blue", pch=3)
  # points(x=269114.1, y=1508975.5, col="blue", pch=3)
  # points(x=268942.1, y=1508926.4, col="blue", pch=3)
  # points(x=269001.1, y=1508980.4, col="blue", pch=3)
  # points(x=269060.1, y=1509034.5, col="blue", pch=3)
  # ## Problem trap points
  # points(x=269046.0, y=1508809, col="orange", pch=3)
  # points(x=269105, y=1508863, col="orange", pch=3)
  # points(x=269163.9, y=1508917, col="orange", pch=3)
  # points(x=268991.9, y=1508868, col="orange", pch=3)
  # points(x=269050.9, y=1508922, col="orange", pch=3)
  # points(x=269109.9, y=1508976, col="orange", pch=3)
  # points(x=268937.9, y=1508927, col="orange", pch=3)
  # points(x=268996.9, y=1508981, col="orange", pch=3)
  # points(x=269055.9, y=1509035, col="orange", pch=3)
  

  # ## Checking. Plot transects ----
  # plot(clipCP)
  # plot(eff_cts[[407]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[419]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[392]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[404]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[377]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[389]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[362]], add=TRUE, col="blue", pch=21)
  # plot(eff_cts[[374]], add=TRUE, col="blue", pch=21)
  # 
  
  dat <- list(tran = vistran, intgrd=intgrd, area=area)
  
  return(dat)
}

