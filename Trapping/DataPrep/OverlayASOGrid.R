## Code to specify grid locations for CP visual surveys

library(rgeos); library(maptools); library(spatialEco)

overlayCP <- function(ASOcaps, cellsize){
  # Section 1. Load coordinates for ASO Ring/Web structures ----
  ## Use csv and convert to spatial
  trapcoords <- read.csv("Data/TrapRingWeb_coords.csv")
  coordinates(trapcoords) =~ UTM.X + UTM.Y
  ## Specify projection
  proj4string(trapcoords) <- CRS("+proj=utm +zone=55 +units=m +datum=WGS84")
  plot(trapcoords, pch=21, cex=0.5)

  
  # Section 2. Create Grid Polygon ----
  ## Desired grid cell to match surveying grid of Closed Population
  cellsize1 = c(16,8)
  ## Overlay a grid of these dimensions across the space of the ASO (but with room for rotation) and then rotate
  bbox <- st_sfc(st_polygon(list(rbind(c(269420,1494410), c(270941.5,1494410), c(270941.5,1495350), c(269420,1495350), c(269420,1494410)))))
  grd <- sf::st_make_grid(bbox, cellsize = cellsize1, square = TRUE)
  rotang = 42.5
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  grd_rot <- (grd - st_centroid(st_union(grd))) * rot(rotang * pi / 180) +
    st_centroid(st_union(grd))
  ## Specify projection again
  st_crs(grd_rot) <- "+proj=utm +zone=55 +units=m +datum=WGS84"
  ## Convert from sfc polygon to Spatial Polygon
  grd_rot <- as_Spatial(grd_rot, cast=TRUE, IDs=paste0("ID", seq_along(grd_rot)))
  plot(grd_rot, add=TRUE, col="red")
  
  plot(grd_rot, col="red", xlim=c(270638.3,270855.6), ylim=c(1494576,1494743))
  plot(trapcoords, add=TRUE, pch=21, cex=0.6)
  
  
  # Section 3. Convert Trap Locations to Grid Cells ----
  ## Find centroid of all grid cells
  eff_cts <- list()
  for(i in 1:length(grd_rot)){
    eff_cts[[i]] <- gCentroid(grd_rot[i])  ## Warnings about geometries, ignore
  }
  
  ## Identify which grid cells have traps in them
  vistran <- setNames(data.frame(matrix(NA, nrow=(27*13), ncol=4)), c("TranID","GridID","x","y"))
  vistran$TranID <- unique(ASOsurv)
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
  bbox2 <- st_sfc(st_polygon(list(rbind(c(xmin(cpSpace)-7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymin(cpSpace)-2), c(xmax(cpSpace)+7,ymax(cpSpace)+2), c(xmin(cpSpace)-7,ymax(cpSpace)+2), c(xmin(cpSpace)-7,ymin(cpSpace)-2)))))
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


