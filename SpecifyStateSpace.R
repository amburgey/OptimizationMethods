### Defining state space for each study area


#### HMU ####

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
## raster of 1 = study area and 0 = outside fence
hmuSpace <- raster("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/Reclass_hmu_71.tif")
hmuSpace <- projectRaster(hmuSpace, crs="+proj=utm +zone=55 +units=m +datum=WGS84")
coordinates(hmuEdgeAll)=~UTME+UTMN
proj4string(hmuEdgeAll) <- CRS("+proj=utm +zone=55 +units=m +datum=WGS84")

## function from http://rstudio-pubs-static.s3.amazonaws.com/273756_4230265fc8484963903aaea122f933f7.html
cellnumbers <- function(x, query, ...) {
  if (inherits(query, "sf")) query <- sf::as(query, "Spatial")
  if (is.na(projection(x)) || is.na(projection(query)) || projection(x) != projection(query)) {
    warning(sprintf("projections not the same \n    x: %s\nquery: %s", projection(x), projection(query)))
  }
  if (inherits(query, "SpatialPolygons")) {
    a <- cellFromPolygon(x, query, ...)
  }
  if (inherits(query, "SpatialLines")) {
    a <- cellFromLine(x, query, ...)
  }
  if (is.matrix(query) | inherits(query, "SpatialPoints")) {
    a <- cellFromXY(x, query)
  }
  d <- dplyr::bind_rows(lapply(a, mat2d_f), .id = "object_")
  if (ncol(d) == 2L) names(d) <- c("object_", "cell_")
  if (ncol(d) == 3L) names(d) <- c("object_", "cell_", "weight_")
  d
  
}

#' @importFrom tibble as_tibble
mat2d_f <- function(x) {
  
  if (is.null(x)) {
    return(NULL)
  }
  tibble::as_tibble(x)
}

tt <- cellnumbers(hmuSpace, hmuEdgeAll[,3:4]) ## get cell identities for each capture location on raster
HMUhabmat <- as.matrix(hmuSpace)
HMUhabmat <- ifelse(HMUhabmat > 0, 1, 0)  ## simplifies pixels that are not fully in state space to 0 or 1
## Find locations of captures (and currently beginning and end of transects) in pixels vs. lat/lon
locsHMUEE <- matrix(NA, nrow=nrow(hmuEdgeAll), ncol=2)
for(i in 1:nrow(tt)){
  locsHMUEE[i,1] <- as.numeric(floor(tt[i,2]/ncol(hmuSpace))) + 1
  locsHMUEE[i,2] <- as.numeric(tt[i,2]) - (as.numeric(floor(tt[i,2]/ncol(hmuSpace)))*ncol(hmuSpace))
}


hmuEEdelta <- 0  ## will need to play with this
XlhmuEE<-min(hmuEdge[,3]) - hmuEEdelta
XuhmuEE<-max(hmuEdge[,3]) + hmuEEdelta
YlhmuEE<-min(hmuEdge[,4]) - hmuEEdelta
YuhmuEE<-max(hmuEdge[,4]) + hmuEEdelta
AhmuEE <- (XuhmuEE-XlhmuEE)*(YuhmuEE-YlhmuEE)



#### NCR ####

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
ncrEEdelta <- 5  ## will need to play with this
XlncrEE<-min(ncrEE[,3]) - ncrEEdelta
XuncrEE<-max(ncrEE[,3]) + ncrEEdelta
YlncrEE<-min(ncrEE[,4]) - ncrEEdelta
YuncrEE<-max(ncrEE[,4]) + ncrEEdelta
AncrEE <- (XuncrEE-XlncrEE)*(YuncrEE-YlncrEE)



#### NWFN ####

cols <- c(0,7.608,seq(15.608,207.608,16),215.608,223.215)
CPlocs <- rbind.data.frame(cbind(cols,rep(0,17),c(paste("PR",c(62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46),sep=""))),
                           cbind(cols,rep(4,17),c("PR63",paste("SWE",0:14,sep=""),"PR45")),
                           cbind(cols,rep(10,17),c("PR64",paste("AA",0:14,sep=""),"PR44")),
                           cbind(cols,rep(18,17),c("PR65",paste("Z",0:14,sep=""),"PR43")),
                           cbind(cols,rep(26,17),c("PR66",paste("Y",0:14,sep=""),"PR42")),
                           cbind(cols,rep(32,17),c("PR67",paste("X",0:14,sep=""),"PR41")),
                           cbind(cols,rep(40,17),c("PR68",paste("W",0:14,sep=""),"PR40")),
                           cbind(cols,rep(47,17),c("PR69",paste("V",0:14,sep=""),"PR39")),
                           cbind(cols,rep(55,17),c("PR70",paste("U",0:14,sep=""),"PR38")),
                           cbind(cols,rep(63,17),c("PR71",paste("T",0:14,sep=""),"PR37")),
                           cbind(cols,rep(71,17),c("PR72",paste("S",0:14,sep=""),"PR36")),
                           cbind(cols,rep(79,17),c("PR73",paste("R",0:14,sep=""),"PR35")),
                           cbind(cols,rep(87,17),c("PR74",paste("Q",0:14,sep=""),"PR34")),
                           cbind(cols,rep(95,17),c("PR75",paste("P",0:14,sep=""),"PR33")),
                           cbind(cols,rep(103,17),c("PR76",paste("O",0:14,sep=""),"PR32")),
                           cbind(cols,rep(110,17),c("PR77",paste("N",0:14,sep=""),"PR31")),
                           cbind(cols,rep(119,17),c("PR78",paste("M",0:14,sep=""),"PR30")),
                           cbind(cols,rep(127,17),c("PR79",paste("L",0:14,sep=""),"PR29")),
                           cbind(cols,rep(135,17),c("PR80",paste("K",0:14,sep=""),"PR28")),
                           cbind(cols,rep(143,17),c("PR81",paste("J",0:14,sep=""),"PR27")),
                           cbind(cols,rep(151,17),c("PR82",paste("I",0:14,sep=""),"PR26")),
                           cbind(cols,rep(159,17),c("PR83",paste("H",0:14,sep=""),"PR25")),
                           cbind(cols,rep(167,17),c("PR84",paste("G",0:14,sep=""),"PR24")),
                           cbind(cols,rep(175,17),c("PR85",paste("F",0:14,sep=""),"PR23")),
                           cbind(cols,rep(183,17),c("PR86",paste("E",0:14,sep=''),"PR22")),
                           cbind(cols,rep(191,17),c("PR87",paste("D",0:14,sep=""),"PR21")),
                           cbind(cols,rep(199,17),c("PR88",paste("C",0:14,sep=""),"PR20")),
                           cbind(cols,rep(207,17),c("PR89",paste("B",0:14,sep=""),"PR19")),
                           cbind(cols,rep(215,17),c("PR90",paste("A",0:14,sep=""),"PR18")),
                           cbind(cols,rep(220,17),c("PR91",paste("NEE",0:14,sep=""),"PR17")),
                           cbind(cols,rep(224,17),paste("PR",0:16,sep="")))


colnames(CPlocs) <- c("x","y","Location")
CPlocs$x <- as.numeric(as.character(CPlocs$x))
CPlocs$y <- as.numeric(as.character(CPlocs$y))
plot(CPlocs[,1:2])
text(CPlocs[,1:2], labels=CPlocs[,3], cex=0.5, font=2)


## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
cpdelta <- 0  ## no buffer as using entire space of CP
Xlcp<-min(CPlocs[,1]) - cpdelta
Xucp<-max(CPlocs[,1]) + cpdelta
Ylcp<-min(CPlocs[,2]) - cpdelta
Yucp<-max(CPlocs[,2]) + cpdelta
Acp <- (Xucp-Xlcp)*(Yucp-Ylcp)  ## 0.16 m over = close enough for me




#### ASOI ####

#### Add coordinates in for captures for ASOI and RTSI
asoi <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_ASOI/TrapRingWeb_coords.csv")

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
asdelta <- 10 
Xlas<-min(asoi[,6]) - asdelta
Xuas<-max(asoi[,6]) + asdelta
Ylas<-min(asoi[,7]) - asdelta
Yuas<-max(asoi[,7]) + asdelta
Aas <- (Xuas-Xlas)*(Yuas-Ylas)




#### RTSI ####
## Similar to HMU, GNWR has a one-way barrier surrounding the perimeter (allowing snakes out but not in)
## Unlike HMU (where a grid was used to define when snakes are in the study area or not), the transects are far enough away from the perimeter fence that we will just use a buffer
bird <- readGPX("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_SSP/SSP Bird Trap Locations.gpx") #AB-EF, 1-16
mouse <- readGPX("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/To Beth_SSP/SSP Mouse Trap Locations.gpx") #A-F, 1-18 transects

rtsi <- rbind(as.data.frame(mouse$waypoints[1:3]), as.data.frame(bird$waypoints[1:3]))
colnames(rtsi) <- c("Lon","Lat","LOCATION")
coordinates(rtsi) =~ Lon + Lat
proj4string(rtsi) <- CRS("+proj=longlat +datum=WGS84")
rallUTM <- spTransform(rtsi, CRS("+proj=utm +zone=55 +ellps=WGS84"))
rall <- as.data.frame(rallUTM)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
rtdelta <- 10 
Xlrt<-min(rall[,3]) - rtdelta
Xurt<-max(rall[,3]) + rtdelta
Ylrt<-min(rall[,2]) - rtdelta
Yurt<-max(rall[,2]) + rtdelta
Art <- (Xurt-Xlrt)*(Yurt-Ylrt)




#### NWFO ####
## Using a standard grid for now as this set of transects seems to be set up like CP (just outside fence)
## 16 m between 13 traps in a transect, 8 m (???) between transects
CPOlocs <- secr::make.grid(nx = 13, ny = 6, spacex = 16, spacey = 8)
rownames(CPOlocs) <- paste(rep(c("BB","CC","DD","EE","FF","GG"), each=13), c(1:13), sep = "")


