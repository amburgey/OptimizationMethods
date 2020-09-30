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

cols <- c(0,3.202,seq(12,204,16),214,217.202)
CPlocs <- rbind.data.frame(cbind(cols,rep(0,17),paste("PR",seq(1:17),sep="")),
                           cbind(cols,rep(3.202,17),c("PR91",paste("NEE",seq(0:14),sep=""),"PR17")),
                           cbind(cols,rep(11,17),c("PR90",paste("A",seq(0:14),sep=""),"PR18")),
                           cbind(cols,rep(19,17),c("PR89",paste("B",seq(0:14),sep=""),"PR19")),
                           cbind(cols,rep(27,17),c("PR88",paste("C",seq(0:14),sep=""),"PR20")),
                           cbind(cols,rep(35,17),c("PR87",paste("D",seq(0:14),sep=""),"PR21")),
                           cbind(cols,rep(43,17),c("PR86",paste("E",seq(0:14),sep=''),"PR22")),
                           cbind(cols,rep(51,17),c("PR85",paste("F",seq(0:14),sep=""),"PR23")),
                           cbind(cols,rep(59,17),c("PR84",paste("G",seq(0:14),sep=""),"PR24")),
                           cbind(cols,rep(67,17),c("PR83",paste("H",seq(0:14),sep=""),"PR25")),
                           cbind(cols,rep(75,17),c("PR82",paste("I",seq(0:14),sep=""),"PR26")),
                           cbind(cols,rep(83,17),c("PR81",paste("J",seq(0:14),sep=""),"PR27")),
                           cbind(cols,rep(91,17),c("PR80",paste("K",seq(0:14),sep=""),"PR28")),
                           cbind(cols,rep(99,17),c("PR79",paste("L",seq(0:14),sep=""),"PR29")),
                           cbind(cols,rep(107,17),c("PR78",paste("M",seq(0:14),sep=""),"PR30")),
                           cbind(cols,rep(115,17),c("PR77",paste("N",seq(0:14),sep=""),"PR31")),
                           cbind(cols,rep(123,17),c("PR76",paste("O",seq(0:14),sep=""),"PR32")),
                           cbind(cols,rep(131,17),c("PR75",paste("P",seq(0:14),sep=""),"PR33")),
                           cbind(cols,rep(139,17),c("PR74",paste("Q",seq(0:14),sep=""),"PR34")),
                           cbind(cols,rep(147,17),c("PR73",paste("R",seq(0:14),sep=""),"PR35")),
                           cbind(cols,rep(155,17),c("PR72",paste("S",seq(0:14),sep=""),"PR36")),
                           cbind(cols,rep(163,17),c("PR71",paste("T",seq(0:14),sep=""),"PR37")),
                           cbind(cols,rep(171,17),c("PR70",paste("U",seq(0:14),sep=""),"PR38")),
                           cbind(cols,rep(179,17),c("PR69",paste("V",seq(0:14),sep=""),"PR39")),
                           cbind(cols,rep(187,17),c("PR68",paste("W",seq(0:14),sep=""),"PR40")),
                           cbind(cols,rep(195,17),c("PR67",paste("X",seq(0:14),sep=""),"PR41")),
                           cbind(cols,rep(203,17),c("PR66",paste("Y",seq(0:14),sep=""),"PR42")),
                           cbind(cols,rep(211,17),c("PR65",paste("Z",seq(0:14),sep=""),"PR43")),
                           cbind(cols,rep(219,17),c("PR64",paste("AA",seq(0:14),sep=""),"PR44")),
                           cbind(cols,rep(227,17),c("PR63",paste("SWE",seq(0:14),sep=""),"PR45")),
                           cbind(cols,rep(230.202,17),c(paste("PR",c(62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46),sep=""))))

colnames(CPlocs) <- c("x","y","Location")
CPlocs$x <- as.numeric(as.character(CPlocs$x))
CPlocs$y <- as.numeric(as.character(CPlocs$y))


## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
cpdelta <- 0  ## no buffer as using entire space of CP
Xlcp<-min(CPlocs[,1]) - cpdelta
Xucp<-max(CPlocs[,1]) + cpdelta
Ylcp<-min(CPlocs[,2]) - cpdelta
Yucp<-max(CPlocs[,2]) + cpdelta
Acp <- (Xucp-Xlcp)*(Yucp-Ylcp)