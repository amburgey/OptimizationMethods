### Import observations
### Overlay 16m x 8m grid over area
### Assign animal observations to grid cell
### Determine where survey started and ended
### Record which grid cells were surveyed per evening


#### HMU ####

hmuSpace <- raster("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Data/KMZ files/Reclass_hmu_71.tif")
hmuSpace <- projectRaster(hmuSpace, crs="+proj=utm +zone=55 +units=m +datum=WGS84")

## Change to grid (res of hmuSpace is 4.11 x 4.12 m but has to be this fine to prevent grid cells hanging over raster edge)
HMUhabmat <- as.matrix(hmuSpace)
HMUhabmat <- ifelse(HMUhabmat > 0, 1, 0) # makes cells that overhang the hmu a little bit part of hmu

## Import capture locations and convert to UTM
hmucaps <- hmuEdge
coordinates(hmucaps)=~CAPLON+CAPLAT
proj4string(hmucaps) <- CRS("+proj=longlat +datum=WGS84")
## Convert points to UTM
hmucaps <- spTransform(hmucaps, CRS("+proj=utm +zone=55"))

## Plot onto raster and figure out cell IDs for all captures
nr <- dim(HMUhabmat)[1]
nc <- dim(HMUhabmat)[2]
hmuEdge$cellIDs <- cellFromXY(hmuSpace, hmucaps)
hmuEdge$row <- NA
hmuEdge$col <- NA
for(i in 1:nrow(hmuEdge)){
  hmuEdge[i,(ncol(hmuEdge)-1)] <- floor(hmuEdge[i,5]/nc) # get row
  hmuEdge[i,ncol(hmuEdge)] <- hmuEdge[i,5] - (floor(hmuEdge[i,5]/nc)*nc) # get col
}

## Where did each survey start and end
hmuEsurv <- tibble(subsurv) %>% 
  filter(SITE == c("HMUI","HMUR") & PROJECTCODE == "EDGE EFFECT VIS")

HMUtran <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/HMU_EElocations.csv")
# coordinates(HMUtran) =~ UTME + UTMN
# proj4string(HMUtran) <- CRS("+proj=utm +zone=55")
# raster::extract(hmuSpace, HMUtran)  ## check that all points lie somewhere on raster within HMU (shifted two points that were slighlty outside)

hmutran <- list()
hmutran[[1]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[1,4],HMUtran[2,4]),c(HMUtran[1,5],HMUtran[2,5]))), ID = "HE1")))
hmutran[[2]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[3,4],HMUtran[4,4]),c(HMUtran[3,5],HMUtran[4,5]))), ID = "HE2")))
hmutran[[3]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[5,4],HMUtran[6,4]),c(HMUtran[5,5],HMUtran[6,5]))), ID = "HE3")))
hmutran[[4]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[7,4],HMUtran[8,4]),c(HMUtran[7,5],HMUtran[8,5]))), ID = "HE4")))
hmutran[[5]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[9,4],HMUtran[10,4]),c(HMUtran[9,5],HMUtran[10,5]))), ID = "HE5")))
hmutran[[6]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[11,4],HMUtran[12,4]),c(HMUtran[11,5],HMUtran[12,5]))), ID = "HE6")))
hmutran[[7]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[13,4],HMUtran[14,4]),c(HMUtran[13,5],HMUtran[14,5]))), ID = "HE7")))
hmutran[[8]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[15,4],HMUtran[16,4]),c(HMUtran[15,5],HMUtran[16,5]))), ID = "HE8")))
hmutran[[9]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[17,4],HMUtran[18,4]),c(HMUtran[17,5],HMUtran[18,5]))), ID = "HE9")))
hmutran[[10]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[19,4],HMUtran[20,4]),c(HMUtran[19,5],HMUtran[20,5]))), ID = "HP01")))
hmutran[[11]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[21,4],HMUtran[22,4]),c(HMUtran[21,5],HMUtran[22,5]))), ID = "HP02")))
hmutran[[12]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[23,4],HMUtran[24,4]),c(HMUtran[23,5],HMUtran[24,5]))), ID = "HP03")))
hmutran[[13]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[25,4],HMUtran[26,4]),c(HMUtran[25,5],HMUtran[26,5]))), ID = "HP04")))
hmutran[[14]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[27,4],HMUtran[28,4]),c(HMUtran[27,5],HMUtran[28,5]))), ID = "HP05")))
hmutran[[15]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[29,4],HMUtran[30,4]),c(HMUtran[29,5],HMUtran[30,5]))), ID = "HP06")))
hmutran[[16]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[31,4],HMUtran[32,4]),c(HMUtran[31,5],HMUtran[32,5]))), ID = "HP07")))
hmutran[[17]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[33,4],HMUtran[34,4]),c(HMUtran[33,5],HMUtran[34,5]))), ID = "HP08")))
hmutran[[18]] <- SpatialLines(list(Lines(Line(cbind(c(HMUtran[35,4],HMUtran[36,4]),c(HMUtran[35,5],HMUtran[36,5]))), ID = "HP09")))


cellIDs <- NA
hmuEsurveff <- as.data.frame(matrix(data=c(rep(unique(HMUtran[,1]), each=9),unique(HMUtran[,2])), nrow=18, ncol=2))

for(i in 1:18){
  cellIDs[i] <- unlist(cellFromLine(hmuSpace, hmutran[[i]]))
}

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

# df <- data.frame(matrix(unlist(cellIDs), nrow=length(cellIDs), byrow=T))






