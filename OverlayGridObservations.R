### Import observations
### Overlay 16m x 8m grid over area
### Assign animal observations to grid cell
### Determine where survey started and ended
### Record which grid cells were surveyed per evening

source("Select&PrepVisualData.R")

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

HMUtran <- read.csv()

