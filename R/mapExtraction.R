## Load needed packages
#library(raster)
#library(rgdal)
#library(dplyr)
#
## Method 3:shapefiles
#library(maptools)
#library(gfcanalysis)
## plotting
#library(ggplot2)
#library(rgeos)
#library(geosphere)
#library(sp)


make_GeodesicBuffer <- function(pts, width) {
  # A) Construct buffers as points at given distance and bearing ---------------
  dg <- seq(from = 0, to = 360, by = 5)
  # Construct equidistant points defining circle shapes (the "buffer points")
  buff.XY <- geosphere::destPoint(p = pts,
                                  b = rep(dg, each = length(pts)),
                                  d = width)
  # Group (split) "buffer points" by id
  buff.XY <- as.data.frame(buff.XY)
  id  <- rep(1:dim(pts)[1], times = length(dg))
  lst <- split(buff.XY, id)
  # Make SpatialPolygons out of the list of coordinates
  poly   <- lapply(lst, sp::Polygon, hole = FALSE)
  polys  <- lapply(list(poly), sp::Polygons, ID = NA)
  spolys <- sp::SpatialPolygons(Srl = polys,
                                proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  # Disaggregate (split in unique polygons)
  spolys <- sp::disaggregate(spolys)
  return(spolys)
}

defoExt<-function(gr){
aoi<-make_GeodesicBuffer(gr,10000)
#library(plotKML)
#kml(aoi, file.name = "pts_buf_100km.kml")
output_folder <- "/Users/alm204/Documents/Cambridge/LandCovDat"
forest_threshold <- 90

tiles <- calc_gfc_tiles(aoi)
print(length(tiles)) # Number of tiles needed to cover AOI
plot(tiles)
plot(aoi, add=TRUE, lty=2, col="#00ff0050")
download_tiles(tiles, output_folder)
gfc_extract <- extract_gfc(aoi, output_folder, filename="NAK_GFC_extract.tif",overwrite=T)
gfc_thresholded <- threshold_gfc(gfc_extract, forest_threshold=forest_threshold,
                                 filename="NAK_GFC_extract_thresholded.tif",overwrite=T)

gs<-SpatialPolygonsDataFrame(aoi,data=gr)
gfc_stats <- gfc_stats(gs, gfc_thresholded)
return(gfc_stats)
}

#def1<-defoExt(gr)
#
##get spatial polygons
#gr<-data.frame(lat=c(35.01051),lon=c(-15.78841))
#
#-15.78841	35.01051#
