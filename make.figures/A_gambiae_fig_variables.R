##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Multipanel figure of environmental variables
# Author(s): Devon A DeRaad, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 04/03/2020 (dd/mm/yyyy)
##################################

# packages
library(raster)
library(maps)

# function for legend
bar_legend <- function (value_range, col, alpha = 1, title = NULL, round = 0,
                        label_x = 0.7, labels_y = c(0.2, 0.85), 
                        legend_coord = c(0.1, 0.2, 0.3, 0.85), 
                        title_coord = c(0.6, 0.525), title_srt = 90,
                        horizontal = FALSE) {
  if (horizontal == TRUE) {
    legend_image <- as.raster(matrix(scales::alpha(rev(col), alpha), nrow = 1))
  } else {
    legend_image <- as.raster(matrix(scales::alpha(rev(col), alpha), ncol = 1))
  }
  text(x = title_coord[1], y = title_coord[2], labels = title, srt = title_srt)
  if (is.numeric(value_range)) {
    vals <- round(value_range, round)
  } else {
    vals <- value_range
  }
  text(x = label_x, y = labels_y, labels = vals, cex = 0.8)
  rasterImage(legend_image, legend_coord[1], legend_coord[2], 
              legend_coord[3], legend_coord[4])
}


# data
## raster layers
pan <- raster("Env_data_extraction_grouping/Data/crop_pt/prec_crop.tif")
tme <- raster("Env_data_extraction_grouping/Data/crop_pt/temp_crop.tif")
han <- raster("merra_2.5m_mean_00s_crop/bio12.tif")
aut <- raster("Env_data_extraction_grouping/Data/Spatial_autocor/spa_autor.tif")


## localities
occ <- read.csv("Env_data_extraction_grouping/phase1_all_vairables_last.csv")
head(occ)

# pre-processing
## cropping raster layers
ext <- extent(-25, 55, -40, 40)

conNA <- pan[] == 786420

pan[conNA] <- NA
tme <- mask(tme, pan)
han <- crop(han, ext)
aut <- crop(aut, ext)


## buffers
occp1 <- unique(occ[, c("longitude", "latitude")])
colnames(occp1)

### make a spatial object from coordinates using all occs
WGS84 <- aut@crs
occ_sp <- SpatialPointsDataFrame(coords = occp1, data = occp1, proj4string = WGS84)

### project the points using their centriods as reference
centroid <- gCentroid(occ_sp, byid = FALSE)

AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", 
                  centroid@coords[1],
                  " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                  sep = ""))

occ_pr <- spTransform(occ_sp, AEQD)

### create a buffer based on a user-defined distance
buff_area <- gBuffer(occ_pr, width = 100000, byid = T) # distance in meters
buff_areap <- spTransform(buff_area, WGS84)


## colors
colp <- rev(viridis::viridis(255))
colbf <- "#E31616"


# plot
## layout
matn <- c(2, 2, 4, 4, 2, 2, 4, 4, 3, 2, 5, 4, 6, 6, 8, 8, 
          6, 6, 8, 8, 7, 6, 9, 8, 1, 1, 1, 1)

jpeg("Env_data_extraction_grouping/map_variables.jpg", width = 16.6, height = 16.6, 
     units = "cm", res = 600)

layout(matrix(matn, ncol = 4, byrow = T), widths = rep(4, 4), 
       heights = c(rep(3, 6), 1))

## legend
par(cex = 0.8, mar = rep(0, 4))
plot.new()
bar_legend(value_range = c("Low", "High"), col = colp, 
           label_x = c(0.055, 0.65), labels_y = 0.3, 
           legend_coord = c(0.075, 0.4, 0.625, 0.8), 
           horizontal = TRUE)
legend(0.7, 1, bty = "n", legend = "Localities (100 km buffer)", pch = 21, pt.cex = 1.5, 
       pt.bg = NA, col = colbf)

## raster and histogram plots
par(cex = 0.7, mar = rep(0, 4))
maps::map(xlim = c(-15, 40), ylim = c(-25, 25), lforce = "n", col = NA)
image(tme, add = TRUE, col = colp)
plot(buff_areap, border = colbf, add = T)
box(which = "figure")

par(cex = 0.5, mar = c(2.5, 4.2, 2, 1))
hist(occ$t_mean / 100, breaks = 14, xlab = "", main = "Mean annual temperature")
box(bty = "l")


par(cex = 0.7, mar = rep(0, 4))
maps::map(xlim = c(-15, 40), ylim = c(-25, 25), lforce = "n", col = NA)
image(pan, add = TRUE, col = colp)
plot(buff_areap, border = colbf, add = T)
box(which = "figure")

par(cex = 0.5, mar = c(2.5, 4.2, 2, 1))
hist(occ$p_annual, breaks = 14, xlab = "", main = "Annual precipitation")
box(bty = "l")


par(cex = 0.7, mar = rep(0, 4))
maps::map(xlim = c(-15, 40), ylim = c(-25, 25), lforce = "n", col = NA)
image(han, add = TRUE, col = colp)
plot(buff_areap, border = colbf, add = T)
box(which = "figure")

par(cex = 0.5, mar = c(2.5, 4.2, 2, 1))
hist(occ$h_annual, breaks = 14, xlab = "", main = "Annual humidity")
box(bty = "l")


par(cex = 0.7, mar = rep(0, 4))
maps::map(xlim = c(-15, 40), ylim = c(-25, 25), lforce = "n", col = NA)
image(aut, add = TRUE, col = colp)
plot(buff_areap, border = colbf, add = T)
box(which = "figure")

par(cex = 0.5, mar = c(2.5, 4.2, 2, 1))
hist(occ$autocorrelated, breaks = 14, xlab = "", main = "Spatially autocorrelated")
box(bty = "l")

dev.off()