##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Phase 2 Environmental data extraction and categorization
# Author(s): Marlon E. Cobos
# Date: 17/05/2020 (dd/mm/yyyy)
##################################

# needed packages
pcakages <- c("raster", "rgdal")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# data
## occurrences
occ <- read.csv("phase2.unique.sampling.locs.csv", row.names = 1) 

## environmental data
vas <- list.files("Env_data_extraction_grouping/Data/Final_variables", 
                  pattern = ".tif$", full.names = TRUE)
hum <- stack(vas[1:2])
vars <- stack(vas[3:6])

# extraction
hvals <- extract(hum, occ[, c("longitude", "latitude")])
tpvals <- extract(vars, occ[, c("longitude", "latitude")])


## complete table
ctable <- data.frame(occ[, c(1, 3, 2)], cbind(tpvals, hvals)) # table

## searching for high and low values of variables
un_ctable <- ctable[!duplicated(paste0(ctable$longitude, "_", ctable$latitude)), ]
categs <- list()
categsall <- list()

for (i in 4:9) {
  low <- sort(un_ctable[, i])[8]
  high <- sort(un_ctable[, i])[22]
  lows <- which(un_ctable[, i] <= low); lows1 <- which(ctable[, i] <= low)
  highs <- which(un_ctable[, i] >= high); highs1 <- which(ctable[, i] >= high)
  len <- length(un_ctable[, i]); len1 <- length(ctable[, i]) 
  vec <- vector("character", len); vec1 <- vector("character", len1)
  vec[lows] <- "LOW"; vec1[lows1] <- "LOW" 
  vec[highs] <- "HIGH"; vec1[highs1] <- "HIGH"
  vec[!1:len %in% c(lows, highs)] <- "MEDIUM"
  vec1[!1:len1 %in% c(lows1, highs1)] <- "MEDIUM"
  categs[[i]] <- vec; categsall[[i]] <- vec1 
}

categs[sapply(categs, is.null)] <- NULL
categsall[sapply(categsall, is.null)] <- NULL

cattable <- data.frame(un_ctable, do.call(cbind, categs)) # categorical table
colnames(cattable)[10:15] <- paste0(colnames(un_ctable)[4:9], "_categ")

cattableall <- data.frame(ctable, do.call(cbind, categsall)) # categorical table
colnames(cattableall)[10:15] <- paste0(colnames(ctable)[4:9], "_categ")

## save
write.csv(ctable, file = "Env_data_extraction_grouping/phase2_localities_all_vairables.csv", 
          row.names = FALSE)
write.csv(cattableall, file = "Env_data_extraction_grouping/phase2_localities_all_vairables_categorical.csv", 
          row.names = FALSE)
write.csv(cattable, file = "Env_data_extraction_grouping/phase2_unique_localities_all_vairables_categorical.csv", 
          row.names = FALSE)

# shapefile
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
occ_sp <- SpatialPointsDataFrame(coords = cattableall[, 2:3], data = cattableall,
                                 proj4string = WGS84)

writeOGR(occ_sp, ".", "Env_data_extraction_grouping/Shapefiles/A_gambiae_variables_phase2", 
         driver = "ESRI Shapefile")

# histograms
par(mfrow = c(6, 1), cex = 0.7, mar = c(4,4,1,1))
hist(cattableall[, 4], breaks = 14)
hist(cattableall[, 5], breaks = 14)
hist(cattableall[, 6], breaks = 14)
hist(cattableall[, 7], breaks = 14)
hist(cattableall[, 8], breaks = 14)
hist(cattableall[, 9], breaks = 14)

## saving figures
png("Env_data_extraction_grouping/A_gambiae_LMH_groups_phase2.png", width = 166, height = 205, units = "mm", res = 300)
par(mfrow = c(6, 1), cex = 0.7, mar = c(4.5,4,1,1))
hist(cattableall[, 6] / 100, breaks = 14, xlab = "Annual mean temperature (°C)", main = "")
hist(cattableall[, 7], breaks = 14, xlab = "Temperature seasonality (SD)", main = "")
hist(cattableall[, 4], breaks = 14, xlab = "Annual precipitation (mm)", main = "")
hist(cattableall[, 5], breaks = 14, xlab = "Precipitation seasonality (CV)", main = "")
hist(cattableall[, 8], breaks = 14, xlab = "Annual humidity (kg water/kg air)", main = "")
hist(cattableall[, 9], breaks = 14, xlab = "Humidity seasonality (CV)", main = "")
dev.off()

