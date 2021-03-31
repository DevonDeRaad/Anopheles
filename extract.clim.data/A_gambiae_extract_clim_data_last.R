##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Environmental data extraction
# Author(s): Marlon E. Cobos
# Date: 26/02/2021 (dd/mm/yyyy)
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
occ <- read.csv("all_samples.csv") 

## environmental data
vas <- list.files("Env_data_extraction_grouping/Data/Final_variables", 
                  pattern = ".tif$", full.names = TRUE)
hum_ext <- stack(vas[c(1:2, 5)])
vars <- stack(vas[c(3:4, 6:7)])

# extraction
## occurrences
occall <- unique(occ[, c("longitude", "latitude")])
colnames(occall)

# make a spatial object from coordinates using all occs
WGS84 <- hum_ext@crs
occ_sp <- SpatialPointsDataFrame(coords = occall, data = occall, proj4string = WGS84)

# project the points using their centriods as reference
centroid <- gCentroid(occ_sp, byid = FALSE)

AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", 
                  centroid@coords[1],
                  " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                  sep = ""))

occ_pr <- spTransform(occ_sp, AEQD)

# create a buffer based on a user-defined distance
buff_area <- gBuffer(occ_pr, width = 100000, byid = T) # distance in meters
buff_areap <- spTransform(buff_area, WGS84)
plot(buff_areap, add = TRUE)

# extract values of each buffer
## extraction
val <- lapply(1:length(buff_areap), function(x) {
  v <- list(
    p_annual = na.omit(mask(vars[[1]], buff_areap[x, ])[]),  
    p_seasonality = na.omit(mask(vars[[2]], buff_areap[x, ])[]),
    t_mean = na.omit(mask(vars[[3]], buff_areap[x, ])[]),  
    t_seasonality = na.omit(mask(vars[[4]], buff_areap[x, ])[]),
    h_annual = na.omit(mask(hum_ext[[1]], buff_areap[x, ])[]),  
    h_seasonality = na.omit(mask(hum_ext[[2]], buff_areap[x, ])[]),
    autocorrelated = na.omit(mask(hum_ext[[3]], buff_areap[x, ])[])
  )
  
  cat(x, "of", length(buff_areap), "\n")
  v
})

## checking quantiles
statv <- lapply(val, function(x) {
  data.frame(
    p_annual = quantile(x[[1]]),
    p_seasonality = quantile(x[[2]]),
    t_mean = quantile(x[[3]]),
    t_seasonality = quantile(x[[4]]),
    h_annual = quantile(x[[5]]),
    h_seasonality = quantile(x[[6]]),
    autocorrelated = quantile(x[[7]])
  )
})

## mean values
vals <- lapply(val, function(x) {
  c(
    p_annual = mean(x[[1]]),
    p_seasonality = mean(x[[2]]),
    t_mean = mean(x[[3]]),
    t_seasonality = mean(x[[4]]),
    h_annual = mean(x[[5]]),
    h_seasonality = mean(x[[6]]),
    autocorrelated = mean(x[[7]])
  )
})

## association of values
allvals <- data.frame(to_match = paste0(occall[, 1], "-", occall[, 2]), 
                      do.call(rbind, vals))

## random values
set.seed(0)
allvals$random <- runif(nrow(allvals))

## complete table
ctab <- data.frame(occ[, c(1:14, 16, 15)], to_match = paste0(occ[, 16], "-", occ[, 15])) # table
ctable <- merge(ctab, allvals, by = "to_match")
ctable$to_match <- NULL

## separating phase 1 and 2
pha2 <- read.csv("phase2.unique.sampling.locs.csv", row.names = 1)
head(pha2)
head(ctable)

pha2_comp <- merge(data.frame(ox_code = pha2$ox_code), ctable, by = "ox_code")
dim(pha2_comp)

pha2ox <- pha2$ox_code
pha1_comp <- ctable[!ctable$ox_code %in% pha2ox, ]
dim(pha1_comp)

### excluding locality with an only individual
pha1_comp <- pha1_comp[pha1_comp$latitude != 6.09449, ]

## unique localities two phases
ctable_un <- ctable[!duplicated(paste0(ctable$longitude, "_", ctable$latitude)), ]
pha1_un <- pha1_comp[!duplicated(paste0(pha1_comp$longitude, "_", pha1_comp$latitude)), ]
pha2_un <- pha2_comp[!duplicated(paste0(pha2_comp$longitude, "_", pha2_comp$latitude)), ]


## save
write.csv(ctable, file = "Env_data_extraction_grouping/localities_all_vairables_last.csv", 
          row.names = FALSE)
write.csv(ctable_un, file = "Env_data_extraction_grouping/unique_localities_all_vairables_last.csv", 
          row.names = FALSE)

write.csv(pha1_comp, file = "Env_data_extraction_grouping/phase1_all_vairables_last.csv", 
          row.names = FALSE)
write.csv(pha1_un, file = "Env_data_extraction_grouping/unique_phase1_all_vairables_last.csv", 
          row.names = FALSE)

write.csv(pha2_comp, file = "Env_data_extraction_grouping/phase2_all_vairables_last.csv", 
          row.names = FALSE)
write.csv(pha2_un, file = "Env_data_extraction_grouping/unique_phase2_all_vairables_last.csv", 
          row.names = FALSE)


## saving figures
png("Env_data_extraction_grouping/A_gambiae_variables_all_localities_last.png", width = 166, 
    height = 190, units = "mm", res = 300)
par(mfrow = c(4, 2), cex = 0.7, mar = c(4.5,4,1,1))
hist(ctable[, 19] / 100, breaks = 14, xlab = "Annual mean temperature (°C)", main = "")
hist(ctable[, 20], breaks = 14, xlab = "Temperature seasonality (SD)", main = "")
hist(ctable[, 17], breaks = 14, xlab = "Annual precipitation (mm)", main = "")
hist(ctable[, 18], breaks = 14, xlab = "Precipitation seasonality (CV)", main = "")
hist(ctable[, 21], breaks = 14, xlab = "Annual humidity (kg water/kg air)", main = "")
hist(ctable[, 22], breaks = 14, xlab = "Humidity seasonality (CV)", main = "")
hist(ctable[, 23], breaks = 14, xlab = "Spatially autocorrelated values", main = "")
hist(ctable[, 24], breaks = 14, xlab = "Random values", main = "")
dev.off()

png("Env_data_extraction_grouping/A_gambiae_variables_phase1_last.png", width = 166, 
    height = 190, units = "mm", res = 300)
par(mfrow = c(4, 2), cex = 0.7, mar = c(4.5,4,1,1))
hist(pha1_comp[, 19] / 100, breaks = 14, xlab = "Annual mean temperature (°C)", main = "")
hist(pha1_comp[, 20], breaks = 14, xlab = "Temperature seasonality (SD)", main = "")
hist(pha1_comp[, 17], breaks = 14, xlab = "Annual precipitation (mm)", main = "")
hist(pha1_comp[, 18], breaks = 14, xlab = "Precipitation seasonality (CV)", main = "")
hist(pha1_comp[, 21], breaks = 14, xlab = "Annual humidity (kg water/kg air)", main = "")
hist(pha1_comp[, 22], breaks = 14, xlab = "Humidity seasonality (CV)", main = "")
hist(pha1_comp[, 23], breaks = 14, xlab = "Spatially autocorrelated values", main = "")
hist(pha1_comp[, 24], breaks = 14, xlab = "Random values", main = "")
dev.off()

png("Env_data_extraction_grouping/A_gambiae_variables_phase2_last.png", width = 166, 
    height = 190, units = "mm", res = 300)
par(mfrow = c(4, 2), cex = 0.7, mar = c(4.5,4,1,1))
hist(pha2_comp[, 19] / 100, breaks = 14, xlab = "Annual mean temperature (°C)", main = "")
hist(pha2_comp[, 20], breaks = 14, xlab = "Temperature seasonality (SD)", main = "")
hist(pha2_comp[, 17], breaks = 14, xlab = "Annual precipitation (mm)", main = "")
hist(pha2_comp[, 18], breaks = 14, xlab = "Precipitation seasonality (CV)", main = "")
hist(pha2_comp[, 21], breaks = 14, xlab = "Annual humidity (kg water/kg air)", main = "")
hist(pha2_comp[, 22], breaks = 14, xlab = "Humidity seasonality (CV)", main = "")
hist(pha2_comp[, 23], breaks = 14, xlab = "Spatially autocorrelated values", main = "")
hist(pha2_comp[, 24], breaks = 14, xlab = "Random values", main = "")
dev.off()


