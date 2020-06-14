##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Detecting levels of spatial autocorrelation of A. gambiae localities
# Author(s): Marlon E. Cobos
# Date:  24/01/2019 (dd/mm/yyyy)
##################################

# packages
pcakages <- c("raster", "rgdal")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# data
## occurrences
occ <- read.csv("15_sample_information.csv")

## raster data processing
### mean temperature
tmax <- list.files(path = "chelsa-climate", pattern = "tmax", full.names = TRUE)
tmin <- list.files(path = "chelsa-climate", pattern = "tmin", full.names = TRUE)

stn <- gregexpr("\\d\\d\\d\\d_\\d\\d", tmin)
stan <- regmatches(tmin, stn)
statn <- unlist(stan)

for (i in 4:length(tmax)) {
  stackr <- stack(tmax[i], tmin[i])
  meanr <- calc(stackr, mean)
  writeRaster(meanr, paste0("chelsa-climate/CHELSA_tmean_", statn[i], "_V1.2.1.tif"),
              format = "GTiff")
  cat("raster", i, "of", length(tmax), "processed\n")
}


### two months processing
rlist <- as.character(occ$tif_date)

stn <- gregexpr(".*-.*", rlist)
stan <- regmatches(rlist, stn)
statn <- unlist(stan)

rasters <- unique(statn)

for (i in 1:length(rasters)) {
  ras <- strsplit(rasters[i], "-")[[1]]
  
  rp1 <- list.files(path = "chelsa-climate", pattern = paste0(".*prec_", ras[1]), 
                    full.names = TRUE)
  rp2 <- list.files(path = "chelsa-climate", pattern = paste0(".*prec_", ras[2]), 
                    full.names = TRUE)
  
  rt1 <- list.files(path = "chelsa-climate", pattern = paste0(".*tmean_", ras[1]), 
                    full.names = TRUE)
  rt2 <- list.files(path = "chelsa-climate", pattern = paste0(".*tmean_", ras[2]), 
                    full.names = TRUE)
  
  meanpr <- calc(stack(rp1, rp2), mean)
  meantr <- calc(stack(rt1, rt2), mean)
  
  writeRaster(meanpr, paste0("chelsa-climate/CHELSA_prec_", rasters[i], "_V1.2.1.tif"),
              format = "GTiff")
  writeRaster(meantr, paste0("chelsa-climate/CHELSA_tmean_", rasters[i], "_V1.2.1.tif"),
              format = "GTiff")
  
  cat("raster", i, "of", length(rasters), "processed\n")
}

## info
names(occ)
occc <- occ[c("ox_code", "tif_date", "longitude", "latitude")]

temp <- vector()
prec <- vector()

rlist <- paste0("_", rlist, "_")

for (i in 1:dim(occc)[1]) {
  p <- list.files(path = "chelsa-climate", pattern = paste0(".*prec", rlist[i]), 
                  full.names = TRUE)
  t <- list.files(path = "chelsa-climate", pattern = paste0(".*tmean", rlist[i]), 
                  full.names = TRUE)
  
  if (i > 1) {
    if (rlist[i - 1] == rlist[i]) {
      
    } else {
      tp <- raster(t)
      pr <- raster(p)
    }
  } else {
    tp <- raster(t)
    pr <- raster(p)
  }
  
  temp[i] <- extract(tp, occc[i, 3:4])
  prec[i] <- extract(pr, occc[i, 3:4])
  
  cat("process", i, "of", dim(occc)[1], "finished\n")
}

temp <- temp / 100
occs <- cbind(temp, prec)

table <- data.frame(occc, temp, prec)

dir.create("Spatial_autocorrelation")
write.csv(table, file = "Spatial_autocorrelation/15_occurrences_values.csv", 
          row.names = FALSE)

# distances
## geographic
gdist <- pointDistance(occc[, 3:4], lonlat = T)
gdistv <- na.omit(c(gdist))[na.omit(c(gdist)) != 0]/1000

## environmental distance
occss <- scale(occs)
edist <- dist(occss)

edistv <- c(edist)

# relationship 
glmd <- glm(edistv ~ gdistv)
glmd
summary(glmd)

## relationship results
sink("Spatial_autocorrelation/glm_results.txt")
cat("GLM results\n")
glmd

cat("\n\nGLM summary\n")
summary(glmd)
sink()

# plot distances
x11()
par(mar = c(4.1, 4, 0.5, 0.5))
plot(gdistv, edistv, pch = 19, col = "grey25", xlab = "Geographic distance (km)", 
     ylab = "Environmental distance")

# saving results
png(paste("Spatial_autocorrelation/A_gambiae_dist_g_e.png", sep = ""), width = 80, 
    height = 80, units = "mm", res = 600)
par(mar = c(4.1, 4, 0.5, 0.5), cex = 0.7)
plot(gdistv, edistv, pch = 19, col = "grey25", xlab = "", ylab = "", cex = 0.8)
title(xlab = "Geographic distance (km)", ylab = "Environmental distance", cex.lab = 1.2)
dev.off()
