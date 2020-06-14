##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Environmental data download and preparation
# Author(s): Marlon E. Cobos
# Date: 17/02/2020 (dd/mm/yyyy)
##################################

# needed packages
pcakages <- c("raster", "rgdal", "rgeos", "maps")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# data
## temperature and prcipitation data 2000 to 2012 (Chelsa-Climate)
dir.create("Env_data_extraction_grouping")
dir.create("Env_data_extraction_grouping/Data") # where to save the data
dir.create("Env_data_extraction_grouping/Data/Precipitation")
dir.create("Env_data_extraction_grouping/Data/Temperature")

years <- c(2000:2012) # years to which the data will be downloaded

## download in loop
for (i in years) {
  for (j in 1:12) {
    if (j < 10) {
      jc <- paste(0, j, sep = "")
    }else {
      jc <- j
    }
    ### precipitation
    namep <- paste0("Data/Precipitation/CHELSA_prec_", # name of files in your computer
                    i, "_", jc, ".tif") 
    urlp <- paste0("https://www.wsl.ch/lud/chelsa/data/timeseries/prec/CHELSA_prec_", # url of files in database
                   i, "_", jc, "_V1.2.1.tif")
    downp <- download.file(url = urlp, destfile = namep, mode = "wb", quiet = TRUE) # donwload the file 
    
    ### temperature
    namet <- paste0("Data/Temperature/CHELSA_tmean_", # name of files in your computer
                    i, "_", jc, ".tif") 
    urlt <- paste0("https://www.wsl.ch/lud/chelsa/data/timeseries/tmean/CHELSA_tmean_", # url of files in database
                   i, "_", jc, "_V1.2.1.tif")
    downt <- download.file(url = urlt, destfile = namet, mode = "wb", quiet = TRUE) # donwload the file
  }
  
  cat("year", i, "of", "2000-2012\n")
}

## humidity data (MERRACLIM)
### downloading the variables (MERRACLIM 2.5 min, mean, 2000s)
download.file(url = "https://datadryad.org/stash/downloads/file_stream/100462", 
              destfile = "2.5m_mean_00s.zip", mode = "wb",
              quiet = FALSE) # donwload the zipped folder 

### subdirectories for organizing information
dir.create("merra_2.5m_mean_00s")

unzip(zipfile = "2.5m_mean_00s.zip", exdir = "merra_2.5m_mean_00s") # unzip
unlink("2.5m_mean_00s.zip") # erasing zipped folder


# processing data
## occurrences
occall <- read.csv("all_samples.csv")
colnames(occall)

# make a spatial object from coordinates using all occs
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
occ_sp <- SpatialPointsDataFrame(coords = occall[, c("longitude", "latitude")], 
                                 data = occall, proj4string = WGS84)

# mapping
map(xlim = c(-35, 85), ylim = c(-28, 25))
points(occ_sp, pch = 16, col = "red")

# project the points using their centriods as reference
centroid <- gCentroid(occ_sp, byid = FALSE)

AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", 
                  centroid@coords[1],
                  " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                  sep = ""))

occ_pr <- spTransform(occ_sp, AEQD)

# create a buffer based on a user-defined distance
buff_area <- gBuffer(occ_pr, width = 100000) # distance in meters
buff_area <- disaggregate(buff_area)
buff_areap <- spTransform(buff_area, WGS84)
plot(buff_areap, add = TRUE)

# crop and mask
## data
temp <- list.files("Env_data_extraction_grouping/Data/Temperature/", 
                   pattern = ".tif$", full.names = TRUE)
tnam <- list.files("Env_data_extraction_grouping/Data/Temperature/", 
                   pattern = ".tif$")
prec <- list.files("Env_data_extraction_grouping/Data/Precipitation/", 
                   pattern = ".tif$", full.names = TRUE)
pnam <- list.files("Env_data_extraction_grouping/Data/Precipitation/", 
                   pattern = ".tif$")
humi <- list.files("merra_2.5m_mean_00s/", pattern = ".tif$", full.names = TRUE)
hnam <- list.files("merra_2.5m_mean_00s/", pattern = ".tif$")


dir.create("Env_data_extraction_grouping/Data/Precipitation_mask")
dir.create("Env_data_extraction_grouping/Data/Temperature_mask")
dir.create("Env_data_extraction_grouping/Data/Humidity_mask")

## masking precipitation and tmperature
for (i in 1:length(temp)) {
  ## read file, mask, write
  ras <- raster(temp[i])
  maskras <- mask(crop(ras, buff_areap), buff_areap)
  writeRaster(maskras, filename = paste0("Env_data_extraction_grouping/Data/Temperature_mask/", 
                                         tnam[i]), format = "GTiff", overwrite = TRUE)
  
  ras <- raster(prec[i])
  maskras <- mask(crop(ras, buff_areap), buff_areap) 
  writeRaster(maskras, filename = paste0("Env_data_extraction_grouping/Data/Precipitation_mask/", 
                                         pnam[i]), format = "GTiff", overwrite = TRUE)
  
  cat("layer", i, "of", length(temp), "\n")
}

## masking humidity
merra_vars <- stack(humi)

### renaming
names(merra_vars) <- gsub("X2_5m_mean_00s_", "", names(merra_vars))

### cropping and masking using bio1 from worldclim to keep only land
bio1_wc <- getData('worldclim', var = 'bio', res = 2.5)[[1]]
bio1_wc <- projectRaster(bio1_wc, merra_vars)

merra_masked <- mask(crop(merra_vars, bio1_wc), bio1_wc)

### saving cropped and masked merra layers
var_names <- paste0("merra_2.5m_mean_00s_crop/", names(merra_masked), ".tif")
saving <- lapply(1:nlayers(merra_masked), function(x) {
  writeRaster(merra_masked[[x]], filename = var_names[x], 
              format = "GTiff", overwrite = TRUE)
})

### masking to buffers
hum <- mask(crop(merra_masked, buff_areap), buff_areap)
var_names <- gsub("merra_2.5m_mean_00s_crop", 
                  "Env_data_extraction_grouping/Data/Humidity_mask", var_names)

saving <- lapply(4:9, function(x) {
  writeRaster(hum[[x]], filename = var_names[x], format = "GTiff", overwrite = TRUE)
})

### placing them in final variables
dir.create("Env_data_extraction_grouping/Data/Final_variables")
writeRaster(hum[[4]], filename = "Env_data_extraction_grouping/Data/Final_variables/hannual.tif", 
            format = "GTiff", overwrite = TRUE)
writeRaster(hum[[7]], filename = "Env_data_extraction_grouping/Data/Final_variables/hseasonality.tif", 
            format = "GTiff", overwrite = TRUE)


# raster calculations to get new variables (precipitation and temperature only)
## mean of temp, standard deviation of temp, sum prec, and coef variation of prec 
temp <- list.files("Env_data_extraction_grouping/Data/Temperature_mask/", 
                   pattern = ".tif$", full.names = TRUE)
prec <- list.files("Env_data_extraction_grouping/Data/Precipitation_mask/", 
                   pattern = ".tif$", full.names = TRUE)

meant <- list()
meanp <- list()

for (j in 1:12) {
  if (j < 10) {
    jc <- paste(0, j, sep = "")
  }else {
    jc <- j
  }
  
  temp <- list.files("Env_data_extraction_grouping/Data/Temperature_mask/", 
                     pattern = paste0(".*", jc, ".tif$"), full.names = TRUE)
  prec <- list.files("Env_data_extraction_grouping/Data/Precipitation_mask/", 
                     pattern = paste0(".*", jc, ".tif$"), full.names = TRUE)
  
  temps <- stack(temp)
  precs <- stack(prec)
  
  meant[[j]] <- calc(temps, mean)
  meanp[[j]] <- calc(precs, mean)
  
  cat(j, "of", length(1:12), "months\n")
}

meantr <- do.call(stack, meant)
meanpr <- do.call(stack, meanp)

tr <- na.omit(values(meantr))
pr <- na.omit(values(meanpr))

mdt <- apply(tr, 1, mean)
mdevt <- meantr[[1]]
mdevt[!is.na(mdevt[])] <- mdt
writeRaster(mdevt, filename = "Env_data_extraction_grouping/Data/Final_variables/tmean.tif", 
            format = "GTiff", overwrite = TRUE)

sdt <- apply(tr, 1, sd)
sdevt <- meantr[[1]]
sdevt[!is.na(sdevt[])] <- sdt
writeRaster(sdevt, filename = "Env_data_extraction_grouping/Data/Final_variables/tseasonality.tif", 
            format = "GTiff", overwrite = TRUE)

sdp <- apply(pr, 1, sum)
sdpr <- meanpr[[1]]
sdpr[!is.na(sdpr[])] <- sdp
writeRaster(sdpr, filename = "Env_data_extraction_grouping/Data/Final_variables/pannual.tif", 
            format = "GTiff", overwrite = TRUE)

cdp <- apply(pr, 1, function(x){sd(x) / mean(x) * 100})
cdpr <- meanpr[[1]]
cdpr[!is.na(cdpr[])] <- cdp
writeRaster(cdpr, filename = "Env_data_extraction_grouping/Data/Final_variables/pseasonality.tif", 
            format = "GTiff", overwrite = TRUE)
