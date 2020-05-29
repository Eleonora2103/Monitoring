R_code_snow.r

setwd("C:/lab/")

# install.packages("ncdf4")
install.packages("ncdf4") 

library(raster)
library(ncdf4)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # import image: snow from may

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) # colorRampPalette good for snow

# Excercise: plot the snow cover with the color ramp palette
plot(snowmay, col=cl)

# Slow manner to import the set
setwd("C:/lab/snow")

snow2000 <- raster("snow2000r.tif")
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

#############

# fast version of import and plot of many data for lazy people!
# lapply: apply the function over a list or vector; we want to apply the function raster (used to import single layer)

rlist <- list.files(pattern="snow") # make the list of all the files with common pattern
rlist

import <- lapply(rlist, raster) # import all the files

snow.multitemp <- stack(import) # stack: put alltogether different layer (snow layers)
snow.multitemp
plot(snow.multitemp, col=cl)

# let's make a prediction
source("prediction.r")



