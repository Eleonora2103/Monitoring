# R code faPAR
# how to look 


library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA))
levelplot(copNDVI)

setwd("C:/lab/")

faPAR10 <- raster("faPAR10.tif")

levelplot(faPAR10)

pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

load("faPAR.RData")

# the original faPAR from Copernicus is 2GB

# let's see how much space is needed for an 8-bit set

library(raster)
library(rasterdiv)

writeRaster(copNDVI, "copNDVI.tif")
# 5.3MB

library(rasterVis)

# faPAR: levelplot this set
levelplot(faPar)


