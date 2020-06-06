# R_code_faPAR.r
# how to look at chemical cycling from sateliites

# levelplot(copNDVI)

library(raster)
library(rasterVis) # used to make levelplot
library(rasterdiv) 

plot(copNDVI) # Copernicus NDVI, present in the rasterdiv library
copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # remove data from 253 to 255 and put no data as value --> remove water from the analysis
levelplot(copNDVI)

setwd("C:/lab/")

faPAR10 <- raster("faPAR10.tif") # import the image; faPAR10 because it is aggregate by a fact=10
levelplot(faPAR10)

pdf("copNDVI.pdf") # make a pdf file; quotes because we are saving the file outside R
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

### Regression model between faPar and NDVI

erosion <- c(12, 14, 16, 24, 26, 40, 55, 67)
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")

model1 <- lm(hm ~ erosion)
summary(model1)
abline(model1)

### faPAR vs. NDVI model
setwd("C:/lab/")

library(raster)
faPAR10 <- raster("faPAR10.tif")

library(rasterdiv)
plot(faPAR10)
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

# make random sample
library(sf) # to call st_* functions

random.points <- function(x,n) # x = raster file; n = number of point that we want select
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}

pts <- random.points(faPAR10,1000)

# make the extraction of the faPAR and copNDVI from a raster; pass values from a maps 
copNDVIp <- extract(copNDVI, pts) # copNDVI points = 1000 points where we put the data
faPAR10p <- extract(faPAR10,pts) # faPAR points

# photosynthesis(y) vs. biomass(x)
model2 <- lm(faPAR10p ~ copNDVIp) # lm = to fit linear models
summary(model2)
## R-squared: 0.4279 → the point are not perfectly in line but they are far from being random, they are related to each other
## p-value: < 2.2e-16 --> the relation is confirmed
plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")




