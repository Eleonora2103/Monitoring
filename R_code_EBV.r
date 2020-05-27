# R_code_EBV.r

#upload image
setwd("C:/lab/") 
library(raster) #to see the image

#raster <- # imports a single layer, yet satellite images are made of more than one layer
snt <- brick("snt_r10.tif") #also named the picture
snt #see the characteristics of the image

#plot the image
plot(snt)

#B1 = blue
#B2 = green
#B3 = red
#B4 = NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin")
plotRGB(snt, 4, 3, 2, stretch="Lin")

pairs(snt)

library(RSToolbox)
#PCA anlysis
sntPCA <- rasterPCA(snt)
sntPCA

summary(sntPCA$model) # info about output of the model
# 70% of the info (component 1)
plot(sntPCA$map)

#plot RGB
plotRGB(sntPCA$map, 1, 2, 3, stretch="Lin")

###calculate standard deviation
#set the moving window
window <- matrix(1, nrow = 5, ncol = 5)

#focal --> calculates values for a neighbourhood of cells
sd_nst <- focal(sntPCA$map$PC1, w = window, fun =sd)   #we use tha map and the frst principal component; w = window; fun = function = mean, mode, max, min, sd (standar deviation)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt,4,3,2, stretch="lin")
plot(sd_snt, col=cl)


####### day 2: cladonia example
# R_code_cladonia_focal.rFile

setwd("C:/lab/")
library(raster)

# import the images: two functions --> raster (one single band); brick (import several layer at time)
clad <- brick("cladonia_stellaris_calaita.JPG")

plotRGB(clad, 1, 2, 3, stretch="lin")

window <- matrix(1, nrow = 3, ncol = 3) # window in which we make the calculation
window

pairs(clad)

library(RStoolbox)

### PCA analysis
cladpca <- rasterPCA(clad)
cladpca

summary(cladpca$model)
# 98% of the info (component 1)
plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd) # do the calculation of the diversity into a certain band of the image

PC1_agg <- aggregate(cladpca$map$PC1, fact=10) # accelerate the calculation of sd
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

# plot the two different set we made
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)


par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)
# plot(sd_clad_agg, col=cl)




 
