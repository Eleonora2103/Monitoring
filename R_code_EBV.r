# R_code_EBV.r

#upload image
setwd("C:/lab/") 
library(raster) #to see the image

#raster <- # imports a single layer, yet satellite images are made of more than one layer
# satellites images (Sentinel) are composed by many layers
snt <- brick("snt_r10.tif") #also named the picture
snt #see the characteristics of the image
# 30 bits image

#plot the image
plot(snt)
## complex system: passing from bare soil, to forest, opening parts related to agriculture, and forest.

#B1 = blue
#B2 = green
#B3 = red
#B4 = NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin") # visible colours
plotRGB(snt, 4, 3, 2, stretch="Lin") # NIR in top of red: vegetation being coloured in red

pairs(snt) # to see if the axis are related to each other

# The analysis to start from multi system and moving to only one layer is the multivariate analysis
library(RSToolbox)
#PCA anlysis (coming out from the multivariate analysis)
sntPCA <- rasterPCA(snt)
sntPCA

# $: linking the model to the analysis we have done (PCA)
summary(sntPCA$model) # info about output of the model 
# 70% of the info (component 1)
plot(sntPCA$map) # PC1 has the highest amount of information

#plot RGB of the first three component
plotRGB(sntPCA$map, 1, 2, 3, stretch="Lin")

### calculate standard deviation on top of PCA
## set the moving window
## matrix of 5x5 rows and columns with value=1 -> stat that there are 5 
## by 5 pixels moving throughout the image (the value 1 not impact the analysis)
window <- matrix(1, nrow = 5, ncol = 5) 
window # see the 5 by 5 window

# focal (raster function) --> calculates values for a neighbourhood of cells to calculate sd
sd_nst <- focal(sntPCA$map$PC1, w = window, fun =sd)   
## we are going to make the analysis on top of the map and we use only the first principal component; 
## w = window (moving window); fun = function = mean, mode, max, min, sd (standard deviation)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl) # plot the final result of measuring diversity from space

# comparison between the original map and the diversity map
par(mfrow=c(1,2))  # two images beside to the other
plotRGB(snt,4,3,2, stretch="lin", main = "Original Image")
plot(sd_snt, col=cl, main = "Diversity")

####### day 2: Cladonia example
# R_code_cladonia_focal.rFile

setwd("C:/lab/")
library(raster)

# import the images: two functions --> raster (one single band); brick (import several layer at time)
# in this case we have directly RGB image (made of three layers) so we made use of the brick function
clad <- brick("cladonia_stellaris_calaita.JPG")

plotRGB(clad, 1, 2, 3, stretch="lin") # RGB bands: red, green, blue

# we make the focal function on top of the image to see where there is the major variability
# set the moving window 
window <- matrix(1, nrow = 3, ncol = 3) # window in which we make the calculation (3x3 pixels); 
# 1: set a value for the pixels which is not impacting further calculation
window

pairs(clad) # to see if the bands are related to each other

### PCA analysis
library(RStoolbox)

cladpca <- rasterPCA(clad)
cladpca # information about pca

summary(cladpca$model) #let's see how much information is explaining the PC1
# 98% of the info (component 1):  this is common since it is an image taken into the visible spectrum 
# (in the visible the bands are really correlated to each other)
plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

# focal (raster function) --> calculates values for a neighbourhood of cells to calculate sd
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd) # do the calculation of the diversity into a certain band of the image

PC1_agg <- aggregate(cladpca$map$PC1, fact=10) # accelerate the calculation of sd
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd) # make the focal function on top of the pc1 agrregated

# plot the two different set we made
par(mfrow=c(1,2)) # two images beside to the other in multiframe
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl) # cladonia set aggregated

# plot the calculation
par(mfrow=c(1,2)) 
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)
# plot(sd_clad_agg, col=cl)
