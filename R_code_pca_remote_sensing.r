# R_code_pca_remote_sensing.r

library(raster)
library(RStoolbox)

setwd("C:/lab/")

p224r63_2011 <- brick("p224r63_2011_masked.grd")  # import the images composed by different layers

## Landsat bands
# b1 blue
# b2 green
# b3 red
# b4 NIR
# b5 SWIR
# b6 thermal infrared
# b7 SWIR
# b8 panchromatic

# RGB   
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# ggplotRGB
library(ggplot2)
ggRGB(p224r63_2011,5,4,3)  # showing the coordinates

# do the same, with the 1988 image

p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

par(mfrow=c(1,2)) 
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# let's see the correlation between band 1 and band 3
names(p224r63_2011) # names: to see the names of the variables
# "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt"  "B7_sre"
# sre: spectrum reflectance; bt: bit transfer -> different kind of reflectance methods

# see how much the bands are correlated each other
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre) # $: link two variables


# PCA
# decrease the resolution: we have huge amount of pixels
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)

# library(RStoolbox) is now needed (to make PCA)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res) 

cl <- colorRampPalette(c('dark grey','grey','light grey'))(100) # let's plot the $map
plot(p224r63_2011_pca$map, col=cl)

 # how much variation is explained by each component
summary(p224r63_2011_pca$model) # summary: to see the information clearly
# PC1 99.83% of the whole variation

pairs(p224r63_2011) # let’s see how the different bands are correlated to each other

plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretc="Lin")) # 1=PC1; 2=PC2; 3=PC3 

# Repeat fot he 1988 image
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)

summary(p224r63_1988_pca$model)
#99.56% by PCI

pairs(p224r63_1988)

# difference between the two pca maps
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)

## since the PC1 is carrying 99% of the whole information, we can plot only this to see the difference
cldif <- colorRampPalette(c('blue','black','yellow'))(100) 
plot(difpca$PC1,col=cldif)
