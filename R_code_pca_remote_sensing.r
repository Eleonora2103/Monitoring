# R_code_pca_remote_sensing.r

library(raster)
library(RStoolbox)

setwd("C:/lab/")

p224r63_2011 <- brick("p224r63_2011_masked.grd")

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
ggRGB(p224r63_2011,5,4,3)

# do the same, with the 1988 image

p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011)
# "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt"  "B7_sre"

plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)

# PCA
# decrease the resolution

p224r63_2011_res <- aggregate(p224r63_2011, fact=10)

# library(RStoolbox) is now needed
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)

cl <- colorRampPalette(c('dark grey','grey','light grey'))(100) # 
plot(p224r63_2011_pca$map, col=cl)

summary(p224r63_2011_pca$model)
# PC1 99.83% of the whole variation

pairs(p224r63_2011)
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretc="Lin"))


p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)

summary(p224r63_1988_pca$model)

 #99.56% by PCI
pairs(p224r63_1988)

#difference
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
cldif <- colorRampPalette(c('blue','black','yellow'))(100) 

plot(difpca$PC1,col=cldif)





