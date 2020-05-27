# R code for remote sensing data analysis

#raster
install.packages("raster")
install.packages("RStoolbox")

setwd("C:/lab/")

#raster <- # imports a single layer, yet satellite images are made of more than one layer
library(raster)

# import images

p224r63_2011 <- brick("p224r63_2011_masked.grd") # brick: import images with several layer at time

plot(p224r63_2011)

cl <- colorRampPalette(c('black','grey','light grey'))(100)    # c= concatenate

# bands of Landstat
# B1: blue
# B2: green
# B3: red
# B4:NIR

# multiframe of different plots       
# c= series of characters

par(mfrow=c(2,2))

# B1: blue
clb <- cl <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)    # $= link different part (band to imagine)

# B2: green 
clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

# B3: red
clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

# B4: NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

# let's change the par
par(mfrow=c(4,1))

clb <- cl <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)    # $= link different part (band to imagine)

clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

#plotRGB
# dev.off to close the previous part
plotRGB(p224r63_2011, r=3, g=2, b= 1, stretch="Lin")         # stretch= stretching the colour to see it better (Lin: linear)

# highlight the vegetation (3 bend at time) respect to the other parts of the image
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") 

# NIR on top of the G component of the RGB
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") 

# Exercise NIR on top of the B component of the RGB 
plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin") 

#############
setwd("C:/lab/")

load("rs.RData")
ls()

p224r63_1988_masked.grd       # masked: no data where there is water

library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd")

plot(p224r63_1988)

# plot in RGB visible 321 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") 

# Plot in false colour RGB 432 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=2, b=2, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") 

# enhance the noise
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist") 

# bands of Landstat
# B1: blue
# B2: green
# B3: red: B3_sre
# B4:NIR: B4_sre

dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi2011, col=cl)

# dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi1988, col=cl)

# difference from one year to the other
diff <- dvi2011 - dvi1988
plot(diff)

# Let's change the grain (dimension of pixel) of our images  # res: resempling

p224r63_2011res <- aggregate(p224r63_2011, fact=10)   # fact: amount of time we want to increase our images (10 time the image)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)

par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")
 
