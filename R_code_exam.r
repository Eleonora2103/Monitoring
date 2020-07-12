# R_code_exam.r

# 1. R first code
# 2. R code multipanel
# 3. R code spatial
# 4. R code point pattern
# 5. R code for multivariate analysis
# 6. R code_rs
# 7. R code_ecosystem functions
# 8. R code_pca_remote_sensing
# 9. R code_faPAR
# 10. R code_radiance
# 11. R code_EBV
# 12. Rcode_snow
# 13. R code_no2
# 14. R code_interpolation
# 15. R code_sdm

# 1. R_code_first

install.packages("sp")

library(sp)

data(meuse)

head(meuse)

attach(meuse)

plot(zinc,copper)
plot(zinc,copper,col="red")
plot(zinc,copper,col="red",pch=19)
plot(zinc,copper,col="red",pch=19,cex=2)

##############################################
##############################################
##############################################

# 2. R_code_mutipanel.r

### Multipanel in R: the second lecture of monitoring Ecosystems

install.packages("sp")
install.packages("GGally")

library(sp)
library(GGally)

data(meuse)

attach(meuse)

head(meuse)

plot(cadmium,zinc)

plot(cadmium,zinc,pch=15)

plot(cadmium,zinc,pch=15,col="red")

pairs(meuse)

pairs(~cadmium+copper+lead+zinc,data=meuse)

pairs(meuse[,3:6])

pairs(meuse[,3:6],pch=19)

pairs(meuse[,3:6],pch=19,col="green")

pairs(meuse[,3:6],pch=19,col="green",cex=2)

# GGally package will prettify the graph
ggpairs(meuse[,3:6])

###################################################
###################################################
###################################################

# 3. R_code_spatial.r

## R code for spatial view points

library(sp)  # you want to use that packages

data(meuse)

head(meuse)

# coordinates
coordinates(meuse) = ~x+y

plot(meuse)

spplot(meuse, "zinc")   

# Plot the spatial amount of copper and change thetitle

spplot(meuse, "copper")

spplot(meuse, "copper", main = "Copper concentration")

bubble(meuse, "zinc")

bubble(meuse, "zinc", main = "Zinc COncentration")

# Exercise: bubble copper in red

bubble(meuse, "copper", main = "Copper concentration", col = "red")

#### Impotring new data

# download covid agg.csv from out teaching site and build a folder called lab into C:
#put the covid_agg.csv file into the folder lab

# setting the working directory: lab

setwd("C:/lab/")             # set WorkingDirectly

covid <- read.table("covid_agg.csv", head=T)       # arrows link a function to an object

head(covid)

attach(covid)
plot(country,cases)

plot(country, cases, las=0)  # parallel label

plot(country, cases, las=1)  # label horizontal

plot(country, cases, las=2)  # label perpendicular

plot(country, cases, las=3)  # vertical labels

plot(country, cases, las=3, cex.axis=0.5)

plot(country, cases, las=3, cex.axis=0.7

plot(country, cases, las=3, cex.axis=0.5)

# ggplot2 package
install.packages("ggplot2")

library(ggplot2)  #require(ggplot2)
     
# load the previously saved .Rdata

setwd("C:/lab/")   

load("ggplot.RData")
     
ls()    # list of objects 

library(ggplot2)

data(mpg)
head(mpg)
     
# key components: data, aes, geometry
     
ggplot(mpg, aes(x=displ,y=hwy)) + geom_point() # geomestry is declared as separate
ggplot(mpg, aes(x=displ,y=hwy)) + geom_line()
ggplot(mpg, aes(x=displ,y=hwy)) + geom_polygon()
     
head(covid)
     
ggplot(covid, aes(x=lon,y=lat, size=cases)) + geom_point()
############################################
############################################
############################################
     
# 4. R_code_point_pattern_analysis.r

## point pattern analysis: density map

install.packages("spatstat")

library(spatstat)

attach(covid)
head(covid)

covids <- ppp(lot, lan, c(-180,180), c(-90,90))     # panel point pattern

d <- density(covids)

plot(d)
points(covids)

setwd("C:/lab/")
load("point_pattern.RData")

ls()

# covids: point pattern
# d: density map
library(spatstat)

plot(d)
points(covids)

install.packages("rgdal")

library(rgdal)

# let's input vector ines (x0y0, x1y1, x2y2..)
coastline <- readOGR("ne_10m_coastline.shp")

plot(coastlines, add=T)         # add:adding the lines to the previous image)

# change of the colour, we are making ramp palette   c=array of colours   (100): all the possible colours from yellow to red
# converting hand-designed `sequential' or `diverging' color schemes into continuous color ramps
cl <- colorRampPalette(c("yellow","orange","red"))(100)
plot(d, col=cl)
points(covids)
plot(coastlines, add=T)

# Exercise: new colour ramp palette
clr <- colorRampPalette(c("light green","yellow","orange","violet"))(100)
plot(d, col=clr, main="densitites of covid-19")
points(covids)
plot(coastlines, add=T)

# export as pdf or png("covid_density.png")
pdf("covid_density.pdf")
clr <- colorRampPalette(c("light green","yellow","orange","violet"))(100)
plot(d, col=clr, main="densitites of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()
     
##############################################
##############################################
##############################################

# 5. R_code_multivariate.r

## R code for multivariate analysis

setwd("C:/lab/")

# vegan = vegetation analysis

library(vegan)

biomes <- read.table("biomes.csv", header=T, sep=",")   # csv is an extension

head(biomes)

# DEtrended CORrespondence ANAlysis 

multivar <- decorana(biomes)
plot(multivar)
multivar

plot(multivar)

# biomes types
biomes_types <- read.table("biomes_types.csv", header=T, sep=",") 

head(biomes_types)

attach(biomes_types)
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3)            # ellipse connecting all the plots of same biomes

# also col=c("green","blue","red","black")

ordispider(multivar, type, col=1:4, label=T)

######################################
######################################
######################################
     
# 6. R_code_rs.r
     
## R code for remote sensing data analysis

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

# plotRGB 
## RGB: to show an image using the visible part of the spectrum we can put every single band into the  
## correspondent component of RGB spectrum. What we are doing is putting the exact bands into each  
## components (of the Landsat band) to see the colours we might see with human eyes.

# dev.off to close the previous part

plotRGB(p224r63_2011, r=3, g=2, b= 1, stretch="Lin")        
# stretch= stretching the colour to see it better (Lin: linear) --> pass from a 
# minimum 50 to a maximum of 100 to a minimum of 0 to a maximum of 100

# highlight the vegetation (3 bend at time) respect to the other parts of the image
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") 

# NIR on top of the G component of the RGB
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") 

# Exercise NIR on top of the B component of the RGB 
plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin") 

############# day 2: data from 1988
setwd("C:/lab/")

load("rs.RData")  # import the data done the previous day
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
# Mode1: stretch again the colour --> hist: we are going to shock the image with a huge amount of 
# different colours immediately from the smallest to the greatest value
# Mode2: multivariate analysis
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
     
###########################################
###########################################
###########################################

# 7. R_code_ecosystem_functions.r

## R code to view biomass over the world and calculate changes in ecosystem functions
## energy
## chemical cycling
## proxies

install.packages("rasterdiv")  # Diversity based Raster Data
install.packages("rasterVis")  # Methods for enhanced visualization and interaction with raster data

library(rasterVis) # to make levelplot
library(rasterdiv)

data(copNDVI) # cop: COpernicus
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA))  # removing water-based colours

levelplot(copNDVI) # Draw Level Plots and Contour plots

copNDVI10 <- aggregate(copNDVI, fact=10)
levelplot(copNDVI10)


copNDVI100 <- aggregate(copNDVI, fact=100)
levelplot(copNDVI100)

###########################
library(ggplot2)

myPalette <- colorRampPalette(c('white','green','dark green'))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))


ggR(copNDVI, geom_raster = TRUE) +
scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
labs(x="Longitude",y="Latitude", fill="")+
# theme(legend.position = "bottom") +
# NULL
# +
# ggtitle("NDVI")


setwd("C:/lab/")
library(raster)

defor1 <- brick("defor1_.jpg")  # function to import images several layers
defor2 <- brick("defor2_.jpg")

# band1: NIR, defor1_.1, defor2_.1
# band2: red,  defor1_.2, defor2_.2
# band3: green

# NIR mounted on the red component
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin") # Red, Green and Blue component
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

dvi1 <- defor1$defor1_.1 - defor1$defor1_.2
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2

cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
par(mfrow=c(1,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

difdvi <- dvi1 - dvi2

dev.off()   # different with difference extents
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld)

hist(difdvi)
     
#######################################
#######################################
#######################################

# 8. R_code_pca_remote_sensing.r

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

#################################
#################################
#################################

# 9. R_code_faPAR.r

## how to look at chemical cycling from sateliites

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

#####################################
#####################################
#####################################

# 1. R_code_radiance.r

## bit example
library(raster)  # format where there are pixels

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13, 1.44, 1.55, 3.4) # put data into iìeach pixel

plot(toy)
text(toy, digits=2) # add the values of toy, using 2 decimal degrees

## we are going to force this data towards lower amount of bits
## Image of 2 bits how it looks like
## assign a code (number) from 0 to 3; with 4 values we are ranging from 0 to 3
toy2bits <- stretch(toy,minv=0,maxv=3) # stretch: function to change the range of values
storage.mode(toy2bits[]) = "integer" # use only integer numbers
plot(toy2bits)
text(toy2bits, digits=2)

# toy with a range of 4 bits (16 combinations)
toy4bits <- stretch(toy,minv=0,maxv=15) # range from 0 to 15
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)

# toy with a range of 4 bits (256 combinations)
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2)

# plot altogether
par(mfrow=c(1,4))
plot(toy)
text(toy, digits=2)

plot(toy2bits)
text(toy2bits, digits=2)

plot(toy4bits)
text(toy4bits, digits=2)

plot(toy8bits)
text(toy8bits, digits=2)

dev.off
library(rasterdiv) 
plot(copNDVI) # see the range of data (bits) of copNDVI
## values ranging from 0 to 255 --> 8 bits images (most used - sometimes 16 bits)

#######################################
#######################################
#######################################

# 11. R_code_EBV.r

## upload image
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

##########################################
##########################################
##########################################

# 12. R_code_snow.r

setwd("C:/lab/")

# install.packages("ncdf4")
install.packages("ncdf4") 

library(raster)
library(ncdf4)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # import image: snow from May

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

########### day 2
setwd("C:/lab/snow")

library(raster)
library(ncdf4)

# Excercise: import all of the snow cover images altogether

rlist <- list.files(pattern="snow")

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

# load("name.RData")
prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# Export the output
# you made the calculation and you want to send the output to a collegue

writeRaster(prediction, "final.tif")  # writing the entire raster object to a file

# final stack
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack, col=cl)

# export the R graph
pdf("my_final_exciting_graph") # create .pdf file
plot(final.stack, col=cl)
dev.off

png("my_final_exciting_graph") # create .png file
plot(final.stack, col=cl)
dev.off

##############
setwd("C:/lab/")

# install.packages("ncdf4")
install.packages("ncdf4") 
library(ncdf4)
library(raster)

snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # import image

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

ext <- c(0, 20, 35, 50)
zoom(snow, ext=ext)

snowitaly <- crop(snow, ext)
plot(snowitaly, col=cl)

zoom(snow, ext=drawExtent())

########################################
########################################
########################################

# 13. R_code_no2.r

 setwd("C:/lab/no2/")

library(raster)

# import the files through the lapply function
rlist <- list.files(pattern="EN")
rlist

import <- lapply(rlist, raster)
EN <- stack(import)
cl <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(EN, col=cl)

# january and march
par(mfrow=c(1,2))
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)

#RGB space
plotRGB(EN, r=1, g=7, b=13, stretch="lin")

# difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(dif, col=cld)

# quantitative estimate: do the box-plot
boxplot(EN)
boxplot(EN, outline=F, horizontal=T)
boxplot(EN,outline=F, horizontal=T, axes=T)

# plot!
plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red") # y=a+bx

setwd("C:/lab/snow/")
# import the snow cover imeages altogether

# fast version of import and plot of many data for lazy people!
rlist <- list.files(pattern="snow20")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)

snow.multitemp
plot(snow.multitemp)

plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1, col="red")

plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r)
abline(0,1,col="red")

######################################
######################################
######################################

# 15. R_code_interpolation.r

## Beach forest - Casentino

setwd("C:/lab/")
 
library(spatstat) # Spatial Point Pattern Analysis
 
# import data - since it is only a table, we will use the read.table function
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T) 
head(imp) # to see the first 6 rows of the dataset

attach(inp)
# plot(inp$X, inp$Y)
plot(X,Y)

summary(inp)
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000)) # Planar Point Pattern

names(inp)
marks(inppp) <- Canopy.cov # label the sigle point

# estimate the canopy cover
canopy <- Smooth(inppp) # interpolate the data (smooth) of the inppp set where they have not been measured
# Warning: adapt the cross validation: chacking the new values with the original one
plot(canopy)
points(inppp, col="green")

marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

# plot the variables together
## output <- stack(canopy,lichs)
## plot(lichs) # don't work due to a lenght difference
par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)

# final plot of the two variables
par(mfrow=c(1,3))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)

dev.off()

############ Psammophilus forest
inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T) 
attach(inp.psam)

head(inp.psam)
plot(E, N)

summaru(inp.psam)
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))

marks(inp.psam.ppp) <- C_org

C <- Smooth(inp.psam.ppp)
 
plot(E,N)

C <- Smooth(inp.psam.ppp)
plot(C)
points(inp.psam.ppp)

####################################
####################################
####################################

# 15. R_code_sdm.r 

install.packages("sdm") # Species Distribution Modelling

library(sdm) 
library(raster) # predictors; make use of the raster set
library(rgdal)# species; manage coordinates and spactial data

# import the file
file <- system.file("external/species.shp", package="sdm") 
## external: there are species data and environmental variables data (inside the sdm package)

# use the graphical part of the file
species <- shapefile(file) 
plot(species) # present/absent plot
species # information about species

plot(species[species$Occurrence == 1,],col='blue',pch=16)

path <- system.file("external", package="sdm")

lst <- list.files(path=path,pattern='asc$',full.names = T) #
lst

preds <- stack(lst)
plot(preds)

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100) # change the colorRampPalette
plot(preds, col=cl)

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# model
d <- sdmData(train=species, predictors=preds) # explain which are the data

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")
p1 <- predict(m1, newdata=preds)

plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)
 
s1 <- stack(preds, p1)
plot(s1, col=cl)
