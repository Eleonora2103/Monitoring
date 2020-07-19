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

# Download and install the packages 
# sp package = to show the points in a map
install.packages("sp") 

# attach an add-on package
library(sp)

# load specific dataset
data(meuse)
# meuse: dataset comprising of four heavy metals (Copper, Zinc, Cadmium, Lead) measures in the top soil along the river Meuse

meuse # let's see how the meuse data set is structured 

# see the first six lines of the whole data set
head(meuse)

# the database is attached in R
attach(meuse)

# plot the database
# Let's plot two variables together to see if they are related
# See if zinc concentration is related to that of copper
plot(zinc,copper)
# col = specify the colours of the point
plot(zinc,copper,col="red") # quotes are used for each object that is outside of R 
# pch = used to specify point shapes
plot(zinc,copper,col="red",pch=19)
# cex = used to make character exageration
plot(zinc,copper,col="red",pch=19,cex=2)

##############################################
##############################################
##############################################

# 2. R_code_mutipanel.r

### Multipanel in R: the second lecture of monitoring Ecosystems

install.packages("sp")
install.packages("GGally")
# GGally is an extension of ggplot2; it adds functions to reduce the complexity of combining geometric objects with transformed data

library(sp)
library(GGally)

data(meuse)
attach(meuse)

# names: function to get the name of all variables in the dataset
names(meuse)
head(meuse)

plot(cadmium,zinc)

plot(cadmium,zinc,pch=15)

plot(cadmium,zinc,pch=15,col="red",cex=2)
# pch = used to specify point shapes
# col = specify the colours of the point
# cex = used to make character exageration

# create a matrix of scatterplot
pairs(meuse)

# error = "the sizze is too large" --> reshape with the mouse the graph window and relaunch the code

pairs(~cadmium+copper+lead+zinc,data=meuse)
# create scatterplots using only Cadmium, Copper, Lead, and Zinc

pairs(meuse[,3:6])
# from coloumn 3 to 6
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

# function to plot the dataset
plot(meuse)

# spplot = plots several layers with a single legend for all maps
spplot(meuse, "zinc") # plot the spatial amount of zinc

# Plot the spatial amount of copper
spplot(meuse, "copper")

spplot(meuse, "copper", main = "Copper concentration")
# main = change the title

# bubble = creates a bubble plot of spatial data
bubble(meuse, "zinc")

bubble(meuse, "zinc", main = "Zinc COncentration")

# Exercise: bubble copper in red

bubble(meuse, "copper", main = "Copper concentration", col = "red")

#### Importing new data

# download covid agg.csv from out teaching site and build a folder called lab into C:
# put the covid_agg.csv file into the folder lab

# setting the working directory: lab

setwd("C:/lab/") # set Working Directory, for Windows

covid <- read.table("covid_agg.csv", head=T) 
# the arrow (<-) links a function to an object
# read.table: reads a file in table format and creates a data frame from it
# head = T/TRUE --> the first line is the header of the table

head(covid) # let's see the first six line of the table

attach(covid)
# let's plot the table
plot(country,cases)

plot(country, cases, las=0)  # parallel label

plot(country, cases, las=1)  # label horizontal

plot(country, cases, las=2)  # label perpendicular

plot(country, cases, las=3)  # vertical labels

# cex.axis = changes the size of the axis annotation
plot(country, cases, las=3, cex.axis=0.5)

plot(country, cases, las=3, cex.axis=0.7

plot(country, cases, las=3, cex.axis=0.5)

# ggplot2 package to create graphics
install.packages("ggplot2")

library(ggplot2)  # require(ggplot2)
     
# load the previously saved .Rdata

setwd("C:/lab/")   

load("ggplot.RData")

# ls = shows what data sets and functions we defined and used in the recalled workingspace
ls() # list of objects 

library(ggplot2)

data(mpg)
# mpg: provides fuel economy data from 1999 to 2008 for 38 popular model of cars
head(mpg)
     
ggplot(mpg, aes(x=displ,y=hwy)) + geom_point() # point geometry 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_line() # line geometry
ggplot(mpg, aes(x=displ,y=hwy)) + geom_polygon() # polygon geometry
# key components: data, aes, geometry
# aes: declare the variables we want to plot
# geometry: declare the geomentry of the plot; they are declared as separate

head(covid)
     
ggplot(covid, aes(x=lon,y=lat, size=cases)) + geom_point() 
# get the size of each element from objects, in this case x
# the size of the points change according to the cases
     
############################################
############################################
############################################
     
# 4. R_code_point_pattern_analysis.r

## point pattern analysis: density map

install.packages("spatstat")
# Spatial Point Pattern Analysis: statystical analysis of spatial point pattern
     
library(spatstat)

# attach the Covid dataset
attach(covid)
head(covid)

covids <- ppp(lot, lan, c(-180,180), c(-90,90)) # planar point pattern, specify the variables and the range
# lot, lan: explain the x and y variables and their range
     
# density map
d <- density(covids)

# plot and show the points within the density map
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

install.packages("rgdal") # rgdal: provides bindings to the 'Geospatial' Data Abstraction Library

library(rgdal)

# let's input vector ines (x0y0, x1y1, x2y2..)
coastline <- readOGR("ne_10m_coastline.shp") 

plot(coastlines, add=T) # add: add the lines to the previous image


cl <- colorRampPalette(c("yellow","orange","red"))(100)
# colorRampPalette: converting hand-designed `sequential' or `diverging' color schemes into continuous color ramps
# change of the colour, we are making ramp palette   
# c: array of colours   
# (100): all the possible colours from yellow to red
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
# shut down the device
dev.off()
     
##############################################
##############################################
##############################################

# 5. R_code_multivariate.r

## R code for multivariate analysis

setwd("C:/lab/")

# vegan = VEGetation ANalysis
library(vegan)

biomes <- read.table("biomes.csv", header=T, sep=",")   # csv is an extension
# sep="," -> separator between the names of the words within the rows
head(biomes) # have a look to the data base

# Multivariate Analysis
# DEtrended CORrespondence ANAlysis 
multivar <- decorana(biomes)
plot(multivar)
multivar

plot(multivar)

# biomes types
biomes_types <- read.table("biomes_types.csv", header=T, sep=",") 
# assign name to a function (<-) + read.table to import from outside (name) + header = true bc first line is header + separator is the comma, so we state it

head(biomes_types)

attach(biomes_types)
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3) 
# ordiellipse: ellipse connecting all the plots of same biomes
# col: specify the colours of the point
# kind: type of grouping
# lwd: dimension of the line of the ellipse

# also col=c("green","blue","red","black")

ordispider(multivar, type, col=1:4, label=T) # draw a 'spider' diagram where each point is connected to the central group with segment

######################################
######################################
######################################
     
# 6. R_code_rs.r
     
## R code for remote sensing data analysis

install.packages("raster")
# raster:  raster is the format with the pixels so matrices of the row, columns and values
install.packages("RStoolbox")
# RStoolbox: package for remote sensing image, processing and analysis such as calculating spectral indeces, principal component transformation

setwd("C:/lab/")

# import images
library(raster)
#raster:  imports a single layer, yet satellite images are made of more than one layer

p224r63_2011 <- brick("p224r63_2011_masked.grd") 
# brick: import images with several layer at time

plot(p224r63_2011)

cl <- colorRampPalette(c('black','grey','light grey'))(100) # c: concatenate

# bands of Landstat
# B1: blue
# B2: green
# B3: red
# B4: NIR

par(mfrow=c(2,2))
# multiframe of different plots       
# c= series of characters
     
# B1: blue
clb <- cl <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb) # $= link different part (band to image)

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
plot(p224r63_2011$B1_sre, col=clb) # $= link different part (band to imagine)

clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)
dev.off()
     
## plotRGB 
## RGB: to show an image using the visible part of the spectrum; we can put every single band into the  
## correspondent component of RGB spectrum. What we are doing is putting the exact bands into each  
## components (of the Landsat band) to see the colours we might see with human eyes.

# dev.off to close the previous part

plotRGB(p224r63_2011, r=3, g=2, b= 1, stretch="Lin")        
# stretch: stretching the colour to see it better (Lin: linear) --> pass from a 
# minimum 50 to a maximum of 100 to a minimum of 0 to a maximum of 100

# NIR on top of the R component of the RGB
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

p224r63_1988_masked.grd # masked: no data where there is water

library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd")

plot(p224r63_1988)

# plot in RGB visible 321 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") 

# Plot in false colour RGB 432 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
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
# B4: NIR: B4_sre
     
# we are going to calculate the DVI
# DVI 2011
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
p224r63_2011res <- aggregate(p224r63_2011, fact=10)   
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)
# fact: amount of time we want to increase our images (10 time the image)
     
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

install.packages("rasterdiv")  # diversity based Raster Data
install.packages("rasterVis")  # methods for enhanced visualization and interaction with raster data

library(rasterVis) # to make levelplot
library(rasterdiv)

data(copNDVI) # cop: COpernicus NDVI
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # removing water-based colours
# reclassify: function that reclassify groups of values to other values

levelplot(copNDVI) # Draw Level Plots and Contour plots

copNDVI10 <- aggregate(copNDVI, fact=10) # fact: aggregating 10 pixel in 1 
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

##########################

# DVI deforestation
setwd("C:/lab/")
library(raster)

defor1 <- brick("defor1_.jpg")  # function to import images with several layers
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

# calculation of DVI
dvi1 <- defor1$defor1_.1 - defor1$defor1_.2
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2

cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
par(mfrow=c(1,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

# Let's see the difference between dvi1 and dvi2
difdvi <- dvi1 - dvi2

dev.off()   
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld)

# make an histogram
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

# image from 2011
# RGB   
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# ggplotRGB
library(ggplot2)
ggRGB(p224r63_2011,5,4,3)  # show the coordinates

# do the same, with the 1988 image

p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

par(mfrow=c(1,2)) 
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

# let's see the correlation between band 1 and band 3
names(p224r63_2011) # names: to see the names of the variables
# "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt" "B7_sre"
# sre: spectrum reflectance; bt: bit transfer (different kind of reflectance methods)

# see how much the bands are correlated each other
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre) # $: link two variables

# PCA: Principal Component Analysis
# decrease the resolution: we have huge amount of pixels
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)

# library(RStoolbox) is now needed (to make PCA)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res) 
# rasterPCA: function to make the PCA

cl <- colorRampPalette(c('dark grey','grey','light grey'))(100) # let's plot the $map
plot(p224r63_2011_pca$map, col=cl)

# $model: how much variation is explained by each component
summary(p224r63_2011_pca$model) # summary: to see the information clearly
# PC1 99.83% of the whole variation

pairs(p224r63_2011) # let’s see how the different bands are correlated to each other

plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretc="Lin")) # 1=PC1; 2=PC2; 3=PC3 

# Repeat for the 1988 image
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

# faPAR: we can see the real power of the forest in keeping carbon
faPAR10 <- raster("faPAR10.tif") # import the image; faPAR10 because it is aggregate by a fact=10
levelplot(faPAR10)

pdf("copNDVI.pdf") # make a pdf file; quotes because we are saving the file outside R
levelplot(copNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

##############################
setwd("C:/lab/")
load("faPAR.RData")

# the original faPAR from Copernicus is 2GB

# let's see how much space is needed for an 8-bit set

library(raster)
library(rasterdiv)

writeRaster(copNDVI, "copNDVI.tif") # write an entire Raster object to a file
# 5.3MB

library(rasterVis)

# faPAR: levelplot this set
levelplot(faPar)

### Regression model between faPar and NDVI

erosion <- c(12, 14, 16, 24, 26, 40, 55, 67)
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")

# fit linear models
model1 <- lm(hm ~ erosion)
summary(model1)
# pvalue: how many times is a random situation
# p < 0.01 -> lower probability (1/100) that is random so the variables are related
abline(model1) # line described by 'a' and 'b'

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
model2 <- lm(faPAR10p ~ copNDVIp) # lm: to fit linear models
summary(model2)
## R-squared: 0.4279 → the point are not perfectly in line but they are far from being random, they are related to each other
## p-value: < 2.2e-16 --> the relation is confirmed
plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")

#####################################
#####################################
#####################################

# 10. R_code_radiance.r

## bit example
library(raster)  # format where there are pixels

# let's create a raster with two rows and two coloumns
toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13, 1.44, 1.55, 3.4) # put data into each pixel

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

# EBVs: Essential Biodiversity Variables

## upload image
setwd("C:/lab/") 
library(raster) #to see the image

# satellites images (Sentinel) are composed by many layers
snt <- brick("snt_r10.tif") #also named the picture
snt # see the characteristics of the image
# 30 bits image

# plot the image
plot(snt)
## complex system: passing from bare soil, to forest, opening parts related to agriculture, and forest

# B1 = blue
# B2 = green
# B3 = red
# B4 = NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin") # visible colours
plotRGB(snt, 4, 3, 2, stretch="Lin") # NIR in top of red: vegetation being coloured in red

pairs(snt) # to see if the axis are related to each other

# The analysis to start from multi system and moving to only one layer is the multivariate analysis
library(RStoolbox)
#PCA anlysis (coming out from the multivariate analysis)
sntPCA <- rasterPCA(snt)
sntPCA

# $: linking the model to the analysis we have done (PCA)
summary(sntPCA$model) # info about output of the model 
# 70% of the info (component 1)
plot(sntPCA$map) # PC1 has the highest amount of information

#plot RGB of the first three component
plotRGB(sntPCA$map, 1, 2, 3, stretch="Lin")

## calculate standard deviation on top of PCA
## set the moving window
## matrix of 5x5 rows and columns with value: 1 --> state that there are 
## 5 by 5 pixels moving throughout the image (the value 1 not impact the analysis)
window <- matrix(1, nrow = 5, ncol = 5) 
window # see the 5 by 5 window

# focal (raster function) --> calculates values for a neighbourhood of cells to calculate sd
sd_nst <- focal(sntPCA$map$PC1, w = window, fun =sd)   
# w = window (moving window); fun = 'function' = mean, mode, max, min, sd (standard deviation)
# we are going to make the analysis on top of the map and we use only the first principal component; 

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
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd) # make the focal function on top of the pc1 aggregated

# plot the two different set we made
par(mfrow=c(1,2)) # two images beside to the other in multiframe
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl) # cladonia set aggregated

# plot the calculation
par(mfrow=c(1,2)) 
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")

# plot(sd_clad_agg, col=cl)
plot(sd_clad, col=cl)


##########################################
##########################################
##########################################

# 12. R_code_snow.r

setwd("C:/lab/")

# install.packages("ncdf4")
install.packages("ncdf4") 
# ncdf4 package to read all the netCDF files. All Copernicus data use this extension

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

rlist <- list.files(pattern="snow") # make the list of all the files with common pattern (snow)
rlist

# lapply: apply the function over a list or vector; we want to apply the function raster (used to import single layer)
import <- lapply(rlist, raster) # import all the files

snow.multitemp <- stack(import) # stack: put altogether different layer (snow layers)
snow.multitemp
plot(snow.multitemp, col=cl)

# let's make a prediction
source("prediction.r")
# causes R to accept its input from the named file or URL or connection or expressions directly

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
# Crop an image
setwd("C:/lab/")

# install.packages("ncdf4")
install.packages("ncdf4") 
library(ncdf4)
library(raster)

snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # import image

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

ext <- c(0, 20, 35, 50)
zoom(snow, ext=ext)
# zoom the image

snowitaly <- crop(snow, ext)
# crop: allows to create a new image (the zoom made previously)
plot(snowitaly, col=cl)

## you can also use drawExtent to create a new image
# zoom(snow, ext=drawExtent())
# snowitaly <- crop(snow, drawExtent())
 
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
plotRGB(EN, r=1, g=7, b=13, stretch="lin") # each colour is associated with the number of the image in the stack

# difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(dif, col=cld)

# quantitative estimate: do the box-plot
boxplot(EN)
boxplot(EN, outline=F, horizontal=T) # remove the outline; move the box-plot horizontally
boxplot(EN, outline=F, horizontal=T, axes=T) # add an axes to make the box-plot easier to read

# plot!
plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red") # y=a+bx
# if there is a decrease in NO2,diff values should lay below the y=x line

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

# 14. R_code_interpolation.r

## Beach forest - Casentino

setwd("C:/lab/")
 
library(spatstat) # Spatial Point Pattern Analysis
 
# import data - since it is only a table, we will use the read.table function
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T) 
head(inp) # to see the first 6 rows of the dataset

attach(inp)
# plot(inp$X, inp$Y) if we did not attach the set
plot(X,Y)

summary(inp) # to know the minum and maxiumum value of x and y
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000)) # Planar Point Pattern

names(inp)
marks(inppp) <- Canopy.cov # label the sigle point

# estimate the canopy cover
canopy <- Smooth(inppp) # interpolate the data (smooth) of the inppp set where they have not been measured
# Warning: adapt the cross validation: chacking the new values with the original one
plot(canopy)
points(inppp, col="green") # add the point to the map

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

summary(inp.psam)
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

# Species Distribution Modelling

install.packages("sdm") # Species Distribution Modelling

library(sdm) 
library(raster) # predictors; make use of the raster set
library(rgdal)# species; manage coordinates and spactial data

# import the file (presence/absence file)
file <- system.file("external/species.shp", package="sdm") 
## external: there are species data and environmental variables data (inside the sdm package)

# use the graphical part of the file
species <- shapefile(file) 
plot(species) # presence/absence plot
species # information about species

# plot the occurance species
plot(species[species$Occurrence == 1,],col='blue',pch=16)

path <- system.file("external", package="sdm")

lst <- list.files(path=path,pattern='asc$',full.names = T) # ASII file: an extension
lst

# make a stack of the data
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

# make the model
d <- sdmData(train=species, predictors=preds) # explain which are the data

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm") # generalised linear model: several predictors altogether
p1 <- predict(m1, newdata=preds)

plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)
 
# we can make the final stack with all the predictors and variables
s1 <- stack(preds, p1)
plot(s1, col=cl)

#########################################
#########################################
#########################################

# R_code_final_project.r

# Maremma Region
# Vegetation analysis on the Maremma region and a zoom on Uccellina Park

library(raster) 
setwd("C:/lab/maremma/")

# Sentinel-2 bands
# b1 = Coastal Aerosol
# B2 = blue
# B3 = green
# B4 = red
# B5 = Vegetation Red Edge
# B6 = Vegetation Red Edge
# B7 = Vegetation Red Edge
# B8 = NIR
# B9 = Water Vapour
# B10 = SWIR-Cirrus
# B11 = SWIR
# B12 = SWIR
# Images from May, 2020

# import the files through the lapply function

rlist <- list.files(pattern="20200523")
rlist
import <- lapply(rlist, raster)
May <- stack(import)
plot(May) # plot the four bands that I picked

# Images from June, 2020
rlist <- list.files(pattern="20200622")
rlist
import1 <- lapply(rlist, raster)
June <- stack(import1)
plot(June)

# Image from July, 2020
rlist <- list.files(pattern="20200707")
rlist
import2 <- lapply(rlist, raster)
July <- stack(import2)
plot(July)

# RGB space

# Plot in the visible
# May
plotRGB(May, r=3, g=2, b=1, stretch="lin") 

# June
plotRGB(June, r=3, g=2, b=1, stretch="lin")

# July
plotRGB(July, r=3, g=2, b=1, stretch="lin")

par(mfrow=c(1,3))
# adjust the parameters so the axes colors are white
par(col.axis = "white", col.lab = "white")
# plot
plotRGB(May,3,2,1,scale= "20000", stretch = "lin", axes = TRUE, main = "23.05.2020")
# set bounding box to white as well; add the axes and the title
box(col = "white")
plotRGB(June,3,2,1,scale= "20000", stretch = "lin", axes =TRUE, main = "22.06.2020")
box(col = "white")
plotRGB(July,3,2,1,scale= "20000", stretch = "lin", axes =TRUE, main = "07.07.2020")
box(col = "white")

# NIR on top of the red component
par(mfrow=c(1,3))
plotRGB(May, 4,3,2, scale= "20000", stretch = "lin", axes = TRUE, main = "23.05.2020")
box(col = "white")
plotRGB(June, 4,3,2, scale= "20000", stretch = "lin", axes =TRUE, main = "22.06.2020")
box(col = "white")
plotRGB(July, 4,3,2, scale= "20000", stretch = "lin", axes =TRUE, main = "07.07.2020")
box(col = "white")

# DVI calculation: NIR - RED

# Names: to see the order of the bands
names(May)
names(June)
names(July)

# DVI
dviMay <- (May$T32TPN_20200523T100559_B08 - May$T32TPN_20200523T100559_B04) 
dviJune <- (June$T32TPN_20200622T100559_B8A_20m - June$T32TPN_20200622T100559_B04_20m) 
dviJuly <- (July$T32TPN_20200707T101031_B8A_20m - July$T32TPN_20200707T101031_B04_20m) 

# plot DVI
clDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(dviMay, col = clNDVI, main = "23.05.2020")
plot(dviJune, col = clNDVI, main = "22.06.2020")
plot(dviJuly, col = clNDVI, main = "07.07.2020")

# Difference between July and June DVI. 
# It was not possible with May because of the difference in resolution
difDVI <- dviJuly - dviMay
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difDVI, col=cld, main = "Difference DVI")

# NDVI
ndviMay <- (May$T32TPN_20200523T100559_B08 - May$T32TPN_20200523T100559_B04) / (May$T32TPN_20200523T100559_B08 + May$T32TPN_20200523T100559_B04)
ndviJune <- (June$T32TPN_20200622T100559_B8A_20m - June$T32TPN_20200622T100559_B04_20m) / (June$T32TPN_20200622T100559_B8A_20m + June$T32TPN_20200622T100559_B04_20m)
ndviJuly <- (July$T32TPN_20200707T101031_B8A_20m - July$T32TPN_20200707T101031_B04_20m) / (July$T32TPN_20200707T101031_B8A_20m + July$T32TPN_20200707T101031_B04_20m)

# Plot NDVI
clNDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(ndviMay, col = clNDVI, main = "23.05.2020")
plot(ndviJune, col = clNDVI, main = "22.06.2020")
plot(ndviJuly, col = clNDVI, main = "07.07.2020")

# Difference NDVI between June and July
diffNDVI <- ndviJuly - ndviJune
plot(diffNDVI, col=clNDVI, main = "Difference NDVI")

# Zoom on Parco dell'Uccellina

# Posso fare dirrettamete così perchè tutte le mie bande hanno una risoluzione di 10m

May2020 <- stack(import)
ext <- c(662000, 680000, 4710000, 4730000) # set the coordinates of the Park
Parcomay <- crop(May2020, ext)
plotRGB(Parcomay, r=3, g=2, b=1, stretch="lin")

------------------------------------------------------

May <- stack(import)
ext <- c(670000, 683000, 4670000, 4705000) # set the coordinates of the Park
argentario_may <- crop(May, ext)
plotRGB(argentario_may, r=3, g=2, b=1, stretch="lin")

July <- stack(import2)
ext <- c(670000, 683000, 4670000, 4705000) # set the coordinates of the Park
argentario_july <- crop(July, ext)
plotRGB(argentario_july, r=3, g=2, b=1, stretch="lin")
par(mfrow=c(1,2))
plotRGB(argentario_may, r=3, g=2, b=1, stretch="lin")
plotRGB(argentario_july, r=3, g=2, b=1, stretch="lin")

#DVI
dvi_may <- (argentario_may$T32TPN_20200523T100559_B08 - argentario_may$T32TPN_20200523T100559_B04)
dvi_july <- (argentario_july$T32TPN_20200707T101031_B8A_20m - argentario_july$T32TPN_20200707T101031_B04_20m) 
clNDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,2))
plot(dvi_may, col = clNDVI, main = "23.05.2020")
plot(dvi_july, col = clNDVI, main = "07.07.2020")
--------------------------------------------
June2020 <- stack(import1)
ext <- c(662000, 680000, 4710000, 4730000) # set the coordinates of the Park
Parcojune <- crop(June2020, ext)
plotRGB(Parcojune, r=3, g=2, b=1, stretch="lin")

July2020 <- stack(import2)
ext <- c(662000, 680000, 4710000, 4730000) # set the coordinates of the Park
Parcojuly <- crop(July2020, ext)
plotRGB(Parcojuly, r=3, g=2, b=1, stretch="lin")

par(mfrow=c(1,3))
plotRGB(Parcomay, r=3, g=2, b=1, stretch="lin")
plotRGB(Parcojune, r=3, g=2, b=1, stretch="lin")
plotRGB(Parcojuly, r=3, g=2, b=1, stretch="lin")


# NDVI on Parco dell'Uccellina
ndvi_park_may <- (Parcomay$T32TPN_20200523T100559_B08 - Parcomay$T32TPN_20200523T100559_B04) / (Parcomay$T32TPN_20200523T100559_B08 + Parcomay$T32TPN_20200523T100559_B04)
ndvi_park_june <- (Parcojune$T32TPN_20200622T100559_B8A_20m - Parcojune$T32TPN_20200622T100559_B04_20m) / (Parcojune$T32TPN_20200622T100559_B8A_20m + Parcojune$T32TPN_20200622T100559_B04_20m)
ndvi_park_july <- (Parcojuly$T32TPN_20200707T101031_B8A_20m - Parcojuly$T32TPN_20200707T101031_B04_20m) / (Parcojuly$T32TPN_20200707T101031_B8A_20m + Parcojuly$T32TPN_20200707T101031_B04_20m)

# plot NDVI
clDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(ndvi_park_may, col = clNDVI, main = "23.05.2020")
plot(ndvi_park_june, col = clNDVI, main = "22.06.2020")
plot(ndvi_park_july, col = clNDVI, main = "07.07.2020")

# DVI
dvi_park_may <- (Parcomay$T32TPN_20200523T100559_B08 - Parcomay$T32TPN_20200523T100559_B04) 
dvi_park_june <- (Parcojune$T32TPN_20200622T100559_B8A_20m - Parcojune$T32TPN_20200622T100559_B04_20m) 
dvi_park_july <- (Parcojuly$T32TPN_20200707T101031_B8A_20m - Parcojuly$T32TPN_20200707T101031_B04_20m) 

# plot DVI
clNDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(dvi_park_may, col = clNDVI, main = "23.05.2020")
plot(dvi_park_june, col = clNDVI, main = "22.06.2020")
plot(dvi_park_july, col = clNDVI, main = "07.07.2020")

# difference DVI
difDVI_park <- dvi_park_july - ndvi_park_june
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difDVI_park, col=cld, main = "Difference DVI")

# quantitative estimate
boxplot(difDVI_park)
boxplot(difDVI_park, outline=F, horizontal=T) # remove the outline; move the box-plot horizontally
boxplot(difDVI_park, outline=F, horizontal=T, axes=T) 

# plot!
plot(dvi_park_july, dvi_park_june)
abline(0,1,col="red") # y=a+bx





setwd("C:/lab/May8/")
band8 <- stack(import)
rlist <- list.files(pattern="20200622")
rlist
import <- lapply(rlist, raster)
T32TPN_20200523T100559_B08 <- resample(T32TPN_20200523T100559_B08,T32TPN_20200523T100559_B8A)

band8 <- brick("T32TPN_20200523T100559_B08")
band8A <- raster("T32TPN_20200523T100559_B8A")


setwd("C:/lab/June60/")

# Images from June, 2020, resolution = 60m
rlist <- list.files(pattern="20200622")
rlist
import <- lapply(rlist, raster)
june60 <- stack(import)
plotRGB(june60, r=1, g=3, b=2, stretch="lin") 
