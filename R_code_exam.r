# R_code_exam.r

# 1. R first code
# 2. R code multipanel
# 3. R spatial
# 4. R code for multivariate analysis
# 5. R_code_rs.r

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
# R code for spatial view points

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
     
##############################################
##############################################
##############################################

# R_code_multivariate.r
# R code for multivariate analysis

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
     
# 5. R_code_rs.r
     
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
