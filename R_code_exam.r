# R_code_exam.r

# R first code
# 2. R spatial

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

###################################################
###################################################
###################################################

# 2. R spatial
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
