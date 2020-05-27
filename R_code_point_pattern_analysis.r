# point pattern analysis: density map

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






