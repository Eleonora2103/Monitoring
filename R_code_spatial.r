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






