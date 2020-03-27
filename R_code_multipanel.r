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
