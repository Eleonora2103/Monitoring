install.packages("sp")

library(sp)

data(meuse)

head(meuse)

attach(meuse)

plot(zinc,copper)
plot(zinc,copper,col="red")
plot(zinc,copper,col="red",pch=19)
plot(zinc,copper,col="red",pch=19,cex=2)
