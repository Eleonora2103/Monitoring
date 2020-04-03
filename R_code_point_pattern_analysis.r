# point pattern analysis: density map

install.packages("spatstat")

library(spatstat)

attach(covid)
head(covid)

covids <- ppp(lot, lan, c(-180,180), c(-90,90))     # panel point pattern

d <- density(covids)

plot(d)
points(covids)




