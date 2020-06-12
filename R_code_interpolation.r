# R_code_interpolation.r
# Beach forest - Casentino

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

