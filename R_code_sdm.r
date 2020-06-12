# R_code_sdm.r 

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
