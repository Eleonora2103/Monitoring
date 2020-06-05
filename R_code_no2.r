# R_code_no2.r

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
plotRGB(EN, r=1, g=7, b=13, stretch="lin")

# difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(dif, col=cld)

# quantitative estimate: do the box-plot
boxplot(EN)
boxplot(EN, outline=F, horizontal=T)
boxplot(EN,outline=F, horizontal=T, axes=T)

# plot!
plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red") # y=a+bx

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

