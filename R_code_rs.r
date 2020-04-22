# R code for remote sensing data analysis

#raster
install.packages("raster")
install.packages("RStoolbox")

setwd("C:/lab/")

library(raster)

# import imagines

p224r63_2011 <- brick("p224r63_2011_masked.grd")

plot(p224r63_2011)

cl <- colorRampPalette(c('black','grey','light grey'))(100)    # c= concatenate

# bands of Landstat
# B1: blue
# B2: green
# B3: red
# B4:NIR

# multiframe of different plots       
# c= series of characters

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

#plotRGB
# dev.off to close the previous part
plotRGB(p224r63_2011, r=3, g=2, b= 1, stretch="Lin")         # stretch= stretching the colour to see it better

# highlight the vegetation (3 bend at time)
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") 

# NIR on top of the G component of the RGB
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") 

plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin") 




