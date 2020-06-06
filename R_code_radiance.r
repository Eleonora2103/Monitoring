# R_code_radiance.r

# bit example
library(raster)  # format where there are pixels

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13, 1.44, 1.55, 3.4) # put data into iìeach pixel

plot(toy)
text(toy, digits=2) # add the values of toy, using 2 decimal degrees

## we are going to force this data towards lower amount of bits
## Image of 2 bits how it looks like
## assign a code (number) from 0 to 3; with 4 values we are ranging from 0 to 3
toy2bits <- stretch(toy,minv=0,maxv=3) # stretch: function to change the range of values
storage.mode(toy2bits[]) = "integer" # use only integer numbers
plot(toy2bits)
text(toy2bits, digits=2)

# toy with a range of 4 bits (16 combinations)
toy4bits <- stretch(toy,minv=0,maxv=15) # range from 0 to 15
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)

# toy with a range of 4 bits (256 combinations)
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2)

# plot altogether
par(mfrow=c(1,4))
plot(toy)
text(toy, digits=2)

plot(toy2bits)
text(toy2bits, digits=2)

plot(toy4bits)
text(toy4bits, digits=2)

plot(toy8bits)
text(toy8bits, digits=2)

dev.off
library(rasterdiv) 
plot(copNDVI) # see the range of data (bits) of copNDVI
## values ranging from 0 to 255 --> 8 bits images (most used - sometimes 16 bits)


