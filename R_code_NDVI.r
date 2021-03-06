# Maremma NDVI 
ndviMay <- (May$T32TPN_20200523T100559_B08 - May$T32TPN_20200523T100559_B04) / (May$T32TPN_20200523T100559_B08 + May$T32TPN_20200523T100559_B04)
ndviJune <- (June$T32TPN_20200622T100559_B8A_20m - June$T32TPN_20200622T100559_B04_20m) / (June$T32TPN_20200622T100559_B8A_20m + June$T32TPN_20200622T100559_B04_20m)
ndviJuly <- (July$T32TPN_20200707T101031_B8A_20m - July$T32TPN_20200707T101031_B04_20m) / (July$T32TPN_20200707T101031_B8A_20m + July$T32TPN_20200707T101031_B04_20m)

# Plot NDVI
clNDVI <- colorRampPalette(c('blue','white','red'))(100)
par(mfrow=c(1,3))
plot(ndviMay, col = clNDVI, main = "23.05.2020")
plot(ndviJune, col = clNDVI, main = "22.06.2020")
plot(ndviJuly, col = clNDVI, main = "07.07.2020")

# Difference NDVI between June and July
diffNDVI <- ndviJuly - ndviJune
plot(diffNDVI, col=clNDVI, main = "Difference NDVI")

# NDVI on Parco dell'Uccellina
ndvi_park_may <- (Parcomay$T32TPN_20200523T100559_B08 - Parcomay$T32TPN_20200523T100559_B04) / (Parcomay$T32TPN_20200523T100559_B08 + Parcomay$T32TPN_20200523T100559_B04)
ndvi_park_june <- (Parcojune$T32TPN_20200622T100559_B8A_20m - Parcojune$T32TPN_20200622T100559_B04_20m) / (Parcojune$T32TPN_20200622T100559_B8A_20m + Parcojune$T32TPN_20200622T100559_B04_20m)
ndvi_park_july <- (Parcojuly$T32TPN_20200707T101031_B8A_20m - Parcojuly$T32TPN_20200707T101031_B04_20m) / (Parcojuly$T32TPN_20200707T101031_B8A_20m + Parcojuly$T32TPN_20200707T101031_B04_20m)

# plot NDVI
clDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(ndvi_park_may, col = clNDVI, main = "23.05.2020")
plot(ndvi_park_june, col = clNDVI, main = "22.06.2020")
plot(ndvi_park_july, col = clNDVI, main = "07.07.2020")

difNDVI_park <- ndvi_park_july$T32TPN_20200707T101031_B8A_20m - ndvi_park_june$T32TPN_20200622T100559_B8A_20m

clNDVI <- colorRampPalette(c("darkblue","yellow","red","black"))(200)
plot(difNDVI_park, col=clNDVI, main = "Difference NDVI")

# DVI
dviMay <- (May$T32TPN_20200523T100559_B08 - May$T32TPN_20200523T100559_B04) 
dviJune <- (June$T32TPN_20200622T100559_B8A_20m - June$T32TPN_20200622T100559_B04_20m) 
dviJuly <- (July$T32TPN_20200707T101031_B8A_20m - July$T32TPN_20200707T101031_B04_20m) 

# plot DVI
clDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(dviMay, col = clNDVI, main = "23.05.2020")
plot(dviJune, col = clNDVI, main = "22.06.2020")
plot(dviJuly, col = clNDVI, main = "07.07.2020")

# DVI
dvi_park_may <- (Parcomay$T32TPN_20200523T100559_B08 - Parcomay$T32TPN_20200523T100559_B04) 
dvi_park_june <- (Parcojune$T32TPN_20200622T100559_B8A_20m - Parcojune$T32TPN_20200622T100559_B04_20m) 
dvi_park_july <- (Parcojuly$T32TPN_20200707T101031_B8A_20m - Parcojuly$T32TPN_20200707T101031_B04_20m) 

# plot DVI
clNDVI = colorRampPalette(c("darkblue","yellow","red","black"))(200)
par(mfrow=c(1,3))
plot(dvi_park_may, col = clNDVI, main = "23.05.2020")
plot(dvi_park_june, col = clNDVI, main = "22.06.2020")
plot(dvi_park_july, col = clNDVI, main = "07.07.2020")
