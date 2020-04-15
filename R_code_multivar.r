# R code for multivariate analysis

setwd("C:/lab/")

# vegan = vegetation analysis

library(vegan)

biomes <- read.table("biomes.csv", header=T, sep=",")   # csv is an extension

head(biomes)

# DEtrended CORrespondence ANAlysis 

multivar <- decorana(biomes)
plot(multivar)
multivar

plot(multivar)

# biomes types
biomes_types <- read.table("biomes_types.csv", header=T, sep=",") 

head(biomes_types)

attach(biomes_types)
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3)            # ellipse connecting all the plots of same biomes

# also col=c("green","blue","red","black")

ordispider(multivar, type, col=1:4, label=T)

