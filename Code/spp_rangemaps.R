## Explore range maps in CAUDATA shapefile

rm(list=ls())

library(rgdal)
library(sp)
library(raster)
library(rgeos)


setwd("~/GitHub/Drought_Synchrony")

allranges<-readOGR("./Data/CaudataShapeFiles/CAUDATA.shp")
states<-readOGR("./Data/statesp020.shp")
conus<-states[!states$STATE %in% c("Alaska","Puerto Rico", "U.S. Virgin Islands", "Hawaii"),]

proj4string(allranges)
proj4string(states)

conus<-spTransform(conus,CRS(proj4string(allranges)))

# let's extract species we're interested in

A_opacum<-allranges[allranges$binomial == "Ambystoma opacum",]
A_maculatum<-allranges[allranges$binomial == "Ambystoma maculatum",]
A_talpoideum<-allranges[allranges$binomial == "Ambystoma talpoideum",]
A_texanum<-allranges[allranges$binomial == "Ambystoma texanum",]
A_cingulatum<-allranges[allranges$binomial == "Ambystoma cingulatum",]
A_bishopi<-allranges[allranges$binomial == "Ambystoma bishopi",]
A_annulatum<-allranges[allranges$binomial == "Ambystoma annulatum",]
C_alleganiensis<-allranges[allranges$binomial == "Cryptobranchus alleganiensis",]

pdf("./Results/rangemaps_rough.pdf", width=11, height=8.5)

par(mar=c(1,1,2,1), mfrow=c(2,4))

plot(A_opacum, col="grey", main="A. opacum")
lines(conus)

plot(A_maculatum, col="grey", main="A. maculatum")
lines(conus)

plot(A_talpoideum, col="grey", main="A. talpoideum")
lines(conus)

plot(A_texanum, col="grey", main="A. texanum")
lines(conus)

plot(A_cingulatum, col="grey", main="A. cingulatum")
lines(conus)

plot(A_bishopi, col="grey", main="A. bishopi")
lines(conus)

plot(A_annulatum, col="grey", main="A. annulatum")
lines(conus)

plot(C_alleganiensis, col="grey", main="C. alleganiensis")
lines(conus)

dev.off()


# writeOGR(A_opacum,"./Data/CaudataShapeFiles/","a_opacum", driver="ESRI Shapefile")
# writeOGR(A_maculatum,"./Data/CaudataShapeFiles/","a_maculatum", driver="ESRI Shapefile")
# writeOGR(A_talpoideum,"./Data/CaudataShapeFiles/","a_talpoideum", driver="ESRI Shapefile")
# writeOGR(A_texanum,"./Data/CaudataShapeFiles/","a_texanum", driver="ESRI Shapefile")
# writeOGR(A_cingulatum,"./Data/CaudataShapeFiles/","a_cingulatum", driver="ESRI Shapefile")
# writeOGR(A_bishopi,"./Data/CaudataShapeFiles/","a_bishopi", driver="ESRI Shapefile")
# writeOGR(A_annulatum,"./Data/CaudataShapeFiles/","a_annulatum", driver="ESRI Shapefile")
# writeOGR(C_alleganiensis,"./Data/CaudataShapeFiles/","c_alleganiensis", driver="ESRI Shapefile")