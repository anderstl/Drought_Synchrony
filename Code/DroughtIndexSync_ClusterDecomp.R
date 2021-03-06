####DROUGHT SYNCHRONY####

rm(list=ls())

#load packages
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(ncf)
library(wsyn)
library(rgdal)
library(RColorBrewer)

#set local working directory
#setwd("C:/Users/Tom/Documents/GitRepos/Drought_Synchrony/")
setwd("~/GitHub/Drought_Synchrony/")

source("./Code/plotClusterMap_20181210.R")

#Load centroids of climate divisions
noaa.cents<-read.csv("Data/noaa.centroids.csv")

#load PHDI and PSDI values for all states
phdi<-read.csv("Data/AllStatesPHDI.csv", na.strings = "-99.99")
pdsi<-read.csv("Data/AllStatesPDSI.csv", na.strings = "-99.99")
stateshp<-readOGR("./Data/statesp020.shp")
stateshp<-stateshp[!stateshp$STATE %in% c("Alaska","U.S. Virgin Islands","Hawaii","Puerto Rico"),]

#create unique ID based on state code division
noaa.cents$Location<-paste(noaa.cents$StateCode,noaa.cents$Division,sep="")

#convert PHDI & PDSI to long format
phdi.long<-phdi%>%
  tidyr::gather(key=Month,value=PHDI,-(ID))

pdsi.long<-pdsi%>% #on hold until i figure out what the ID in the PDSI data represents
  tidyr::gather(key=Month,value=PDSI,-(ID))

#droughtmerge<-merge(phdi.long, pdsi.long, by="ID")

#create columns for states, divisions, years and locations (combo of state and division)
phdi.long$StateCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),1,1)),as.numeric(substr(as.character(phdi.long$ID),1,2)))
phdi.long$DivisionCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),2,3)),as.numeric(substr(as.character(phdi.long$ID),3,4)))
phdi.long$Year<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),6,9)),as.numeric(substr(as.character(phdi.long$ID),7,10)))
phdi.long$Location<-paste(phdi.long$StateCode,phdi.long$DivisionCode,sep="")

pdsi.long$StateCode<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),1,1)),as.numeric(substr(as.character(pdsi.long$ID),1,2)))
pdsi.long$DivisionCode<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),2,3)),as.numeric(substr(as.character(pdsi.long$ID),3,4)))
pdsi.long$Year<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),6,9)),as.numeric(substr(as.character(pdsi.long$ID),7,10)))
pdsi.long$Location<-paste(pdsi.long$StateCode,pdsi.long$DivisionCode,sep="")

#combine with centroid data
phdi.long<-merge(phdi.long,noaa.cents[,-1],by=c("Location","StateCode"))
pdsi.long<-merge(pdsi.long,noaa.cents[,-1],by=c("Location","StateCode"))

#convert back to wide format, using the annual average and dropping 2019 from the data
phdi.wide<-phdi.long%>%
  dplyr::filter(Year<2019)%>%
  group_by(Location,Year)%>%
  dplyr::summarise(PHDI=mean(PHDI))%>%
  spread(key = Year,value = PHDI)

pdsi.wide<-pdsi.long%>%
  dplyr::filter(Year<2019)%>%
  group_by(Location,Year)%>%
  dplyr::summarise(PDSI=mean(PDSI))%>%
  spread(key = Year,value = PDSI)

# All years: cross-correlations, wavelet mean field and clustering
phdi.clean<-cleandat(as.matrix(phdi.wide[,-1]),clev=5,times=1895:2018)
pdsi.clean<-cleandat(as.matrix(pdsi.wide[,-1]),clev=5,times=1895:2018)

phdi.clust<-clust(phdi.clean$cdat, times=1895:2018, coords=data.frame(lon=unique(phdi.long$lon),lat=unique(phdi.long$lat)), method="pearson")
phdi.clust<-addwmfs(phdi.clust)
phdi.clust<-addwpmfs(phdi.clust)

pdsi.clust<-clust(pdsi.clean$cdat, times=1895:2018, coords=data.frame(lon=unique(pdsi.long$lon),lat=unique(pdsi.long$lat)), method="pearson")
pdsi.clust<-addwmfs(pdsi.clust)
pdsi.clust<-addwpmfs(pdsi.clust)

## do cluster map plotting for phdi -- all years, all of the US
pdf("./Results/clustermap_phdi_CONUS_allyrs.pdf", width=6.5, height=4.5)
plotmap(phdi.clust)
dev.off()

pdf("./Results/modulewmfs_phdi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(phdi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,1.6))
}
dev.off()

pdf("./Results/modulewpmfs_phdi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(phdi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii), zlim=c(0,1))
}
dev.off()

print(phdi.clust$clusters)

## do cluster map plotting for pdsi -- all years, all of the US
pdf("./Results/clustermap_pdsi_CONUS_allyrs.pdf", width=6.5, height=4.5)
plotmap(pdsi.clust)
dev.off()

pdf("./Results/modulewmfs_pdsi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(pdsi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,1.7))
}
dev.off()

pdf("./Results/modulewpmfs_pdsi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(pdsi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii), zlim=c(0,1))
}
dev.off()

print(pdsi.clust$clusters)


################################################
## Now just for aman range in breeding season ##

aman.phdi<-filter(phdi.long,Location%in%c(232,233,234,235, #missouri
                                          31,32,33,34,35,37, #arkansas
                                          343,345,346,348,349) & Year<2019) #oklahoma

amanshp<-stateshp[stateshp$STATE %in% c("Missouri","Arkansas","Oklahoma"),]

#filter PHDI for the breeding season only(Sept and Oct), and take the average by division
aman.breed.phdi<-aman.phdi%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(phdi=mean(PHDI))

#convert to wide
aman.breed.phdi.wide<-aman.breed.phdi%>%
  spread(Year,phdi)

#clean and standardize data
aman.breed.phdi.clean<-cleandat(as.matrix(aman.breed.phdi.wide[,-1]),clev=5,times=1895:2018)

aman.breed.phdi.clust<-clust(aman.breed.phdi.clean$cdat, times=1895:2018, 
                             coords=data.frame(lon=unique(aman.phdi$lon),lat=unique(aman.phdi$lat)), method="pearson")
aman.breed.phdi.clust<-addwmfs(aman.breed.phdi.clust)
aman.breed.phdi.clust<-addwpmfs(aman.breed.phdi.clust)

## do cluster map plotting for phdi -- AMAN, breeding season
pdf("./Results/clustermap_phdi_aman_breed.pdf", width=6.5, height=4.5)
plotClusterMap(aman.breed.phdi.clust, basemap=amanshp)
dev.off()

pdf("./Results/modulewmfs_phdi_aman_breed.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(aman.breed.phdi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,2.5))
}
dev.off()

pdf("./Results/modulewpmfs_phdi_aman_breed.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(aman.breed.phdi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii),zlim=c(0,1))
}
dev.off()

print(aman.breed.phdi.clust$clusters)

######################################################
## Try out k-means clustering as an alternative 
######################################################

aman.phdi.kmeans<-kmeans(aman.breed.phdi.clean$cdat, 2)
print(aman.breed.phdi.clust$clusters[[2]])
print(aman.phdi.kmeans$cluster) #clusters are virtually identical

phdi.kmeans<-kmeans(phdi.wide,4)
print(phdi.clust$clusters[[4]])
print(phdi.kmeans$cluster)


pdf("./Results/kmeansmap_phdi_CONUS_alltime.pdf", width=6.5, height=4.5)
plot(phdi.long$lon, phdi.long$lat, pch=16, col=phdi.kmeans$cluster)
dev.off()
