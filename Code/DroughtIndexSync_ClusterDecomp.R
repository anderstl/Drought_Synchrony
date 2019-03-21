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

# pdsi.long<-pdsi%>% #on hold until i figure out what the ID in the PDSI data represents
#   tidyr::gather(key=Month,value=PDSI,-(ID))
# 
# droughtmerge<-merge(phdi.long, pdsi.long, by="ID")

#create columns for states, divisions, years and locations (combo of state and division)
phdi.long$StateCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),1,1)),as.numeric(substr(as.character(phdi.long$ID),1,2)))
phdi.long$DivisionCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),2,3)),as.numeric(substr(as.character(phdi.long$ID),3,4)))
phdi.long$Year<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),6,9)),as.numeric(substr(as.character(phdi.long$ID),7,10)))
phdi.long$Location<-paste(phdi.long$StateCode,phdi.long$DivisionCode,sep="")

#combine with centroid data
phdi.long<-merge(phdi.long,noaa.cents[,-1],by=c("Location","StateCode"))

#convert back to wide format, using the annual average and dropping 2019 from the data
phdi.wide<-phdi.long%>%
  dplyr::filter(Year<2019)%>%
  group_by(Location,Year)%>%
  dplyr::summarise(PHDI=mean(PHDI))%>%
  spread(key = Year,value = PHDI)

# All years: cross-correlations, wavelet mean field and clustering
phdi.clean<-cleandat(as.matrix(phdi.wide[,-1]),clev=5,times=1895:2018)

phdi.clust<-clust(phdi.clean$cdat, times=1895:2018, coords=data.frame(lon=unique(phdi.long$lon),lat=unique(phdi.long$lat)), method="pearson")
phdi.clust<-addwmfs(phdi.clust)
phdi.clust<-addwpmfs(phdi.clust)

## do cluster map plotting for phdi -- all years, all of the US
pdf("./Results/clustermap_phdi_CONUS_allyrs.pdf", width=6.5, height=4.5)
plotmap(phdi.clust)
dev.off()

pdf("./Results/modulewmfs_phdi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(phdi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii))
}
dev.off()

pdf("./Results/modulewpmfs_phdi_CONUS_allyrs.pdf", width=6.5, height=9)
par(mfrow=c(4,1))
for(ii in 1:4){
  plotmag(phdi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii))
}
dev.off()

print(phdi.clust$clusters)
