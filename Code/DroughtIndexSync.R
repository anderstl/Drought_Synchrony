####DROUGHT SYNCHRONY####

#load packages
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(ncf)
library(wsyn)

#set local working directory
setwd("C:/Users/Tom/Documents/GitRepos/Drought_Synchrony/")

#Load centroids of climate divisions
noaa.cents<-read.csv("Data/noaa.centroids.csv")

#load PHDI and PSDI values for all states
phdi<-read.csv("Data/AllStatesPHDI.csv")
phsi<-read.csv("Data/AllStatesPHDI.csv")

#create unique ID based on state code division
noaa.cents$Location<-paste(noaa.cents$StateCode,noaa.cents$Division,sep="")

#convert PHDI to long format
phdi.long<-phdi%>%
  tidyr::gather(key=Month,value=PHDI,-(ID))

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
phdi.ncf<-Sncf(x=unique(phdi.long$lon),y=unique(phdi.long$lat),z=phdi.clean$cdat,latlon=T)
phdi.wmf<-wmf(phdi.clean$cdat,times=1895:2018)
phdi.clust<-clust(phdi.clean$cdat,times=1895:2018,coords=data.frame(lon=unique(phdi.long$lon),lat=unique(phdi.long$lat)),method="spearman")

#first half (pre-1956: cross-correlations, wavelet mean field and clustering
phdiPre56.clean<-cleandat(as.matrix(phdi.wide[,2:63]),clev=5,times=1895:1956)
phdiPre56.clust<-clust(phdiPre56.clean$cdat,times=1895:1956,coords=data.frame(lon=unique(phdi.long$lon),lat=unique(phdi.long$lat)),method="spearman")
phdiPre56.ncf<-Sncf(x=unique(phdi.long$lon),y=unique(phdi.long$lat),z=phdiPre56.clean$cdat,latlon=T)
phdiPre56.wmf<-wmf(phdiPre56.clean$cdat,times=1895:1956)

#second half (post-1956): cross-correlations, wavelet mean field and clustering
phdiPost56.clean<-cleandat(as.matrix(phdi.wide[,64:125]),clev=5,times=1957:2018)
phdiPost56.wmf<-wmf(phdiPost56.clean$cdat,times=1957:2018)
phdiPost56.ncf<-Sncf(x=unique(phdi.long$lon),y=unique(phdi.long$lat),z=phdiPost56.clean$cdat,latlon=T)
phdiPost56.clust<-clust(phdiPost56.clean$cdat,times=1957:2018,coords=data.frame(lon=unique(phdi.long$lon),lat=unique(phdi.long$lat)),method="spearman")

#make plots
#pdf(file="Results/PHDI.pdf",height=9,width=9)
par(mfrow=c(3,3))
plot(phdi.ncf)
title("1895-2018")
plot(phdiPre56.ncf)
title("1895-1956")
plot(phdiPost56.ncf)
title("1957-2018")
plotmap(phdi.clust)
plotmap(phdiPre56.clust)
plotmap(phdiPost56.clust)
plotmag(phdi.wmf)
plotmag(phdiPre56.wmf)
plotmag(phdiPost56.wmf)
#dev.off()


####Drought synchchrony for the Ringed Salamander (Amybstoma annulatum)####

aman.phdi<-filter(phdi.long,Location%in%c(232,233,234,235, #missouri
                          31,32,33,34,35,37, #arkansas
                          343,345,346,348,349) & Year<2019) #oklahoma

#filter PHDI for the breeding season only(Sept and Oct), and take the average by division
aman.breed.phdi<-aman.phdi%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(phdi=mean(PHDI))

#plot raw time series by Location
ggplot(aman.breed.phdi,aes(Year,phdi,color=as.factor(Location),group=Location))+geom_point()+geom_line()

#convert to wide
aman.breed.phdi.wide<-aman.breed.phdi%>%
  spread(Year,phdi)

#clean and standardize data
aman.breed.phdi.clean<-cleandat(as.matrix(aman.breed.phdi.wide[,-1]),clev=5,times=1895:2018)

#run wmf, clustering and cross-correlation tests
aman.phdi.wmf<-wmf(aman.breed.phdi.clean$cdat,times=1895:2018)
aman.phdi.clust<-clust(aman.breed.phdi.clean$cdat,times=1895:2018,coords=data.frame(lon=unique(aman.phdi$lon),lat=unique(aman.phdi$lat)),method="spearman")
phdi.ncf<-Sncf(x=unique(aman.phdi$lon),y=unique(aman.phdi$lat),z=aman.breed.phdi.clean$cdat,latlon=T)

#plots
par(mfrow=c(1,3))
plot(phdi.ncf)
plotmap(aman.phdi.clust)
plotmag(aman.phdi.wmf)
