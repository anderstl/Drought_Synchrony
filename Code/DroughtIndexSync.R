####DROUGHT SYNCHRONY####

#load packages
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(ncf)
library(wsyn)

#set local working directory
#setwd("~/GitHub/Drought_Synchrony/")
setwd("~/GitRepos/Drought_Synchrony/")

#Load centroids of climate divisions
noaa.cents<-read.csv("Data/noaa.centroids.csv")

#load PHDI and PSDI values for all states
phdi<-read.csv("Data/AllStatesPHDI.csv")
pdsi<-read.csv("Data/AllStatesPDSI.csv")
enso<-read.csv("Data/ENSO_SST3.4_1981_2010meanremoved.csv")
nao<-read.csv("Data/NCAR_NAOdata.csv")
names(nao)[1]<-"Year"
enso<-enso[-((dim(enso)[1]-6):dim(enso)[1]),]

#create unique ID based on state code division
noaa.cents$Location<-paste(noaa.cents$StateCode,noaa.cents$Division,sep="")

#convert all data to long format
phdi.long<-phdi%>%
  tidyr::gather(key=Month,value=PHDI,-(ID))
pdsi.long<-pdsi%>%
  tidyr::gather(key=Month,value=PDSI,-(ID))
nao.long<-nao%>%
  tidyr::gather(key=Month,value=NAO,-(Year))%>%
  filter(Year%in%c(1895:2018))
enso.long<-enso%>%
  tidyr::gather(key=Month,value=ENSO,-(Year))%>%
  filter(Year%in%c(1895:2018))

#create columns for states, divisions, years and locations (combo of state and division)
phdi.long$StateCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),1,1)),as.numeric(substr(as.character(phdi.long$ID),1,2)))
phdi.long$DivisionCode<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),2,3)),as.numeric(substr(as.character(phdi.long$ID),3,4)))
phdi.long$Year<-ifelse(nchar(trunc(phdi.long$ID))==9,as.numeric(substr(as.character(phdi.long$ID),6,9)),as.numeric(substr(as.character(phdi.long$ID),7,10)))
phdi.long$Location<-paste(phdi.long$StateCode,phdi.long$DivisionCode,sep="")

pdsi.long$StateCode<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),1,1)),as.numeric(substr(as.character(pdsi.long$ID),1,2)))
pdsi.long$DivisionCode<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),2,3)),as.numeric(substr(as.character(pdsi.long$ID),3,4)))
pdsi.long$Year<-ifelse(nchar(trunc(pdsi.long$ID))==9,as.numeric(substr(as.character(pdsi.long$ID),6,9)),as.numeric(substr(as.character(pdsi.long$ID),7,10)))
pdsi.long$Location<-paste(pdsi.long$StateCode,pdsi.long$DivisionCode,sep="")
pdsi.long<-pdsi.long[pdsi.long$StateCode<50,]

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

nao.wide<-nao.long%>%
  group_by(Year)%>%
  dplyr::summarise(NAO=mean(NAO))
nao.mat<-matrix(nao.wide$NAO,nrow=dim(phdi.wide)[1],ncol=length(nao.wide$NAO),byrow=T)
             
enso.wide<-enso.long%>%
  group_by(Year)%>%
  dplyr::summarise(ENSO=as.numeric(mean(as.numeric(ENSO),na.rm=T)))
enso.mat<-matrix(enso.wide$ENSO,nrow=dim(phdi.wide)[1],ncol=length(enso.wide$ENSO),byrow=T)

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

#clean data
phdi.clean<-cleandat(as.matrix(phdi.wide[,-1]),clev=5,times=1895:2018)
pdsi.clean<-cleandat(as.matrix(pdsi.wide[,-1]),clev=5,times=1895:2018)
nao.clean<-cleandat(as.matrix(nao.mat),clev=5,times=1895:2018)
enso.clean<-cleandat(as.matrix(enso.mat),clev=5,times=1895:2018)

#Do spatial coherence between phdi and psdi with NAO and ENSO
bands<-rbind(c(2,4),c(4,50))

nao.phdi<-coh(dat1=phdi.clean$cdat,dat2=nao.clean$cdat,times=1895:2018,norm="powall",
               sigmethod="fast",nrand=1000,f0=1)
nao.pdsi<-coh(dat1=pdsi.clean$cdat,dat2=nao.clean$cdat,times=1895:2018,norm="powall",
               sigmethod="fast",nrand=1000,f0=1)
enso.phdi<-coh(dat1=phdi.clean$cdat,dat2=enso.clean$cdat,times=1895:2018,norm="powall",
    sigmethod="fast",nrand=1000,f0=1)
enso.pdsi<-coh(dat1=pdsi.clean$cdat,dat2=enso.clean$cdat,times=1895:2018,norm="powall",
               sigmethod="fast",nrand=1000,f0=1)

nao.phdi.test<-bandtest(nao.phdi,bands[1,])
nao.phdi.test$bandp
nao.phdi.test<-bandtest(nao.phdi,bands[2,])
nao.phdi.test$bandp

nao.pdsi.test<-bandtest(nao.pdsi,bands[1,])
nao.pdsi.test$bandp
nao.pdsi.test<-bandtest(nao.pdsi,bands[2,])
nao.pdsi.test$bandp

enso.phdi.test<-bandtest(enso.phdi,bands[1,])
enso.phdi.test$bandp
enso.phdi.test<-bandtest(enso.phdi,bands[2,])
enso.phdi.test$bandp

enso.pdsi.test<-bandtest(enso.pdsi,bands[1,])
enso.pdsi.test$bandp
enso.pdsi.test<-bandtest(enso.pdsi,bands[2,])
enso.pdsi.test$bandp
