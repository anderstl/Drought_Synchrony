#set working directory
setwd("C:/Users/Tom/Documents/Projects/ClimateSynchrony/")
rm(list=ls())

#load packages
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(ncf)
library(wsyn)

#load drought files
filenames <- list.files("DroughtData", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)

#merge files together
drought.dat<-do.call("rbind",ldf)

#drop "County" from county names
drought.dat$County<-as.factor(gsub(' County', '' , drought.dat$County))

#load centroid data
centroids<-read.csv("County_Centroids_Ambystoma_annulatum_range_Jan-2019.csv")

#rename counties and states so they match between data frames
centroids$state<-revalue(centroids$state, c("Arkansas"="AR", "Missouri"="MO","Oklahoma"="OK"))
levels(centroids$county)[levels(centroids$county)=="Dekalb"] <- "DeKalb"
levels(centroids$county)[levels(centroids$county)=="St. Francis"] <- "Saint Francis"
levels(centroids$county)[levels(centroids$county)=="St Francois"] <- "Saint Francois"
levels(centroids$county)[levels(centroids$county)=="St Charles"] <- "Saint Charles"
levels(centroids$county)[levels(centroids$county)=="St Clair"] <- "Saint Clair"
levels(centroids$county)[levels(centroids$county)=="Ste Genevieve"] <- "Sainte Genevieve"
levels(centroids$county)[levels(centroids$county)=="St Louis City"] <- "Saint Louis City"
levels(centroids$county)[levels(centroids$county)=="St Louis"] <- "Saint Louis"

#combine data sets
all.dat<-merge(drought.dat,centroids,by.x=c("State","County"),by.y=c("state","county"))

#fiter to counties where AMAN has known presences
aman.dat<-all.dat%>%filter(Aa.Present=="Present")

#Make function to change dates
extrdate<-function(dat,name){
  dat$Year<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%Y"))
  dat$Month<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%m"))
  dat$Day<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%d"))
  dat$Jdate<-as.numeric(format(as.Date(as.character(dat[,name]),"%m/%d/%Y"),"%j"))
  return(dat)
}

#extract month, day and year info from start date of each week
aman.dat<-extrdate(aman.dat,"ValidStart")

ggplot(aman.dat,aes(ValidStart,D0,group=County))+
  geom_point(show.legend=F)+
  geom_line(aes(color=County),show.legend=F)+
  facet_grid(~State)

aman.breeding<-aman.dat%>%
  filter(Month%in%9:10)%>%
  group_by(State,County,Year)%>%
  dplyr::summarise(D0mean=mean(D0),D1mean=mean(D1),D2mean=mean(D2))
ggplot(aman.breeding,aes(Year,D2mean,group=County))+
  geom_point(show.legend=F)+
  geom_line(aes(color=County),show.legend=F)

aman.meta<-aman.dat%>%
  filter(Month%in%4:6)%>%
  group_by(State,County,Year)%>%
  dplyr::summarise(D0mean=mean(D0),D1mean=mean(D1),D2mean=mean(D2))
ggplot(aman.meta,aes(Year,D2mean,group=County))+
  geom_point(show.legend=F)+
  geom_line(aes(color=County),show.legend=F)

#breeding synchrony
breed.d0wide<-spread(aman.breeding[,c("State","County","Year","D0mean")],Year,D0mean)
breed.d0wideClean<-cleandat(as.matrix(breed.d0wide[,3:dim(breed.d0wide)[2]]),times=2000:2018,clev=3)
breed.sync<-synmat(breed.d0wideClean$cdat,times = 2000:2018,method="spearman")
aman.dist<-gcdist(x=unique(aman.dat$lon),y=unique(aman.dat$lat))
breed.ncf<-Sncf(x=unique(aman.dat$lon),y=unique(aman.dat$lat),z=breed.d0wideClean$cdat,latlon = T)
plot(breed.ncf)
breed.clust<-clust(breed.d0wideClean$cdat,times=2000:2018,coords=data.frame(lon=unique(aman.dat$lon),lat=unique(aman.dat$lat)),method="spearman")
plotmap(breed.clust)
plotClusterMap(breed.clust)
# library(maps)
# map("state",c('missouri','arkansas','oklahoma'))
# points(breed.clust$coords$lon,breed.clust$coords$lat)


#metamorphosis synchrony
meta.d0wide<-spread(aman.meta[,c("State","County","Year","D0mean")],Year,D0mean)
meta.d0wideClean<-cleandat(as.matrix(meta.d0wide[,3:dim(meta.d0wide)[2]]),times=2000:2018,clev=3)
meta.ncf<-Sncf(x=unique(aman.dat$lon),y=unique(aman.dat$lat),z=meta.d0wideClean$cdat,latlon = T)
plot(meta.ncf)
meta.clust<-clust(meta.d0wideClean$cdat,times=2000:2018,coords=data.frame(lon=unique(aman.dat$lon),lat=unique(aman.dat$lat)),method="spearman")
plotmap(meta.clust)


