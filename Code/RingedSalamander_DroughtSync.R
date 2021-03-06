

####Drought synchrony for the Ringed Salamander (Amybstoma annulatum)####
source(Code/DroughtIndexSync.R)

aman.phdi<-filter(phdi.long,Location%in%c(232,233,234,235, #missouri
                                          31,32,33,34,35,37, #arkansas
                                          343,345,346,348,349) & Year<2019) #oklahoma
aman.pdsi<-filter(pdsi.long,Location%in%c(232,233,234,235, #missouri
                                          31,32,33,34,35,37, #arkansas
                                          343,345,346,348,349) & Year<2019) #oklahoma

#filter PHDI and PDSI for the breeding season only(Sept and Oct), and take the average by division
aman.breed.phdi<-aman.phdi%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(phdi=mean(PHDI))

aman.breed.pdsi<-aman.pdsi%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(pdsi=mean(PDSI))

aman.meta.phdi<-aman.phdi%>%
  filter(Month%in%c(c("Apr","May")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(phdi=mean(PHDI))

aman.meta.pdsi<-aman.pdsi%>%
  filter(Month%in%c(c("Apr","May")))%>%
  group_by(Location,Year)%>%
  dplyr::summarise(pdsi=mean(PDSI))

#plot raw time series by Location
ggplot(aman.breed.phdi,aes(Year,phdi,color=as.factor(Location),group=Location))+
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  labs(x="Year",y="Breeding PHDI")

ggplot(aman.breed.pdsi,aes(Year,pdsi,color=as.factor(Location),group=Location))+
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  labs(x="Year",y="Breeding PDSI")

ggplot(aman.meta.phdi,aes(Year,phdi,color=as.factor(Location),group=Location))+
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  labs(x="Year",y="Metamorphosis PHDI")

ggplot(aman.meta.pdsi,aes(Year,pdsi,color=as.factor(Location),group=Location))+
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  labs(x="Year",y="Metamorphosis PDSI")

##convert to wide format
#breeding data
aman.breed.phdi.wide<-aman.breed.phdi%>%
  spread(Year,phdi)
aman.breed.pdsi.wide<-aman.breed.pdsi%>%
  spread(Year,pdsi)
#metamorphosis data
aman.meta.phdi.wide<-aman.meta.phdi%>%
  spread(Year,phdi)
aman.meta.pdsi.wide<-aman.meta.pdsi%>%
  spread(Year,pdsi)

#clean and standardize data
aman.breed.phdi.clean<-cleandat(as.matrix(aman.breed.phdi.wide[,-1]),clev=5,times=1895:2018)
aman.breed.pdsi.clean<-cleandat(as.matrix(aman.breed.pdsi.wide[,-1]),clev=5,times=1895:2018)
aman.meta.phdi.clean<-cleandat(as.matrix(aman.meta.phdi.wide[,-1]),clev=5,times=1895:2018)
aman.meta.pdsi.clean<-cleandat(as.matrix(aman.meta.pdsi.wide[,-1]),clev=5,times=1895:2018)

#run cluster analysis
#phdi-breeding
aman.breed.phdi.clust<-clust(aman.breed.phdi.clean$cdat, times=1895:2018, 
                             coords=data.frame(lon=unique(aman.phdi$lon),lat=unique(aman.phdi$lat)), method="pearson")
aman.breed.phdi.clust<-addwmfs(aman.breed.phdi.clust)
aman.breed.phdi.clust<-addwpmfs(aman.breed.phdi.clust)

#pdsi-breeding
aman.breed.pdsi.clust<-clust(aman.breed.pdsi.clean$cdat, times=1895:2018, 
                             coords=data.frame(lon=unique(aman.pdsi$lon),lat=unique(aman.pdsi$lat)), method="pearson")
aman.breed.pdsi.clust<-addwmfs(aman.breed.pdsi.clust)
aman.breed.pdsi.clust<-addwpmfs(aman.breed.pdsi.clust)

#phdi-metamorphosis
aman.meta.phdi.clust<-clust(aman.meta.phdi.clean$cdat, times=1895:2018, 
                             coords=data.frame(lon=unique(aman.phdi$lon),lat=unique(aman.phdi$lat)), method="pearson")
aman.meta.phdi.clust<-addwmfs(aman.meta.phdi.clust)
aman.meta.phdi.clust<-addwpmfs(aman.meta.phdi.clust)

#pdsi-metamorphosis
aman.meta.pdsi.clust<-clust(aman.meta.pdsi.clean$cdat, times=1895:2018, 
                             coords=data.frame(lon=unique(aman.pdsi$lon),lat=unique(aman.pdsi$lat)), method="pearson")
aman.meta.pdsi.clust<-addwmfs(aman.meta.pdsi.clust)
aman.meta.pdsi.clust<-addwpmfs(aman.meta.pdsi.clust)

## do cluster map plotting for phdi -- all years, all of the US
pdf("./Results/clustermap_aman.pdf", width=6.5, height=4.5)
par(mfrow=c(2,2))
plotClusterMap(aman.breed.phdi.clust, basemap=amanshp)
title("PHDI-AMAN breeding")
plotClusterMap(aman.breed.pdsi.clust, basemap=amanshp)
title("PDSI-AMAN breeding")
plotClusterMap(aman.meta.phdi.clust, basemap=amanshp)
title("PHDI-AMAN metamorphosis")
plotClusterMap(aman.meta.pdsi.clust, basemap=amanshp)
title("PDSI-AMAN metamorphosis")
dev.off()

pdf("./Results/modulewmfs_aman.pdf", width=6.5, height=9)
par(mfrow=c(2,1))
for(ii in 1:2){
  plotmag(aman.breed.phdi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,2.5))
}
par(mfrow=c(2,1))
for(ii in 1:2){
  plotmag(aman.breed.pdsi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,3))
}
par(mfrow=c(2,1))
for(ii in 1:2){
  plotmag(aman.meta.phdi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,2.5))
}
par(mfrow=c(2,1))
for(ii in 1:2){
  plotmag(aman.meta.pdsi.clust$wmfs[[ii]][[ii]],title=paste0("wmf module ",ii),zlim=c(0,2.5))
}
dev.off()

pdf("./Results/modulewpmfs__aman.pdf", width=6.5, height=9)
par(mfrow=c(2,1))
for(ii in 1:2){
  plotmag(aman.breed.phdi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii),zlim=c(0,1))
}
for(ii in 1:2){
  plotmag(aman.breed.pdsi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii),zlim=c(0,1))
}
for(ii in 1:2){
  plotmag(aman.meta.phdi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii),zlim=c(0,1))
}
for(ii in 1:2){
  plotmag(aman.meta.pdsi.clust$wpmfs[[ii]][[ii]],title=paste0("wpmf module ",ii),zlim=c(0,1))
}
dev.off()

print(aman.breed.phdi.clust$clusters)

#compare synchrony to climate indices
enso.meta<-enso.long%>%
  filter(Month%in%c(c("Apr","May")))%>%
  group_by(Year)%>%
  dplyr::summarise(ENSO=as.numeric(mean(as.numeric(ENSO),na.rm=T)))
enso.meta.mat<-matrix(enso.meta$ENSO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)
pdo.meta<-pdo.long%>%
  filter(Month%in%c(c("Apr","May")))%>%
  group_by(Year)%>%
  dplyr::summarise(PDO=mean(PDO))
pdo.meta.mat<-matrix(pdo.meta$PDO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)
nao.meta<-nao.long%>%
  filter(Month%in%c(c("Apr","May")))%>%
  group_by(Year)%>%
  dplyr::summarise(NAO=mean(NAO))
nao.meta.mat<-matrix(nao.meta$NAO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)

enso.breed<-enso.long%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Year)%>%
  dplyr::summarise(ENSO=as.numeric(mean(as.numeric(ENSO),na.rm=T)))
enso.breed.mat<-matrix(enso.breed$ENSO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)
pdo.breed<-pdo.long%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Year)%>%
  dplyr::summarise(PDO=mean(PDO))
pdo.breed.mat<-matrix(pdo.breed$PDO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)
nao.breed<-nao.long%>%
  filter(Month%in%c(c("Sep","Oct")))%>%
  group_by(Year)%>%
  dplyr::summarise(NAO=mean(NAO))
nao.breed.mat<-matrix(nao.breed$NAO,nrow=length(unique(aman.breed.pdsi$Location)),ncol=length(unique(aman.breed.pdsi$Year)),byrow=T)

nao.breed.clean<-cleandat(nao.breed.mat[,-124],clev=5,times=1895:2018) 
enso.breed.clean<-cleandat(enso.breed.mat,clev=5,times=1895:2018)
pdo.breed.clean<-cleandat(pdo.breed.mat,clev=5,times=1895:2018)
nao.meta.clean<-cleandat(nao.meta.mat[,-124],clev=5,times=1895:2018) 
enso.meta.clean<-cleandat(enso.meta.mat,clev=5,times=1895:2018)
pdo.meta.clean<-cleandat(pdo.meta.mat,clev=5,times=1895:2018)

#test for coherence between indices during metaing season
nao.phdi.breed<-coh(dat1=aman.breed.phdi.clean$cdat,dat2=nao.clean$cdat,times=1895:2018,norm="powall",
              sigmethod="fast",nrand=1000,f0=1)
pdo.phdi.breed<-coh(dat1=aman.breed.phdi.clean$cdat,dat2=pdo.breed.clean$cdat,times=1895:2018,norm="powall",
              sigmethod="fast",nrand=1000,f0=1)
enso.phdi.breed<-coh(dat1=aman.breed.phdi.clean$cdat,dat2=enso.breed.clean$cdat,times=1895:2018,norm="powall",
               sigmethod="fast",nrand=1000,f0=1)
nao.pdsi.meta<-coh(dat1=aman.meta.pdsi.clean$cdat,dat2=nao.clean$cdat,times=1895:2018,norm="powall",
                    sigmethod="fast",nrand=1000,f0=1)
pdo.pdsi.meta<-coh(dat1=aman.meta.pdsi.clean$cdat,dat2=pdo.meta.clean$cdat,times=1895:2018,norm="powall",
                    sigmethod="fast",nrand=1000,f0=1)
enso.pdsi.meta<-coh(dat1=aman.meta.pdsi.clean$cdat,dat2=enso.meta.clean$cdat,times=1895:2018,norm="powall",
                     sigmethod="fast",nrand=1000,f0=1)


#pick some arbitrary timescale bands to start with for significance tests
bands<-rbind(c(2,4),c(4,50))

#test coherence between NAO and phdi during breeding
nao.phdi.test<-bandtest(nao.phdi.breed,bands[1,])
nao.phdi.test$bandp
nao.phdi.test<-bandtest(nao.phdi.breed,bands[2,])
nao.phdi.test$bandp

#test coherence between NAO and pdsi during breeding
nao.pdsi.test<-bandtest(nao.pdsi.breed,bands[1,])
nao.pdsi.test$bandp
nao.pdsi.test<-bandtest(nao.pdsi.breed,bands[2,])
nao.pdsi.test$bandp

#test coherence between ENSO and phdi during breeding
enso.phdi.test<-bandtest(enso.phdi.breed,bands[1,])
enso.phdi.test$bandp
enso.phdi.test<-bandtest(enso.phdi.breed,bands[2,])
enso.phdi.test$bandp

#test coherence between ENSO and pdsi during breeding
enso.pdsi.test<-bandtest(enso.pdsi.breed,bands[1,])
enso.pdsi.test$bandp
enso.pdsi.test<-bandtest(enso.pdsi.breed,bands[2,])
enso.pdsi.test$bandp

#test coherence between PDO and phdi during breeding
pdo.phdi.test<-bandtest(pdo.phdi.breed,bands[1,])
pdo.phdi.test$bandp
pdo.phdi.test<-bandtest(pdo.phdi.breed,bands[2,])
pdo.phdi.test$bandp

#test coherence between pdo and pdsi during breeding
pdo.pdsi.test<-bandtest(pdo.pdsi.breed,bands[1,])
pdo.pdsi.test$bandp
pdo.pdsi.test<-bandtest(pdo.pdsi.breed,bands[2,])
pdo.pdsi.test$bandp

#test for effects during metamorphosis

#test coherence between ENSO and phdi during metamorphosis
enso.phdi.test<-bandtest(enso.phdi.meta,bands[1,])
enso.phdi.test$bandp
enso.phdi.test<-bandtest(enso.phdi.meta,bands[2,])
enso.phdi.test$bandp

#test coherence between ENSO and pdsi during metamorphosis
enso.pdsi.test<-bandtest(enso.pdsi.meta,bands[1,])
enso.pdsi.test$bandp
enso.pdsi.test<-bandtest(enso.pdsi.meta,bands[2,])
enso.pdsi.test$bandp

#test coherence between PDO and phdi during metamorphosis
pdo.phdi.test<-bandtest(pdo.phdi.meta,bands[1,])
pdo.phdi.test$bandp
pdo.phdi.test<-bandtest(pdo.phdi.meta,bands[2,])
pdo.phdi.test$bandp

#test coherence between pdo and pdsi during metamorphosis
pdo.pdsi.test<-bandtest(pdo.pdsi.meta,bands[1,])
pdo.pdsi.test$bandp
pdo.pdsi.test<-bandtest(pdo.pdsi.meta,bands[2,])
pdo.pdsi.test$bandp
