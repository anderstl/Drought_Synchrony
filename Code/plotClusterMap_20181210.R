#' Map clusters for a \code{clust} object
#' 
#' Produces a map of the locations of sampling for a \code{clust} object, with models
#' indicated with colors.
#' 
#' @param inclust A \code{clust} object, as created with \code{clust}
#' @param spltlvl The split level in the clustering to use. This is the index of inclust$clusters.
#' Default the final split.
#' @param basemap If provided, a \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' to plot beneath the node locations. See Details.
#' @param nodewgt Weights controlling the sizes of nodes. If \code{NULL} (the default),
#' node size does not vary. If \code{mod.decomp}, the (scaled) decomposition of the node's 
#' contribution to the total modularity is used as the weight. Optionally, the user may instead 
#' provide a vector with length equal to the number of locations.
#' @param nodesize A length = 2 vector giving the minimum and maximum node size for plotting. Defaults to c(1,3).
#' If \code{is.null(nodewgt)}, the first entry (the minumum) controls node size.
#' @param nodecol A vector of colors to identify modules. Defaults to \code{"Set1"} from 'RColorBrewer'.
#' @param edgewgt Weights controlling the width of edges (links) between nodes. If \code{NULL} 
#' (the default), no edges are drawn. If a single number q between 0 and 1 is provided, plot edges 
#' with synchrony values greater than the qth quantile of \code{inclust$adj}. 
#' If \code{"adj"}, use the adjacency matrix to produce the weights; these may be binary. If a matrix having 
#' dimensions identical to \code{inclust$adj}, the matrix entries give the widths of the plotted lines.
#' @param edgesize A length = 2 vector giving the minimum and maximum edge width for plotting. Defaults to c(0.1,1).
#' @param addlegend \code{TRUE} (the default) or \code{FALSE}, add a legend in the right-hand margin.
#' @param filename a filename, possibly including path info, but without a file extension. If present,
#' exports the plot as a .pdf using the specified filename.
#' 
#' @return \code{plotClusterMap} produces a cluster plot map.
#' 
#' @details \code{basemap} must be in the same coordinate system as the location coordinates.
#' Objects of class \code{SpatialPolygons} and \code{SpatialPolygonsDataFrame} are produced by the
#' 'sp' package. 
#' 
#' @export

plotClusterMap<-function(inclust,spltlvl=length(inclust$clusters),basemap=NULL,nodewgt=NULL,
                         nodesize=c(1,3),nodecol=brewer.pal(9,"Set1"), edgewgt=NULL, edgesize=c(0.1,1),
                         addlegend=TRUE, filename=NA)
{
  
  #some checking of validity of inputs
  if(!is.null(nodecol)){
    if(length(unique(nodecol)) < max(unlist(inclust$clusters[spltlvl]))){
      stop("more clusters than colors in nodecol")
    }
  }
  
  #convert inclust$coords to common format
  if(all(c("X","Y") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$X,
                                                                     Y=inclust$coords$Y)}
  if(all(c("lat","lon") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$lon,
                                                                        Y=inclust$coords$lat)}
  if(all(c("latitude","longitude") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$longitude,
                                                                                   Y=inclust$coords$latitude)}
  
  if(!is.na(filename)){
    pdf(paste0(filename,".pdf"))
  }
  
  #expand right side margin if a legend is being added
  if(addlegend){
    par.mar<-par("mar")
    mar.new<-par.mar
    mar.new[4]<-6.1
    par(mar=mar.new,xpd=T)
  }
  
  if(!is.null(basemap)){
    #give a warning if basemap extent does not contain all points
    ext<-bbox(basemap)
    if(any(c(min(coords$X)<ext[1,1],
             max(coords$X)>ext[1,2],
             min(coords$Y)<ext[2,1],
             max(coords$Y)>ext[2,2]))
    ){
      warning("basemap extent does not contain all nodes; possible coordinate system mismatch?")
    }
    plot(basemap, axes=T)
  }
  
  if(is.null(basemap)){
    plot(coords$X,coords$Y,pch=16,col=NA,xlab="X",ylab="Y")
  }
  
  x0<-rep(coords$X, each=nrow(coords))
  y0<-rep(coords$Y, each=nrow(coords))
  x1<-rep(coords$X, times=nrow(coords))
  y1<-rep(coords$Y, times=nrow(coords))
  
  if(!is.null(edgewgt)){
    if(length(edgewgt)!=1 & length(edgewgt)!=length(inclust$adj)){
      stop("invalid 'edgewgt' argument")
    }
    if(all(dim(edgewgt)==dim(inclust$adjmat))){
      edgewgts<-edgewgt
    }
    if(edgewgt=="adj"){
      ewgt<-inclust$adj-min(inclust$adj,na.rm=T)
      ewgt<-ewgt/max(ewgt,na.rm=T)
      edgewgts<-ewgt*(edgesize[2]-edgesize[1]) + edgesize[1]
    }
    if(length(edgewgt)==1 & is.numeric(edgewgt)){
      wgtthresh<-quantile(inclust$adj,edgewgt,na.rm=T)
      edgewgts<-rep(0,prod(dim(inclust$adj)))
      edgewgts[c(inclust$adj)>wgtthresh]<-1
    }

    segments(x0,y0,x1,y1,lwd=edgewgts)
  }
  
  if(is.null(nodewgt)){
    points(coords[,1], coords[,2], pch=16, cex=nodesize[1], col=nodecol[unlist(inclust$clusters[spltlvl])])
  }
  else{ 
    if(nodewgt=="mod.decomp"){
      membwgt<-inclust$modres[[spltlvl]]$nodeQ-min(inclust$modres[[spltlvl]]$nodeQ)
      membwgt<-membwgt/max(membwgt)
      nodecex<-membwgt*(nodesize[2]-nodesize[1]) + nodesize[1]
    }
    if(length(nodewgt)==nrow(coords)){
      membwgt<-nodewgt-min(nodewgt)
      membwgt<-membwgt/max(membwgt)
      nodecex<-membwgt*(nodesize[2]-nodesize[1]) + nodesize[1]
    }
    points(coords[,1], coords[,2], pch=16, cex=nodecex, col=nodecol[unlist(inclust$clusters[spltlvl])])
  }

  if(addlegend){
    legx<-par('usr')[2] + 0.01*abs(diff(par('usr')[1:2]))
    legy1<-par('usr')[4]
    leg1<-legend(legx,legy1,legend=paste0("module ",1:max(unlist(inclust$clusters[spltlvl]))), pch=16, col=
             nodecol[1:max(unlist(inclust$clusters[spltlvl]))],title="Membership",bty="n")
    if(!is.null(nodewgt)){
      legy2<-legy1 - leg1$rect$h
      labs=round(c(min(inclust$modres[[spltlvl]]$nodeQ,na.rm=T),
                   mean(inclust$modres[[spltlvl]]$nodeQ,na.rm=T),
                   max(inclust$modres[[spltlvl]]$nodeQ,na.rm=T)),digits=3)
      sizes=c(min(nodecex),mean(nodecex),max(nodecex))
      leg2<-legend(legx,legy2,legend=labs,pt.cex=sizes,pch=1,title="Node weight",bty="n")
    }
    if(!is.null(edgewgt)){
      if(edgewgt=="adj" | is.matrix(edgewgt)){
        if(is.null(nodewgt)){legy3<-legy1-leg1$rect$h}
        else{legy3<-legy1-leg1$rect$h-leg2$rect$h}
        labs2=round(c(min(inclust$adj,na.rm=T),mean(inclust$adj,na.rm=T),
                     max(inclust$adj,na.rm=T)),digits=3)
        wdths=c(min(edgewgts,na.rm=T),mean(edgewgts,na.rm=T),max(edgewgts,na.rm=T))
        legend(legx,legy3,legend=labs2,lty=1,lwd=wdths,title="Edge weight",bty="n")
      }
    }
    par(mar=par.mar) #reset 'mar' graphics parameter
  }
  
  if(!is.na(filename)){dev.off()}

}