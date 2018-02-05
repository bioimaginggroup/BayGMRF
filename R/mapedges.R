map.edges<-function(map,col,add=FALSE)
{
  require(parallel)
  N<-length(map)
  if(!add)
    {
    .range<-function(mappart)
  {
    return(apply(mappart,2,range))
  }
  all.range<-unlist(mclapply(map,.range))
  all.range<-array(all.range,c(4,N))
  minX<-min(all.range[1,],na.rm=TRUE)
  maxX<-max(all.range[2,],na.rm=TRUE)
  minY<-min(all.range[3,],na.rm=TRUE)
  maxY<-max(all.range[4,],na.rm=TRUE)
  plot(c(minX,maxX),c(minY,maxY),col="white",axes=FALSE,xlab="",ylab="")
  }
  
  for (i in 1:N)
    lines(map[[i]][,2:3],col=col[i],lwd=3)

}