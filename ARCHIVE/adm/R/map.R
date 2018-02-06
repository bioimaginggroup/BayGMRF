map<-
function(data, map, cutpoints=seq(min(data),max(data),length=256), autoscale=FALSE, legend=TRUE, append=FALSE)
{
  if (autoscale)
  {
    data = (data-min(data))/(max(data)-min(data)+1e-8)
  }
  #cutpoints = c(-1e9,cutpoints, 1e9)
  
  farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))
  
  
  xmin <- 1:length(map)
  xmax <- 1:length(map)
  ymin <- 1:length(map)
  ymax <- 1:length(map)
  
  for(i in 1:length(map))
  {
    xmin[i] <- min(map[[i]][,2],na.rm=T)
    xmax[i] <- max(map[[i]][,2],na.rm=T)
    ymin[i] <- min(map[[i]][,3],na.rm=T)
    ymax[i] <- max(map[[i]][,3],na.rm=T)
  }
  
  breite <- c(min(xmin),max(xmax))
  hoehe <- c(min(ymin),max(ymax))
  
  if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")
  
  
  for(k in length(map):1)
  {
    polygon(map[[k]][,2],map[[k]][,3],col=farben[map[[k]][1,1]])
  }
  
  if (legend)
  {
    for(i in 1:8)
    {
      polygon(c(5800,5800,6100,6100),c(500+350*(i-1),500+350*i,500+350*i,500+350*(i-1)),col=gray(i/9))
    }
    for(i in 2:8)
    {
      text(6375,500+350*(i-1),as.character(round(cutpoints[i/8*length(cutpoints)], 2)),cex=.7,col=rgb(0,0,0))
    }
  }
}