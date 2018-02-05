scotland.data<-read.table("data/Scotland.txt",header=TRUE)
scotland.graph<-scan("data/scotland.graph")

scotland.coords<-c()
l<-length(scotland.graph)
counter<-1
scotland.non<-c()
while(counter<l)
{
  counter<-counter+1
  i<-scotland.graph[counter]
  counter<-counter+1
  if (scotland.graph[counter]==0)
    {
    scotland.non=c(scotland.non,i)
  }
  else
  {
  for (j in 1:scotland.graph[counter])
  {
    counter<-counter+1
    k<-scotland.graph[counter]
    if(i<k)
      scotland.coords<-rbind(scotland.coords,c(i,k))
  }
  }
}

for (j in rev(scotland.non))
{
  scotland.coords[scotland.coords>j]<-scotland.coords[scotland.coords>j]-1
}
scotland.data$Region<-scotland.data$Region+1
scotland.data<-scotland.data[-scotland.non,]
for (j in rev(scotland.non))
{
  scotland.data$Region[scotland.data$Region>j]<-scotland.data$Region[scotland.data$Region>j]-1
}

library("shapefiles")
scotland.shp<-read.shp("data/scot.shp")

scotland.shape<-vector(length=52,mode="list")
counter<-0
counter2<-0
for (i in 1:56)
{
  if (!(any(as.integer(i)==scotland.non)))
  {
    counter<-counter+1   
    counter2<-counter2+1   
    grenzen<-scotland.shp$shp[[i]]$points
    grenzen<-cbind(counter,grenzen)
    parts<-c(scotland.shp$shp[[i]]$parts,dim(grenzen)[1])
    for (j in 1:scotland.shp$shp[[i]]$num.parts)
    scotland.shape[[counter2]]<-grenzen[(parts[j]+1):parts[j+1],]
  }
}

map(1:52,scotland.shape)
