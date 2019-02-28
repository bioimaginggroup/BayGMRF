#' Prepare data from regular array with mask
#'
#' @param data data
#' @param mask mask
#'
#' @return list with data as vector, matrix of coordinates, matrix of neighbour pairs
#' @export
#' @import parallel
prepare.regular <- function(data, mask) {
  threeD<-FALSE
  if (length(dim(mask))==3)threeD<-TRUE
    
  find.nn<-function(i,coords)
  {
    nn<-c()
    nc<-coord.index(coords,coords[i,]+c(1,0))
    if (length(nc)==1)nn<-c(i,nc)
    nc<-coord.index(coords,coords[i,]+c(0,1))
    if (length(nc)==1)nn<-c(nn,c(i,nc))
    return(nn)
  }
  
  find.nn.3d<-function(i,coords)
  {
    nn<-c()
    nc<-coord.index.3d(coords,coords[i,]+c(1,0,0))
    if (length(nc)==1)nn<-c(i,nc)
    nc<-coord.index.3d(coords,coords[i,]+c(0,1,0))
    if (length(nc)==1)nn<-c(nn,c(i,nc))
    nc<-coord.index.3d(coords,coords[i,]+c(0,0,1))
    if (length(nc)==1)nn<-c(nn,c(i,nc))
    return(nn)
  }

  coord.index<-function(coords,coord)
  {
    w1<-(coords[,1]==coord[1])
    w2<-(coords[,2]==coord[2])
    return(which(w1&w2))
  }
  
  coord.index.3d<-function(coords,coord)
  {
    w1<-(coords[,1]==coord[1])
    w2<-(coords[,2]==coord[2])
    w3<-(coords[,3]==coord[3])
    return(which(w1&w2&w3))
  }
  
X<-dim(data)[1]
Y<-dim(data)[2]
T<-dim(data)[3]
if (threeD)Z<-dim(data)[3]
if (threeD)T<-dim(data)[4]
N<-sum(mask)
newdata<-array(NA,c(N,T))
coords<-array(NA,c(N,2+threeD))
counter<-0
if (!threeD)
for (i in 1:X)
  for (j in 1:Y)
      if (mask[i,j])
      {
        counter<-counter+1
        coords[counter,]<-c(i,j)
        newdata[counter,]<-data[i,j,]
      }
if (threeD)
  for (i in 1:X)
    for (j in 1:Y)
      for (k in 1:Z)
        if (mask[i,j,k])
      {
        counter<-counter+1
        coords[counter,]<-c(i,j,k)
        newdata[counter,]<-data[i,j,k,]
      }

if(!threeD)nn <- matrix(unlist(parallel::mclapply(1:N,find.nn,coords)),nrow=2)
if(threeD)nn <- matrix(unlist(parallel::mclapply(1:N,find.nn.3d,coords)),nrow=2)
return(list("data"=newdata,"coords"=coords,"nn"=t(nn)))
}


