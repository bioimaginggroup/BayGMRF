#' Prepare data from regular array with mask
#'
#' @param data data
#' @param mask mask
#'
#' @return list with data as vector, matrix of coordinates, matrix of neighbour pairs
#' @export
#' @import parallel
prepare.regular <- function(data, mask) {
  find.nn<-function(i,coords)
  {
    nn<-c()
    nc<-coord.index(coords,coords[i,]+c(1,0))
    if (length(nc)==1)nn<-c(i,nc)
    nc<-coord.index(coords,coords[i,]+c(0,1))
    if (length(nc)==1)nn<-c(nn,c(i,nc))
    return(nn)
  }
  
  coord.index<-function(coords,coord)
  {
    w1<-(coords[,1]==coord[1])
    w2<-(coords[,2]==coord[2])
    return(which(w1&w2))
  }
  

X<-dim(data)[1]
Y<-dim(data)[2]
T<-dim(data)[3]
N<-sum(mask)
newdata<-array(NA,c(N,T))
coords<-array(NA,c(N,2))
counter<-0
for (i in 1:X)
  for (j in 1:Y)
      if (mask[i,j])
      {
        counter<-counter+1
        coords[counter,]<-c(i,j)
        newdata[counter,]<-data[i,j,]
      }
nn <- matrix(unlist(parallel::mclapply(1:N,find.nn,coords)),nrow=2)
return(list("data"=newdata,"coords"=coords,"nn"=t(nn)))
}


