weightedQ<-function(weights,coords)
{
  require(Matrix)
  require(parallel)
  n<-max(coords)
  diags<-unlist(mclapply(1:n,function(i,weights,coords){sum(weights[coords[,1]==i])+sum(weights[coords[,2]==i])},weights,coords))
  return(sparseMatrix(i=c(1:n,coords[,1]),j=c(1:n,coords[,2]),x=c(diags,-weights),symmetric=TRUE))
}
