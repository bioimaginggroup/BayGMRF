#' make GMRF precision matrix with weights
#'
#' @param weights weights
#' @param nn matrix with pair of neighbours per row.
#'
#' @return sparseMatrix of dimension max(nn)
#' @export
#'
#' @import parallel Matrix
#' @examples
#' Q<-weightedQ(1:4,matrix(c(1,1,2,3,2,3,3,4),ncol=2))
#' print(Q)
weightedQ<-function(weights,nn)
{
  n<-max(nn)
  if (length(weights)==1)weights<-rep(weights,dim(nn)[1])
  diags<-unlist(parallel::mclapply(1:n,function(i,weights,nn){sum(weights[nn[,1]==i])+sum(weights[nn[,2]==i])},weights,nn))
  return(Matrix::sparseMatrix(i=c(1:n,nn[,1]),j=c(1:n,nn[,2]),x=c(diags,-weights),symmetric=TRUE))
}
