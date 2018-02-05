#' Make a 3D array from a list of coodinates
#'
#' @param x list of values
#' @param coords list of coordinates
#' @param dims first two dimensions
#'
#' @import oro.nifti Matrix
#' @return
#' @export
make.3D.from.coords<-function(x,coords,dims)
{
  n<-length(x)
  m<-array(NA,c(dims,n))
  for (i in 1:n)
  {
    s<-dim(coords[[i]])[1]
    for (j in 1:s)
      m[coords[[i]][j,1],coords[[i]][j,2],i]<-x[[i]][j]
  }
  return(oro.nifti::as.nifti(m))
}
