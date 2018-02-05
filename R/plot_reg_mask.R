#' Plot 2D image from coordinate data
#'
#' @param x values
#' @param coords coordinates (2D)
#' @param dims dimension of image
#'
#' @export
#' @import Matrix
#' @export
img.reg.mask<-function(x,coords,dims=apply(coords,2,max))
{
  m<-Matrix::sparseMatrix(i=coords[,1],j=coords[,2],x = x,dims=dims)
  image(m,xlab="",ylab="")
}