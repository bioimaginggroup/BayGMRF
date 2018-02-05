#' Make precision matrix
#'
#' @param n dimension of Q
#' @param type "rw1" or "I" or "rw1+I"
#'
#' @return sparseMatrix (Matrix package)
#' @export
#' @import Matrix
#'
#' @examples
makeQ<-function(n, type, params=NULL)
{
  if (type=="rw1")
  {
    Q<-Matrix::sparseMatrix(i = c(1:n, 2:n), j=c(1:n, 1:(n-1)),
                            x = c(1,rep(2,n-2),1,rep(-1,n-1)), symmetric = TRUE)
  }
  if (type=="I")
  {
    Q<-Matrix::sparseMatrix(i = c(1:n), j=c(1:n), x = rep(1,n), symmetric = TRUE)
  }
  if (type=="rw1+I")
  {
    Q<-Matrix::sparseMatrix(i = c(1:n, 2:n), j=c(1:n, 1:(n-1)),
                            x = c(1+params[1],rep(2+params[1],n-2),1+params[1],rep(-1,n-1)), symmetric = TRUE)
  }
  return(Q)
}
