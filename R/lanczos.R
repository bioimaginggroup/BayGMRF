#' lanczos algorithm for sparse Q
#'
#' @param Q
#'
#' @return
#' @export
#'
#' @examples
lanczos<-function(Q,r)
{
  n=dim(Q)[1]
  #V<-matrix(0,nrow=n,ncol=r+1)
  v.old=rep(0,n)
  beta=rep(0,r)
  alpha=rep(0,r)
  v<-rnorm(n,0,.1)
  v<-v/sqrt(sum(v^2))
  #V[,2]<-v
  for (j in 1:(r-1))
  {
    #w = (Q%*%V[,j+1])+beta[j]*V[,j]
    w = Q%*%v-beta[j]*v.old
    #w = Q%*%v #WIKI
    #alpha[j]=V[,j+1]%*%w
    alpha[j]=v%*%w
    #w=w-alpha[j]*V[,j+1]
    w=w[,1]-alpha[j]%*%v
    #w=w[,1]-alpha[j]%*%v-beta[j]*v.old #WIKI
    beta[j+1] = sqrt(sum(w*w))
    #V[,j+1] = w[,1]/beta[j+1]
    v.old=v
    v=w[1,]/beta[j+1]
    #V[,j+1] = w[,1]/beta[j+1]
  }
  w = Q%*%v+beta[r]*v.old
  #alpha[j]=V[,j+1]%*%w
  alpha[r]=v%*%w
  T <- Matrix::sparseMatrix(i=c(1:r,2:r),j=c(1:r,1:(r-1)),
                    x=c(alpha,beta[-1]), symmetric=TRUE)
  return(T)
  #return(list("V"=V,"beta"=beta,"alpha"=alpha))
}
