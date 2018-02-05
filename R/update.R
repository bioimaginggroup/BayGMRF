#' Update tau
#'
#' @param i 
#' @param tau 
#' @param beta 
#' @param Q 
#' @param Q.X 
#' @param d 
#' @param hyper 
#' @param options 
#'
#' @return
#' @export
#'
#' @examples
update.tau<-function(i,tau,beta,Q,Q.X,d,hyper,options)
{
  if(options$tau[i]){
    if (!is.null(options$share.tau))return(update.tau.share(i,tau,beta,Q,Q.X,d,hyper,options))
    aa <- hyper$tau.a[i] + (d[Q.X[i]]-1)/2
    bb <- hyper$tau.b[i] + as.vector(0.5*(t(beta[[i]])%*%Q[[Q.X[i]]]%*%beta[[i]]))
    return(rgamma(n = 1,  shape = aa, rate = bb))
  }
  else
  {
    return(tau[i])
  }
}

update.tau.share<-function(i,tau,beta,Q,Q.X,d,hyper,options)
{
  aa <- hyper$tau.a[i]
  bb <- hyper$tau.b[i]
  w<-which(options$share.tau==i)
  for (j in w)
  {
    aa <- aa + (d[Q.X[j]]-1)/2
    bb <- bb + as.vector(0.5*(t(beta[[j]])%*%Q[[Q.X[j]]]%*%beta[[j]]))
  }
  temp <- rgamma(n = 1,  shape = aa, rate = bb)
  if(is.na(temp))temp<-hyper$tau.a[i]/hyper$tau.b[i]
  return(temp)
}
