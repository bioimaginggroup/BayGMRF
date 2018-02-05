det2<-function(Q)
{
  e<-eigen(Q)$values
  e<-e[e>1e-10]
  loge<-sum(log(e))
  rank<-length(e)
  return(list("logdet"=loge,"rank"=rank))
}
