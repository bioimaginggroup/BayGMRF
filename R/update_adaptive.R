update.adaptive<-function(...)
{
  return(switch(method,
     mh = update.adaptive.mh(...) 
  ))
}


sum.for.diag<-function(weights,nn)
{
  p<-max(nn)
  w<-parallel::mclapply(1:p,function(i,nn)return(c(which(nn[,1]==i),which(nn[,2]==i))),nn)
  d<-unlist(parallel::mclapply(w,function(i,weights)return(sum(weights[i])),weights))
  return(d)
}

update.adaptive.mh<-function(weights,coords,tau,beta,nu.a,nu.b)
{
    n<-length(weights)
    p <- max(Q.list)
    diags<-sum.for.diag(weights,coords)
    Q <- Matrix::sparseMatrix(i=c(1:p,coords[,1]),j=c(1:p,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
    ev<-log(diag(chol(Q[-p,-p])))
    for (i in n:1)
    {
      b.temp <- (tau*(beta[coords[i,1]]-beta[coords[i,2]])^2/2 + nu.b)
      a.temp <- nu.a+1/2 
      weights.proposed <- weights
      weights.proposed[i] <-rgamma(1,a.temp,b.temp) 
      diags.proposed<-diags
      diags.proposed[coords[i,1]]<-diags.proposed[coords[i,1]]-weights[i]+weights.proposed[i]
      diags.proposed[coords[i,2]]<-diags.proposed[coords[i,2]]-weights[i]+weights.proposed[i]
      Q.proposed<-Matrix::sparseMatrix(i=c(1:p,coords[,1]),j=c(1:p,coords[,2]),x=c(diags.proposed,-weights.proposed),symmetric=TRUE)
        
      start<-min(coords[i,])
      new.ev<-log(diag(chol(Q.proposed[start:(p-1),start:(p-1)])))
      a<-sum(ev[start:(p-1)]-new.ev)
      a<-exp(a)*sqrt(weights.proposed[i]/weights[i])
      if(runif(1)<a)
      {
        diags<-diags.proposed
        weights<-weights.proposed
        Q<-Q.proposed
        ev[start:(p-1)]<-new.ev
      }
      }
    }
