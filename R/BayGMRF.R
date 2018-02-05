#' Bayesian latent GMRF main function
#'
#' @param Y X
#' @param X X
#' @param E E
#' @param nn list of neighbours
#' @param mcmc mcmc options
#' @param model model
#' @param options options
#' @param hyper hyperparameters
#'
#' @import Matrix
#' @return list of samples
#' @export
#'
BayGMRF <- function(Y,X,Q.list,mcmc,options,hyper,model="Gaussian",E=1,Q.X=1)
{
  if (is.null(dim(Y)))Y <- array(Y,c(length(Y),1))
  N <- dim(Y)[1]  
  T <- dim(Y)[2]

  if (!is.list(X))X<-list(X)
  s <- length(X)
  p <- lapply(X,length)
  Q.s<-length(Q.list)
  n <- d <- rep(NA,Q.s)
  Q<-vector("mode"="list","length"=Q.s)
  for (i in 1:Q.s)
    {
    if (is.character(Q.list[[i]]))
    {
      Q[[i]]<-BayGMRF::makeQ(as.numeric(Q.list[[i]][2]),Q.list[[i]][1])
    }
    else
    {
    if (dim(Q.list[[i]])[2]==2)
      {
        Q[[i]]<-BayGMRF::weightedQ(rep(1,max(dim(Q.list[[i]]))),Q.list[[i]])
      }
      n[i]<-dim(Q.list[[i]])[1]
    }
    d[i] <- dim(Q[[i]])[1]
  }
  
  if (is.null(options$partial))
  {
    partial=vector("list",s)
    for (i in 1:s)partial[[i]]<-1:N
  }
  else
  {
    partial=options$partial
  }
  if (is.null(options$X.true))
  {
    X.t<-rep(TRUE,s)
  }
  else
  {
    X.t<-options$X.true
  }
  if (is.null(options$adaptive))
  {
    adaptive<-rep(FALSE,Q.s)
  }
  else
  {
    adaptive<-options$adaptive
    if(any(adaptive))
      {
      weights<-parallel::mclapply(1:Q.s, 
        function(i,adaptive,options,n){
          if(adaptive[i]){
            return(rep(options$nu.a[i]/options$nu.b[i],n[i]))
            }else{
            return(NULL)}},
      adaptive,options,n)
    }
  }
  
  
  nr.samples<-floor((mcmc$iterations-mcmc$burnin)/mcmc$thin)
  
  beta.sample <- vector("list",s)
  for (i in 1:s)beta.sample[[i]]<-array(NA,c(nr.samples,d[Q.X[i]]))
  beta<-parallel::mclapply(1:s,function(i,d,Q.X){rep(0,d[Q.X[i]])},d=d,Q.X=Q.X)
  sigma2.sample <- rep(NA,nr.samples)
  s.tau<-s
  if (!is.null(options$share.tau))
    {
    share.tau<-options$share.tau
    s.tau<-max(share.tau)
  }
  else
  {
    share.tau<-1:s
  }
  tau.sample <- array(NA,c(nr.samples,s.tau))

  tau <- hyper$tau.a/hyper$tau.b
  if (length(tau)!=s.tau)tau<-rep(tau,s.tau)
  sigma2 <- hyper$sigma2.b/hyper$sigma2.a
  XX <- lapply(X,function(X)sum(X^2))
  if (model=="Gaussian")resid <- Y
  if (model=="Poisson")eta <- resid<-log(Y+1)  
  if(options$print.time)start.time<-proc.time()
  cat("Starting iterations.\n")
  for (iteration in 1:mcmc$iterations){
    for (i in 1:s)
    {
      P <- tau[share.tau[i]]*Q[[Q.X[i]]] + Matrix::Diagonal(d[Q.X[i]],XX[[i]]/sigma2)
      if (X.t[i])
        {
        resid[partial[[i]],] <- resid[partial[[i]],]+beta[[i]]*X[[i]]
        }
      else
        {
          resid[partial[[i]],] <- resid[partial[[i]],]+beta[[i]]%*%t(X[[i]])
        }
      if (model=="Gaussian")b <- resid[partial[[i]],]%*%X[[i]]/sigma2
      if (model=="Poisson")b <- hyper$eta*eta[partial[[i]],]
      beta[[i]] <- as.vector(BayGMRF::rmvnormcanon(n = 1, b = b, P = P))
      resid[partial[[i]],] <- resid[partial[[i]],]-beta[[i]]%*%t(X[[i]])
    }  
    
    tau<-unlist(lapply(1:s.tau,update.tau,tau,beta,Q,Q.X,d,hyper,options))
    
    if(options$sigma2){
      aa <- hyper$sigma2.a + N*T/2
      bb <- hyper$sigma2.b + as.vector(0.5*sum(resid^2))
      sigma2 <- 1/rgamma(n = 1,  shape = aa, rate = bb)
    }
    
    if (model=="Poisson")
      {
        m.temp<-log(E)*hyper$eta + Y
        for (i in 1:s)
          {
          m.temp<- beta
        }
        eta<-sample.eta(eta, m.temp, hyper$eta, E)
      }

    sample<-(iteration-mcmc$burnin)/mcmc$thin
    if ((iteration%%mcmc$thin) == 0 && sample > 0 && sample <= nr.samples){ 
      for (i in 1:s)
        {
        beta.sample[[i]][sample,] <- beta[[i]]
      }
      tau.sample[sample,] <- tau
      sigma2.sample[sample] <- sigma2
    }    
    
    if(any(adaptive))
    {
      Q<-parallel::mclapply(1:Q.s,update.adaptive)
    }
    
    if(options$print.time){
      time<-(proc.time()-start.time)[3]
      estimatedtime<-time/iteration*mcmc$iterations-time
      for (i in 1:200)cat("\b")
      cat(paste(roundtime(time),"elapsed -",roundtime(estimatedtime),"to go."))
    }
    
    
    
  }#end iteration
  
  returnlist<-list("beta"=beta.sample,"tau"=tau.sample,"sigma2"=sigma2.sample)
  return(returnlist)
}
