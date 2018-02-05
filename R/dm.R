dm <- function(y, e, coords, nr.it, burnin, thin, nu=2, do.kappa=TRUE,  kappa=1, ab.kappa=c(1,1), tau=1000, print.time=FALSE){
  
  
  d <- length(y)  
  nu.a<-nu.b<-nu/2
  if(length(nu)==2)
  {
    nu.a<-nu[1]
    nu.b<-nu[2]
  }
  p<-max(dim(coords))
  
  gamma <- vector(length = d) 
  eta <- vector(length = d) 
  
  gamma <- rep(x=0, time=d)  
  
  Q<-weightedQ(rep(1,p),coords)
  
  
  diagM <- array(0,c(p,d))
  for (i in 1:p)
  {
    diagM[i,as.vector(coords[i,])]<-1
  }
  

  for(u in 1:length(y)){
    if(y[u] == 0){
      eta[u] <- 0
    }
    else
    {
      eta[u] <- log(y[u]/e[u])
    }
  }  
  
  nr.samples<-floor((nr.it-burnin)/thin)
  gamma.sample <- array(NA,c(nr.samples,d))
  alpha.sample <- array(NA,c(nr.samples,d))
  eta.sample <- array(NA,c(nr.samples,d))
  
  kappa.sample <- rep(NA,nr.samples)
  
  
  for (iteration in 1:nr.it){
    if(print.time)start.time<-proc.time()
    if(iteration%%10==0){cat(".")}
    P <- kappa*Q + diag(tau,d)
    b <- tau*eta
    gamma <- as.vector(rmvnorm(n = 1, b = b, P = P))
    
    # kappa
    if(do.kappa){
      aa <- ab.kappa[1] + (d-1)/2
      bb <- ab.kappa[2] + as.vector(0.5*(t(gamma)%*%Q%*%gamma))
      kappa <- rgamma(n = 1,  shape = aa, rate = bb)
    }
    
    # eta
    m.temp <- (gamma +log(e))*tau + y
    eta<-sample.eta(eta, m.temp, tau, e)
    
    
    
      


    sample<-(iteration-burnin)/thin
    if ((iteration%%thin) == 0 && sample > 0 && sample <= nr.samples){ 
      gamma.sample[sample,] <- gamma
      kappa.sample[sample] <- kappa
      eta.sample[sample,] <- eta
      alpha.sample[sample,] <- eta-gamma

    }    
    
    if(print.time)print((proc.time()-start.time)[3])
  }#end iteration
  
  returnlist<-list("gamma"=gamma.sample,"kappa"=kappa.sample,"eta"=eta.sample,
                   "alpha"=alpha.sample)
  return(returnlist)
  
} 