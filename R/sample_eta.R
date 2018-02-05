sample.eta<-function(eta.sample.aktuell, m.temp, tau.sample.aktuell,e)
{
  require(truncnorm)
  n<-length(eta.sample.aktuell)
  x <- eta.sample.aktuell+log(e)
  u <- rexp(n)+exp(x)
  eta<-rtruncnorm(n,b=log(u),mean=(m.temp+log(e))/tau.sample.aktuell,sd=sqrt(1/tau.sample.aktuell))
  return(eta-log(e))
}
