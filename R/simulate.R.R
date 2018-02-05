simulate<-function(dim, mu, c, x, beta){

  X<-array(0,dim)
  X<-mu+c*x*beta
  return(X)
}
