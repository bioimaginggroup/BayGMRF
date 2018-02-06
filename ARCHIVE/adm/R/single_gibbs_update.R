single.gibbs.update<-function(i,coords,d,Q1,weights,kappa,gamma,nu.b)
{
  xy<-as.vector(coords[i,])
  xy<-xy[xy!=d]
  
  Q0<-Q1[xy,xy]
  
  if (length(xy)==2)
  {
    zaehler = diff(Q0)
    zaehler = weights[i]*t(zaehler)%*%zaehler
    qtilde0=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
    Q0=Q0+zaehler/(1-weights[i]*qtilde0)
    qtilde=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
  }
  else
  { 
    zaehler = weights[i]*t(Q0)%*%Q0
    qtilde0=Q0
    qtilde=Q0+zaehler/(1-weights[i]*Q0)
  }
  if (length(xy)==2){
    b.temp <- (kappa*(gamma[xy[1]]-gamma[xy[2]])^2/2)
  }
  else{
    b.temp <- (kappa*(gamma[xy]-gamma[d])^2/2)
  }
  a.temp <- 1.5
  #print(c(i,b.temp))

#write(c(qtilde,weights[i]),file="qt.txt",append=TRUE)
  test<-rstgamma(n=1,shape=a.temp,rate=b.temp,shift=-1/qtilde,a=0,b=1)
  return(test)

}
