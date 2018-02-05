library(Matrix)

#Funktion zum Ziehen aus Normalverteilung mit Präzision Q und Erwartungswert Q^-1 m
rcannorm<-function (m, Q) 
{
  L = chol(Q)
  w = solve(t(L), m)
  mu = solve(L, w)
  z = rnorm(length(m))
  v = solve(L, z)
  return = mu + v
}

T<-100
nriter<-600 # Anzahl Iterationen
burnin<-100 # Anzahl Iterationen Burn-in, damit effektiv 500 Iterationen

# Simuliere I Zeitreihen
I<-5
t<-seq(0,100,length=T)
somedata<-t*cos(t)
plot(t,somedata,type="l")
points(t,somedata)

data<-array(NA,c(T,I))
for (i in 1:I)data[,i]<-rnorm(T,somedata,10)
plot(t,data[,1])

# Vorbesetzungen der Parameter
lambda<-1+(1:(T-1)/T)
Q<-sparseMatrix(i=c(1:T,1:(T-1)),j=c(1:T,2:T),x=c(c(lambda,0)+c(0,lambda),-lambda),symmetric=TRUE)
tau<-1
sigma2<-1
logdetQ<-determinant(Q,logarithm=TRUE)$modulus

# Arrays zum Abspeichern
theta.save<-array(NA,c(nriter-burnin,T))
lambda.save<-array(NA,c(nriter-burnin,T-1))
tau.save<-rep(NA,nriter-burnin)
sigma2.save<-rep(NA,nriter-burnin)

data.sum<-apply(data,1,sum)

system.time({ # Laufzeit messen
  for (iter in 1:nriter)
  {
    # Gibbs-Sampler für theta (MVNorm)
    V<-tau*Q+diag(I/sigma2,T)
    m<-data.sum/sigma2
    theta<-rcannorm(m,V)
    
    # Gibbs-Sampler für tau (Gamma)
    tau<-rgamma(1,1+(T-1)/2,1+0.5*(t(theta)%*%Q%*%theta)[1,1])
    
    # Gibbs-Sampler für sigma^2 (InvGamma)
    bb<-0.001
    for (i in 1:I)bb<-bb+0.5*sum((data[,i]-theta)^2)
    sigma2<- 1/rgamma(1,1+I*T/2,bb)
    
    # M-H-Sampler für Gewichte
    for (i in 1:(T-1))
    {
      lprop<-lambda
      lprop[i]<-rgamma(1,1,1+0.5*(theta[i+1]-theta[i])^2)
      Qprop<-sparseMatrix(i=c(1:T,1:(T-1)),j=c(1:T,2:T),x=c(c(lprop,0)+c(0,lprop),-lprop),symmetric=TRUE)
      logdetQprop<-determinant(Qprop,logarithm=TRUE)$modulus
      alpha<-sqrt(exp(logdetQprop-logdetQ))
      if(is.na(alpha))alpha<-0
      if(runif(1)<alpha)
      {
        lambda<-lprop
        logdetQ<-logdetQprop
      }
    }
    # Berechne Q aus Gewichten
    Q<-sparseMatrix(i=c(1:T,1:(T-1)),j=c(1:T,2:T),x=c(c(lambda,0)+c(0,lambda),-lambda),symmetric=TRUE)
    
    # Nach dem burnin Werte speichern
    if (iter>burnin)
    {
      theta.save[iter-burnin,]<-theta[,1]
      lambda.save[iter-burnin,]<-lambda
      tau.save[iter-burnin]<-tau
      sigma2.save[iter-burnin]<-sigma2
    }
  }
})

plot(theta.save[,12]) # Prüfe Mixing
plot(tau.save) # Prüfe Mixing
plot(lambda.save[,12]) # Prüfe Mixing


plot(somedata,type="l",lwd=2) # Wahre Kurve
for (i in 1:I)points(data[,i],cex=0.5,col=i) # Daten
lines(apply(theta.save,2,mean),col="green",lwd=2) # Punktschätzer
lines(apply(theta.save,2,median),col="blue",lwd=2) # Punktschätzer 
lines(apply(theta.save,2,quantile,.025),col="blue",lwd=1,lty=2) # Punktweise CI
lines(apply(theta.save,2,quantile,.975),col="blue",lwd=1,lty=2)

plot(apply(lambda.save,2,median),type="l") # Punktschätzer Gewichte
