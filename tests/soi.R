library(BayGMRF)
data("soi-sim",package="BayGMRF")
z<-FixEf
zz<-Matrix(t(FixEf)%*%FixEf)
x<-Matrix(X)
y <- Matrix(Y,ncol=nPat)

n <- apply(delta,1,function(x)sum(x!=0))
i <- rep(1:nVox,n)
j <- as.vector(t(delta))
j <- j[j!=0]

skip <- j>i
i <- i[!skip]
j <- j[!skip]

n1<-length(i)

i<-c(1:nVox,i)
j<-c(1:nVox,j)
ein <- c(n,rep(-1,n1))
Q<-sparseMatrix(i,j,x = ein,symmetric = TRUE)

test<-soi(Y,X,FixEf,Q,all=TRUE,goldsmith=FALSE,mcmc.controls=list("nr.iterations"=1000,"burnin"=100))
  

index<-function(i,j,k)
  return((k-1)*50*50+(j-1)*50+i)

selection1<-c()
for (i in seq(1,50,by=2))
  for (j in seq(1,50,by=2))
    for (k in seq(1,5,by=2))
    selection1<-c(selection1,index(i,j,k))
for (i in seq(2,50,by=2))
  for (j in seq(2,50,by=2))
    for (k in seq(1,5,by=2))
      selection1<-c(selection1,index(i,j,k))
for (i in seq(2,50,by=2))
  for (j in seq(1,50,by=2))
    for (k in seq(2,5,by=2))
      selection1<-c(selection1,index(i,j,k))
for (i in seq(1,50,by=2))
  for (j in seq(2,50,by=2))
    for (k in seq(2,5,by=2))
      selection1<-c(selection1,index(i,j,k))


selection2<-c()
for (i in seq(1,50,by=2))
  for (j in seq(2,50,by=2))
    for (k in seq(1,5,by=2))
      selection2<-c(selection2,index(i,j,k))
for (i in seq(2,50,by=2))
  for (j in seq(1,50,by=2))
    for (k in seq(1,5,by=2))
      selection2<-c(selection2,index(i,j,k))
for (i in seq(2,50,by=2))
  for (j in seq(2,50,by=2))
    for (k in seq(2,5,by=2))
      selection2<-c(selection2,index(i,j,k))
for (i in seq(1,50,by=2))
  for (j in seq(1,50,by=2))
    for (k in seq(2,5,by=2))
      selection2<-c(selection2,index(i,j,k))

gamma <- gammasum <- rep(1,nVox)
beta <- betasum <- beta2sum <- beta3sum <- rep(0,nVox)
alpha <- alphasum <- alpha2sum <- alpha3sum <- rep(0,nPat)

burnin = 100
niter = 1000

all = FALSE
goldsmith=FALSE

tau = 1/100
sigma2 = 1/1000
a=-4
b=.8

par(mfrow=c(2,2))
for (iter in 1:niter)
{
  cat("\b\b\b\b\b\b\b")
  cat(iter)
  
  m <- (y - beta%*%x)%*%z
  m <- t(m)/sigma2
  alpha <- rmvnorm(1,m,zz/sigma2)
  
  if(all)
    {
    m <- x%*%t(y-alpha%*%t(z))
  m <- m/sigma2
  Q0 <- xx/sigma2 + Q/tau
  w<-gamma==1
  Q0 <- Q0[w,w]
  m<-m[w]
  beta[w] <- rmvnorm(1,m,Q0)
  beta[-w] <- 0
  }
  
  if(!all&!goldsmith)
    {
  Q0 <- xx/sigma2 + Q/tau
  m <- x%*%t(y-alpha%*%t(z))
  m <- m/sigma2
  
  for (i in sample(nVox))
  {
    if (gamma[i]==0)
      {
      beta[i]=0.0
      }
    else
      {
        Qaa <- Q0[i,i]
        Qab <- Q0[-i,i]%*%beta[-i]
        beta[i] = rnorm(1,(m[i]-Qab)/Qaa,sqrt(1/Qaa))
      }
  }
  
  res <- (y-alpha%*%t(z))
  res2 <- res%*%t(res)
  
  for (i in sample(nVox))
  {
    nn <- delta[i,]
    nn <- nn[nn!=0]
    g=exp(-a+b*sum(1-2*gamma[nn]))
    pi = 1/(1+g)
    gammastar <- runif(1,g)
    if ((gamma==0)&(gammastar==1))
      {
        wstar<-w<-(gamma==1)
        Xw<-Matrix(X[w,])
        S <- res2 - res%*%t(Xw)%*%solve(Xw%*%t(Xw))%*%Xw%*%t(res)
        S1 <- 
        logL1 <- -.5*nVox*(y-alpha%*%t(z)) 
        L0 <- 1 #??
      }
  }
  } #if all
  
  resid=as.vector(y-alpha%*%t(z)-beta%*%X)
  sumgamma = sum(gamma)
  
  if (goldsmith)
  {
    for (i in sample(nVox))
    {
      nn <- delta[i,]
    nn <- nn[nn!=0]
    mu <- (y-alpha%*%t(z)-beta[-i]%*%X[-i,])/sigma2 #time!
    sumbeta<-sum(beta[nn])
    lnn<-length(nn)
    mu<- mu%*%X[i,] + sumbeta/tau
    sig <- t(X[i,])%*%X[i,]/sigma2 + lnn/tau
    sig<-1/sig[1,1]
    mu<-as.vector(mu*sig)
    betastar <-rnorm(1,mu,sqrt(sig))

    res1<-resid+beta[i]*X[i,]
    res2<-Matrix(res1-betastar*X[i,])
    res1<-Matrix(res1)

    loggl <- -.5*(res1%*%t(res1)-res2%*%t(res2))/sigma2
      
      loggl<- loggl+length(nn)*.5*(betastar-sumbeta/lnn)^2/tau
      loggl <- loggl-a+b*(nVox-2*sumgamma)

      gl <-exp(loggl)*sqrt(2*pi*tau/length(nn))
      
      gl=gl[1,1]
      gl=min(gl,1)
      gl=max(gl,0)
      if(is.na(gl))gl=.5
      sumgamma<-sumgamma-gamma[i]
      gamma[i]<-rbinom(1,1,gl)
      if (gamma[i]==0){
        beta[i]<-0
        resid<-res1
      }
      if (gamma[i]==1)
        {
        beta[i]=betastar
        resid=res2
        sumgamma=sumgamma+1
      }
    }
  }
    if (iter > burnin)
  {
    alphasum<-alphasum+alpha
    alpha2sum<-alpha2sum+alpha^2
    alpha3sum<-alpha3sum+alpha^3
    betasum<-betasum+beta
    beta2sum<-beta2sum+beta^2
    beta3sum<-beta3sum+beta^3
  }
  
  if (iter%%1==0){
  betaimg<-array(beta,c(50,50,5))
  for (i in 2:5)
    bioimagetools::img(betaimg[,,i])
  
  gammaimg<-array(gamma,c(50,50,5))
  for (i in 2:5)
    bioimagetools::img(gammaimg[,,i])
  }
}
