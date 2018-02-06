rm(.Random.seed)
library(adm)
data(scotland)
#par(mfrow=c(1,3))
true<-rep(0,52)
true[18:40]<-1
true[23:25]<-true[29:32]<-true[35:37]<-2
#map(true,scotland.shape)

scotland.sim<-scotland.data[,-3]
scotland.sim$E<-rpois(52,scotland.data$E-1)+1
globalgamma<-sample(c(-1,0,1))
gamma.true<-rnorm(52,globalgamma[true+1],.01)
#map(gamma.true,scotland.shape)

scotland.sim$Counts<-rpois(52,scotland.sim$E*exp(gamma.true))
#map(scotland.sim$Counts,scotland.shape)
#map(scotland.sim$Counts/scotland.sim$E,scotland.shape)

options(mc.cores=4)
scoti<-dm(scotland.sim$Counts, scotland.sim$E, scotland.coords, nr.it=20000, burnin=10000, thin=5, nu=1, do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)


gammai<-apply(scoti$gamma,2,mean,na.rm=TRUE)
#cuts<-seq(-1.7,1.7,length=32)
#map(gammai,scotland.shape,cutpoints=cuts)
#map(gamma.true,scotland.shape,cutpoints=cuts)
scota<-adm(scotland.sim$Counts, scotland.sim$E, scotland.coords, nr.it=20000, burnin=10000, thin=5, nu=1, adaptive=TRUE, method="gibbs2beta", do.kappa=FALSE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)
gammaa<-apply(scota$gamma,2,mean,na.rm=TRUE)
#map(gammaa,scotland.shape,cutpoints=cuts)

#par(mfrow=c(1,1))
#plot(gamma.true,gammai)
#points(gamma.true,gammaa,col="blue")
error.i<-(gamma.true-gammai)^2
error.a<-(gamma.true-gammaa)^2

#plot(gammaa,gammai)
#lines(c(-10,10),c(-10,10))

code=paste0(sample(c(LETTERS,letters,0:9),6,rep=TRUE),collapse="")
write(c(sum(error.i),sum(error.a)),file="scot_error.txt",append=TRUE)
save.image(file=paste0("scot/",code,".Rdata"))