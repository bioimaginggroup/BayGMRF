library(adm)
library(parallel)
data(Germanydata)
data(germanycoords)
data(germanySimBl)

options(mc.cores=16)

goforit<-function(i)
{
return(adm(Y.sim[[i]], Germany$E, coords, nr.it=4000, burnin=1000, thin=5, nu=1, adaptive=TRUE, method="gibbs3beta", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE))
}

sim.all<-mclapply(1:10,goforit, mc.cores=10)
save.image(file="d-sim-bl-gibbs.Rdata")


