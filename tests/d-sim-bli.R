library(adm)
library(parallel)
data(Germanydata)
data(germanycoords)
data(germanySimBl)

options(mc.cores=32)

goforit<-function(i)
{
return(adm(Y.sim[[i]], Germany$E, coords, nr.it=1500, burnin=500, thin=1, nu=c(1,1), adaptive=FALSE, method="", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE))
}

sim.alli<-mclapply(1:10,goforit)
save.image(file="d-sim-bli.Rdata")


