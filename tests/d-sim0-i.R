library(adm)
data(Germanydata)
data(germanycoords)
data(germanySim)

options(mc.cores=16)
sim0i<-adm(Y.sim0, Germany$E, coords, nr.it=1500, burnin=500, thin=1, nu=1, adaptive=FALSE, method="", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)
save.image(file="d-sim0-i.Rdata")
