library(adm)
data(Germanydata)
data(germanycoords)
data(germanySim)

options(mc.cores=16)
sim2i<-adm(Y.sim2, Germany$E, coords, nr.it=1500, burnin=500, thin=1, nu=1, adaptive=FALSE, do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)
save.image(file="d-sim2-i.Rdata")
