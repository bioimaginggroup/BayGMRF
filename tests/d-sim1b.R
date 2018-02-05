library(adm)
data(Germanydata)
data(germanycoords)
data(germanySim)

options(mc.cores=24)
sim1b<-adm(Y.sim1, Germany$E, coords, nr.it=1500, burnin=500, thin=1, nu=c(1,1e-6), adaptive=TRUE, method="gibbs3", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,1), tau=1000, print.time=FALSE)
save.image(file="d-sim1-b.Rdata")