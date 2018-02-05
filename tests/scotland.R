options(mc.cores=2)
system.time({scoti<-adm(scotland.data$Counts, scotland.data$E, scotland.coords, nr.it=1500, burnin=500, thin=1, nu=1, adaptive=FALSE, do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)})
system.time({scot<-adm(scotland.data$Counts, scotland.data$E, scotland.coords, nr.it=1500, burnin=500, thin=1, nu=1, adaptive=TRUE, method="gibbs3beta", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, print.time=FALSE)})


gammai<-apply(scoti$gamma,2,mean)
gamma<-apply(scot$gamma,2,mean)

par(mfrow=c(1,3))
map(gammai,scotland.shape)
title(main="gmrf")
map(gamma,scotland.shape)
title(main="agmrf")
map(scotland.data$Counts/scotland.data$E,scotland.shape)
title(main="rates")

par(mfrow=c(1,1))
map(gamma-gammai,scotland.shape)

