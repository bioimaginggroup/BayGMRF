remove(list=ls())
library(adm)
data("germanySim")

gamma.sample<-sim1a$gamma
gamma1a<-apply(gamma.sample,2,median)
gamma.sample<-sim1b$gamma
gamma1b<-apply(gamma.sample,2,median)
gamma.sample<-sim1c$gamma
gamma1c<-apply(gamma.sample,2,median)

gamma.sample<-sim1i$gamma
gamma1i<-apply(gamma.sample,2,median)

data("germanymap")
par(mfrow=c(2,2))
cutpoints<-seq(0,.9,length=128)
map(abs(truegamma1-gamma1a),germany,cutpoints=cutpoints)
title("a")
map(abs(truegamma1-gamma1b),germany,cutpoints=cutpoints)
title("b")
map(abs(truegamma1-gamma1c),germany,cutpoints=cutpoints)
title("c")
map(abs(truegamma1-gamma1i),germany,cutpoints=cutpoints)
title("gmrf")

w1a<-apply(sim1a$w,2,median)
w1b<-apply(sim1b$w,2,median)
w1c<-apply(sim1c$w,2,median)

gamma.sample<-sim0i$gamma
gamma0i<-apply(gamma.sample,2,median)
gamma.sample<-sim0a$gamma
gamma0a<-apply(gamma.sample,2,median)

data("germanymap")
par(mfrow=c(1,2))
map(truegamma0,germany)
map(gamma0i,germany)
map(gamma0a,germany)

w0a<-apply(sim0a$w,2,median)

gamma.sample<-sim2a$gamma
gamma2a<-apply(gamma.sample,2,median)
gamma.sample<-sim2i$gamma
gamma2i<-apply(gamma.sample,2,median)

truegamma2<-truegamma2-mean(truegamma2)
gamma2i<-gamma2i-mean(gamma2i)
gamma2a<-gamma2a-mean(gamma2a)

data("germanymap")
par(mfrow=c(1,2))
map(truegamma2,germany,cutpoints=seq(-1,1,length=64))
map(gamma2i,germany,cutpoints=seq(-1,1,length=64))
map(gamma2a,germany,cutpoints=seq(-1,1,length=64))

map(abs(truegamma2-gamma2i),germany,cutpoints=cutpoints)
title("gmrf")
map(abs(truegamma2-gamma2a),germany,cutpoints=cutpoints)
title("adaptive")

w2a<-apply(sim2a$w,2,median)

gamma.sample<-sim3i$gamma
gamma3i<-apply(gamma.sample,2,median)
gamma.sample<-sim3a$gamma
gamma3a<-apply(gamma.sample,2,median)

data("germanymap")
par(mfrow=c(1,2))
map(truegamma3,germany)
map(gamma3i,germany)
map(gamma3a,germany)

cutpoints<-seq(0,1,length=128)
map(abs(truegamma3-gamma3i),germany,cutpoints=cutpoints)
title("gmrf")
map(abs(truegamma3-gamma3a),germany,cutpoints=cutpoints)
title("adaptive")

par(mfrow=c(1,1))
plot(abs(truegamma3-gamma3i),abs(truegamma3-gamma3a))
lines(c(-1,2),c(-1,2))
