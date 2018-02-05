gamma.i<-gamma.a<-array(NA,c(10,544))
for (i in 1:10)
{
  gamma.i[i,]<-apply(sim.alli[[i]]$gamma,2,median)
  gamma.a[i,]<-apply(sim.all[[i]]$gamma,2,median)
}

for (i in 1:10)
{
  print(c(sum((truegamma[[i]]-gamma.i[i,])^2),sum((truegamma[[i]]-gamma.a[i,])^2)),4)
}

par(mfrow=c(3,3))
for (i in 1:9)
  map(abs(truegamma[[i]]-gamma.i[i,]),germany)
for (i in 1:9)
  map(abs(truegamma[[i]]-gamma.a[i,]),germany)

error.i<-error.a<-array(NA,c(10,544))
for (i in 1:10)
{
  error.i[i,]<-truegamma[[i]]-gamma.i[i,]
  error.a[i,]<-truegamma[[i]]-gamma.a[i,]
}
error.a<-abs(error.a)
error.i<-abs(error.i)

par(mfrow=c(1,2))
map(apply(error.i,2,mean),germany,cutpoints=seq(0,.4,length=128))
title("global gmrf")
map(apply(error.a,2,mean),germany,cutpoints=seq(0,.4,length=128))
title("adaptive gmrf")

par(mfrow=c(1,1))
map(apply(error.i,2,mean)-apply(error.a,2,mean),germany,cutpoints=seq(-.1,.1,length=5))
   
