remove(list=ls())
library(adm)
data(germanymap)
data(Germanydata)
Y.sim<-ground.bl<-truegamma<-vector(length=10,mode="list")
for (i in 1:10)
{  
ground.bl[[i]] <-rnorm(16,1,.5)
bl <- 1+
  ifelse((1:544)<16,1,0)+
  ifelse((1:544)<17,1,0)+
  ifelse((1:544)<64,1,0)+
  ifelse((1:544)<66,1,0)+
  ifelse((1:544)<120,1,0)+
  ifelse((1:544)<146,1,0)+
  ifelse((1:544)<182,1,0)+
  ifelse((1:544)<226,1,0)+
  ifelse((1:544)<322,1,0)+
  ifelse((1:544)<328,1,0)+
  ifelse((1:544)<330,1,0)+
  ifelse((1:544)<374,1,0)+
  ifelse((1:544)<411,1,0)+
  ifelse((1:544)<465,1,0)+
  ifelse((1:544)<505,1,0)
#map(bl,germany)
ground<-ground.bl[[i]][bl]
#map(ground,germany)
truegamma[[i]]<-rnorm(544,ground,.01)
#map(truegamma,germany)
truelambda<-Germany$E*exp(truegamma[[i]])
Y.sim[[i]]<-rpois(544,truelambda)
#map(log(Y.sim/Germany$E),germany)
}
save.image(file="~/Dropbox/projects/approxGibbs/adm/data/germanySimBl.RData")
