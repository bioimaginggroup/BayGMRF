remove(list=ls())
library(adm)
data(germanymap)
data(Germanydata)
data(germanycoords)

ground1<-1+.5*ifelse((1:544)<328,1,0)
truegamma1<-rnorm(544,ground1,.1)
truegamma1<-1+truegamma1-mean(truegamma1)
map(truegamma1,germany)

truelambda1<-Germany$E*exp(truegamma1)
Y.sim1<-rpois(544,truelambda1)

ground2<-1+.5*ifelse((1:544)<328,1,0)+.5*ifelse((1:544)<120,1,0)
truegamma2<-rnorm(544,ground2,.1)
truegamma2<-1+truegamma2-mean(truegamma2)
map(truegamma2,germany)

truelambda2<-Germany$E*exp(truegamma2)
Y.sim2<-rpois(544,truelambda2)

ground0<-0.5
truegamma0<-rnorm(544,ground0,.1)
map(truegamma0,germany)

truelambda0<-Germany$E*exp(truegamma0)
Y.sim0<-rpois(544,truelambda0)

ground3<-1+.5*(ifelse((1:544)<328,1,0)+ifelse((1:544)<120,1,0)+ifelse((1:544)<226,1,0))
map(ground3,germany)
truegamma3<-rnorm(544,ground3,.1)
truegamma3<-1+truegamma3-mean(truegamma3)
map(truegamma3,germany)

truelambda3<-Germany$E*exp(truegamma3)
Y.sim3<-rpois(544,truelambda3)

save.image(file="~/Dropbox/projects/approxGibbs/adm/data/germanySim.RData")
