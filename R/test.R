X<-20
Y<-20

mu<-rnorm(c(X,Y),0,.1)
c<-array(0,c(X,Y))
c[,floor(Y/2):Y]<-1
beta<-array(0,c(X,Y))
for (i in 1:Y)beta[,i]<-i/Y
