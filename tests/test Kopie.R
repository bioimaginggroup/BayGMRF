X<-20
Y<-15

mu<-array(rnorm(X*Y,0,1),c(X,Y))
c<-array(0,c(X,Y))
c[,floor(Y/2):Y]<-1
beta<-array(0,c(X,Y))
for (i in 1:Y)beta[,i]<-i/Y
x<-10+array(rpois(X*Y,10),c(X,Y))

data<-simulate(c(X,Y), mu, c, x, beta)

image(data)
image(beta)
image(x)
image(mu)
image(c)

