n=10^4
Q<-makeQ(n,"rw1+I",.5)
T<-lanczos(Q,n)
plot(eigen(T)$values,type="l")
lines(eigen(Q)$values,col="blue",lty=3)
