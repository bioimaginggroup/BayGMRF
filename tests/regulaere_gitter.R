library(Matrix) #für spärlich besetzte Matrizen
N<-50

# Random Walk erster Ordnung
rw1<-function(n)
{
  i<-c(1:n,2:n)
  j<-c(1:n,(1:(n-1)))
  x<-c(1,rep(2,n-2),1,rep(-1,n-1))
  return(sparseMatrix(i=i,j=j,x=x,symmetric=TRUE))
}

R1<-rw1(N)
image(R1)
rankMatrix(R1)
determinant(R1)

I<-sparseMatrix(1:N,1:N,x=1)
rankMatrix(R1+I)
determinant(R1+I)
determinant(R1+0.00001*I)

plot(eigen(R1)$values)

#Random Walk zweiter Ordnung
rw2<-function(n)
{
  i<-c(1:n,2:n,3:n)
  j<-c(1:n,1:(n-1),1:(n-2))
  x<-c(1,5,rep(6,n-4),5,1,-2,rep(-4,n-3),-2,rep(1,n-2))
  return(sparseMatrix(i=i,j=j,x=x,symmetric=TRUE))
}

R2<-rw2(N)
image(R2)
rankMatrix(R2)

# 2D, erste Nachbarn
Q1 = kronecker(R1, I) + kronecker(I, R1) # FUN="+" funktioniert nicht mit sparseMatrix
image(Q1)
image(Q1[1001:1150,1001:1150])

rankMatrix(Q1, method="qr")

# 2D, zweite Nachbarn
Q2 = kronecker(R2, I) + kronecker(I, R2) 
image(Q2)
image(Q2[1001:1150,1001:1150])

rankMatrix(Q2, method="qr") #?
eigen(Q2)$values

# 2D, diagonale Nachbarn
Q = kronecker(R1, R1) 
image(Q)
image(Q[1001:1150,1001:1150])
print(Q[1001:1070,1001:1070])

print((N-1)*(N-1))
rankMatrix(Q, method="qr") 

# 2D, 12 nächste Nachbarn
Q12 = Q2 - Q
image(Q12)
image(Q12[1001:1150,1001:1150])
print(Q12[1001:1070,1001:1070])

rankMatrix(Q12, method="qr") #?
eigen(Q12)$values #?

## Rue-Block-Algorithmus
sample.rue<-function(m,Q)
{
  L = chol(Q)
  w = solve(L,m)
  u = solve(t(L),w)
  z <- rnorm(length(m))
  v = solve(t(L),z)
  return(u+v)
}

m <- seq(1,5,length=100)
Q <- rw2(100)+sparseMatrix(1:100,1:100,x=1)
L = chol(Q)
image(L)
w = solve(L,m)
plot(w)
u = solve(t(L),w)
plot(u)
z <- rnorm(length(m))
plot(z)
v = solve(t(L),z)
plot(v)
plot(u+v)

#Kovarianzmatrix
Sigma<-solve(Q)
print(Sigma)
image(Sigma)
#Erwartungswert
mu<-Sigma%*%m
plot(mu)

plot(sample.rue(1:50, rw2(50)+sparseMatrix(1:50,1:50,x=1)))
plot(sample.rue(1:50, rw2(50)+100*sparseMatrix(1:50,1:50,x=1)))
plot(sample.rue(1:50, 100*sparseMatrix(1:50,1:50,x=1)))
plot(sample.rue(1:50, 100*rw2(50)+100*sparseMatrix(1:50,1:50,x=1)))



