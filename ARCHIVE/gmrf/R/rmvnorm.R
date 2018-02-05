rmvnorm <- function(n, b, P){
  L <- chol(P) 
  u <- solve(t(L), b) 
  v <- solve(L, u) 
  z <- rnorm(n = length(b), mean = 0, sd = 1)
  m <- solve(L, z) 
  Gamma <- v + m
  return(as.vector(Gamma))
}
