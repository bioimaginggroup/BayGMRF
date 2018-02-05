rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}