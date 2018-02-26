qsgamma <- function(p, shape, rate = 1, shift = 0, scale = 1/rate, lower.tail = TRUE)
{
  return(qgamma(p, shape, rate=1/scale, lower.tail = lower.tail, log.p = FALSE) + shift)
}