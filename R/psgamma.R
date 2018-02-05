psgamma <- function(q, shape, rate = 1, shift = 0, scale = 1/rate, lower.tail = TRUE)
{
  return(pgamma(q-shift, shape, rate=1/scale, lower.tail = lower.tail, log.p = FALSE))
}