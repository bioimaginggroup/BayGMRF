# rstgamma<-function(a,b,shift,trunc=0,trunc2=1)
# {
#   x<-trunc-1
#   counter<-0
#   while((x<trunc)|(x>trunc2)|(x==Inf))
#   {
#     x=rgamma(1,a,b)+shift
#     counter=counter+1
#     if (counter>3)if(x<trunc)x <- rtrunc(1,shape=a,rate=b,shift=shift,a=trunc,spec="sgamma")
#     if(counter>10)return(1)
#     if(counter>10)cat(".")
#     
#   }
#   return(x)
# }

rsgamma<-function(n,shape,rate,shift)
{
return(rgamma(n,shape,rate)+shift)
}

  
rstgamma<-function(n,shape,rate,shift,a,b)
{
  r<-rtrunc(n=n,spec="sgamma",a=a,b=b,shape=shape,rate=rate,shift=shift)
  if(r==-Inf)r=a
  if(r==Inf)r=b
  return(r)
}
