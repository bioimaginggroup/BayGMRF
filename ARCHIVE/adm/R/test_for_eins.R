test.for.eins<-function(cc,coords)
  return((sum(as.vector(coords)==cc[1])==1)|(sum(as.vector(coords)==cc[2])==1))
