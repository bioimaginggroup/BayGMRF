data.germany <- read.table("red_Germany.dat", skip=1)
names(data.germany) <- c("Region", "Expected", "Observed", "Kov.Smoke")
data.germany$Region <- data.germany$Region + 1
x.germany <- data.germany$Kov.Smoke
y.germany <- data.germany$Observed
e.germany <- data.germany$Expected
d<-length(y.germany)


size <- 1100 # Samplingumfang
thin <- 1 # jeder thin-te Wert  
burn <- 100 # erst ab burn-ten Wert

coords<-read.table("d_coords.txt")
# for (i in 1:(d-1))
#   for (j in i:d)
#     if (Q.germanyknapp[i,j]==-1)
#       coords<-rbind(coords,c(i,j))
coords<-as.matrix(coords)

setwd("~/Dropbox/projects/approxGibbs")
#Rprof()
options(mc.cores=24)
t.gibbs3<-system.time({  # Laufzeit messen
  adm(size, y.germany, e.germany, coords, 1, burn, thin, adaptive=TRUE, method="gibbs3", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,.001), tau=1000, path="Gibbs3",print.time=FALSE)
})
print(t.gibbs3)
#Rprof()
#summaryRprof(tmp)

