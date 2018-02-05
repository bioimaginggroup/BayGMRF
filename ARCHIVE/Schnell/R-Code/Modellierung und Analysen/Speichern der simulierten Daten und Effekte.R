redscotQ <- as.matrix(read.table("redscotQ.txt"))
colnames(redscotQ) <- c()

# Entferne Inseln
dim.scotQ <- dim(redscotQ)[1]
redscotQ <- redscotQ[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ), -c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)] 
  
  
#######################################
## Abspeichern der simulierten Daten ##
#######################################

# Setting 1

y.sim.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sim.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sim.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
  
for(i in 1:10){
  load(paste("Simulierte Daten/y.sim", i,".setting1.RData", sep=""))
  load(paste("Simulierte Daten/Gamma.sim", i,".setting1.RData", sep=""))
  load(paste("Simulierte Daten/alpha.sim", i,".setting1.RData", sep=""))
    
  y.sim.setting1.mx[,i] <- y.sim
  Gamma.sim.setting1.mx[,i] <- Gamma.standardisiert
  alpha.sim.setting1.mx[,i] <- alpha.standardisiert
} 

sum.Gamma.alpha.sim.setting1.mx <- Gamma.sim.setting1.mx + alpha.sim.setting1.mx

save(y.sim.setting1.mx, file="Simulierte Daten/y.sim.setting1.mx.RData")
save(Gamma.sim.setting1.mx, file="Simulierte Daten/Gamma.sim.setting1.mx.RData")
save(alpha.sim.setting1.mx, file="Simulierte Daten/alpha.sim.setting1.mx.RData")
save(sum.Gamma.alpha.sim.setting1.mx, file="Simulierte Daten/sum.Gamma.alpha.sim.setting1.mx.RData")


# Setting 2

y.sim.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sim.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sim.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])

for(i in 1:10){
  load(paste("Simulierte Daten/y.sim", i,".setting2.RData", sep=""))
  load(paste("Simulierte Daten/Gamma.sim", i,".setting2.RData", sep=""))
  load(paste("Simulierte Daten/alpha.sim", i,".setting2.RData", sep=""))
  
  y.sim.setting2.mx[,i] <- y.sim
  Gamma.sim.setting2.mx[,i] <- Gamma.standardisiert
  alpha.sim.setting2.mx[,i] <- alpha.standardisiert
} 

sum.Gamma.alpha.sim.setting2.mx <- Gamma.sim.setting2.mx + alpha.sim.setting2.mx

save(y.sim.setting2.mx, file="Simulierte Daten/y.sim.setting2.mx.RData")
save(Gamma.sim.setting2.mx, file="Simulierte Daten/Gamma.sim.setting2.mx.RData")
save(alpha.sim.setting2.mx, file="Simulierte Daten/alpha.sim.setting2.mx.RData")
save(sum.Gamma.alpha.sim.setting2.mx, file="Simulierte Daten/sum.Gamma.alpha.sim.setting2.mx.RData")


# Setting 3

y.sim.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sim.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sim.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])

for(i in 1:10){
  load(paste("Simulierte Daten/y.sim", i,".setting3.RData", sep=""))
  load(paste("Simulierte Daten/Gamma.sim", i,".setting3.RData", sep=""))
  load(paste("Simulierte Daten/alpha.sim", i,".setting3.RData", sep=""))
  
  y.sim.setting3.mx[,i] <- y.sim
  Gamma.sim.setting3.mx[,i] <- Gamma.standardisiert
  alpha.sim.setting3.mx[,i] <- alpha.standardisiert
} 

sum.Gamma.alpha.sim.setting3.mx <- Gamma.sim.setting3.mx + alpha.sim.setting3.mx

save(y.sim.setting3.mx, file="Simulierte Daten/y.sim.setting3.mx.RData")
save(Gamma.sim.setting3.mx, file="Simulierte Daten/Gamma.sim.setting3.mx.RData")
save(alpha.sim.setting3.mx, file="Simulierte Daten/alpha.sim.setting3.mx.RData")
save(sum.Gamma.alpha.sim.setting3.mx, file="Simulierte Daten/sum.Gamma.alpha.sim.setting3.mx.RData")


# Setting 4

y.sim.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sim.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sim.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])

for(i in 1:10){
  load(paste("Simulierte Daten/y.sim", i,".setting4.RData", sep=""))
  load(paste("Simulierte Daten/Gamma.sim", i,".setting4.RData", sep=""))
  load(paste("Simulierte Daten/alpha.sim", i,".setting4.RData", sep=""))
  
  y.sim.setting4.mx[,i] <- y.sim
  Gamma.sim.setting4.mx[,i] <- Gamma.standardisiert
  alpha.sim.setting4.mx[,i] <- alpha.standardisiert
} 

sum.Gamma.alpha.sim.setting4.mx <- Gamma.sim.setting4.mx + alpha.sim.setting4.mx

save(y.sim.setting4.mx, file="Simulierte Daten/y.sim.setting4.mx.RData")
save(Gamma.sim.setting4.mx, file="Simulierte Daten/Gamma.sim.setting4.mx.RData")
save(alpha.sim.setting4.mx, file="Simulierte Daten/alpha.sim.setting4.mx.RData")
save(sum.Gamma.alpha.sim.setting4.mx, file="Simulierte Daten/sum.Gamma.alpha.sim.setting4.mx.RData")



########################################################################
## Abspeichern der simulierten Verteilungscharakteristika nach Sample ##
########################################################################

# Setting 1

eta.sample.mean.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sample.mean.zent.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sample.mean.zent.setting1.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
weight_sr.sample.mean.setting1.mx <- matrix(data=NA, ncol=10, nrow=117)


for(i in 1:10){

  load(paste("Simulierte Samples/Setting1/y.sim", i,"/Verteilungscharakteristika/eta.sample.mean.RData", sep=""))  
  load(paste("Simulierte Samples/Setting1/y.sim", i,"/Verteilungscharakteristika/Ga.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting1/y.sim", i,"/Verteilungscharakteristika/alpha.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting1/y.sim", i,"/Verteilungscharakteristika/weight_sr.sample.mean.RData", sep=""))  
  
  
  eta.sample.mean.setting1.mx[,i] <- eta.sample.mean
  Gamma.sample.mean.zent.setting1.mx[,i] <- Ga.sample.mean.zent
  alpha.sample.mean.zent.setting1.mx[,i] <- alpha.sample.mean.zent
  weight_sr.sample.mean.setting1.mx[,i] <- weight_sr.sample.mean
} 

sum.Gamma.alpha.sample.mean.setting1.mx <- Gamma.sample.mean.zent.setting1.mx + alpha.sample.mean.zent.setting1.mx


save(eta.sample.mean.setting1.mx, file="Simulierte Samples/eta.sample.mean.setting1.mx.RData")
save(Gamma.sample.mean.zent.setting1.mx, file="Simulierte Samples/Gamma.sample.mean.zent.setting1.mx.RData")
save(alpha.sample.mean.zent.setting1.mx, file="Simulierte Samples/alpha.sample.mean.zent.setting1.mx.RData")
save(sum.Gamma.alpha.sample.mean.setting1.mx, file="Simulierte Samples/sum.Gamma.alpha.sample.mean.setting1.mx.RData")
save(weight_sr.sample.mean.setting1.mx, file="Simulierte Samples/weight_sr.sample.mean.setting1.mx.RData")



# Setting 2

eta.sample.mean.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sample.mean.zent.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sample.mean.zent.setting2.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
weight_sr.sample.mean.setting2.mx <- matrix(data=NA, ncol=10, nrow=117)


for(i in 1:10){
  
  load(paste("Simulierte Samples/Setting2/y.sim", i,"/Verteilungscharakteristika/eta.sample.mean.RData", sep=""))  
  load(paste("Simulierte Samples/Setting2/y.sim", i,"/Verteilungscharakteristika/Ga.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting2/y.sim", i,"/Verteilungscharakteristika/alpha.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting2/y.sim", i,"/Verteilungscharakteristika/weight_sr.sample.mean.RData", sep=""))  
  
  
  eta.sample.mean.setting2.mx[,i] <- eta.sample.mean
  Gamma.sample.mean.zent.setting2.mx[,i] <- Ga.sample.mean.zent
  alpha.sample.mean.zent.setting2.mx[,i] <- alpha.sample.mean.zent
  weight_sr.sample.mean.setting2.mx[,i] <- weight_sr.sample.mean
} 

sum.Gamma.alpha.sample.mean.setting2.mx <- Gamma.sample.mean.zent.setting2.mx + alpha.sample.mean.zent.setting2.mx


save(eta.sample.mean.setting2.mx, file="Simulierte Samples/eta.sample.mean.setting2.mx.RData")
save(Gamma.sample.mean.zent.setting2.mx, file="Simulierte Samples/Gamma.sample.mean.zent.setting2.mx.RData")
save(alpha.sample.mean.zent.setting2.mx, file="Simulierte Samples/alpha.sample.mean.zent.setting2.mx.RData")
save(sum.Gamma.alpha.sample.mean.setting2.mx, file="Simulierte Samples/sum.Gamma.alpha.sample.mean.setting2.mx.RData")
save(weight_sr.sample.mean.setting2.mx, file="Simulierte Samples/weight_sr.sample.mean.setting2.mx.RData")



# Setting 3

eta.sample.mean.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sample.mean.zent.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sample.mean.zent.setting3.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
weight_sr.sample.mean.setting3.mx <- matrix(data=NA, ncol=10, nrow=117)


for(i in 1:10){
  
  load(paste("Simulierte Samples/Setting3/y.sim", i,"/Verteilungscharakteristika/eta.sample.mean.RData", sep=""))  
  load(paste("Simulierte Samples/Setting3/y.sim", i,"/Verteilungscharakteristika/Ga.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting3/y.sim", i,"/Verteilungscharakteristika/alpha.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting3/y.sim", i,"/Verteilungscharakteristika/weight_sr.sample.mean.RData", sep=""))  
  
  
  eta.sample.mean.setting3.mx[,i] <- eta.sample.mean
  Gamma.sample.mean.zent.setting3.mx[,i] <- Ga.sample.mean.zent
  alpha.sample.mean.zent.setting3.mx[,i] <- alpha.sample.mean.zent
  weight_sr.sample.mean.setting3.mx[,i] <- weight_sr.sample.mean
} 

sum.Gamma.alpha.sample.mean.setting3.mx <- Gamma.sample.mean.zent.setting3.mx + alpha.sample.mean.zent.setting3.mx


save(eta.sample.mean.setting3.mx, file="Simulierte Samples/eta.sample.mean.setting3.mx.RData")
save(Gamma.sample.mean.zent.setting3.mx, file="Simulierte Samples/Gamma.sample.mean.zent.setting3.mx.RData")
save(alpha.sample.mean.zent.setting3.mx, file="Simulierte Samples/alpha.sample.mean.zent.setting3.mx.RData")
save(sum.Gamma.alpha.sample.mean.setting3.mx, file="Simulierte Samples/sum.Gamma.alpha.sample.mean.setting3.mx.RData")
save(weight_sr.sample.mean.setting3.mx, file="Simulierte Samples/weight_sr.sample.mean.setting3.mx.RData")



# Setting 4

eta.sample.mean.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
Gamma.sample.mean.zent.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
alpha.sample.mean.zent.setting4.mx <- matrix(data=NA, ncol=10, nrow=dim(redscotQ)[1])
weight_sr.sample.mean.setting4.mx <- matrix(data=NA, ncol=10, nrow=117)


for(i in 1:10){
  
  load(paste("Simulierte Samples/Setting4/y.sim", i,"/Verteilungscharakteristika/eta.sample.mean.RData", sep=""))  
  load(paste("Simulierte Samples/Setting4/y.sim", i,"/Verteilungscharakteristika/Ga.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting4/y.sim", i,"/Verteilungscharakteristika/alpha.sample.mean.zent.RData", sep=""))
  load(paste("Simulierte Samples/Setting4/y.sim", i,"/Verteilungscharakteristika/weight_sr.sample.mean.RData", sep=""))  
  
  
  eta.sample.mean.setting4.mx[,i] <- eta.sample.mean
  Gamma.sample.mean.zent.setting4.mx[,i] <- Ga.sample.mean.zent
  alpha.sample.mean.zent.setting4.mx[,i] <- alpha.sample.mean.zent
  weight_sr.sample.mean.setting4.mx[,i] <- weight_sr.sample.mean
} 

sum.Gamma.alpha.sample.mean.setting4.mx <- Gamma.sample.mean.zent.setting4.mx + alpha.sample.mean.zent.setting4.mx


save(eta.sample.mean.setting4.mx, file="Simulierte Samples/eta.sample.mean.setting4.mx.RData")
save(Gamma.sample.mean.zent.setting4.mx, file="Simulierte Samples/Gamma.sample.mean.zent.setting4.mx.RData")
save(alpha.sample.mean.zent.setting4.mx, file="Simulierte Samples/alpha.sample.mean.zent.setting4.mx.RData")
save(sum.Gamma.alpha.sample.mean.setting4.mx, file="Simulierte Samples/sum.Gamma.alpha.sample.mean.setting4.mx.RData")
save(weight_sr.sample.mean.setting4.mx, file="Simulierte Samples/weight_sr.sample.mean.setting4.mx.RData")

