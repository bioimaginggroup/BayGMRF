###########################################################################################
###########################################################################################
###########################################################################################
#################################### Konvergenzanalyse ####################################
###########################################################################################
###########################################################################################
###########################################################################################



##################
## Lade Package ##
##################

library(coda)



#######################
## Laden der Samples ##
#######################

# Test für 25.000 Iterationen bei Deutschland
# burn <- 3000
# limit <- 1250
# limit <- 2500
limit <- 10001
burn <- 1

# weights

load("Samples/1. Lauf/weight_sr.sample.RData")
weight_sr.sample1 <- t(weight_sr.sample[,burn:limit])

load("Samples/2. Lauf/weight_sr.sample.RData")
weight_sr.sample2 <- t(weight_sr.sample[,burn:limit])

load("Samples/3. Lauf/weight_sr.sample.RData")
weight_sr.sample3 <- t(weight_sr.sample[,burn:limit])

load("Samples/4. Lauf/weight_sr.sample.RData")
weight_sr.sample4 <- t(weight_sr.sample[,burn:limit])

load("Samples/5. Lauf/weight_sr.sample.RData")
weight_sr.sample5 <- t(weight_sr.sample[,burn:limit])


# Gamma

load("Samples/1. Lauf/Ga.sample.RData")
Ga.sample1 <- t(Ga.sample[,burn:limit])

load("Samples/2. Lauf/Ga.sample.RData")
Ga.sample2 <- t(Ga.sample[,burn:limit])

load("Samples/3. Lauf/Ga.sample.RData")
Ga.sample3 <- t(Ga.sample[,burn:limit])

load("Samples/4. Lauf/Ga.sample.RData")
Ga.sample4 <- t(Ga.sample[,burn:limit])

load("Samples/5. Lauf/Ga.sample.RData")
Ga.sample5 <- t(Ga.sample[,burn:limit])


# alpha

load("Samples/1. Lauf/alpha.sample.RData")
alpha.sample1 <- t(alpha.sample[,burn:limit])

load("Samples/2. Lauf/alpha.sample.RData")
alpha.sample2 <- t(alpha.sample[,burn:limit])

load("Samples/3. Lauf/alpha.sample.RData")
alpha.sample3 <- t(alpha.sample[,burn:limit])

load("Samples/4. Lauf/alpha.sample.RData")
alpha.sample4 <- t(alpha.sample[,burn:limit])

load("Samples/5. Lauf/alpha.sample.RData")
alpha.sample5 <- t(alpha.sample[,burn:limit])



# weights Simulation 

load("Samples/weight_sr.sample.setting1.RData")
weight_sr.sample1 <- t(weight_sr.sample[,burn:limit])

load("Samples/weight_sr.sample.setting2.RData")
weight_sr.sample2 <- t(weight_sr.sample[,burn:limit])

load("Samples/weight_sr.sample.setting3.RData")
weight_sr.sample3 <- t(weight_sr.sample[,burn:limit])

load("Samples/weight_sr.sample.setting4.RData")
weight_sr.sample4 <- t(weight_sr.sample[,burn:limit])




###############################
## Umwandeln in MCMC Objekte ##
###############################


mcmc.weight_sr.sample1 <- mcmc(weight_sr.sample1)
mcmc.weight_sr.sample2 <- mcmc(weight_sr.sample2)
mcmc.weight_sr.sample3 <- mcmc(weight_sr.sample3)
mcmc.weight_sr.sample4 <- mcmc(weight_sr.sample4)
mcmc.weight_sr.sample5 <- mcmc(weight_sr.sample5)

mcmc.Ga.sample1 <- mcmc(Ga.sample1)
mcmc.Ga.sample2 <- mcmc(Ga.sample2)
mcmc.Ga.sample3 <- mcmc(Ga.sample3)
mcmc.Ga.sample4 <- mcmc(Ga.sample4)
mcmc.Ga.sample5 <- mcmc(Ga.sample5)

mcmc.alpha.sample1 <- mcmc(alpha.sample1)
mcmc.alpha.sample2 <- mcmc(alpha.sample2)
mcmc.alpha.sample3 <- mcmc(alpha.sample3)
mcmc.alpha.sample4 <- mcmc(alpha.sample4)
mcmc.alpha.sample5 <- mcmc(alpha.sample5)



#####################
## Zusammenfassung ##
#####################

#summary(mcmc.weight_sr.sample1)
#summary(mcmc.weight_sr.sample2)
#summary(mcmc.weight_sr.sample3)
#summary(mcmc.weight_sr.sample4)
#summary(mcmc.weight_sr.sample5)


weight_sr.list <- mcmc.list(list(mcmc.weight_sr.sample1, mcmc.weight_sr.sample2, mcmc.weight_sr.sample3, mcmc.weight_sr.sample4, mcmc.weight_sr.sample5))

Ga.list <- mcmc.list(list(mcmc.Ga.sample1, mcmc.Ga.sample2, mcmc.Ga.sample3, mcmc.Ga.sample4, mcmc.Ga.sample5))

alpha.list <- mcmc.list(list(mcmc.alpha.sample1, mcmc.alpha.sample2, mcmc.alpha.sample3, mcmc.alpha.sample4, mcmc.alpha.sample5))



##########################################
## Verwende Gelman and Rubin Diagnostik ##
##########################################

gelman.diag(weight_sr.list)

pdf("Konvergenzanalyse über Potential Scale Reduction Factor - Gewichte.pdf")
gelman.plot(weight_sr.list, xlab="Ausgedünnte Iterationen", ylab="R")
dev.off()


gelman.diag(Ga.list)

pdf("Konvergenzanalyse über Potential Scale Reduction Factor - strukturierte Effekte.pdf")
gelman.plot(Ga.list, xlab="Ausgedünnte Iterationen", ylab="R")
dev.off()


gelman.diag(alpha.list)

pdf("Konvergenzanalyse über Potential Scale Reduction Factor - unstrukturierte Effekte.pdf")
gelman.plot(alpha.list, xlab="Ausgedünnte Iterationen", ylab="R")
dev.off()
