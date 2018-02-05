########################
## Einlesen der Daten ##
########################

# Schottland Daten laden
data.scotland <- read.table("red.Scotland.dat")
# Übergebe Variablennamen
names(data.scotland) <- c("District", "Observed", "Expected", "PcAFF", "Latitude", "Longitude")

# Laden der reduzierten schottischen Nachbarschaftsmatrix
redscotQ <- as.matrix(read.table("redscotQ.txt"))
colnames(redscotQ) <- c()
redscotQknapp <- knappeMatrix(redscotQ)
# Entferne Inseln
dim.scotQ <- dim(redscotQ)[1]
redscotQ <- redscotQ[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ), -c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)] 
redscotQknapp <- redscotQknapp[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ), -c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)]  


# Deutschland Daten laden
data.germany <- read.table("red.Germany.dat", skip=1)
# Übergebe Variablennamen
names(data.germany) <- c("Region", "Expected", "Observed", "Kov.Smoke")
#Zähler der Region um 1 rauf gesetzt, da Start bei 0
data.germany$Region <- data.germany$Region + 1

# Laden der reduzierten deutschen Nachbarschaftsmatrix
redgermanyQ <- as.matrix(read.table("redgermanyQ.txt"))
colnames(redgermanyQ) <- c()
redgermanyQknapp <- knappeMatrix(redgermanyQ)


# Aufbereitung für den Hybrid-Sampler

# Festlegung der Parameter des hybrid.sampler für SChottlanddaten ohne Inseln
x.scot <- data.scotland$PcAFF[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)]
y.scot <- data.scotland$Observed[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)]
e.scot <- data.scotland$Expected[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)]
Q.scot <- redscotQ
Q.scotknapp <- redscotQknapp

# Festlegung der Parameter des hybrid.sampler für Deutschlanddaten
x.germany <- data.germany$Kov.Smoke
y.germany <- data.germany$Observed
e.germany <- data.germany$Expected
Q.germany <- redgermanyQ
Q.germanyknapp <- redgermanyQknapp


# Laden der Samples
load("Samples/mu.sample.RData")
load("Samples/beta.sample.RData")
load("Samples/Ga.sample.RData")
load("Samples/kappa.sample.RData")
load("Samples/eta.sample.RData")
load("Samples/tau.sample.RData")
load("Samples/alpha.sample.RData")
load("Samples/weight_sr.sample.RData")
load("Samples/Q_weight.sample.aktuell.RData")
load("Samples/alpha.weight.proz.sample.RData")
load("Samples/alpha.eta.proz.sample.RData")
load("Samples/eta.zufall.sample.RData")

# size <- length(mu.sample)
# burn <- 1000
# thin <- 10



#################################
## Grafische Veranschaulichung ##
#################################


plot(as.vector(eta.zufall.sample[1,]))

plot(mu.sample[burn:size], type="l", xlab="Iteration", ylab="mu", main="mu Samplingpfad")
plot(alpha.eta.proz.sample[burn:size], type="l", ylim=c(0,1), main="Akzeptanzwahrscheinlichkeit eta", ylab="Wahrscheinlichkeit", xlab="Iteration")
plot(alpha.weight.proz.sample[-(1:burn)], type="l", ylim=c(0,1), main="Akzeptanzwahrscheinlichkeit weight", ylab="Wahrscheinlichkeit", xlab="Iteration")

plot(beta.sample[burn:size], type="l", xlab="Iteration", ylab="beta", main="beta Samplingpfad")
plot(kappa.sample[burn:size], type="l", xlab="Iteration", ylab="kappa", main="kappa Samplingpfad")
plot(tau.sample[burn:size], type="l", xlab="Iteration", ylab="tau", main="tau Samplingpfad")

plot(Ga.sample[53,burn:size], type="l", xlab="Iteration", ylab="gamma", main="gamma[1] Samplingpfad")
plot(density(Ga.sample[1,burn:size]), xlab="gamma", ylab="Dichte", main="gamma Sampler", lwd=4)
acf(Ga.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="gamma[1] Autokorrelation")

plot(alpha.sample[47,burn:size], type="l", xlab="Iteration", ylab="alpha", main="alpha[1] Samplingpfad")
plot(density(alpha.sample[1,burn:size]), xlab="alpha", ylab="Dichte", main="alpha Sampler", lwd=4)
acf(alpha.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="alpha[1] Autokorrelation")

plot(eta.sample[1,1:size], type="l", xlab="Iteration", ylab="eta", main="eta[1] Samplingpfad")
plot(density(eta.sample[1,burn:size]), xlab="eta", ylab="Dichte", main="eta Sampler", lwd=4)
acf(eta.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="eta[1] Autokorrelation")

plot(weight_sr.sample[2,burn:size], type="l", xlab="Iteration", ylab="weight", main="weigth[1] Samplingpfad")
plot(density(weight_sr.sample[1,burn:size]), xlab="weight", ylab="Dichte", main="weight Sampler", lwd=4)
acf(weight_sr.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="weigth[1] Autokorrelation")





# Bestimmung der Verteilungscharakteristika
mean(mu.sample)
median(mu.sample)

mean(beta.sample[burn:size])
median(beta.sample[burn:size])

mean(kappa.sample[burn:size])
median(kappa.sample[burn:size])

mean(tau.sample[burn:size])
median(tau.sample[burn:size])

apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean) # Punktschätzer
range(apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean)) # Punktschätzer

apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean) # Punktschätzer
range(apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean))

apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=mean) # Punktschätzer
range(apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=mean))

apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean) # Punktschätzer
range(apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean))

mean(y.germany)
mean(exp(apply(X=eta.sample[,(burn+1):(size)], MARGIN=1, FUN=mean)))



################################################
## Abspeichern der Verteilungscharakteristika ##
################################################

#Berechne Mittelwerte
Ga.sample.mean <- apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean)
alpha.sample.mean <- apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean)
eta.sample.mean <- apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean)
weight_sr.sample.mean <- apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=mean)

#Berechne Standardabweichung
Ga.sample.sd <- apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=sd)
alpha.sample.sd <- apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=sd)
eta.sample.sd <- apply(X=eta.sample[,burn:size], MARGIN=1, FUN=sd)
weight_sr.sample.sd <- apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=sd)

# Schätzung für Intercept
Ga.sample.mean.zent <- Ga.sample.mean - mean(Ga.sample.mean)
alpha.sample.mean.zent <- alpha.sample.mean - mean(alpha.sample.mean)
mu.sample.schaetz <- mean(Ga.sample.mean) + mean(alpha.sample.mean)

# Kredibilitätsintervall Gamma
Ga.sample.zent <- Ga.sample - mean(Ga.sample.mean) # Datenzentrierung
Ga.sample.025 <- apply(X=Ga.sample.zent[,burn:size], MARGIN=1, FUN=quantile,.025) # für Intervallschätzer
Ga.sample.075 <- apply(X=Ga.sample.zent[,burn:size], MARGIN=1, FUN=quantile,.975) # für Intervallschätzer

Ga.sample.sign <- vector(length = length(Ga.sample.mean))

for(i in 1:length(Ga.sample.mean)){
  if(Ga.sample.025[i] > 0 && Ga.sample.075[i] > 0){ # falls signifikant größer 0
    Ga.sample.sign[i] <- 1
  } else if(Ga.sample.025[i] < 0 && Ga.sample.075[i] < 0){ # falls signifikant kleiner 0
    Ga.sample.sign[i] <- -1
  } else{ # falls nicht signifikant größer anders 0
    Ga.sample.sign[i] <- 0
  }
}
#Ga.sample.sign

# Kredibilitätsintervall alpha
alpha.sample.zent <- alpha.sample - mean(alpha.sample.mean) # Datenzentrierung
alpha.sample.025 <- apply(X=alpha.sample.zent[,burn:size], MARGIN=1, FUN=quantile,.025) # für Intervallschätzer
alpha.sample.075 <- apply(X=alpha.sample.zent[,burn:size], MARGIN=1, FUN=quantile,.975) # für Intervallschätzer

alpha.sample.sign <- vector(length = length(alpha.sample.mean))

for(i in 1:length(alpha.sample.mean)){
  if(alpha.sample.025[i] > 0 && alpha.sample.075[i] > 0){
    alpha.sample.sign[i] <- 1
  } else if(alpha.sample.025[i] < 0 && alpha.sample.075[i] < 0){
    alpha.sample.sign[i] <- -1
  } else{
    alpha.sample.sign[i] <- 0
  }
}
#alpha.sample.sign


save(Ga.sample.mean, file="Verteilungscharakteristika/Ga.sample.mean.RData")
save(Ga.sample.sd, file="Verteilungscharakteristika/Ga.sample.sd.RData")
save(Ga.sample.mean.zent, file="Verteilungscharakteristika/Ga.sample.mean.zent.RData")
save(Ga.sample.sign, file="Verteilungscharakteristika/Ga.sample.sign.RData")
save(alpha.sample.mean, file="Verteilungscharakteristika/alpha.sample.mean.RData")
save(alpha.sample.sd, file="Verteilungscharakteristika/alpha.sample.sd.RData")
save(alpha.sample.mean.zent, file="Verteilungscharakteristika/alpha.sample.mean.zent.RData")
save(alpha.sample.sign, file="Verteilungscharakteristika/alpha.sample.sign.RData")
save(eta.sample.mean, file="Verteilungscharakteristika/eta.sample.mean.RData")
save(eta.sample.sd, file="Verteilungscharakteristika/eta.sample.sd.RData")
save(weight_sr.sample.mean, file="Verteilungscharakteristika/weight_sr.sample.mean.RData")
save(weight_sr.sample.sd, file="Verteilungscharakteristika/weight_sr.sample.sd.RData")

einzelparameter <- data.frame(data=c(burn*thin, mean(mu.sample[burn:size]), median(mu.sample[burn:size]), 
                                     mean(Ga.sample.mean), mean(alpha.sample.mean), mu.sample.schaetz,
                                     mean(beta.sample[burn:size]), median(beta.sample[burn:size]),
                                     mean(kappa.sample[burn:size]), median(kappa.sample[burn:size]),
                                     mean(tau.sample[burn:size]), median(tau.sample[burn:size]),
                                     mean(y.scot), mean(exp(apply(X=eta.sample[,(burn):(size)], MARGIN=1, FUN=mean))),
                                     range(apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean))[1],
                                     range(apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean))[2],
                                     range(apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean))[1],
                                     range(apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean))[2],
                                     range(apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean))[1],
                                     range(apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean))[2],
                                     alpha.eta.proz.sample[size]))
colnames(einzelparameter) <- "Einzelparameter:"
rownames(einzelparameter) <- c("Burn-In Phase", "mu.sample.mean","mu.sample.median", "Gamma.mu.anteil", "alpha.mu.anteil",
                               "Geschätztes mu", "beta.sample.mean", "beta.sample.median",
                               "kappa.sample.mean", "kappa.sample.median", "tau.sample.mean", "tau.sample.median",
                               "y.scot.mean", "exp.eta.sample.mean", "range.Ga.sample.mean.unten", 
                               "range.Ga.sample.mean.oben", "range.alpha.sample.mean.unten", "range.alpha.sample.mean.oben",
                               "range.eta.sample.mean.unten", "range.eta.sample.mean.oben", "Akzeptanzwahrscheinlichkeit eta")

einzelparameter <- round(einzelparameter, 4)

write.table(einzelparameter, "Verteilungscharakteristika/Einzelparameter.txt", row.names = TRUE, col.names = TRUE, quote=FALSE)
