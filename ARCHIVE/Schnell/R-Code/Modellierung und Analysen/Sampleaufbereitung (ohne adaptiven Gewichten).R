# Laden der Samples
load("Samples/mu.sample.RData")
load("Samples/beta.sample.RData")
load("Samples/Ga.sample.RData")
load("Samples/kappa.sample.RData")
load("Samples/eta.sample.RData")
load("Samples/tau.sample.RData")
load("Samples/alpha.sample.RData")
load("Samples/alpha.eta.proz.sample.RData")
load("Samples/eta.zufall.sample.RData")


# size <- length(mu.sample)
# burn <- 5000



#################################
## Grafische Veranschaulichung ##
#################################

plot(as.vector(eta.zufall.sample[48,-(1:burn)]))

plot(alpha.eta.proz.sample[-(1:burn)], type="l", ylim=c(0,1), main="Akzeptanzwahrscheinlichkeit eta", ylab="Wahrscheinlichkeit", xlab="Iteration")

plot(mu.sample[burn:size], type="l", xlab="Iteration", ylab="mu", main="mu Samplingpfad")
plot(beta.sample[burn:size], type="l", xlab="Iteration", ylab="beta", main="beta Samplingpfad")
plot(density(beta.sample[burn:size]), xlab="beta", ylab="Dichte", main="beta Sampler", lwd=4)
acf(beta.sample[burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="beta Autokorrelation")


plot(kappa.sample[burn:size], type="l", xlab="Iteration", ylab="kappa", main="kappa Samplingpfad")
plot(density(kappa.sample[burn:size]), xlab="kappa", ylab="Dichte", main="kappa Sampler", lwd=4)
acf(kappa.sample[burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="kappa Autokorrelation")


plot(tau.sample[burn:size], type="l", xlab="Iteration", ylab="tau", main="tau Samplingpfad")
plot(density(tau.sample[burn:size]), xlab="tau", ylab="Dichte", main="tau Sampler", lwd=4)
acf(tau.sample[burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="tau Autokorrelation")


plot(Ga.sample[3,burn:size], type="l", xlab="Iteration", ylab="gamma", main="gamma[1] Samplingpfad")
plot(density(Ga.sample[1,burn:size]), xlab="gamma", ylab="Dichte", main="gamma Sampler", lwd=4)
acf(Ga.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="gamma[1] Autokorrelation")


plot(alpha.sample[2,burn:size], type="l", xlab="Iteration", ylab="alpha", main="alpha[1] Samplingpfad")
plot(density(alpha.sample[1,burn:size]), xlab="alpha", ylab="Dichte", main="alpha Sampler", lwd=4)
acf(alpha.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="alpha[1] Autokorrelation")


plot(eta.sample[53,1:size], type="l", xlab="Iteration", ylab="eta", main="eta[1] Samplingpfad")
plot(density(eta.sample[1,burn:size]), xlab="eta", ylab="Dichte", main="eta Sampler", lwd=4)
acf(eta.sample[1,burn:size], lag.max=size, xlab="Iteration", ylab="empirische ACF", main="eta[1] Autokorrelation")



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

apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean) # Punktschätzer
range(apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean))

mean(y.scot)
mean(exp(apply(X=eta.sample[,(burn:size)], MARGIN=1, FUN=mean)))






################################################
## Abspeichern der Verteilungscharakteristika ##
################################################

#Berechne Mittelwerte
Ga.sample.mean <- apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean)
alpha.sample.mean <- apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean)
eta.sample.mean <- apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean)

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
Ga.sample.sign

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
alpha.sample.sign

save(Ga.sample.mean, file="Verteilungscharakteristika/Ga.sample.mean.RData")
save(Ga.sample.mean.zent, file="Verteilungscharakteristika/Ga.sample.mean.zent.RData")
save(Ga.sample.sign, file="Verteilungscharakteristika/Ga.sample.sign.RData")
save(alpha.sample.mean, file="Verteilungscharakteristika/alpha.sample.mean.RData")
save(alpha.sample.mean.zent, file="Verteilungscharakteristika/alpha.sample.mean.zent.RData")
save(alpha.sample.sign, file="Verteilungscharakteristika/alpha.sample.sign.RData")
save(eta.sample.mean, file="Verteilungscharakteristika/eta.sample.mean.RData")

#thin <- 10

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
