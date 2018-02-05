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
load("Samples/weight_sr.sample.RData")
load("Samples/alpha.weight.proz.sample.RData")


size <- length(mu.sample)
burn <- 2000


# Berechne gleitenden Mittelwert
gleit.Mittelwert <- function(parameter.sample){
  rm <- 0
    for(i in 1:length(parameter.sample[(burn):size])){
    rm[i] <- (sum(parameter.sample[(burn):(burn+i-1)]))/i
  }
  return(rm)
}


#country <- "scot"
country <- "german"



#################################
## Grafische Veranschaulichung ##
#################################


# beta

pdf(paste("Plots/", country, "mG_Samplingpfad_beta.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(beta.sample[burn:size], xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
     ylab="", cex.lab=1.8, cex.axis=1.8, type="l")

dev.off()


pdf(paste("Plots/", country, "mG_Dichte_beta.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(density(beta.sample[burn:size]), xlab=expression(paste(beta)), 
     lwd=4, main="", ylab="", cex.lab=1.8, cex.axis=1.8)

dev.off()



pdf(paste("Plots/", country, "mG_ACF_beta.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

acf(beta.sample[burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
    cex.lab=1.8, cex.axis=1.8 , main="", ylab="")

dev.off()


pdf(paste("Plots/", country, "mG_MA_beta.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(gleit.Mittelwert(beta.sample), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
     ylab="", cex.lab=1.8, cex.axis=1.8, lwd=2)

dev.off()



# kappa

pdf(paste("Plots/", country, "mG_Samplingpfad_kappa.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(kappa.sample[burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
     ylab="", cex.lab=1.8, cex.axis=1.8)

dev.off()


pdf(paste("Plots/", country, "mG_Dichte_kappa.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(density(kappa.sample[burn:size]), xlab=expression(paste(kappa)), 
     lwd=4, cex.lab=1.8, cex.axis=1.8 , main="", ylab="")

dev.off()


pdf(paste("Plots/", country, "mG_ACF_kappa.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

acf(kappa.sample[burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
    cex.lab=1.8, cex.axis=1.8 , main="", ylab="")

dev.off()



pdf(paste("Plots/", country, "mG_MA_kappa.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(gleit.Mittelwert(kappa.sample), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
     ylab="", cex.lab=1.8, cex.axis=1.8, lwd=2)

dev.off()




# Strukturierter Effekt

apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean)


pdf(paste("Plots/", country, "mG_", length(Ga.sample[,size]), "Samplingpfade_gamma.pdf", sep=""), width=8.27, height=11.69)
par(mfrow=c(5,2))
for(i in 1:length(Ga.sample[,size])){
    plot(Ga.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
     cex.lab=1, cex.axis=1, main=paste("Samplingpfad strukturierter Effekt (Region s=", i, ")", sep=""))
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(Ga.sample[,size]), "Dichten_gamma.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(Ga.sample[,size])){
plot(density(Ga.sample[i,burn:size]), xlab=expression(paste(gamma, " der Region s")), 
     lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte strukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(Ga.sample[,size]), "ACF_gamma.pdf", sep=""), width=9, height=11.69) 
par(mfrow=c(5,3))
for(i in 1:length(Ga.sample[,size])){
acf(Ga.sample[i,burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
    cex.lab=1, cex.axis=1 , main=paste("ACF strukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
}
dev.off()



pdf(paste("Plots/", country, "mG_", length(Ga.sample[,size]), "MA_gamma.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(Ga.sample[,size])){
plot(gleit.Mittelwert(Ga.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
     ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. strukturierter Effekt (Region s=", i, ")", sep=""))
}
dev.off()





# Unstrukturierter Effekt

apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean)


pdf(paste("Plots/", country, "mG_", length(alpha.sample[,size]), "Samplingpfade_alpha.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(alpha.sample[,size])){
  plot(alpha.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
       cex.lab=1, cex.axis=1, main=paste("Samplingpfad unstrukturierter Effekt (Region s=", i, ")", sep=""))
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(alpha.sample[,size]), "Dichten_alpha.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(alpha.sample[,size])){
  plot(density(alpha.sample[i,burn:size]), xlab=expression(paste(alpha, " der Region s")), 
       lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte unstrukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(alpha.sample[,size]), "ACF_alpha.pdf", sep=""), width=9, height=11.69) 
par(mfrow=c(5,3))
for(i in 1:length(alpha.sample[,size])){
  acf(alpha.sample[i,burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
      cex.lab=1, cex.axis=1 , main=paste("ACF unstrukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
}
dev.off()



pdf(paste("Plots/", country, "mG_", length(alpha.sample[,size]), "MA_alpha.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(alpha.sample[,size])){
  plot(gleit.Mittelwert(alpha.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
       ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. unstrukturierter Effekt (Region s=", i, ")", sep=""))
}
dev.off()






# Prädiktor eta

apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean)


#alpha.eta.proz.sample[1] <- 1

pdf(paste("Plots/", country, "mG_Akzeptanzwahrscheinlichkeit_eta.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))
plot(alpha.eta.proz.sample[-1], type="l", ylim=c(0,1), cex.lab=1.8, cex.axis=1.8, 
     ylab="", xlab="Ausgedünnte Iterationen", lwd=4)
dev.off()


pdf(paste("Plots/", country, "mG_", length(eta.sample[,size]), "Samplingpfade_eta.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(eta.sample[,size])){
  plot(eta.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
       cex.lab=1, cex.axis=1, main=paste("Samplingpfad Prädiktor (Region s=", i, ")", sep=""))
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(eta.sample[,size]), "Dichten_eta.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(eta.sample[,size])){
  plot(density(eta.sample[i,burn:size]), xlab=expression(paste(eta, " der Region s")), 
       lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte Prädiktor (Region s=", i, ")", sep=""), ylab="")
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(eta.sample[,size]), "ACF_eta.pdf", sep=""), width=9, height=11.69) 
par(mfrow=c(5,3))
for(i in 1:length(eta.sample[,size])){
  if(eta.sample[i,burn] != eta.sample[i,size]){
    acf(eta.sample[i,burn:size], lag.max=100, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
        cex.lab=1, cex.axis=1 , main=paste("ACF Prädiktor (Region s=", i, ")", sep=""), ylab="")
  }
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(eta.sample[,size]), "MA_eta.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(eta.sample[,size])){
  plot(gleit.Mittelwert(eta.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
       ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. Prädiktor (Region s=", i, ")", sep=""))
}
dev.off()





# Gewichte

apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=mean)


#alpha.weight.proz.sample[1] <- 1

pdf(paste("Plots/", country, "mG_Akzeptanzwahrscheinlichkeit_weight_sr.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))
plot(alpha.weight.proz.sample[-1], type="l", ylim=c(0,1), cex.lab=1.8, cex.axis=1.8, 
     ylab="", xlab="Ausgedünnte Iterationen", lwd=4)
dev.off()



pdf(paste("Plots/", country, "mG_", length(weight_sr.sample[,size]), "Samplingpfade_weight_sr.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(weight_sr.sample[,size])){
  plot(weight_sr.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
       cex.lab=1, cex.axis=1, main=paste("Samplingpfad Gewicht (Nachbarschaft g=", i, ")", sep=""))
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(weight_sr.sample[,size]), "Dichten_weight_sr.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(weight_sr.sample[,size])){
  plot(density(weight_sr.sample[i,burn:size]), xlab="Gewicht der Nachbarschaft g", 
       lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte Gewicht (Nachbarschaft g=", i, ")", sep=""), ylab="")
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(weight_sr.sample[,size]), "ACF_weight_sr.pdf", sep=""), width=9, height=11.69) 
par(mfrow=c(5,3))
for(i in 1:length(weight_sr.sample[,size])){
  if(weight_sr.sample[i,burn] != weight_sr.sample[i,size]){
    acf(weight_sr.sample[i,burn:size], lag.max=100, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
        cex.lab=1, cex.axis=1 , main=paste("ACF Gewicht (Nachbarschaft g=", i, ")", sep=""), ylab="")
  }
}
dev.off()


pdf(paste("Plots/", country, "mG_", length(weight_sr.sample[,size]), "MA_weight_sr.pdf", sep=""), width=8.27, height=11.69) 
par(mfrow=c(5,2))
for(i in 1:length(weight_sr.sample[,size])){
  plot(gleit.Mittelwert(weight_sr.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
       ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. Gewicht (Nachbarschaft g=", i, ")", sep=""))
}
dev.off()




# Beispiel "Bad Mixing" stufig

bad.mixing <- c(rnorm(n=1500, mean=5, sd=1.5), rnorm(n=800, mean=10, sd=0.8), rnorm(n=1700, mean=15, sd=1.2), rnorm(n=3000, mean=8, sd=1), rnorm(n=3000, mean=15, sd=1.2))

pdf("Plots/bad_mixing_stufig.pdf", width=8, height=8) 
par(mfrow=c(1,1))

plot(bad.mixing[-1], type="l", xlab="Iterationen", ylab=expression(paste(theta)), 
     cex.lab=1.5, cex.axis=1.5)

dev.off()


# Beispiel "Bad Mixing" slow

pdf("Plots/bad_mixing_slow.pdf", width=8, height=8) 
par(mfrow=c(1,1))

plot(eta.sample[47,-1], type="l", xlab="Iterationen", ylab=expression(paste(theta)), 
     cex.lab=1.5, cex.axis=1.5)

dev.off()


# Beispiel "Good Mixing"

pdf("Plots/good_mixing.pdf", width=8, height=8) 
par(mfrow=c(1,1))

plot(beta.sample[-1], xlab="Iterationen", ylab=expression(paste(theta)), 
     cex.lab=1.5, cex.axis=1.5, type="l")

dev.off()


# Berechne gleitenden Mittelwert "Bad Mixing" für Beispiel
gleit.Mittelwert.bad <- function(parameter.sample){
  rm <- 0
  for(i in 1:length(parameter.sample[(burn+1):size])){
    if(i < 3000){
      rm[i] <- sum(parameter.sample[(burn+1):(burn+i)])/i
    }else if(i >=3000){
      if(i < 5000){
        rm[i] <- sum(parameter.sample[(burn+1):(burn+i)])/(i+(i-3000)) 
      }else if(i >= 5000){
        rm[i] <- sum(parameter.sample[(burn+1):(burn+i)])/(i+(i-3000)-(i-5000))
      }
    }  
  }
  return(rm)
}


# Mittelwertsentwicklung beu schlechtem Mixing

pdf("Plots/scotoG_badMA_kappa.pdf", width=10, height=8) 
par(mfrow=c(1,1))

plot(gleit.Mittelwert.bad(kappa.sample), type="l", xlab="Iterationen",
     ylab="", cex.lab=1.5, cex.axis=1.5, lwd=2)

dev.off()


pdf("Plots/scotoG_goodMA_kappa.pdf", width=10, height=8) 
par(mfrow=c(1,1))

plot(gleit.Mittelwert(kappa.sample), type="l", xlab="Iterationen",
     ylab="", cex.lab=1.5, cex.axis=1.5, lwd=2)

dev.off()





# Beispiel ACF Bad Mixing

a <- 1

for(i in 2:10000){
  a[i] <- a[i-1]*0.99
}
a

pdf("Plots/ACF_bad.pdf", width=10, height=8) 
par(mfrow=c(1,1))

acf(a,lag.max=50, xlab="Verzögerte Iterationen", 
    cex.lab=1.8, cex.axis=1.8 , main="", ylab="")

dev.off()


# Beispiel ACF Good Mixing

b <- 1

for(i in 2:1000){
  b[i] <- b[i-1]*0.9
}
b

pdf("Plots/ACF_good.pdf", width=10, height=8) 
par(mfrow=c(1,1))

acf(b,lag.max=50, xlab="Verzögerte Iterationen", 
    cex.lab=1.8, cex.axis=1.8 , main="", ylab="")

dev.off()

