###########################################################################################
###########################################################################################
###########################################################################################
##################################### MSE Berechnung ######################################
###########################################################################################
###########################################################################################
###########################################################################################



#################################
## Lade Package und Funktionen ##
#################################

library(Matrix)
library(lattice)
library(ggplot2)
library(reshape2)
library(car)
library(RColorBrewer)

source("extrafunktionen.R") # Laden der benötigten Funktionen
# beinhaltet die Funktionen:
# 1) "readandchangegraph" zum Einlesen einer .graph Datei (Nachbarschaftsmatrix)
# 2) "knappeMatrix" zur Umformung einer Matrix zu einer sparse Matrix
# 3) "MatrizenIndizes" zum Ausgeben der besetzten Zeilen- und Spaltenindizes einer Matrix
# 4) "rmvnorm_RBA" zum Ziehen aus einer multivariaten Normalverteilung
# 5) "sortierte.indizes" Funktion für Bestandteile einer "sparseMatrix" in geordneter Form



########################
## Einlesen der Daten ##
########################

# Laden der reduzierten schottischen Nachbarschaftsmatrix
redscotQ <- as.matrix(read.table("redscotQ.txt"))
colnames(redscotQ) <- c()

# Entferne Inseln
dim.scotQ <- dim(redscotQ)[1]
redscotQ <- redscotQ[-c(dim.scotQ-2, dim.scotQ-1, dim.scotQ), -c(dim.scotQ-2, dim.scotQ-1, dim.scotQ)] 

# Festlegung der Parameter des hybrid.sampler für SChottlanddaten ohne Inseln
Q.scot <- redscotQ

# Bestimmung der Gewichtseinträge der Nachbarschaftsmatrix
oD_Q <- matrix(data=0, nrow=dim(Q.scot)[1], ncol=dim(Q.scot)[2])
for(r in 1:dim(Q.scot)[1]){
  for(c in 1:dim(Q.scot)[2]){
    if(r < c){
      oD_Q[r,c] <- Q.scot[r,c]
    }
  }
}

# setze Diagonalelemente 0, da dort die negative Summe der Gewichte der jeweiligen Zeile eingesetzt wird ohne extra Sampling
diag(oD_Q) <- 0 
# Abspeichern der Zeilenindizes ohne Diagonale, für die gesampelt werden muss!
zeilenindizes.oD_Q <- MatrizenIndizes(oD_Q)[[1]] 
# Abspeichern der Spaltenindizes ohne Diagonale, für die gesampelt werden muss!
spaltenindizes.oD_Q <- MatrizenIndizes(oD_Q)[[2]]


# Lege Farbpalette fest

my.blue <- brewer.pal(9, "Blues")
my.red <- brewer.pal(9, "Reds")
my.col <- c(my.blue[4], my.red[4])

my.orange <- brewer.pal(9, "Oranges")
my.green <- brewer.pal(9, "Greens")
my.col2 <- c(my.green[4], my.orange[4]) 

my.col3 <- brewer.pal(12, "Paired")


####################################
## Laden der simulierten Gewichte ##
####################################

# Setting 1

Q.setting1 <- as.matrix(read.table("redscotQ.stark.setting1.txt"))
colnames(Q.setting1) <- c()
Q.setting1 <- Q.setting1[-(54:56),-(54:56)] # ohne Inseln

weight_sr.setting1 <- vector(length = length(zeilenindizes.oD_Q))

for(p in 1:length(zeilenindizes.oD_Q)){ 
  weight_sr.setting1[p] <- (Q.setting1[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]])*(-1)
}

table(weight_sr.setting1)


# Setting 2

Q.setting2 <- as.matrix(read.table("redscotQ.stark.setting2.txt"))
colnames(Q.setting2) <- c()
Q.setting2 <- Q.setting2[-(54:56),-(54:56)] # ohne Inseln

weight_sr.setting2 <- vector(length = length(zeilenindizes.oD_Q))

for(p in 1:length(zeilenindizes.oD_Q)){ 
  weight_sr.setting2[p] <- (Q.setting2[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]])*(-1)
}

table(weight_sr.setting2)


# Setting 3

Q.setting3 <- as.matrix(read.table("redscotQ.stark.setting3.txt"))
colnames(Q.setting3) <- c()
Q.setting3 <- Q.setting3[-(54:56),-(54:56)] # ohne Inseln

weight_sr.setting3 <- vector(length = length(zeilenindizes.oD_Q))

for(p in 1:length(zeilenindizes.oD_Q)){ 
  weight_sr.setting3[p] <- (Q.setting3[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]])*(-1)
}

table(weight_sr.setting3)


# Setting 4

Q.setting4 <- as.matrix(read.table("redscotQ.stark.setting4.txt"))
colnames(Q.setting4) <- c()
Q.setting4 <- Q.setting4[-(54:56),-(54:56)] # ohne Inseln

weight_sr.setting4 <- vector(length = length(zeilenindizes.oD_Q))

for(p in 1:length(zeilenindizes.oD_Q)){ 
  weight_sr.setting4[p] <- (Q.setting4[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]])*(-1)
}

table(weight_sr.setting4)


########################################
## Lade simulierten Daten und Effekte ##
########################################

# Setting 1

load("Simulierte Daten/y.sim.setting1.mx.RData")
load("Simulierte Daten/Gamma.sim.setting1.mx.RData")
load("Simulierte Daten/alpha.sim.setting1.mx.RData")
load("Simulierte Daten/sum.Gamma.alpha.sim.setting1.mx.RData")


# Setting 2

load("Simulierte Daten/y.sim.setting2.mx.RData")
load("Simulierte Daten/Gamma.sim.setting2.mx.RData")
load("Simulierte Daten/alpha.sim.setting2.mx.RData")
load("Simulierte Daten/sum.Gamma.alpha.sim.setting2.mx.RData")


# Setting 3

load("Simulierte Daten/y.sim.setting3.mx.RData")
load("Simulierte Daten/Gamma.sim.setting3.mx.RData")
load("Simulierte Daten/alpha.sim.setting3.mx.RData")
load("Simulierte Daten/sum.Gamma.alpha.sim.setting3.mx.RData")


# Setting 4

load("Simulierte Daten/y.sim.setting4.mx.RData")
load("Simulierte Daten/Gamma.sim.setting4.mx.RData")
load("Simulierte Daten/alpha.sim.setting4.mx.RData")
load("Simulierte Daten/sum.Gamma.alpha.sim.setting4.mx.RData")



#####################################
## Lade Verteilungscharakteristika ##
#####################################

for(i in 1:4){
  load(paste("Simulierte Samples/Gamma.sample.mean.zent.setting", i, ".mx.RData", sep=""))  
  load(paste("Simulierte Samples/alpha.sample.mean.zent.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Samples/sum.Gamma.alpha.sample.mean.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Samples/eta.sample.mean.setting", i, ".mx.RData", sep=""))  
  load(paste("Simulierte Samples/weight_sr.sample.mean.setting", i, ".mx.RData", sep=""))  
}


##################
## Berechne MSE ##
##################

# Berechne MSE Setting 1 - summierte räumliche Effekte
MSE.sGa.setting1 <- (apply((sum.Gamma.alpha.sample.mean.setting1.mx - sum.Gamma.alpha.sim.setting1.mx)^2, MARGIN=2, FUN=sum))/53
summary(MSE.sGa.setting1)

# Berechne MSE Setting 2 - summierte räumliche Effekte
MSE.sGa.setting2 <- (apply((sum.Gamma.alpha.sample.mean.setting2.mx - sum.Gamma.alpha.sim.setting2.mx)^2, MARGIN=2, FUN=sum))/53
summary(MSE.sGa.setting2)

# Berechne MSE Setting 3 - summierte räumliche Effekte
MSE.sGa.setting3 <- (apply((sum.Gamma.alpha.sample.mean.setting3.mx - sum.Gamma.alpha.sim.setting3.mx)^2, MARGIN=2, FUN=sum))/53
summary(MSE.sGa.setting3)

# Berechne MSE Setting 4 - summierte räumliche Effekte
MSE.sGa.setting4 <- (apply((sum.Gamma.alpha.sample.mean.setting4.mx - sum.Gamma.alpha.sim.setting4.mx)^2, MARGIN=2, FUN=sum))/53
summary(MSE.sGa.setting4)

wilcox.test(MSE.sGa.setting1, MSE.sGa.setting2, exact=TRUE, alternative="less")
wilcox.test(MSE.sGa.setting1, MSE.sGa.setting3, exact=TRUE, alternative="greater")
wilcox.test(MSE.sGa.setting1, MSE.sGa.setting4, exact=TRUE, alternative="less")
wilcox.test(MSE.sGa.setting2, MSE.sGa.setting3, exact=TRUE, alternative="greater")
wilcox.test(MSE.sGa.setting2, MSE.sGa.setting4, exact=TRUE, alternative="greater")
wilcox.test(MSE.sGa.setting3, MSE.sGa.setting4, exact=TRUE, alternative="less")



# Berechne MSE Setting 1 - Gewichte
MSE.weights.setting1 <- (apply((weight_sr.sample.mean.setting1.mx - weight_sr.setting1)^2, MARGIN=2, FUN=sum))/117
summary(MSE.weights.setting1)

# Berechne MSE Setting 2 - Gewichte
MSE.weights.setting2 <- (apply((weight_sr.sample.mean.setting2.mx - weight_sr.setting2)^2, MARGIN=2, FUN=sum))/117
summary(MSE.weights.setting2)

# Berechne MSE Setting 3 - Gewichte
MSE.weights.setting3 <- (apply((weight_sr.sample.mean.setting3.mx - weight_sr.setting3)^2, MARGIN=2, FUN=sum))/117
summary(MSE.weights.setting3)

# Berechne MSE Setting 4 - Gewichte
MSE.weights.setting4 <- (apply((weight_sr.sample.mean.setting4.mx - weight_sr.setting4)^2, MARGIN=2, FUN=sum))/117
summary(MSE.weights.setting4)

wilcox.test(MSE.weights.setting1, MSE.weights.setting2, exact=TRUE, alternative="less")
wilcox.test(MSE.weights.setting1, MSE.weights.setting3, exact=TRUE, alternative="less")
wilcox.test(MSE.weights.setting1, MSE.weights.setting4, exact=TRUE, alternative="less")
wilcox.test(MSE.weights.setting2, MSE.weights.setting3, exact=TRUE, alternative="less")
wilcox.test(MSE.weights.setting2, MSE.weights.setting4, exact=TRUE, alternative="less")
wilcox.test(MSE.weights.setting3, MSE.weights.setting4, exact=TRUE, alternative="less")


################
## Plotte MSE ##
################


# geschätzte summierte mittlere Effekte vs. simulierte summierte standardisierte Effekte 

pdf(paste("Plots/MSE_sGa_setting1.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.sGa.setting1, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.sGa.setting1), col="blue", lty=2, lwd=2)
points(MSE.sGa.setting1, pch=16, col="grey40")
legend("bottomleft",col="blue", paste("mittlerer MSE = ", round(mean(MSE.sGa.setting1),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_sGa_setting2.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.sGa.setting2, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.sGa.setting2), col="blue", lty=2, lwd=2)
points(MSE.sGa.setting2, pch=16, col="grey40")
legend("topright",col="blue", paste("mittlerer MSE = ", round(mean(MSE.sGa.setting2),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_sGa_setting3.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.sGa.setting3, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.sGa.setting3), col="blue", lty=2, lwd=2)
points(MSE.sGa.setting3, pch=16, col="grey40")
legend("topright",col="blue", paste("mittlerer MSE = ", round(mean(MSE.sGa.setting3),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_sGa_setting4.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.sGa.setting4, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.sGa.setting4), col="blue", lty=2, lwd=2)
points(MSE.sGa.setting4, pch=16, col="grey40")
legend("bottomright",col="blue", paste("mittlerer MSE = ", round(mean(MSE.sGa.setting4),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_sGa_allsettings.pdf", sep=""), width=10, height=7) 
par(mfrow=c(1,1))

plot(MSE.sGa.setting2, pch=16, col=my.red[6], ylab="Mittlerer quadratischer Fehler", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.sGa.setting2), col=my.red[6], lty=2, lwd=2)
abline(h=mean(MSE.sGa.setting1), col=my.blue[6], lty=2, lwd=2)
abline(h=mean(MSE.sGa.setting3), col=my.green[6], lty=2, lwd=2)
abline(h=mean(MSE.sGa.setting4), col=my.orange[6], lty=2, lwd=2)
lines(MSE.sGa.setting2, pch=16, col=my.red[6], lty=3, lwd=2)
lines(MSE.sGa.setting1, pch=16, col=my.blue[6], lty=3, lwd=1.5)
lines(MSE.sGa.setting3, pch=16, col=my.green[6], lty=3, lwd=2)
lines(MSE.sGa.setting4, pch=16, col=my.orange[6], lty=3, lwd=2)
points(MSE.sGa.setting2, pch=16, col=my.red[6])
points(MSE.sGa.setting1, pch=16, col=my.blue[6])
points(MSE.sGa.setting3, pch=16, col=my.green[6])
points(MSE.sGa.setting4, pch=16, col=my.orange[6])
legend("topright", col=c("white", my.blue[6], my.red[6], my.green[6], my.orange[6]), 
       c("Mittlerer MSE:", paste("Setting 1 = ", round(mean(MSE.sGa.setting1),2), sep=""), 
         paste("Setting 2 = ", round(mean(MSE.sGa.setting2),2), sep=""),
         paste("Setting 3 = ", round(mean(MSE.sGa.setting3),2), sep=""),
         paste("Setting 4 = ", round(mean(MSE.sGa.setting4),2), sep="")),
         lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()




MSE.sGa.allsettings <- c(MSE.sGa.setting1, MSE.sGa.setting2, MSE.sGa.setting3, MSE.sGa.setting4)
setting.nr <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10))
MSE.sGa.df <- as.data.frame(cbind(MSE.sGa.allsettings, setting.nr))


# einfarbig

plot.MSE.sGa <- ggplot(MSE.sGa.df, aes(factor(setting.nr), MSE.sGa.allsettings))
plot.MSE.sGa + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2,outlier.colour = NULL) + theme_bw() + labs(x="Setting", y="Mittlerer quadrierter Fehler") +
           theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
           stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") 


# mehrfarbig

plot.MSE.sGa <- ggplot(MSE.sGa.df, aes(factor(setting.nr), MSE.sGa.allsettings, fill = factor(setting.nr)))
plot.MSE.sGa + geom_boxplot(size=1.1, fatten=1.8, outlier.colour = NULL) + theme_bw() + labs(x="Setting", y="Mittlerer quadratischer Fehler - Räumliche Effekte") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") +
  scale_fill_manual(name = "Setting", values = c(my.col, my.col2))




# geschätzte mittlere Gewichte vs. festgelegte Gewichte 

pdf(paste("Plots/MSE_weights_setting1.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.weights.setting1, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.weights.setting1), col="blue", lty=2, lwd=2)
points(MSE.weights.setting1, pch=16, col="grey40")
legend("bottomleft",col="blue", paste("mittlerer MSE = ", round(mean(MSE.weights.setting1),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_weights_setting2.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.weights.setting2, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.weights.setting2), col="blue", lty=2, lwd=2)
points(MSE.weights.setting2, pch=16, col="grey40")
legend("bottomright",col="blue", paste("mittlerer MSE = ", round(mean(MSE.weights.setting2),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_weights_setting3.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.weights.setting3, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.weights.setting3), col="blue", lty=2, lwd=2)
points(MSE.weights.setting3, pch=16, col="grey40")
legend("bottomleft",col="blue", paste("mittlerer MSE = ", round(mean(MSE.weights.setting3),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()


pdf(paste("Plots/MSE_weights_setting4.pdf", sep=""), width=10, height=8) 
par(mfrow=c(1,1))

plot(MSE.weights.setting4, pch=16, col="grey40", ylab="Mittlerer quadrierter Fehler (MSE)", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.weights.setting4), col="blue", lty=2, lwd=2)
points(MSE.weights.setting4, pch=16, col="grey40")
legend("topright",col="blue", paste("mittlerer MSE = ", round(mean(MSE.weights.setting4),2), sep="") , lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()



MSE.weights.allsettings <- c(MSE.weights.setting1, MSE.weights.setting2, MSE.weights.setting3, MSE.weights.setting4)
setting.nr <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10))
MSE.weights.df <- as.data.frame(cbind(MSE.weights.allsettings, setting.nr))

# einfarbig

plot.MSE.weights <- ggplot(MSE.weights.df, aes(factor(setting.nr), MSE.weights.allsettings))
plot.MSE.weights + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2, outlier.colour = NULL) + theme_bw() + labs(x="Setting", y="Mittlerer quadrierter Fehler") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40")


# mehrfarbig

plot.MSE.weights <- ggplot(MSE.weights.df, aes(factor(setting.nr), MSE.weights.allsettings, fill = factor(setting.nr)))
plot.MSE.weights + geom_boxplot(size=1.1, fatten=1.2,outlier.colour = NULL) + theme_bw() + labs(x="Setting", y="Mittlerer quadratischer Fehler") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") +
  scale_fill_manual(name = "Setting", values = c(my.col, my.col2))



pdf(paste("Plots/MSE_weights_allsettings.pdf", sep=""), width=10, height=7) 
par(mfrow=c(1,1))

plot(MSE.weights.setting2, pch=16, col=my.red[6], ylim=c(0.13, 0.38), ylab="Mittlerer quadratischer Fehler", xlab="Simulationsdurchlauf", xaxt="n", panel.first=grid(col="gray80"))
axis(1, at=1:10, labels=1:10, las=1, cex.axis=1)
abline(v=1:10, col="gray80", lty=3)
abline(h=mean(MSE.weights.setting2), col=my.red[6], lty=2, lwd=2)
abline(h=mean(MSE.weights.setting1), col=my.blue[6], lty=2, lwd=2)
abline(h=mean(MSE.weights.setting3), col=my.green[6], lty=2, lwd=2)
abline(h=mean(MSE.weights.setting4), col=my.orange[6], lty=2, lwd=2)
lines(MSE.weights.setting2, pch=16, col=my.red[6], lty=3, lwd=2)
lines(MSE.weights.setting1, pch=16, col=my.blue[6], lty=3, lwd=1.5)
lines(MSE.weights.setting3, pch=16, col=my.green[6], lty=3, lwd=2)
lines(MSE.weights.setting4, pch=16, col=my.orange[6], lty=3, lwd=2)
points(MSE.weights.setting2, pch=16, col=my.red[6])
points(MSE.weights.setting1, pch=16, col=my.blue[6])
points(MSE.weights.setting3, pch=16, col=my.green[6])
points(MSE.weights.setting4, pch=16, col=my.orange[6])
legend("bottomleft", col=c("white", my.blue[6], my.red[6], my.green[6], my.orange[6]), 
       c("Mittlerer MSE:", paste("Setting 1 = ", round(mean(MSE.weights.setting1),2), sep=""), 
         paste("Setting 2 = ", round(mean(MSE.weights.setting2),2), sep=""),
         paste("Setting 3 = ", round(mean(MSE.weights.setting3),2), sep=""),
         paste("Setting 4 = ", round(mean(MSE.weights.setting4),2), sep="")),
       lty=2, lwd=1.5, inset=0.025, bg="white")

dev.off()






# Geschätzte mittlere Gewichte Setting 1

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting1 <- weight_sr.sample.mean.setting1.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting1 <- c(allweights.mean.sample.setting1, weight_sr.sample.mean.setting1.mx[,i])
}

allweights.mean.sample.setting1.df <- as.data.frame(cbind(allweights.mean.sample.setting1, simlauf))

plot.mean.sample.setting1 <- ggplot(allweights.mean.sample.setting1.df, aes(factor(simlauf), allweights.mean.sample.setting1))
plot.mean.sample.setting1 + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2,outlier.colour = NULL) + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") 



# Geschätzte mittlere Gewichte Setting 2

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting2 <- weight_sr.sample.mean.setting2.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting2 <- c(allweights.mean.sample.setting2, weight_sr.sample.mean.setting2.mx[,i])
}

allweights.mean.sample.setting2.df <- as.data.frame(cbind(allweights.mean.sample.setting2, simlauf))

plot.mean.sample.setting2 <- ggplot(allweights.mean.sample.setting2.df, aes(factor(simlauf), allweights.mean.sample.setting2))
plot.mean.sample.setting2 + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2,outlier.colour = NULL) + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") 


# Geschätzte mittlere Gewichte Setting 3

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting3 <- weight_sr.sample.mean.setting3.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting3 <- c(allweights.mean.sample.setting3, weight_sr.sample.mean.setting3.mx[,i])
}

allweights.mean.sample.setting3.df <- as.data.frame(cbind(allweights.mean.sample.setting3, simlauf))

plot.mean.sample.setting3 <- ggplot(allweights.mean.sample.setting3.df, aes(factor(simlauf), allweights.mean.sample.setting3))
plot.mean.sample.setting3 + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2,outlier.colour = NULL) + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") 



# Geschätzte mittlere Gewichte Setting 4

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting4 <- weight_sr.sample.mean.setting4.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting4 <- c(allweights.mean.sample.setting4, weight_sr.sample.mean.setting4.mx[,i])
}

allweights.mean.sample.setting4.df <- as.data.frame(cbind(allweights.mean.sample.setting4, simlauf))

plot.mean.sample.setting4 <- ggplot(allweights.mean.sample.setting4.df, aes(factor(simlauf), allweights.mean.sample.setting4))
plot.mean.sample.setting4 + geom_boxplot(fill = "white", colour = "#3366FF", size=0.8, fatten=2, outlier.colour = NULL) + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") +
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  stat_summary(fun.y=mean, geom="point", color="gray40") + stat_summary(fun.y=mean, geom="line", aes(group=1), color="gray40") 







# Alternative mit Splitting nach Faktor simulierte Gewichte

# Geschätzte mittlere Gewichte Setting 1

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting1 <- weight_sr.sample.mean.setting1.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting1 <- c(allweights.mean.sample.setting1, weight_sr.sample.mean.setting1.mx[,i])
}

allweights.mean.sample.setting1.df <- as.data.frame(cbind(allweights.mean.sample.setting1, simlauf, rep(weight_sr.setting1, 10)))
colnames(allweights.mean.sample.setting1.df) <- c("weights.mean.set1", "simlauf", "weights.sim.set1")
allweights.mean.sample.setting1.df$simlauf <- as.factor(allweights.mean.sample.setting1.df$simlauf)
allweights.mean.sample.setting1.df$weights.sim.set1 <- as.factor(allweights.mean.sample.setting1.df$weights.sim.set1)

plot.mean.sample.setting1 <- ggplot(allweights.mean.sample.setting1.df, aes(simlauf, weights.mean.set1, fill = factor(weights.sim.set1)))
plot.mean.sample.setting1 + geom_boxplot() + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) + scale_fill_manual(name = "sim. Gewichte", values = my.red[4])



# Geschätzte mittlere Gewichte Setting 2

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting2 <- weight_sr.sample.mean.setting2.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting2 <- c(allweights.mean.sample.setting2, weight_sr.sample.mean.setting2.mx[,i])
}

allweights.mean.sample.setting2.df <- as.data.frame(cbind(allweights.mean.sample.setting2, simlauf, rep(weight_sr.setting2, 10)))
colnames(allweights.mean.sample.setting2.df) <- c("weights.mean.set2", "simlauf", "weights.sim.set2")
allweights.mean.sample.setting2.df$simlauf <- as.factor(allweights.mean.sample.setting2.df$simlauf)
allweights.mean.sample.setting2.df$weights.sim.set2 <- as.factor(allweights.mean.sample.setting2.df$weights.sim.set2)

plot.mean.sample.setting2 <- ggplot(allweights.mean.sample.setting2.df, aes(simlauf, weights.mean.set2, fill = factor(weights.sim.set2)))
plot.mean.sample.setting2 + geom_boxplot() + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) + scale_fill_manual(name = "sim. Gewichte", values = my.col)



# Geschätzte mittlere Gewichte Setting 3

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting3 <- weight_sr.sample.mean.setting3.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting3 <- c(allweights.mean.sample.setting3, weight_sr.sample.mean.setting3.mx[,i])
}

allweights.mean.sample.setting3.df <- as.data.frame(cbind(allweights.mean.sample.setting3, simlauf, rep(weight_sr.setting3, 10)))
colnames(allweights.mean.sample.setting3.df) <- c("weights.mean.set3", "simlauf", "weights.sim.set3")
allweights.mean.sample.setting3.df$simlauf <- as.factor(allweights.mean.sample.setting3.df$simlauf)
allweights.mean.sample.setting3.df$weights.sim.set3 <- as.factor(allweights.mean.sample.setting3.df$weights.sim.set3)

plot.mean.sample.setting3 <- ggplot(allweights.mean.sample.setting3.df, aes(simlauf, weights.mean.set3, fill = factor(weights.sim.set3)))
plot.mean.sample.setting3 + geom_boxplot() + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) + scale_fill_manual(name = "sim. Gewichte", values = my.col)



# Geschätzte mittlere Gewichte Setting 4

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

allweights.mean.sample.setting4 <- weight_sr.sample.mean.setting4.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting4 <- c(allweights.mean.sample.setting4, weight_sr.sample.mean.setting4.mx[,i])
}

allweights.mean.sample.setting4.df <- as.data.frame(cbind(allweights.mean.sample.setting4, simlauf, rep(weight_sr.setting4, 10)))
colnames(allweights.mean.sample.setting4.df) <- c("weights.mean.set4", "simlauf", "weights.sim.set4")
allweights.mean.sample.setting4.df$simlauf <- as.factor(allweights.mean.sample.setting4.df$simlauf)
allweights.mean.sample.setting4.df$weights.sim.set4 <- as.factor(allweights.mean.sample.setting4.df$weights.sim.set4)

plot.mean.sample.setting4 <- ggplot(allweights.mean.sample.setting4.df, aes(simlauf, weights.mean.set4, fill = factor(weights.sim.set4)))
plot.mean.sample.setting4 + geom_boxplot() + theme_bw() + labs(x="Simulationsdurchlauf", y="Geschätzte mittlere Gewichte") + 
    theme(panel.grid.major=element_line(color="gray", linetype="dotted")) + scale_fill_manual(name = "sim. Gewichte", values = my.col)



# Plotte geschätzte mittlere Gewichte vs. Settings nach Faktor simulierte Gewichte

simlauf <- c(rep(1, 117), rep(2, 117), rep(3, 117), rep(4, 117), rep(5, 117), 
             rep(6, 117), rep(7, 117), rep(8, 117), rep(9, 117), rep(10, 117))

# Gewichte Setting 1

allweights.mean.sample.setting1 <- weight_sr.sample.mean.setting1.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting1 <- c(allweights.mean.sample.setting1, weight_sr.sample.mean.setting1.mx[,i])
}

for.correct.levels <- data.frame(1, 1, 0.1, 1) 

allweights.mean.sample.setting1.df <- as.data.frame(cbind(allweights.mean.sample.setting1, simlauf, rep(weight_sr.setting1, 10), rep(1, 1170)))
colnames(for.correct.levels) <- colnames(allweights.mean.sample.setting1.df)
allweights.mean.sample.setting1.df <- rbind(allweights.mean.sample.setting1.df, for.correct.levels)
colnames(allweights.mean.sample.setting1.df) <- c("weights.mean.set", "simlauf", "weights.sim.set", "Setting")
allweights.mean.sample.setting1.df$simlauf <- as.factor(allweights.mean.sample.setting1.df$simlauf)
allweights.mean.sample.setting1.df$weights.sim.set <- as.factor(allweights.mean.sample.setting1.df$weights.sim.set)
allweights.mean.sample.setting1.df$Setting <- as.factor(allweights.mean.sample.setting1.df$Setting)
allweights.mean.sample.setting1.df <- allweights.mean.sample.setting1.df[-1171,]


# Gewichte Setting 2

allweights.mean.sample.setting2 <- weight_sr.sample.mean.setting2.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting2 <- c(allweights.mean.sample.setting2, weight_sr.sample.mean.setting2.mx[,i])
}

allweights.mean.sample.setting2.df <- as.data.frame(cbind(allweights.mean.sample.setting2, simlauf, rep(weight_sr.setting2, 10), rep(2, 1170)))
colnames(allweights.mean.sample.setting2.df) <- c("weights.mean.set", "simlauf", "weights.sim.set", "Setting")
allweights.mean.sample.setting2.df$simlauf <- as.factor(allweights.mean.sample.setting2.df$simlauf)
allweights.mean.sample.setting2.df$weights.sim.set <- as.factor(allweights.mean.sample.setting2.df$weights.sim.set)
allweights.mean.sample.setting2.df$Setting <- as.factor(allweights.mean.sample.setting2.df$Setting)


# Gewichte Setting 3

allweights.mean.sample.setting3 <- weight_sr.sample.mean.setting3.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting3 <- c(allweights.mean.sample.setting3, weight_sr.sample.mean.setting3.mx[,i])
}

allweights.mean.sample.setting3.df <- as.data.frame(cbind(allweights.mean.sample.setting3, simlauf, rep(weight_sr.setting3, 10), rep(3, 1170)))
colnames(allweights.mean.sample.setting3.df) <- c("weights.mean.set", "simlauf", "weights.sim.set", "Setting")
allweights.mean.sample.setting3.df$simlauf <- as.factor(allweights.mean.sample.setting3.df$simlauf)
allweights.mean.sample.setting3.df$weights.sim.set <- as.factor(allweights.mean.sample.setting3.df$weights.sim.set)
allweights.mean.sample.setting3.df$Setting <- as.factor(allweights.mean.sample.setting3.df$Setting)



# Gewichte Setting 4

allweights.mean.sample.setting4 <- weight_sr.sample.mean.setting4.mx[,1]
for(i in 2:10){
  allweights.mean.sample.setting4 <- c(allweights.mean.sample.setting4, weight_sr.sample.mean.setting4.mx[,i])
}

allweights.mean.sample.setting4.df <- as.data.frame(cbind(allweights.mean.sample.setting4, simlauf, rep(weight_sr.setting4, 10), rep(4, 1170)))
colnames(allweights.mean.sample.setting4.df) <- c("weights.mean.set", "simlauf", "weights.sim.set", "Setting")
allweights.mean.sample.setting4.df$simlauf <- as.factor(allweights.mean.sample.setting4.df$simlauf)
allweights.mean.sample.setting4.df$weights.sim.set <- as.factor(allweights.mean.sample.setting4.df$weights.sim.set)
allweights.mean.sample.setting4.df$Setting <- as.factor(allweights.mean.sample.setting4.df$Setting)


allweights.sample.df <- rbind(allweights.mean.sample.setting1.df, allweights.mean.sample.setting2.df, 
                              allweights.mean.sample.setting3.df, allweights.mean.sample.setting4.df)


plot.mean.sample <- ggplot(allweights.sample.df, aes(Setting, weights.mean.set, fill = factor(weights.sim.set)))
plot.mean.sample + geom_boxplot(size=0.8, fatten=1.8) + theme_bw() + labs(x="Setting", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) + 
  scale_fill_manual(name = "sim. Gewichte", values = my.col3[1:2]) 




#################################################################
## Plotte die durchschnittliche Anzahl an Nachbarn pro Gewicht ##
#################################################################

avN.perweight <- 0

for(i in 1:length(zeilenindizes.oD_Q)){
  avN.perweight[i] <- (diag(Q.scot)[zeilenindizes.oD_Q[i]] + diag(Q.scot)[spaltenindizes.oD_Q[i]])/2
}

avN.perweight <- as.factor(avN.perweight)
avN.perweight  <- recode(avN.perweight , "c(1.5, 2.5,3,3.5)=1; c(4,4.5,5)=2; c(5.5,6,6.5)=3; c(7,7.5,8,8.5,10)=4") #Umcodiert
levels(avN.perweight) <- c("1 - 3.5", "4 - 5", "5.5 - 6.5", "7 - 10")


table(avN.perweight)

allweights.sample.df.erw <- cbind(allweights.sample.df, rep(avN.perweight, 40))
colnames(allweights.sample.df.erw) <- c("weights.mean.set", "simlauf", "weights.sim.set", "Setting", "avN.perweight")



# nach Faktor avN.perweight
plot.mean.sample <- ggplot(allweights.sample.df.erw, aes(Setting, weights.mean.set, fill = factor(avN.perweight)))
plot.mean.sample + geom_boxplot(size=1.1, fatten=1.8) + theme_bw() + labs(x="Setting", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) +
  scale_fill_manual(name = "Nacharn pro Gewicht", values = c(my.col, my.col2))


# nach Faktor Setting
plot.mean.sample <- ggplot(allweights.sample.df.erw, aes(avN.perweight, weights.mean.set, fill = factor(Setting)))
plot.mean.sample + geom_boxplot(size=0.8, fatten=1.8) + theme_bw() + labs(x="Mittlere Anzahl an Nachbarn pro Gewicht", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted")) +
  scale_fill_manual(name = "Setting", values = c(my.col, my.col2))


# nur nach Setting
plot.mean.sample <- ggplot(allweights.sample.df.erw, aes(Setting, weights.mean.set, fill = factor(Setting)))
plot.mean.sample + geom_boxplot(size=1.1, fatten=1.8) + theme_bw() + labs(x="Setting", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  scale_fill_manual(name = "Setting", values = c(my.col, my.col2))




#############################
## Wilcoxon-Rangsummentest ##
#############################


# Für Settings

allweights.sampling.setting <- data.frame(allweights.mean.sample.setting1, allweights.mean.sample.setting2,
                                         allweights.mean.sample.setting3, allweights.mean.sample.setting4)

colnames(allweights.sampling.setting) <- c("Setting1", "Setting2", "Setting3", "Setting4")


wilcox.test(allweights.sampling.setting$Setting1, allweights.sampling.setting$Setting2, exact=FALSE, alternative="two.sided")
wilcox.test(allweights.sampling.setting$Setting1, allweights.sampling.setting$Setting3, exact=FALSE, alternative="two.sided")
wilcox.test(allweights.sampling.setting$Setting1, allweights.sampling.setting$Setting4, exact=FALSE, alternative="two.sided")
wilcox.test(allweights.sampling.setting$Setting2, allweights.sampling.setting$Setting3, exact=FALSE, alternative="two.sided")
wilcox.test(allweights.sampling.setting$Setting2, allweights.sampling.setting$Setting4, exact=FALSE, alternative="two.sided")
wilcox.test(allweights.sampling.setting$Setting3, allweights.sampling.setting$Setting4, exact=FALSE, alternative="two.sided")



# Für Setting unterteilt nach sim. Gewichte

weights.sampling.setting2.simweights <- data.frame(allweights.mean.sample.setting2, rep(weight_sr.setting2, 10))
colnames(weights.sampling.setting2.simweights) <- c("weights", "simweights")

weights.sampling.setting3.simweights <- data.frame(allweights.mean.sample.setting3, rep(weight_sr.setting3, 10))
colnames(weights.sampling.setting3.simweights) <- c("weights", "simweights")

weights.sampling.setting4.simweights <- data.frame(allweights.mean.sample.setting4, rep(weight_sr.setting4, 10))
colnames(weights.sampling.setting4.simweights) <- c("weights", "simweights")

wilcox.test(subset(weights.sampling.setting2.simweights, simweights==1.5)[,1], subset(weights.sampling.setting2.simweights, simweights==0.1)[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights, simweights==1.5)[,1], subset(weights.sampling.setting3.simweights, simweights==0.1)[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights, simweights==1.5)[,1], subset(weights.sampling.setting4.simweights, simweights==0.1)[,1], exact=FALSE, alternative="greater")



# Für Nachbarn/Gewicht nach Setting

weights.sampling.setting1.simweights.nperw <- data.frame(allweights.mean.sample.setting1, rep(weight_sr.setting1, 10), rep(avN.perweight, 10))
colnames(weights.sampling.setting1.simweights.nperw) <- c("weights", "simweights", "avN.perweight")

weights.sampling.setting2.simweights.nperw <- data.frame(allweights.mean.sample.setting2, rep(weight_sr.setting2, 10), rep(avN.perweight, 10))
colnames(weights.sampling.setting2.simweights.nperw) <- c("weights", "simweights", "avN.perweight")

weights.sampling.setting3.simweights.nperw <- data.frame(allweights.mean.sample.setting3, rep(weight_sr.setting3, 10), rep(avN.perweight, 10))
colnames(weights.sampling.setting3.simweights.nperw) <- c("weights", "simweights", "avN.perweight")

weights.sampling.setting4.simweights.nperw <- data.frame(allweights.mean.sample.setting4, rep(weight_sr.setting4, 10), rep(avN.perweight, 10))
colnames(weights.sampling.setting4.simweights.nperw) <- c("weights", "simweights", "avN.perweight")



table(weights.sampling.setting1.simweights.nperw$simweights, weights.sampling.setting1.simweights.nperw$avN.perweight)
table(weights.sampling.setting2.simweights.nperw$simweights, weights.sampling.setting2.simweights.nperw$avN.perweight)
table(weights.sampling.setting3.simweights.nperw$simweights, weights.sampling.setting3.simweights.nperw$avN.perweight)
table(weights.sampling.setting4.simweights.nperw$simweights, weights.sampling.setting4.simweights.nperw$avN.perweight)



wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="4 - 5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], subset(weights.sampling.setting1.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")

wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="4 - 5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], subset(weights.sampling.setting2.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")

wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="4 - 5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], subset(weights.sampling.setting3.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")

wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="4 - 5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="1 - 3.5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="4 - 5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")
wilcox.test(subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="5.5 - 6.5")[,1], subset(weights.sampling.setting4.simweights.nperw, avN.perweight=="7 - 10")[,1], exact=FALSE, alternative="greater")




######################
## Empirische Daten ##
######################

load("Gewichte empirische Daten/Schottland/weight_sr.sample.scot.RData")
weight_sr.sample.scot <- weight_sr.sample

load("Gewichte empirische Daten/Schottland/weight_sr.sample.mean.scot.RData")
weight_sr.sample.mean.scot <- weight_sr.sample.mean


load("Gewichte empirische Daten/Deutschland/weight_sr.sample.german.RData")
weight_sr.sample.german <- weight_sr.sample

load("Gewichte empirische Daten/Deutschland/weight_sr.sample.mean.german.RData")
weight_sr.sample.mean.german <- weight_sr.sample.mean



plot(weight_sr.sample.mean.scot, ylim=c(1, 1.62), cex.lab=0.9, cex.axis=0.9, pch=16, cex=0.7, col=1, ylab="Geschätzte mittlere Gewichte", xlab="", xaxt="n")
grid()
points(weight_sr.sample.mean.scot, pch=16, cex=0.7, col=1)
points(x=which(weight_sr.sample.mean.scot > 1.5), 
       y=subset(weight_sr.sample.mean.scot, weight_sr.sample.mean.scot > 1.5),
       pch=16, cex=0.7, col="firebrick1")
# axis(1, at=1:117, labels=FALSE)


plot(weight_sr.sample.mean.german, ylim=c(0.37, 0.66), cex.lab=1, cex.axis=1, pch=16, cex=0.7, col=1, ylab="Geschätzte mittlere Gewichte", xlab="", xaxt="n")     
grid()
points(weight_sr.sample.mean.german, pch=16, cex=0.7, col=1)
points(x=which(weight_sr.sample.mean.german > 0.55), 
       y=subset(weight_sr.sample.mean.german, weight_sr.sample.mean.german > 0.55),
       pch=16, cex=0.7, col="firebrick1")
# axis(1, at=1:1416, labels=FALSE)
dev.off()


# Durchschnittliche Anzahl an Nachbarn pro Gewicht für Schottland

weight_sr.nperw.scot <- data.frame(weight_sr.sample.mean.scot, avN.perweight)

plot.mean.sample <- ggplot(weight_sr.nperw.scot, aes(avN.perweight, weight_sr.sample.mean.scot, fill = factor(avN.perweight)))
plot.mean.sample + geom_boxplot(size=0.8, fatten=1.8) + theme_bw() + labs(x="Mittlere Anzahl an Nachbarn pro Gewicht", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  scale_fill_manual(name = "Setting", values = c(my.blue[2], my.blue[4], my.blue[6], my.blue[8]))


# Rangsummen-Test
sub.weight_sr.nperw.scot.kat1 <- subset(weight_sr.nperw.scot, avN.perweight == "1 - 3.5")
sub.weight_sr.nperw.scot.kat2 <- subset(weight_sr.nperw.scot, avN.perweight == "4 - 5")
sub.weight_sr.nperw.scot.kat3 <- subset(weight_sr.nperw.scot, avN.perweight == "5.5 - 6.5")
sub.weight_sr.nperw.scot.kat4 <- subset(weight_sr.nperw.scot, avN.perweight == "7 - 10")

wilcox.test(sub.weight_sr.nperw.scot.kat1[,1], sub.weight_sr.nperw.scot.kat2[,1], exact=TRUE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.scot.kat1[,1], sub.weight_sr.nperw.scot.kat3[,1], exact=TRUE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.scot.kat1[,1], sub.weight_sr.nperw.scot.kat4[,1], exact=TRUE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.scot.kat2[,1], sub.weight_sr.nperw.scot.kat3[,1], exact=TRUE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.scot.kat2[,1], sub.weight_sr.nperw.scot.kat4[,1], exact=TRUE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.scot.kat3[,1], sub.weight_sr.nperw.scot.kat4[,1], exact=TRUE, alternative="greater")



which(weight_sr.sample.mean.scot > 1.5)
zeilenindizes.oD_Q[c(116,117)]
spaltenindizes.oD_Q[c(116,117)]


# Berechne durchschnittliche Anzahl an Nachbarn pro Gewicht für Deutschland 

# Laden der reduzierten deutschen Nachbarschaftsmatrix
redgermanQ <- as.matrix(read.table("redgermanyQ.txt"))
colnames(redgermanQ) <- c()

# Festlegung der Parameter des hybrid.sampler 
Q.german <- redgermanQ

# Bestimmung der Gewichtseinträge der Nachbarschaftsmatrix
oD_Qg <- matrix(data=0, nrow=dim(Q.german)[1], ncol=dim(Q.german)[2])
for(r in 1:dim(Q.german)[1]){
  for(c in 1:dim(Q.german)[2]){
    if(r < c){
      oD_Qg[r,c] <- Q.german[r,c]
    }
  }
}

# setze Diagonalelemente 0, da dort die negative Summe der Gewichte der jeweiligen Zeile eingesetzt wird ohne extra Sampling
diag(oD_Qg) <- 0 
# Abspeichern der Zeilenindizes ohne Diagonale, für die gesampelt werden muss!
zeilenindizes.oD_Qg <- MatrizenIndizes(oD_Qg)[[1]] 
# Abspeichern der Spaltenindizes ohne Diagonale, für die gesampelt werden muss!
spaltenindizes.oD_Qg <- MatrizenIndizes(oD_Qg)[[2]]


avN.perweightg <- 0

for(i in 1:length(zeilenindizes.oD_Qg)){
  avN.perweightg[i] <- (diag(Q.german)[zeilenindizes.oD_Qg[i]] + diag(Q.german)[spaltenindizes.oD_Qg[i]])/2
}

avN.perweightg <- as.factor(avN.perweightg)
avN.perweightg  <- recode(avN.perweightg , "c(2, 2.5,3,3.5)=1; c(4,4.5,5)=2; c(5.5,6,6.5)=3; c(7,7.5,8,8.5,9,9.5,10,10.5)=4") #Umcodiert
levels(avN.perweightg) <- c("1 - 3.5", "4 - 5", "5.5 - 6.5", "7 - 10.5")


table(avN.perweightg)

weight_sr.nperw.german <- data.frame(weight_sr.sample.mean.german, avN.perweightg)

plot.mean.sample <- ggplot(weight_sr.nperw.german, aes(avN.perweightg, weight_sr.sample.mean.german, fill = factor(avN.perweightg)))
plot.mean.sample + geom_boxplot(size=0.8, fatten=1.8) + theme_bw() + labs(x="Mittlere Anzahl an Nachbarn pro Gewicht", y="Geschätzte mittlere Gewichte") + 
  theme(panel.grid.major=element_line(color="gray", linetype="dotted"), legend.position="none") +
  scale_fill_manual(name = "Setting", values = c(my.green[2], my.green[4], my.green[6], my.green[8]))

sum(table(subset(weight_sr.nperw.german, weight_sr.sample.mean.german > 0.55)$avN.perweightg))


# Rangsummen-Test
sub.weight_sr.nperw.german.kat1 <- subset(weight_sr.nperw.german, avN.perweightg == "1 - 3.5")
sub.weight_sr.nperw.german.kat2 <- subset(weight_sr.nperw.german, avN.perweightg == "4 - 5")
sub.weight_sr.nperw.german.kat3 <- subset(weight_sr.nperw.german, avN.perweightg == "5.5 - 6.5")
sub.weight_sr.nperw.german.kat4 <- subset(weight_sr.nperw.german, avN.perweightg == "7 - 10.5")

wilcox.test(sub.weight_sr.nperw.german.kat1[,1], sub.weight_sr.nperw.german.kat2[,1], exact=FALSE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.german.kat1[,1], sub.weight_sr.nperw.german.kat3[,1], exact=FALSE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.german.kat1[,1], sub.weight_sr.nperw.german.kat4[,1], exact=FALSE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.german.kat2[,1], sub.weight_sr.nperw.german.kat3[,1], exact=FALSE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.german.kat2[,1], sub.weight_sr.nperw.german.kat4[,1], exact=FALSE, alternative="greater")
wilcox.test(sub.weight_sr.nperw.german.kat3[,1], sub.weight_sr.nperw.german.kat4[,1], exact=FALSE, alternative="greater")
