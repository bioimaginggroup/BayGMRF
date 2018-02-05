###########################################################################################
###########################################################################################
###########################################################################################
################################## Schottlanddaten Plots ##################################
###########################################################################################
###########################################################################################
###########################################################################################



##############################################################
## Einlesen der schottischen Karten- und Dateninformationen ##
##############################################################

# Für die Aufbereitung der schottischen Daten 
# schlägt Bivand (Quelle: bivand2008 S. 90 f.) vor nach dem Import 
# der Karteninformationen als CRS, diese als "British National Grid"
# zu transformieren. Dazu stellt Bivand folgenden Code bereit:

# Laden des Pakets rgdal
library(rgdal)

# Einlesen der Daten über .shp, .shx und .dbf File.
scot_LL <- readOGR(".", "scot")

# Importieren der Landkreisgrenzen und übergeben der korrekten CRS
proj4string(scot_LL) <- CRS("+proj=longlat ellps=WGS84")

# Transformieren zu "British National Grid"
EPSG <- make_EPSG()
EPSG[grep("British National Grid", EPSG$note), 1:2]
scot_BNG0 <- spTransform(scot_LL, CRS("+init=epsg:27700"))


# Die ursprünglichen Daten von Clayton and Kaldor (1987)
# verwenden die gleichen Landkreis ID's, aber in einer
# anderen Reihenfolge. Diese müssen deshalb mit den Datensatz
# gematcht werden, um bei den Landkreisgrenzpolygonen die
# gleiche Reihenfolge herzustellen.


# Laden des Pakets maptools
library(maptools)

# Laden der schottischen Daten
scot_dat <- read.table("scotland.dat", skip=1)
names(scot_dat) <- c("District", "Observed", "Expected", "PcAFF", "Latitude", "Longitude")
row.names(scot_dat) <- formatC(scot_dat$District, width=2, flag="0")

summary(scot_dat)

# Vorbereitung und matching
ID <- formatC(scot_BNG0$ID, width=2, flag="0")
scot_BNG1 <- spChFIDs(scot_BNG0, ID)
scot_BNG <- spCbind(scot_BNG1, scot_dat[match(ID, row.names(scot_dat)),])

scot_BNG



########################################
## Deskriptive Aufbereitung der Daten ##
########################################

# Plot der schottischen Landkreisgrenzen
oopar <- par(mar=c(0,0,0,0))
plot(scot_BNG, border="black")

# Die Originaldaten beinhalten die erwartete Anzahl an Fällen, die auf den
# Alterseffekten basieren nach den der Methodik von Mantel und Stark (1968).
# Diese sind proportional zu einer Risikopopulation nachdem solche Effekte
# in die Beschreibung der Daten mit aufgenommen wurden.

# Festlegung einer Farbpalette
bluepal <- colorRampPalette(c("white", "midnightblue"))

# Festlegung eines Intervalls für scot_BNG$Observed
brks <- c()
for(i in 2:15){
  brks[1] <- 0
  brks[i] <- max(scot_BNG$Observed)/16 * (i-1)
  brks[16] <- max(scot_BNG$Observed)          
}

# Lade das Paket classInt um Daten in Intervalle aufzuteilen
library(classInt)
O_CI <- classIntervals(scot_BNG$Observed, style="fixed", fixedBreaks=brks)
E_CI <- classIntervals(scot_BNG$Expected, style="fixed", fixedBreaks=brks)



# Bestimmung der Farben
O_CIc <- findColours(O_CI, bluepal(2))
E_CIc <- findColours(E_CI, bluepal(2))


# Plot der beobachteten Daten
plot(scot_BNG, col=O_CIc)
legend("topleft", fill=attr(O_CIc, "palette"), legend=names(attr(O_CIc, "table")), bty="n")

# Plot der erwarteten Daten
plot(scot_BNG, col=E_CIc)
legend("topleft", fill=attr(E_CIc, "palette"), legend=names(attr(E_CIc, "table")), bty="n")


###################################################

# Berechne standardisierte Sterblichkeitsrate
scot_BNG$SMR <- scot_BNG$Observed/scot_BNG$Expected 


# Berechne logarithmierte Variable Observed. Problem Nuller im Datensatz! Deswegen For-Schleife
scot_BNG$logObserved <- 0

for(i in 1:length(scot_BNG$Observed)){
  if(scot_BNG$Observed[i] == 0){
    scot_BNG$logObserved[i] <- 0
  }else{
    scot_BNG$logObserved[i] <- log(scot_BNG$Observed[i])
  }
}


# Berechne logarithmierte Variable Expected
scot_BNG$logExpected <- log(scot_BNG$Expected)




# Alternative Plots über spplot (normal und logarithmiert). Jeweils einmal mit und einmal ohne Rahmenlinien. Muss dann evtl. zusammenkopiert werden.
# spplot(scot_BNG, c("SMR"), col.regions=bluepal(16), par.settings = list(axis.line = list(col = 'transparent')))

# Bestimmung des Plotbereichs
xlim.scot <- c(7101.358-30000, 468284.2+30000)
ylim.scot <- c(529569.050-30000, 1218409.5+30000)

library(lattice)

tps <- list(fontsize=list(text=16))
trellis.par.set(tps)

spplot(scot_BNG, c("Expected"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)
spplot(scot_BNG, c("logExpected"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)

spplot(scot_BNG, c("Observed"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)
spplot(scot_BNG, c("logObserved"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)

spplot(scot_BNG, c("SMR"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)
spplot(scot_BNG, c("PcAFF"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)



# Einfärbung einzelner Distrikte (Inseln und wichtige Städte)
scot_BNG$InselStadt <- rep(0, length(scot_BNG$Observed))
scot_BNG$InselStadt[c(47,48,52)] <- 100 # Inseln: Shetland, Orkney und Äußere Hybriden
scot_BNG$InselStadt[c(8,25,29,46)] <- 50 # Städte: 8=Dundee, 25=Glasgow, 29=Edinburgh, 46=Aberdeen 

grayredpal <- colorRampPalette(c("white", "firebrick1", "gray80"))
brks.InselStadt <- c(1,50,100)

InselStadt_CI <- classIntervals(scot_BNG$InselStadt, style="fixed", fixedBreaks=brks.InselStadt)
InselStadt_CIc <- findColours(InselStadt_CI, grayredpal(3))

plot(scot_BNG, col=InselStadt_CIc)




################################################
## Kartenaufbereitung für Simulationssettings ##
################################################

# Setting 1 der Simulation
scot_BNG$Setting1 <- rep(100, length(scot_BNG$Observed))
scot_BNG$Setting1[c(47,48,52)] <- 0 # Inseln: Shetland, Orkney und Äußere Hybriden

graypal <- colorRampPalette(c("white", "gray40"))
brks.Setting1 <- c(0,100)

Setting1_CI <- classIntervals(scot_BNG$Setting1, style="fixed", fixedBreaks=brks.Setting1)
Setting1_CIc <- findColours(Setting1_CI, graypal(2))

plot(scot_BNG, col=Setting1_CIc)


# Festlegung der Grenzen für Setting 2 der Simulation
scot_BNG$Setting2 <- rep(0, length(scot_BNG$Observed))
scot_BNG$Setting2[c(1:8,12,14,19,45:46,49:51,53:55)] <- 50 # Distrikte oberhalb der Grenze
scot_BNG$Setting2[c(9:11,13,15:18,20:44,56)] <- 100 # Distrikte unterhalb der Grenze

graypal <- colorRampPalette(c("white", "gray40"))
brks.Setting2 <- c(50,100)

Setting2_CI <- classIntervals(scot_BNG$Setting2, style="fixed", fixedBreaks=brks.Setting2)
Setting2_CIc <- findColours(Setting2_CI, graypal(2))

plot(scot_BNG, col=Setting2_CIc)


# Festlegung der Grenzen für Setting 3 der Simulation
scot_BNG$Setting3 <- rep(0, length(scot_BNG$Observed))
scot_BNG$Setting3[c(1:8,12,14,19,45,46,49:51,53:55)] <- 67 # Distrikte oberhalb der Grenze
scot_BNG$Setting3[c(9:11,13,15:18,20:31,56)] <- 100 # Distrikte in der Mitte
scot_BNG$Setting3[c(32:44)] <- 33 # Distrikte unterhalb der Grenze

graypal <- colorRampPalette(c("white", "gray40"))
brks.Setting3 <- c(50,100)

Setting3_CI <- classIntervals(scot_BNG$Setting3, style="fixed", fixedBreaks=brks.Setting3)
Setting3_CIc <- findColours(Setting3_CI, graypal(2))

plot(scot_BNG, col=Setting3_CIc)


# Festlegung der Grenzen für Setting 4 der Simulation
scot_BNG$Setting4 <- rep(0, length(scot_BNG$Observed))
scot_BNG$Setting4[c(1:3,14,19,49:51,53,55)] <- 100 # Distrikte links oberhalb der Grenze
scot_BNG$Setting4[c(4:8,12,45:46,54)] <- 50 # Distrikte rechts oberhalb der Grenze
scot_BNG$Setting4[c(9:11,13,15:18,20:31,56)] <- 25 # Distrikte in der Mitte
scot_BNG$Setting4[c(32:44)] <- 75 # Distrikte unterhalb der Grenze

graypal <- colorRampPalette(c("white", "gray40"))
brks.Setting4 <- c(50,100)

Setting4_CI <- classIntervals(scot_BNG$Setting4, style="fixed", fixedBreaks=brks.Setting4)
Setting4_CIc <- findColours(Setting4_CI, graypal(2))

plot(scot_BNG, col=Setting4_CIc)



par(mfrow=c(2,2))
plot(scot_BNG, border="black")
plot(scot_BNG, col=Setting2_CIc)
plot(scot_BNG, col=Setting3_CIc)
plot(scot_BNG, col=Setting4_CIc)

par(mfrow=c(1,1))
plot(scot_BNG, border="black")
plot(scot_BNG, col=Setting2_CIc)
plot(scot_BNG, col=Setting3_CIc)
plot(scot_BNG, col=Setting4_CIc)




##############################################################################
## Kartenaufbereitung für Schätzergebnisse - empirische Daten ohne Gewichte ##
##############################################################################

########## Achtung!!! ##########
# Beim plotten der Ergebnisse aufpassen auf ID Zuordnung!!!

# Rücktransformation der Reihenfolge der Ergebnisse

# Reihenfolge der Originalindizes
scotRedIndizes <- c(19,3,31,13,34,26,17,55,32,27,30,36,35,45,23,54,44,46,40,43,52,47,29,37,48,53,41,
                    51,50,39,21,9,1,42,20,25,24,49,14,33,38,15,6,28,22,16,12,8,18,0,4,11,2,10,7,5)

# alle Werte +1, da bei Index 0 gestartet wird
scotRedIndizes <- scotRedIndizes + 1


load("Ohne Gewichte/Ga.sample.mean.RData")
load("Ohne Gewichte/Ga.sample.mean.zent.RData")
load("Ohne Gewichte/Ga.sample.sign.RData")
load("Ohne Gewichte/alpha.sample.mean.RData")
load("Ohne Gewichte/alpha.sample.mean.zent.RData")
load("Ohne Gewichte/alpha.sample.sign.RData")
load("Ohne Gewichte/eta.sample.mean.RData")



# mit Inseln auf 1000, sind dann weiß. da außerhalb des Intervals
Ga.sample.mean <- c(Ga.sample.mean,1000,1000,1000)
Ga.sample.mean.zent <- c(Ga.sample.mean.zent,1000,1000,1000)
Ga.sample.sign <- c(Ga.sample.sign,0,0,0)
alpha.sample.mean <- c(alpha.sample.mean,1000,1000,1000)
alpha.sample.mean.zent <- c(alpha.sample.mean.zent,1000,1000,1000)
alpha.sample.sign <- c(alpha.sample.sign,0,0,0)
eta.sample.mean <- c(eta.sample.mean,1000,1000,1000)


# Stelle alte Reihenfolge vor reduzierter Form wieder her
ergebnisse <- data.frame(scotRedIndizes, Ga.sample.mean, Ga.sample.mean.zent, Ga.sample.sign,
                         alpha.sample.mean, alpha.sample.mean.zent, alpha.sample.sign,
                         eta.sample.mean)

ergebnisse <- ergebnisse[order(ergebnisse$scotRedIndizes),]

# Berücksichtige Sortierung nach der Variable ID
wandler <- 0
for(i in 1:length(scot_BNG$ID)){
  wandler[i] <- which(scot_BNG$ID == i)
}
wandler

ergebnisse <- cbind(wandler, ergebnisse)
ergebnisse <- ergebnisse[order(ergebnisse$wandler),]

ergebnisse

scot_BNG$Ga.sample.mean <- ergebnisse$Ga.sample.mean
scot_BNG$Ga.sample.mean.zent <- ergebnisse$Ga.sample.mean.zent
scot_BNG$Ga.sample.sign <- ergebnisse$Ga.sample.sign
scot_BNG$alpha.sample.mean <- ergebnisse$alpha.sample.mean
scot_BNG$alpha.sample.mean.zent <- ergebnisse$alpha.sample.mean.zent
scot_BNG$alpha.sample.sign <- ergebnisse$alpha.sample.sign
scot_BNG$eta.sample.mean <- ergebnisse$eta.sample.mean
scot_BNG$exp.eta.sample.mean <- exp(ergebnisse$eta.sample.mean)
scot_BNG$logObserved2 <- scot_BNG$logObserved
scot_BNG$logObserved2[c(47,48,52)] <- 0
scot_BNG$Observed2 <- scot_BNG$Observed
scot_BNG$Observed2[c(47,48,52)] <- 0
scot_BNG$sum.Ga.alpha.sample.mean.zent <- scot_BNG$Ga.sample.mean.zent + scot_BNG$alpha.sample.mean.zent


summary(scot_BNG)


# Plots ohne Gewichtsadaption

tps <- list(fontsize=list(text=16))
trellis.par.set(tps)

redbluepal <- colorRampPalette(c("firebrick", "white", "midnightblue"))


# ursprüngliche strukturierte Effekte (nicht zentriert)
ogrenze <- max((ceiling(sort(abs(Ga.sample.mean), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("Ga.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)


# ursprüngliche unstrukturierte Effekte (nicht zentriert)
ogrenze <- max((ceiling(sort(abs(alpha.sample.mean), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("alpha.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)


# zentrierte strukturierte Effekte
ogrenze <- max((ceiling(sort(abs(Ga.sample.mean.zent), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("Ga.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

# zentrierte unstrukturierte Effekte
ogrenze <- max((ceiling(sort(abs(alpha.sample.mean.zent), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("alpha.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

# zentrierte summierte räumliche Effekte
ogrenze <- max((ceiling(sort(abs(scot_BNG$sum.Ga.alpha.sample.mean.zent), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("sum.Ga.alpha.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)




# eta mean
ogrenze <- max((ceiling(sort(abs(scot_BNG$eta.sample.mean), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("eta.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

spplot(scot_BNG, c("logObserved2"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)

# exp eta mean
spplot(scot_BNG, c("exp.eta.sample.mean"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)
spplot(scot_BNG, c("Observed2"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)




# Signifikanzen
scot_BNG$Ga.sample.sign
scot_BNG$Ga.sample.sign[1] <- -1
scot_BNG$Ga.sample.sign <- factor(scot_BNG$Ga.sample.sign, labels=c("neg", "0", "pos"))
scot_BNG$Ga.sample.sign[1] <- 0
scot_BNG$Ga.sample.sign[c(47,48,52)] <- 1000 # Inseln raus


spplot(scot_BNG, "Ga.sample.sign", col="black", col.regions=c("firebrick", "gray90", "midnightblue"), xlim = xlim.scot, ylim = ylim.scot)



scot_BNG$alpha.sample.sign <- factor(scot_BNG$alpha.sample.sign, labels=c("neg", "0", "pos"))
scot_BNG$alpha.sample.sign[c(47,48,52)] <- 1000 # Inseln raus

spplot(scot_BNG, "alpha.sample.sign", col="black", col.regions=c("firebrick", "gray90", "midnightblue"), xlim = xlim.scot, ylim = ylim.scot)





################################################################
## Kartenaufbereitung für Schätzergebnisse - simulierte Daten ##
################################################################

########## Achtung!!! ##########
# Beim plotten der Ergebnisse aufpassen auf ID Zuordnung!!!

# Rücktransformation der Reihenfolge der Ergebnisse



# Reihenfolge der Originalindizes
scotRedIndizes <- c(19,3,31,13,34,26,17,55,32,27,30,36,35,45,23,54,44,46,40,43,52,47,29,37,48,53,41,
                    51,50,39,21,9,1,42,20,25,24,49,14,33,38,15,6,28,22,16,12,8,18,0,4,11,2,10,7,5)

# alle Werte +1, da bei Index 0 gestartet wird
scotRedIndizes <- scotRedIndizes + 1

for(i in 1:4){

  load(paste("Simulierte Daten/Simuliert/y.sim.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Simuliert/Gamma.sim.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Simuliert/alpha.sim.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Simuliert/sum.Gamma.alpha.sim.setting", i, ".mx.RData", sep=""))
  
  
  
  load(paste("Simulierte Daten/Gesamplet/eta.sample.mean.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Gesamplet/Gamma.sample.mean.zent.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Gesamplet/alpha.sample.mean.zent.setting", i, ".mx.RData", sep=""))
  load(paste("Simulierte Daten/Gesamplet/sum.Gamma.alpha.sample.mean.setting", i, ".mx.RData", sep=""))
  
  
}


# Vorbereitung der Matrizen
# mit Inseln auf 1000, sind dann weiß. da außerhalb des Intervals
Insel_tausender <- matrix(data=1000, nrow=3, ncol=10)
Insel_nuller <- matrix(data=0, nrow=3, ncol=10)

# Setting 1

y.sim.setting1.mx <- rbind(y.sim.setting1.mx, Insel_nuller)
colnames(y.sim.setting1.mx) <- c("y.sim1.setting1", "y.sim2.setting1", "y.sim3.setting1", "y.sim4.setting1", 
                                 "y.sim5.setting1","y.sim6.setting1","y.sim7.setting1", "y.sim8.setting1", 
                                 "y.sim9.setting1", "y.sim10.setting1")

Gamma.sim.setting1.mx <- rbind(Gamma.sim.setting1.mx, Insel_tausender)
colnames(Gamma.sim.setting1.mx) <- c("Gamma.sim1.setting1", "Gamma.sim2.setting1", "Gamma.sim3.setting1", "Gamma.sim4.setting1", 
                                 "Gamma.sim5.setting1","Gamma.sim6.setting1","Gamma.sim7.setting1", "Gamma.sim8.setting1", 
                                 "Gamma.sim9.setting1", "Gamma.sim10.setting1")

alpha.sim.setting1.mx <- rbind(alpha.sim.setting1.mx, Insel_tausender)
colnames(alpha.sim.setting1.mx) <- c("alpha.sim1.setting1", "alpha.sim2.setting1", "alpha.sim3.setting1", "alpha.sim4.setting1", 
                                     "alpha.sim5.setting1","alpha.sim6.setting1","alpha.sim7.setting1", "alpha.sim8.setting1", 
                                     "alpha.sim9.setting1", "alpha.sim10.setting1")

sum.Gamma.alpha.sim.setting1.mx <- rbind(sum.Gamma.alpha.sim.setting1.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sim.setting1.mx) <- c("sGa.sim1.setting1", "sGa.sim2.setting1", "sGa.sim3.setting1", "sGa.sim4.setting1", 
                                     "sGa.sim5.setting1","sGa.sim6.setting1","sGa.sim7.setting1", "sGa.sim8.setting1", 
                                     "sGa.sim9.setting1", "sGa.sim10.setting1")

eta.sample.mean.setting1.mx <- rbind(eta.sample.mean.setting1.mx, Insel_tausender)
colnames(eta.sample.mean.setting1.mx) <- c("eta.sample1.setting1", "eta.sample2.setting1", "eta.sample3.setting1", "eta.sample4.setting1", 
                                               "eta.sample5.setting1","eta.sample6.setting1","eta.sample7.setting1", "eta.sample8.setting1", 
                                               "eta.sample9.setting1", "eta.sample10.setting1")

Gamma.sample.mean.zent.setting1.mx <- rbind(Gamma.sample.mean.zent.setting1.mx, Insel_tausender)
colnames(Gamma.sample.mean.zent.setting1.mx) <- c("Gamma.sample1.setting1", "Gamma.sample2.setting1", "Gamma.sample3.setting1", "Gamma.sample4.setting1", 
                                           "Gamma.sample5.setting1","Gamma.sample6.setting1","Gamma.sample7.setting1", "Gamma.sample8.setting1", 
                                           "Gamma.sample9.setting1", "Gamma.sample10.setting1")

alpha.sample.mean.zent.setting1.mx <- rbind(alpha.sample.mean.zent.setting1.mx, Insel_tausender)
colnames(alpha.sample.mean.zent.setting1.mx) <- c("alpha.sample1.setting1", "alpha.sample2.setting1", "alpha.sample3.setting1", "alpha.sample4.setting1", 
                                                  "alpha.sample5.setting1","alpha.sample6.setting1","alpha.sample7.setting1", "alpha.sample8.setting1", 
                                                  "alpha.sample9.setting1", "alpha.sample10.setting1")

sum.Gamma.alpha.sample.mean.setting1.mx <- rbind(sum.Gamma.alpha.sample.mean.setting1.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sample.mean.setting1.mx) <- c("sGa.sample1.setting1", "sGa.sample2.setting1", "sGa.sample3.setting1", "sGa.sample4.setting1", 
                                                  "sGa.sample5.setting1","sGa.sample6.setting1","sGa.sample7.setting1", "sGa.sample8.setting1", 
                                                  "sGa.sample9.setting1", "sGa.sample10.setting1")



# Setting 2

y.sim.setting2.mx <- rbind(y.sim.setting2.mx, Insel_nuller)
colnames(y.sim.setting2.mx) <- c("y.sim1.setting2", "y.sim2.setting2", "y.sim3.setting2", "y.sim4.setting2", 
                                 "y.sim5.setting2","y.sim6.setting2","y.sim7.setting2", "y.sim8.setting2", 
                                 "y.sim9.setting2", "y.sim10.setting2")

Gamma.sim.setting2.mx <- rbind(Gamma.sim.setting2.mx, Insel_tausender)
colnames(Gamma.sim.setting2.mx) <- c("Gamma.sim1.setting2", "Gamma.sim2.setting2", "Gamma.sim3.setting2", "Gamma.sim4.setting2", 
                                     "Gamma.sim5.setting2","Gamma.sim6.setting2","Gamma.sim7.setting2", "Gamma.sim8.setting2", 
                                     "Gamma.sim9.setting2", "Gamma.sim10.setting2")

alpha.sim.setting2.mx <- rbind(alpha.sim.setting2.mx, Insel_tausender)
colnames(alpha.sim.setting2.mx) <- c("alpha.sim1.setting2", "alpha.sim2.setting2", "alpha.sim3.setting2", "alpha.sim4.setting2", 
                                     "alpha.sim5.setting2","alpha.sim6.setting2","alpha.sim7.setting2", "alpha.sim8.setting2", 
                                     "alpha.sim9.setting2", "alpha.sim10.setting2")

sum.Gamma.alpha.sim.setting2.mx <- rbind(sum.Gamma.alpha.sim.setting2.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sim.setting2.mx) <- c("sGa.sim1.setting2", "sGa.sim2.setting2", "sGa.sim3.setting2", "sGa.sim4.setting2", 
                                               "sGa.sim5.setting2","sGa.sim6.setting2","sGa.sim7.setting2", "sGa.sim8.setting2", 
                                               "sGa.sim9.setting2", "sGa.sim10.setting2")

eta.sample.mean.setting2.mx <- rbind(eta.sample.mean.setting2.mx, Insel_tausender)
colnames(eta.sample.mean.setting2.mx) <- c("eta.sample1.setting2", "eta.sample2.setting2", "eta.sample3.setting2", "eta.sample4.setting2", 
                                           "eta.sample5.setting2","eta.sample6.setting2","eta.sample7.setting2", "eta.sample8.setting2", 
                                           "eta.sample9.setting2", "eta.sample10.setting2")

Gamma.sample.mean.zent.setting2.mx <- rbind(Gamma.sample.mean.zent.setting2.mx, Insel_tausender)
colnames(Gamma.sample.mean.zent.setting2.mx) <- c("Gamma.sample1.setting2", "Gamma.sample2.setting2", "Gamma.sample3.setting2", "Gamma.sample4.setting2", 
                                                  "Gamma.sample5.setting2","Gamma.sample6.setting2","Gamma.sample7.setting2", "Gamma.sample8.setting2", 
                                                  "Gamma.sample9.setting2", "Gamma.sample10.setting2")

alpha.sample.mean.zent.setting2.mx <- rbind(alpha.sample.mean.zent.setting2.mx, Insel_tausender)
colnames(alpha.sample.mean.zent.setting2.mx) <- c("alpha.sample1.setting2", "alpha.sample2.setting2", "alpha.sample3.setting2", "alpha.sample4.setting2", 
                                                  "alpha.sample5.setting2","alpha.sample6.setting2","alpha.sample7.setting2", "alpha.sample8.setting2", 
                                                  "alpha.sample9.setting2", "alpha.sample10.setting2")

sum.Gamma.alpha.sample.mean.setting2.mx <- rbind(sum.Gamma.alpha.sample.mean.setting2.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sample.mean.setting2.mx) <- c("sGa.sample1.setting2", "sGa.sample2.setting2", "sGa.sample3.setting2", "sGa.sample4.setting2", 
                                                            "sGa.sample5.setting2","sGa.sample6.setting2","sGa.sample7.setting2", "sGa.sample8.setting2", 
                                                            "sGa.sample9.setting2", "sGa.sample10.setting2")


# Setting 3

y.sim.setting3.mx <- rbind(y.sim.setting3.mx, Insel_nuller)
colnames(y.sim.setting3.mx) <- c("y.sim1.setting3", "y.sim2.setting3", "y.sim3.setting3", "y.sim4.setting3", 
                                 "y.sim5.setting3","y.sim6.setting3","y.sim7.setting3", "y.sim8.setting3", 
                                 "y.sim9.setting3", "y.sim10.setting3")

Gamma.sim.setting3.mx <- rbind(Gamma.sim.setting3.mx, Insel_tausender)
colnames(Gamma.sim.setting3.mx) <- c("Gamma.sim1.setting3", "Gamma.sim2.setting3", "Gamma.sim3.setting3", "Gamma.sim4.setting3", 
                                     "Gamma.sim5.setting3","Gamma.sim6.setting3","Gamma.sim7.setting3", "Gamma.sim8.setting3", 
                                     "Gamma.sim9.setting3", "Gamma.sim10.setting3")

alpha.sim.setting3.mx <- rbind(alpha.sim.setting3.mx, Insel_tausender)
colnames(alpha.sim.setting3.mx) <- c("alpha.sim1.setting3", "alpha.sim2.setting3", "alpha.sim3.setting3", "alpha.sim4.setting3", 
                                     "alpha.sim5.setting3","alpha.sim6.setting3","alpha.sim7.setting3", "alpha.sim8.setting3", 
                                     "alpha.sim9.setting3", "alpha.sim10.setting3")

sum.Gamma.alpha.sim.setting3.mx <- rbind(sum.Gamma.alpha.sim.setting3.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sim.setting3.mx) <- c("sGa.sim1.setting3", "sGa.sim2.setting3", "sGa.sim3.setting3", "sGa.sim4.setting3", 
                                               "sGa.sim5.setting3","sGa.sim6.setting3","sGa.sim7.setting3", "sGa.sim8.setting3", 
                                               "sGa.sim9.setting3", "sGa.sim10.setting3")

eta.sample.mean.setting3.mx <- rbind(eta.sample.mean.setting3.mx, Insel_tausender)
colnames(eta.sample.mean.setting3.mx) <- c("eta.sample1.setting3", "eta.sample2.setting3", "eta.sample3.setting3", "eta.sample4.setting3", 
                                           "eta.sample5.setting3","eta.sample6.setting3","eta.sample7.setting3", "eta.sample8.setting3", 
                                           "eta.sample9.setting3", "eta.sample10.setting3")

Gamma.sample.mean.zent.setting3.mx <- rbind(Gamma.sample.mean.zent.setting3.mx, Insel_tausender)
colnames(Gamma.sample.mean.zent.setting3.mx) <- c("Gamma.sample1.setting3", "Gamma.sample2.setting3", "Gamma.sample3.setting3", "Gamma.sample4.setting3", 
                                                  "Gamma.sample5.setting3","Gamma.sample6.setting3","Gamma.sample7.setting3", "Gamma.sample8.setting3", 
                                                  "Gamma.sample9.setting3", "Gamma.sample10.setting3")

alpha.sample.mean.zent.setting3.mx <- rbind(alpha.sample.mean.zent.setting3.mx, Insel_tausender)
colnames(alpha.sample.mean.zent.setting3.mx) <- c("alpha.sample1.setting3", "alpha.sample2.setting3", "alpha.sample3.setting3", "alpha.sample4.setting3", 
                                                  "alpha.sample5.setting3","alpha.sample6.setting3","alpha.sample7.setting3", "alpha.sample8.setting3", 
                                                  "alpha.sample9.setting3", "alpha.sample10.setting3")

sum.Gamma.alpha.sample.mean.setting3.mx <- rbind(sum.Gamma.alpha.sample.mean.setting3.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sample.mean.setting3.mx) <- c("sGa.sample1.setting3", "sGa.sample2.setting3", "sGa.sample3.setting3", "sGa.sample4.setting3", 
                                                            "sGa.sample5.setting3","sGa.sample6.setting3","sGa.sample7.setting3", "sGa.sample8.setting3", 
                                                            "sGa.sample9.setting3", "sGa.sample10.setting3")


# Setting 4

y.sim.setting4.mx <- rbind(y.sim.setting4.mx, Insel_nuller)
colnames(y.sim.setting4.mx) <- c("y.sim1.setting4", "y.sim2.setting4", "y.sim3.setting4", "y.sim4.setting4", 
                                 "y.sim5.setting4","y.sim6.setting4","y.sim7.setting4", "y.sim8.setting4", 
                                 "y.sim9.setting4", "y.sim10.setting4")

Gamma.sim.setting4.mx <- rbind(Gamma.sim.setting4.mx, Insel_tausender)
colnames(Gamma.sim.setting4.mx) <- c("Gamma.sim1.setting4", "Gamma.sim2.setting4", "Gamma.sim3.setting4", "Gamma.sim4.setting4", 
                                     "Gamma.sim5.setting4","Gamma.sim6.setting4","Gamma.sim7.setting4", "Gamma.sim8.setting4", 
                                     "Gamma.sim9.setting4", "Gamma.sim10.setting4")

alpha.sim.setting4.mx <- rbind(alpha.sim.setting4.mx, Insel_tausender)
colnames(alpha.sim.setting4.mx) <- c("alpha.sim1.setting4", "alpha.sim2.setting4", "alpha.sim3.setting4", "alpha.sim4.setting4", 
                                     "alpha.sim5.setting4","alpha.sim6.setting4","alpha.sim7.setting4", "alpha.sim8.setting4", 
                                     "alpha.sim9.setting4", "alpha.sim10.setting4")

sum.Gamma.alpha.sim.setting4.mx <- rbind(sum.Gamma.alpha.sim.setting4.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sim.setting4.mx) <- c("sGa.sim1.setting4", "sGa.sim2.setting4", "sGa.sim3.setting4", "sGa.sim4.setting4", 
                                               "sGa.sim5.setting4","sGa.sim6.setting4","sGa.sim7.setting4", "sGa.sim8.setting4", 
                                               "sGa.sim9.setting4", "sGa.sim10.setting4")

eta.sample.mean.setting4.mx <- rbind(eta.sample.mean.setting4.mx, Insel_tausender)
colnames(eta.sample.mean.setting4.mx) <- c("eta.sample1.setting4", "eta.sample2.setting4", "eta.sample3.setting4", "eta.sample4.setting4", 
                                           "eta.sample5.setting4","eta.sample6.setting4","eta.sample7.setting4", "eta.sample8.setting4", 
                                           "eta.sample9.setting4", "eta.sample10.setting4")

Gamma.sample.mean.zent.setting4.mx <- rbind(Gamma.sample.mean.zent.setting4.mx, Insel_tausender)
colnames(Gamma.sample.mean.zent.setting4.mx) <- c("Gamma.sample1.setting4", "Gamma.sample2.setting4", "Gamma.sample3.setting4", "Gamma.sample4.setting4", 
                                                  "Gamma.sample5.setting4","Gamma.sample6.setting4","Gamma.sample7.setting4", "Gamma.sample8.setting4", 
                                                  "Gamma.sample9.setting4", "Gamma.sample10.setting4")

alpha.sample.mean.zent.setting4.mx <- rbind(alpha.sample.mean.zent.setting4.mx, Insel_tausender)
colnames(alpha.sample.mean.zent.setting4.mx) <- c("alpha.sample1.setting4", "alpha.sample2.setting4", "alpha.sample3.setting4", "alpha.sample4.setting4", 
                                                  "alpha.sample5.setting4","alpha.sample6.setting4","alpha.sample7.setting4", "alpha.sample8.setting4", 
                                                  "alpha.sample9.setting4", "alpha.sample10.setting4")

sum.Gamma.alpha.sample.mean.setting4.mx <- rbind(sum.Gamma.alpha.sample.mean.setting4.mx, Insel_tausender)
colnames(sum.Gamma.alpha.sample.mean.setting4.mx) <- c("sGa.sample1.setting4", "sGa.sample2.setting4", "sGa.sample3.setting4", "sGa.sample4.setting4", 
                                                            "sGa.sample5.setting4","sGa.sample6.setting4","sGa.sample7.setting4", "sGa.sample8.setting4", 
                                                            "sGa.sample9.setting4", "sGa.sample10.setting4")


# Zusätzlich direkte Abspeicherung für Zugriff über "get"


# Setting 1

for(j in 1:10){
  assign(paste("y.sim", j,".setting1", sep=""), y.sim.setting1.mx[,j])
  assign(paste("Gamma.sim", j,".setting1", sep=""), Gamma.sim.setting1.mx[,j])
  assign(paste("alpha.sim", j,".setting1", sep=""), alpha.sim.setting1.mx[,j])
  assign(paste("sGa.sim", j,".setting1", sep=""), sum.Gamma.alpha.sim.setting1.mx[,j])
  
  assign(paste("eta.sample", j,".setting1", sep=""), eta.sample.mean.setting1.mx[,j])
  assign(paste("Gamma.sample", j,".setting1", sep=""), Gamma.sample.mean.zent.setting1.mx[,j])
  assign(paste("alpha.sample", j,".setting1", sep=""), alpha.sample.mean.zent.setting1.mx[,j])
  assign(paste("sGa.sample", j,".setting1", sep=""), sum.Gamma.alpha.sample.mean.setting1.mx[,j])
}


# Setting 2

for(j in 1:10){
  assign(paste("y.sim", j,".setting2", sep=""), y.sim.setting2.mx[,j])
  assign(paste("Gamma.sim", j,".setting2", sep=""), Gamma.sim.setting2.mx[,j])
  assign(paste("alpha.sim", j,".setting2", sep=""), alpha.sim.setting2.mx[,j])
  assign(paste("sGa.sim", j,".setting2", sep=""), sum.Gamma.alpha.sim.setting2.mx[,j])
  
  assign(paste("eta.sample", j,".setting2", sep=""), eta.sample.mean.setting2.mx[,j])
  assign(paste("Gamma.sample", j,".setting2", sep=""), Gamma.sample.mean.zent.setting2.mx[,j])
  assign(paste("alpha.sample", j,".setting2", sep=""), alpha.sample.mean.zent.setting2.mx[,j])
  assign(paste("sGa.sample", j,".setting2", sep=""), sum.Gamma.alpha.sample.mean.setting2.mx[,j])
}


# Setting 3

for(j in 1:10){
  assign(paste("y.sim", j,".setting3", sep=""), y.sim.setting3.mx[,j])
  assign(paste("Gamma.sim", j,".setting3", sep=""), Gamma.sim.setting3.mx[,j])
  assign(paste("alpha.sim", j,".setting3", sep=""), alpha.sim.setting3.mx[,j])
  assign(paste("sGa.sim", j,".setting3", sep=""), sum.Gamma.alpha.sim.setting3.mx[,j])
  
  assign(paste("eta.sample", j,".setting3", sep=""), eta.sample.mean.setting3.mx[,j])
  assign(paste("Gamma.sample", j,".setting3", sep=""), Gamma.sample.mean.zent.setting3.mx[,j])
  assign(paste("alpha.sample", j,".setting3", sep=""), alpha.sample.mean.zent.setting3.mx[,j])
  assign(paste("sGa.sample", j,".setting3", sep=""), sum.Gamma.alpha.sample.mean.setting3.mx[,j])
}


# Setting 4

for(j in 1:10){
  assign(paste("y.sim", j,".setting4", sep=""), y.sim.setting4.mx[,j])
  assign(paste("Gamma.sim", j,".setting4", sep=""), Gamma.sim.setting4.mx[,j])
  assign(paste("alpha.sim", j,".setting4", sep=""), alpha.sim.setting4.mx[,j])
  assign(paste("sGa.sim", j,".setting4", sep=""), sum.Gamma.alpha.sim.setting4.mx[,j])
  
  assign(paste("eta.sample", j,".setting4", sep=""), eta.sample.mean.setting4.mx[,j])
  assign(paste("Gamma.sample", j,".setting4", sep=""), Gamma.sample.mean.zent.setting4.mx[,j])
  assign(paste("alpha.sample", j,".setting4", sep=""), alpha.sample.mean.zent.setting4.mx[,j])
  assign(paste("sGa.sample", j,".setting4", sep=""), sum.Gamma.alpha.sample.mean.setting4.mx[,j])
}



# Stelle alte Reihenfolge vor reduzierter Form wieder her
ergebnisse <- data.frame(scotRedIndizes, y.sim.setting1.mx, y.sim.setting2.mx, y.sim.setting3.mx, y.sim.setting4.mx,
                         Gamma.sim.setting1.mx, Gamma.sim.setting2.mx, Gamma.sim.setting3.mx, Gamma.sim.setting4.mx,
                         alpha.sim.setting1.mx, alpha.sim.setting2.mx, alpha.sim.setting3.mx, alpha.sim.setting4.mx,
                         sum.Gamma.alpha.sim.setting1.mx, sum.Gamma.alpha.sim.setting2.mx, 
                         sum.Gamma.alpha.sim.setting3.mx, sum.Gamma.alpha.sim.setting4.mx,
                         eta.sample.mean.setting1.mx, eta.sample.mean.setting2.mx, 
                         eta.sample.mean.setting3.mx,eta.sample.mean.setting4.mx,
                         Gamma.sample.mean.zent.setting1.mx, Gamma.sample.mean.zent.setting2.mx,
                         Gamma.sample.mean.zent.setting3.mx, Gamma.sample.mean.zent.setting4.mx,
                         alpha.sample.mean.zent.setting1.mx, alpha.sample.mean.zent.setting2.mx,
                         alpha.sample.mean.zent.setting3.mx, alpha.sample.mean.zent.setting4.mx,
                         sum.Gamma.alpha.sample.mean.setting1.mx, sum.Gamma.alpha.sample.mean.setting2.mx,
                         sum.Gamma.alpha.sample.mean.setting3.mx, sum.Gamma.alpha.sample.mean.setting4.mx)

ergebnisse <- ergebnisse[order(ergebnisse$scotRedIndizes),]


# Berücksichtige Sortierung nach der Variable ID
wandler <- 0
for(i in 1:length(scot_BNG$ID)){
  wandler[i] <- which(scot_BNG$ID == i)
}
wandler

ergebnisse <- cbind(wandler, ergebnisse)
ergebnisse <- ergebnisse[order(ergebnisse$wandler),]

ergebnisse


# Anfügen in Datensatz

#Setting 1

scot_BNG$y.sim1.setting1 <- ergebnisse$y.sim1.setting1
scot_BNG$y.sim2.setting1 <- ergebnisse$y.sim2.setting1
scot_BNG$y.sim3.setting1 <- ergebnisse$y.sim3.setting1
scot_BNG$y.sim4.setting1 <- ergebnisse$y.sim4.setting1
scot_BNG$y.sim5.setting1 <- ergebnisse$y.sim5.setting1
scot_BNG$y.sim6.setting1 <- ergebnisse$y.sim6.setting1
scot_BNG$y.sim7.setting1 <- ergebnisse$y.sim7.setting1
scot_BNG$y.sim8.setting1 <- ergebnisse$y.sim8.setting1
scot_BNG$y.sim9.setting1 <- ergebnisse$y.sim9.setting1
scot_BNG$y.sim10.setting1 <- ergebnisse$y.sim10.setting1

scot_BNG$Gamma.sim1.setting1 <- ergebnisse$Gamma.sim1.setting1
scot_BNG$Gamma.sim2.setting1 <- ergebnisse$Gamma.sim2.setting1
scot_BNG$Gamma.sim3.setting1 <- ergebnisse$Gamma.sim3.setting1
scot_BNG$Gamma.sim4.setting1 <- ergebnisse$Gamma.sim4.setting1
scot_BNG$Gamma.sim5.setting1 <- ergebnisse$Gamma.sim5.setting1
scot_BNG$Gamma.sim6.setting1 <- ergebnisse$Gamma.sim6.setting1
scot_BNG$Gamma.sim7.setting1 <- ergebnisse$Gamma.sim7.setting1
scot_BNG$Gamma.sim8.setting1 <- ergebnisse$Gamma.sim8.setting1
scot_BNG$Gamma.sim9.setting1 <- ergebnisse$Gamma.sim9.setting1
scot_BNG$Gamma.sim10.setting1 <- ergebnisse$Gamma.sim10.setting1

scot_BNG$alpha.sim1.setting1 <- ergebnisse$alpha.sim1.setting1
scot_BNG$alpha.sim2.setting1 <- ergebnisse$alpha.sim2.setting1
scot_BNG$alpha.sim3.setting1 <- ergebnisse$alpha.sim3.setting1
scot_BNG$alpha.sim4.setting1 <- ergebnisse$alpha.sim4.setting1
scot_BNG$alpha.sim5.setting1 <- ergebnisse$alpha.sim5.setting1
scot_BNG$alpha.sim6.setting1 <- ergebnisse$alpha.sim6.setting1
scot_BNG$alpha.sim7.setting1 <- ergebnisse$alpha.sim7.setting1
scot_BNG$alpha.sim8.setting1 <- ergebnisse$alpha.sim8.setting1
scot_BNG$alpha.sim9.setting1 <- ergebnisse$alpha.sim9.setting1
scot_BNG$alpha.sim10.setting1 <- ergebnisse$alpha.sim10.setting1

scot_BNG$sGa.sim1.setting1 <- ergebnisse$sGa.sim1.setting1
scot_BNG$sGa.sim2.setting1 <- ergebnisse$sGa.sim2.setting1
scot_BNG$sGa.sim3.setting1 <- ergebnisse$sGa.sim3.setting1
scot_BNG$sGa.sim4.setting1 <- ergebnisse$sGa.sim4.setting1
scot_BNG$sGa.sim5.setting1 <- ergebnisse$sGa.sim5.setting1
scot_BNG$sGa.sim6.setting1 <- ergebnisse$sGa.sim6.setting1
scot_BNG$sGa.sim7.setting1 <- ergebnisse$sGa.sim7.setting1
scot_BNG$sGa.sim8.setting1 <- ergebnisse$sGa.sim8.setting1
scot_BNG$sGa.sim9.setting1 <- ergebnisse$sGa.sim9.setting1
scot_BNG$sGa.sim10.setting1 <- ergebnisse$sGa.sim10.setting1


scot_BNG$eta.sample1.setting1 <- ergebnisse$eta.sample1.setting1
scot_BNG$eta.sample2.setting1 <- ergebnisse$eta.sample2.setting1
scot_BNG$eta.sample3.setting1 <- ergebnisse$eta.sample3.setting1
scot_BNG$eta.sample4.setting1 <- ergebnisse$eta.sample4.setting1
scot_BNG$eta.sample5.setting1 <- ergebnisse$eta.sample5.setting1
scot_BNG$eta.sample6.setting1 <- ergebnisse$eta.sample6.setting1
scot_BNG$eta.sample7.setting1 <- ergebnisse$eta.sample7.setting1
scot_BNG$eta.sample8.setting1 <- ergebnisse$eta.sample8.setting1
scot_BNG$eta.sample9.setting1 <- ergebnisse$eta.sample9.setting1
scot_BNG$eta.sample10.setting1 <- ergebnisse$eta.sample10.setting1

scot_BNG$Gamma.sample1.setting1 <- ergebnisse$Gamma.sample1.setting1
scot_BNG$Gamma.sample2.setting1 <- ergebnisse$Gamma.sample2.setting1
scot_BNG$Gamma.sample3.setting1 <- ergebnisse$Gamma.sample3.setting1
scot_BNG$Gamma.sample4.setting1 <- ergebnisse$Gamma.sample4.setting1
scot_BNG$Gamma.sample5.setting1 <- ergebnisse$Gamma.sample5.setting1
scot_BNG$Gamma.sample6.setting1 <- ergebnisse$Gamma.sample6.setting1
scot_BNG$Gamma.sample7.setting1 <- ergebnisse$Gamma.sample7.setting1
scot_BNG$Gamma.sample8.setting1 <- ergebnisse$Gamma.sample8.setting1
scot_BNG$Gamma.sample9.setting1 <- ergebnisse$Gamma.sample9.setting1
scot_BNG$Gamma.sample10.setting1 <- ergebnisse$Gamma.sample10.setting1

scot_BNG$alpha.sample1.setting1 <- ergebnisse$alpha.sample1.setting1
scot_BNG$alpha.sample2.setting1 <- ergebnisse$alpha.sample2.setting1
scot_BNG$alpha.sample3.setting1 <- ergebnisse$alpha.sample3.setting1
scot_BNG$alpha.sample4.setting1 <- ergebnisse$alpha.sample4.setting1
scot_BNG$alpha.sample5.setting1 <- ergebnisse$alpha.sample5.setting1
scot_BNG$alpha.sample6.setting1 <- ergebnisse$alpha.sample6.setting1
scot_BNG$alpha.sample7.setting1 <- ergebnisse$alpha.sample7.setting1
scot_BNG$alpha.sample8.setting1 <- ergebnisse$alpha.sample8.setting1
scot_BNG$alpha.sample9.setting1 <- ergebnisse$alpha.sample9.setting1
scot_BNG$alpha.sample10.setting1 <- ergebnisse$alpha.sample10.setting1

scot_BNG$sGa.sample1.setting1 <- ergebnisse$sGa.sample1.setting1
scot_BNG$sGa.sample2.setting1 <- ergebnisse$sGa.sample2.setting1
scot_BNG$sGa.sample3.setting1 <- ergebnisse$sGa.sample3.setting1
scot_BNG$sGa.sample4.setting1 <- ergebnisse$sGa.sample4.setting1
scot_BNG$sGa.sample5.setting1 <- ergebnisse$sGa.sample5.setting1
scot_BNG$sGa.sample6.setting1 <- ergebnisse$sGa.sample6.setting1
scot_BNG$sGa.sample7.setting1 <- ergebnisse$sGa.sample7.setting1
scot_BNG$sGa.sample8.setting1 <- ergebnisse$sGa.sample8.setting1
scot_BNG$sGa.sample9.setting1 <- ergebnisse$sGa.sample9.setting1
scot_BNG$sGa.sample10.setting1 <- ergebnisse$sGa.sample10.setting1


#Setting 2

scot_BNG$y.sim1.setting2 <- ergebnisse$y.sim1.setting2
scot_BNG$y.sim2.setting2 <- ergebnisse$y.sim2.setting2
scot_BNG$y.sim3.setting2 <- ergebnisse$y.sim3.setting2
scot_BNG$y.sim4.setting2 <- ergebnisse$y.sim4.setting2
scot_BNG$y.sim5.setting2 <- ergebnisse$y.sim5.setting2
scot_BNG$y.sim6.setting2 <- ergebnisse$y.sim6.setting2
scot_BNG$y.sim7.setting2 <- ergebnisse$y.sim7.setting2
scot_BNG$y.sim8.setting2 <- ergebnisse$y.sim8.setting2
scot_BNG$y.sim9.setting2 <- ergebnisse$y.sim9.setting2
scot_BNG$y.sim10.setting2 <- ergebnisse$y.sim10.setting2

scot_BNG$Gamma.sim1.setting2 <- ergebnisse$Gamma.sim1.setting2
scot_BNG$Gamma.sim2.setting2 <- ergebnisse$Gamma.sim2.setting2
scot_BNG$Gamma.sim3.setting2 <- ergebnisse$Gamma.sim3.setting2
scot_BNG$Gamma.sim4.setting2 <- ergebnisse$Gamma.sim4.setting2
scot_BNG$Gamma.sim5.setting2 <- ergebnisse$Gamma.sim5.setting2
scot_BNG$Gamma.sim6.setting2 <- ergebnisse$Gamma.sim6.setting2
scot_BNG$Gamma.sim7.setting2 <- ergebnisse$Gamma.sim7.setting2
scot_BNG$Gamma.sim8.setting2 <- ergebnisse$Gamma.sim8.setting2
scot_BNG$Gamma.sim9.setting2 <- ergebnisse$Gamma.sim9.setting2
scot_BNG$Gamma.sim10.setting2 <- ergebnisse$Gamma.sim10.setting2

scot_BNG$alpha.sim1.setting2 <- ergebnisse$alpha.sim1.setting2
scot_BNG$alpha.sim2.setting2 <- ergebnisse$alpha.sim2.setting2
scot_BNG$alpha.sim3.setting2 <- ergebnisse$alpha.sim3.setting2
scot_BNG$alpha.sim4.setting2 <- ergebnisse$alpha.sim4.setting2
scot_BNG$alpha.sim5.setting2 <- ergebnisse$alpha.sim5.setting2
scot_BNG$alpha.sim6.setting2 <- ergebnisse$alpha.sim6.setting2
scot_BNG$alpha.sim7.setting2 <- ergebnisse$alpha.sim7.setting2
scot_BNG$alpha.sim8.setting2 <- ergebnisse$alpha.sim8.setting2
scot_BNG$alpha.sim9.setting2 <- ergebnisse$alpha.sim9.setting2
scot_BNG$alpha.sim10.setting2 <- ergebnisse$alpha.sim10.setting2

scot_BNG$sGa.sim1.setting2 <- ergebnisse$sGa.sim1.setting2
scot_BNG$sGa.sim2.setting2 <- ergebnisse$sGa.sim2.setting2
scot_BNG$sGa.sim3.setting2 <- ergebnisse$sGa.sim3.setting2
scot_BNG$sGa.sim4.setting2 <- ergebnisse$sGa.sim4.setting2
scot_BNG$sGa.sim5.setting2 <- ergebnisse$sGa.sim5.setting2
scot_BNG$sGa.sim6.setting2 <- ergebnisse$sGa.sim6.setting2
scot_BNG$sGa.sim7.setting2 <- ergebnisse$sGa.sim7.setting2
scot_BNG$sGa.sim8.setting2 <- ergebnisse$sGa.sim8.setting2
scot_BNG$sGa.sim9.setting2 <- ergebnisse$sGa.sim9.setting2
scot_BNG$sGa.sim10.setting2 <- ergebnisse$sGa.sim10.setting2


scot_BNG$eta.sample1.setting2 <- ergebnisse$eta.sample1.setting2
scot_BNG$eta.sample2.setting2 <- ergebnisse$eta.sample2.setting2
scot_BNG$eta.sample3.setting2 <- ergebnisse$eta.sample3.setting2
scot_BNG$eta.sample4.setting2 <- ergebnisse$eta.sample4.setting2
scot_BNG$eta.sample5.setting2 <- ergebnisse$eta.sample5.setting2
scot_BNG$eta.sample6.setting2 <- ergebnisse$eta.sample6.setting2
scot_BNG$eta.sample7.setting2 <- ergebnisse$eta.sample7.setting2
scot_BNG$eta.sample8.setting2 <- ergebnisse$eta.sample8.setting2
scot_BNG$eta.sample9.setting2 <- ergebnisse$eta.sample9.setting2
scot_BNG$eta.sample10.setting2 <- ergebnisse$eta.sample10.setting2

scot_BNG$Gamma.sample1.setting2 <- ergebnisse$Gamma.sample1.setting2
scot_BNG$Gamma.sample2.setting2 <- ergebnisse$Gamma.sample2.setting2
scot_BNG$Gamma.sample3.setting2 <- ergebnisse$Gamma.sample3.setting2
scot_BNG$Gamma.sample4.setting2 <- ergebnisse$Gamma.sample4.setting2
scot_BNG$Gamma.sample5.setting2 <- ergebnisse$Gamma.sample5.setting2
scot_BNG$Gamma.sample6.setting2 <- ergebnisse$Gamma.sample6.setting2
scot_BNG$Gamma.sample7.setting2 <- ergebnisse$Gamma.sample7.setting2
scot_BNG$Gamma.sample8.setting2 <- ergebnisse$Gamma.sample8.setting2
scot_BNG$Gamma.sample9.setting2 <- ergebnisse$Gamma.sample9.setting2
scot_BNG$Gamma.sample10.setting2 <- ergebnisse$Gamma.sample10.setting2

scot_BNG$alpha.sample1.setting2 <- ergebnisse$alpha.sample1.setting2
scot_BNG$alpha.sample2.setting2 <- ergebnisse$alpha.sample2.setting2
scot_BNG$alpha.sample3.setting2 <- ergebnisse$alpha.sample3.setting2
scot_BNG$alpha.sample4.setting2 <- ergebnisse$alpha.sample4.setting2
scot_BNG$alpha.sample5.setting2 <- ergebnisse$alpha.sample5.setting2
scot_BNG$alpha.sample6.setting2 <- ergebnisse$alpha.sample6.setting2
scot_BNG$alpha.sample7.setting2 <- ergebnisse$alpha.sample7.setting2
scot_BNG$alpha.sample8.setting2 <- ergebnisse$alpha.sample8.setting2
scot_BNG$alpha.sample9.setting2 <- ergebnisse$alpha.sample9.setting2
scot_BNG$alpha.sample10.setting2 <- ergebnisse$alpha.sample10.setting2

scot_BNG$sGa.sample1.setting2 <- ergebnisse$sGa.sample1.setting2
scot_BNG$sGa.sample2.setting2 <- ergebnisse$sGa.sample2.setting2
scot_BNG$sGa.sample3.setting2 <- ergebnisse$sGa.sample3.setting2
scot_BNG$sGa.sample4.setting2 <- ergebnisse$sGa.sample4.setting2
scot_BNG$sGa.sample5.setting2 <- ergebnisse$sGa.sample5.setting2
scot_BNG$sGa.sample6.setting2 <- ergebnisse$sGa.sample6.setting2
scot_BNG$sGa.sample7.setting2 <- ergebnisse$sGa.sample7.setting2
scot_BNG$sGa.sample8.setting2 <- ergebnisse$sGa.sample8.setting2
scot_BNG$sGa.sample9.setting2 <- ergebnisse$sGa.sample9.setting2
scot_BNG$sGa.sample10.setting2 <- ergebnisse$sGa.sample10.setting2


#Setting 3

scot_BNG$y.sim1.setting3 <- ergebnisse$y.sim1.setting3
scot_BNG$y.sim2.setting3 <- ergebnisse$y.sim2.setting3
scot_BNG$y.sim3.setting3 <- ergebnisse$y.sim3.setting3
scot_BNG$y.sim4.setting3 <- ergebnisse$y.sim4.setting3
scot_BNG$y.sim5.setting3 <- ergebnisse$y.sim5.setting3
scot_BNG$y.sim6.setting3 <- ergebnisse$y.sim6.setting3
scot_BNG$y.sim7.setting3 <- ergebnisse$y.sim7.setting3
scot_BNG$y.sim8.setting3 <- ergebnisse$y.sim8.setting3
scot_BNG$y.sim9.setting3 <- ergebnisse$y.sim9.setting3
scot_BNG$y.sim10.setting3 <- ergebnisse$y.sim10.setting3

scot_BNG$Gamma.sim1.setting3 <- ergebnisse$Gamma.sim1.setting3
scot_BNG$Gamma.sim2.setting3 <- ergebnisse$Gamma.sim2.setting3
scot_BNG$Gamma.sim3.setting3 <- ergebnisse$Gamma.sim3.setting3
scot_BNG$Gamma.sim4.setting3 <- ergebnisse$Gamma.sim4.setting3
scot_BNG$Gamma.sim5.setting3 <- ergebnisse$Gamma.sim5.setting3
scot_BNG$Gamma.sim6.setting3 <- ergebnisse$Gamma.sim6.setting3
scot_BNG$Gamma.sim7.setting3 <- ergebnisse$Gamma.sim7.setting3
scot_BNG$Gamma.sim8.setting3 <- ergebnisse$Gamma.sim8.setting3
scot_BNG$Gamma.sim9.setting3 <- ergebnisse$Gamma.sim9.setting3
scot_BNG$Gamma.sim10.setting3 <- ergebnisse$Gamma.sim10.setting3

scot_BNG$alpha.sim1.setting3 <- ergebnisse$alpha.sim1.setting3
scot_BNG$alpha.sim2.setting3 <- ergebnisse$alpha.sim2.setting3
scot_BNG$alpha.sim3.setting3 <- ergebnisse$alpha.sim3.setting3
scot_BNG$alpha.sim4.setting3 <- ergebnisse$alpha.sim4.setting3
scot_BNG$alpha.sim5.setting3 <- ergebnisse$alpha.sim5.setting3
scot_BNG$alpha.sim6.setting3 <- ergebnisse$alpha.sim6.setting3
scot_BNG$alpha.sim7.setting3 <- ergebnisse$alpha.sim7.setting3
scot_BNG$alpha.sim8.setting3 <- ergebnisse$alpha.sim8.setting3
scot_BNG$alpha.sim9.setting3 <- ergebnisse$alpha.sim9.setting3
scot_BNG$alpha.sim10.setting3 <- ergebnisse$alpha.sim10.setting3

scot_BNG$sGa.sim1.setting3 <- ergebnisse$sGa.sim1.setting3
scot_BNG$sGa.sim2.setting3 <- ergebnisse$sGa.sim2.setting3
scot_BNG$sGa.sim3.setting3 <- ergebnisse$sGa.sim3.setting3
scot_BNG$sGa.sim4.setting3 <- ergebnisse$sGa.sim4.setting3
scot_BNG$sGa.sim5.setting3 <- ergebnisse$sGa.sim5.setting3
scot_BNG$sGa.sim6.setting3 <- ergebnisse$sGa.sim6.setting3
scot_BNG$sGa.sim7.setting3 <- ergebnisse$sGa.sim7.setting3
scot_BNG$sGa.sim8.setting3 <- ergebnisse$sGa.sim8.setting3
scot_BNG$sGa.sim9.setting3 <- ergebnisse$sGa.sim9.setting3
scot_BNG$sGa.sim10.setting3 <- ergebnisse$sGa.sim10.setting3


scot_BNG$eta.sample1.setting3 <- ergebnisse$eta.sample1.setting3
scot_BNG$eta.sample2.setting3 <- ergebnisse$eta.sample2.setting3
scot_BNG$eta.sample3.setting3 <- ergebnisse$eta.sample3.setting3
scot_BNG$eta.sample4.setting3 <- ergebnisse$eta.sample4.setting3
scot_BNG$eta.sample5.setting3 <- ergebnisse$eta.sample5.setting3
scot_BNG$eta.sample6.setting3 <- ergebnisse$eta.sample6.setting3
scot_BNG$eta.sample7.setting3 <- ergebnisse$eta.sample7.setting3
scot_BNG$eta.sample8.setting3 <- ergebnisse$eta.sample8.setting3
scot_BNG$eta.sample9.setting3 <- ergebnisse$eta.sample9.setting3
scot_BNG$eta.sample10.setting3 <- ergebnisse$eta.sample10.setting3

scot_BNG$Gamma.sample1.setting3 <- ergebnisse$Gamma.sample1.setting3
scot_BNG$Gamma.sample2.setting3 <- ergebnisse$Gamma.sample2.setting3
scot_BNG$Gamma.sample3.setting3 <- ergebnisse$Gamma.sample3.setting3
scot_BNG$Gamma.sample4.setting3 <- ergebnisse$Gamma.sample4.setting3
scot_BNG$Gamma.sample5.setting3 <- ergebnisse$Gamma.sample5.setting3
scot_BNG$Gamma.sample6.setting3 <- ergebnisse$Gamma.sample6.setting3
scot_BNG$Gamma.sample7.setting3 <- ergebnisse$Gamma.sample7.setting3
scot_BNG$Gamma.sample8.setting3 <- ergebnisse$Gamma.sample8.setting3
scot_BNG$Gamma.sample9.setting3 <- ergebnisse$Gamma.sample9.setting3
scot_BNG$Gamma.sample10.setting3 <- ergebnisse$Gamma.sample10.setting3

scot_BNG$alpha.sample1.setting3 <- ergebnisse$alpha.sample1.setting3
scot_BNG$alpha.sample2.setting3 <- ergebnisse$alpha.sample2.setting3
scot_BNG$alpha.sample3.setting3 <- ergebnisse$alpha.sample3.setting3
scot_BNG$alpha.sample4.setting3 <- ergebnisse$alpha.sample4.setting3
scot_BNG$alpha.sample5.setting3 <- ergebnisse$alpha.sample5.setting3
scot_BNG$alpha.sample6.setting3 <- ergebnisse$alpha.sample6.setting3
scot_BNG$alpha.sample7.setting3 <- ergebnisse$alpha.sample7.setting3
scot_BNG$alpha.sample8.setting3 <- ergebnisse$alpha.sample8.setting3
scot_BNG$alpha.sample9.setting3 <- ergebnisse$alpha.sample9.setting3
scot_BNG$alpha.sample10.setting3 <- ergebnisse$alpha.sample10.setting3

scot_BNG$sGa.sample1.setting3 <- ergebnisse$sGa.sample1.setting3
scot_BNG$sGa.sample2.setting3 <- ergebnisse$sGa.sample2.setting3
scot_BNG$sGa.sample3.setting3 <- ergebnisse$sGa.sample3.setting3
scot_BNG$sGa.sample4.setting3 <- ergebnisse$sGa.sample4.setting3
scot_BNG$sGa.sample5.setting3 <- ergebnisse$sGa.sample5.setting3
scot_BNG$sGa.sample6.setting3 <- ergebnisse$sGa.sample6.setting3
scot_BNG$sGa.sample7.setting3 <- ergebnisse$sGa.sample7.setting3
scot_BNG$sGa.sample8.setting3 <- ergebnisse$sGa.sample8.setting3
scot_BNG$sGa.sample9.setting3 <- ergebnisse$sGa.sample9.setting3
scot_BNG$sGa.sample10.setting3 <- ergebnisse$sGa.sample10.setting3


#Setting 4

scot_BNG$y.sim1.setting4 <- ergebnisse$y.sim1.setting4
scot_BNG$y.sim2.setting4 <- ergebnisse$y.sim2.setting4
scot_BNG$y.sim3.setting4 <- ergebnisse$y.sim3.setting4
scot_BNG$y.sim4.setting4 <- ergebnisse$y.sim4.setting4
scot_BNG$y.sim5.setting4 <- ergebnisse$y.sim5.setting4
scot_BNG$y.sim6.setting4 <- ergebnisse$y.sim6.setting4
scot_BNG$y.sim7.setting4 <- ergebnisse$y.sim7.setting4
scot_BNG$y.sim8.setting4 <- ergebnisse$y.sim8.setting4
scot_BNG$y.sim9.setting4 <- ergebnisse$y.sim9.setting4
scot_BNG$y.sim10.setting4 <- ergebnisse$y.sim10.setting4

scot_BNG$Gamma.sim1.setting4 <- ergebnisse$Gamma.sim1.setting4
scot_BNG$Gamma.sim2.setting4 <- ergebnisse$Gamma.sim2.setting4
scot_BNG$Gamma.sim3.setting4 <- ergebnisse$Gamma.sim3.setting4
scot_BNG$Gamma.sim4.setting4 <- ergebnisse$Gamma.sim4.setting4
scot_BNG$Gamma.sim5.setting4 <- ergebnisse$Gamma.sim5.setting4
scot_BNG$Gamma.sim6.setting4 <- ergebnisse$Gamma.sim6.setting4
scot_BNG$Gamma.sim7.setting4 <- ergebnisse$Gamma.sim7.setting4
scot_BNG$Gamma.sim8.setting4 <- ergebnisse$Gamma.sim8.setting4
scot_BNG$Gamma.sim9.setting4 <- ergebnisse$Gamma.sim9.setting4
scot_BNG$Gamma.sim10.setting4 <- ergebnisse$Gamma.sim10.setting4

scot_BNG$alpha.sim1.setting4 <- ergebnisse$alpha.sim1.setting4
scot_BNG$alpha.sim2.setting4 <- ergebnisse$alpha.sim2.setting4
scot_BNG$alpha.sim3.setting4 <- ergebnisse$alpha.sim3.setting4
scot_BNG$alpha.sim4.setting4 <- ergebnisse$alpha.sim4.setting4
scot_BNG$alpha.sim5.setting4 <- ergebnisse$alpha.sim5.setting4
scot_BNG$alpha.sim6.setting4 <- ergebnisse$alpha.sim6.setting4
scot_BNG$alpha.sim7.setting4 <- ergebnisse$alpha.sim7.setting4
scot_BNG$alpha.sim8.setting4 <- ergebnisse$alpha.sim8.setting4
scot_BNG$alpha.sim9.setting4 <- ergebnisse$alpha.sim9.setting4
scot_BNG$alpha.sim10.setting4 <- ergebnisse$alpha.sim10.setting4

scot_BNG$sGa.sim1.setting4 <- ergebnisse$sGa.sim1.setting4
scot_BNG$sGa.sim2.setting4 <- ergebnisse$sGa.sim2.setting4
scot_BNG$sGa.sim3.setting4 <- ergebnisse$sGa.sim3.setting4
scot_BNG$sGa.sim4.setting4 <- ergebnisse$sGa.sim4.setting4
scot_BNG$sGa.sim5.setting4 <- ergebnisse$sGa.sim5.setting4
scot_BNG$sGa.sim6.setting4 <- ergebnisse$sGa.sim6.setting4
scot_BNG$sGa.sim7.setting4 <- ergebnisse$sGa.sim7.setting4
scot_BNG$sGa.sim8.setting4 <- ergebnisse$sGa.sim8.setting4
scot_BNG$sGa.sim9.setting4 <- ergebnisse$sGa.sim9.setting4
scot_BNG$sGa.sim10.setting4 <- ergebnisse$sGa.sim10.setting4


scot_BNG$eta.sample1.setting4 <- ergebnisse$eta.sample1.setting4
scot_BNG$eta.sample2.setting4 <- ergebnisse$eta.sample2.setting4
scot_BNG$eta.sample3.setting4 <- ergebnisse$eta.sample3.setting4
scot_BNG$eta.sample4.setting4 <- ergebnisse$eta.sample4.setting4
scot_BNG$eta.sample5.setting4 <- ergebnisse$eta.sample5.setting4
scot_BNG$eta.sample6.setting4 <- ergebnisse$eta.sample6.setting4
scot_BNG$eta.sample7.setting4 <- ergebnisse$eta.sample7.setting4
scot_BNG$eta.sample8.setting4 <- ergebnisse$eta.sample8.setting4
scot_BNG$eta.sample9.setting4 <- ergebnisse$eta.sample9.setting4
scot_BNG$eta.sample10.setting4 <- ergebnisse$eta.sample10.setting4

scot_BNG$Gamma.sample1.setting4 <- ergebnisse$Gamma.sample1.setting4
scot_BNG$Gamma.sample2.setting4 <- ergebnisse$Gamma.sample2.setting4
scot_BNG$Gamma.sample3.setting4 <- ergebnisse$Gamma.sample3.setting4
scot_BNG$Gamma.sample4.setting4 <- ergebnisse$Gamma.sample4.setting4
scot_BNG$Gamma.sample5.setting4 <- ergebnisse$Gamma.sample5.setting4
scot_BNG$Gamma.sample6.setting4 <- ergebnisse$Gamma.sample6.setting4
scot_BNG$Gamma.sample7.setting4 <- ergebnisse$Gamma.sample7.setting4
scot_BNG$Gamma.sample8.setting4 <- ergebnisse$Gamma.sample8.setting4
scot_BNG$Gamma.sample9.setting4 <- ergebnisse$Gamma.sample9.setting4
scot_BNG$Gamma.sample10.setting4 <- ergebnisse$Gamma.sample10.setting4

scot_BNG$alpha.sample1.setting4 <- ergebnisse$alpha.sample1.setting4
scot_BNG$alpha.sample2.setting4 <- ergebnisse$alpha.sample2.setting4
scot_BNG$alpha.sample3.setting4 <- ergebnisse$alpha.sample3.setting4
scot_BNG$alpha.sample4.setting4 <- ergebnisse$alpha.sample4.setting4
scot_BNG$alpha.sample5.setting4 <- ergebnisse$alpha.sample5.setting4
scot_BNG$alpha.sample6.setting4 <- ergebnisse$alpha.sample6.setting4
scot_BNG$alpha.sample7.setting4 <- ergebnisse$alpha.sample7.setting4
scot_BNG$alpha.sample8.setting4 <- ergebnisse$alpha.sample8.setting4
scot_BNG$alpha.sample9.setting4 <- ergebnisse$alpha.sample9.setting4
scot_BNG$alpha.sample10.setting4 <- ergebnisse$alpha.sample10.setting4

scot_BNG$sGa.sample1.setting4 <- ergebnisse$sGa.sample1.setting4
scot_BNG$sGa.sample2.setting4 <- ergebnisse$sGa.sample2.setting4
scot_BNG$sGa.sample3.setting4 <- ergebnisse$sGa.sample3.setting4
scot_BNG$sGa.sample4.setting4 <- ergebnisse$sGa.sample4.setting4
scot_BNG$sGa.sample5.setting4 <- ergebnisse$sGa.sample5.setting4
scot_BNG$sGa.sample6.setting4 <- ergebnisse$sGa.sample6.setting4
scot_BNG$sGa.sample7.setting4 <- ergebnisse$sGa.sample7.setting4
scot_BNG$sGa.sample8.setting4 <- ergebnisse$sGa.sample8.setting4
scot_BNG$sGa.sample9.setting4 <- ergebnisse$sGa.sample9.setting4
scot_BNG$sGa.sample10.setting4 <- ergebnisse$sGa.sample10.setting4

summary(scot_BNG)



# Effekte plotten

tps <- list(fontsize=list(text=18))
trellis.par.set(tps)


redbluepal <- colorRampPalette(c("firebrick", "white", "midnightblue"))


for(j in 1:4){
  for(i in 1:10){

    # simulierte Daten
    pdf(paste("Plots/y_sim", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    print(spplot(scot_BNG, c(paste("y.sim", i,".setting", j, sep="")), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off()

    
    # standardisierte strukturierte Effekte
    pdf(paste("Plots/Gamma_sim", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("Gamma.sim", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,2)
    print(spplot(scot_BNG, c(paste("Gamma.sim", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off()
    
    # zentrierte strukturierte mittlere Effekte
    pdf(paste("Plots/Gamma_sample", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("Gamma.sample", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,2)
    print(spplot(scot_BNG, c(paste("Gamma.sample", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off() 
        
  
    # standardisierte unstrukturierte Effekte
    pdf(paste("Plots/alpha_sim", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("alpha.sim", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,3)
    print(spplot(scot_BNG, c(paste("alpha.sim", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off()
        
    # zentrierte unstrukturierte mittlere Effekte
    pdf(paste("Plots/alpha_sample", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("alpha.sample", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,3)
    print(spplot(scot_BNG, c(paste("alpha.sample", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off() 
    
        
    # summierte standardisierte Effekte
    pdf(paste("Plots/sGa_sim", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("sGa.sim", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,2)
    print(spplot(scot_BNG, c(paste("sGa.sim", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off() 
    
    # summierte zentrierte mittlere Effekte
    pdf(paste("Plots/sGa_sample", i, "_setting", j, ".pdf", sep=""), width=8.5, height=11) 
    ogrenze <- max((ceiling(sort(abs(get(paste("sGa.sample", i, ".setting", j, sep=""))), decreasing = TRUE)[4]*10))/10)
    cuts <- seq(-ogrenze, ogrenze, length.out=17)
    cuts <- round(cuts,2)
    print(spplot(scot_BNG, c(paste("sGa.sample", i, ".setting", j, sep="")), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot, auto.key=T, par.settings=list(fontsize=list(text=23))))
    dev.off() 
  
  }
}






#############################################################################
## Kartenaufbereitung für Schätzergebnisse - empirische Daten mit Gewichte ##
#############################################################################

########## Achtung!!! ##########
# Beim plotten der Ergebnisse aufpassen auf ID Zuordnung!!!

# Rücktransformation der Reihenfolge der Ergebnisse

# Reihenfolge der Originalindizes
scotRedIndizes <- c(19,3,31,13,34,26,17,55,32,27,30,36,35,45,23,54,44,46,40,43,52,47,29,37,48,53,41,
                    51,50,39,21,9,1,42,20,25,24,49,14,33,38,15,6,28,22,16,12,8,18,0,4,11,2,10,7,5)

# alle Werte +1, da bei Index 0 gestartet wird
scotRedIndizes <- scotRedIndizes + 1


load("Mit Gewichte/Ga.sample.mean.RData")
load("Mit Gewichte/Ga.sample.mean.zent.RData")
load("Mit Gewichte/Ga.sample.sign.RData")
load("Mit Gewichte/alpha.sample.mean.RData")
load("Mit Gewichte/alpha.sample.mean.zent.RData")
load("Mit Gewichte/alpha.sample.sign.RData")
load("Mit Gewichte/eta.sample.mean.RData")



# mit Inseln auf 1000, sind dann weiß. da außerhalb des Intervals
Ga.sample.mean <- c(Ga.sample.mean,1000,1000,1000)
Ga.sample.mean.zent <- c(Ga.sample.mean.zent,1000,1000,1000)
Ga.sample.sign <- c(Ga.sample.sign,0,0,0)
alpha.sample.mean <- c(alpha.sample.mean,1000,1000,1000)
alpha.sample.mean.zent <- c(alpha.sample.mean.zent,1000,1000,1000)
alpha.sample.sign <- c(alpha.sample.sign,0,0,0)
eta.sample.mean <- c(eta.sample.mean,1000,1000,1000)


# Stelle alte Reihenfolge vor reduzierter Form wieder her
ergebnisse <- data.frame(scotRedIndizes, Ga.sample.mean, Ga.sample.mean.zent, Ga.sample.sign,
                         alpha.sample.mean, alpha.sample.mean.zent, alpha.sample.sign,
                         eta.sample.mean)

ergebnisse <- ergebnisse[order(ergebnisse$scotRedIndizes),]

# Berücksichtige Sortierung nach der Variable ID
wandler <- 0
for(i in 1:length(scot_BNG$ID)){
  wandler[i] <- which(scot_BNG$ID == i)
}
wandler

ergebnisse <- cbind(wandler, ergebnisse)
ergebnisse <- ergebnisse[order(ergebnisse$wandler),]

ergebnisse

scot_BNG$Ga.sample.mean <- ergebnisse$Ga.sample.mean
scot_BNG$Ga.sample.mean.zent <- ergebnisse$Ga.sample.mean.zent
scot_BNG$Ga.sample.sign <- ergebnisse$Ga.sample.sign
scot_BNG$alpha.sample.mean <- ergebnisse$alpha.sample.mean
scot_BNG$alpha.sample.mean.zent <- ergebnisse$alpha.sample.mean.zent
scot_BNG$alpha.sample.sign <- ergebnisse$alpha.sample.sign
scot_BNG$eta.sample.mean <- ergebnisse$eta.sample.mean
scot_BNG$exp.eta.sample.mean <- exp(ergebnisse$eta.sample.mean)
scot_BNG$logObserved2 <- scot_BNG$logObserved
scot_BNG$logObserved2[c(47,48,52)] <- 0
scot_BNG$sum.Ga.alpha.sample.mean.zent <- scot_BNG$Ga.sample.mean.zent + scot_BNG$alpha.sample.mean.zent


summary(scot_BNG)


# Plots ohne Gewichtsadaption

tps <- list(fontsize=list(text=16))
trellis.par.set(tps)

redbluepal <- colorRampPalette(c("firebrick", "white", "midnightblue"))


# ursprüngliche strukturierte Effekte (nicht zentriert)
ogrenze <- max((ceiling(sort(abs(Ga.sample.mean), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("Ga.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)


# ursprüngliche unstrukturierte Effekte (nicht zentriert)
ogrenze <- max((ceiling(sort(abs(alpha.sample.mean), decreasing = TRUE)[4]*100))/100)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,3)
spplot(scot_BNG, c("alpha.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)


# zentrierte strukturierte Effekte
ogrenze <- max((ceiling(sort(abs(Ga.sample.mean.zent), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("Ga.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

# zentrierte unstrukturierte Effekte
ogrenze <- max((ceiling(sort(abs(alpha.sample.mean.zent), decreasing = TRUE)[4]*100))/100)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,3)
spplot(scot_BNG, c("alpha.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

# zentrierte summierte räumliche Effekte
ogrenze <- max((ceiling(sort(abs(scot_BNG$sum.Ga.alpha.sample.mean.zent), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("sum.Ga.alpha.sample.mean.zent"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)




# eta mean
ogrenze <- max((ceiling(sort(abs(scot_BNG$eta.sample.mean), decreasing = TRUE)[4]*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
spplot(scot_BNG, c("eta.sample.mean"), col.regions=redbluepal(16), colorkey = list(labels = list(at = cuts)), at=cuts, xlim = xlim.scot, ylim = ylim.scot)

spplot(scot_BNG, c("logObserved2"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)

# exp eta mean
spplot(scot_BNG, c("exp.eta.sample.mean"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)
spplot(scot_BNG, c("Observed"), col.regions=bluepal(16), xlim = xlim.scot, ylim = ylim.scot)




# Signifikanzen
scot_BNG$Ga.sample.sign
scot_BNG$Ga.sample.sign[1] <- -1
scot_BNG$Ga.sample.sign <- factor(scot_BNG$Ga.sample.sign, labels=c("neg", "0", "pos"))
scot_BNG$Ga.sample.sign[1] <- 0
scot_BNG$Ga.sample.sign[c(47,48,52)] <- 1000 # Inseln raus

spplot(scot_BNG, "Ga.sample.sign", col="black", col.regions=c("firebrick", "gray90", "midnightblue"), xlim = xlim.scot, ylim = ylim.scot)


scot_BNG$alpha.sample.sign
scot_BNG$alpha.sample.sign[1] <- -1
scot_BNG$alpha.sample.sign[2] <- 1
scot_BNG$alpha.sample.sign <- factor(scot_BNG$alpha.sample.sign, labels=c("neg", "0", "pos"))
scot_BNG$alpha.sample.sign[1:2] <- 0
scot_BNG$alpha.sample.sign[c(47,48,52)] <- 1000 # Inseln raus


spplot(scot_BNG, "alpha.sample.sign", col="black", col.regions=c("firebrick", "gray90", "midnightblue"), xlim = xlim.scot, ylim = ylim.scot)

