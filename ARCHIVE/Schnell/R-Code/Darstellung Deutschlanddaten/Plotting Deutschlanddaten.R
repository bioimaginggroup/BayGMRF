############################################################################################
############################################################################################
############################################################################################
################################## Deutschlanddaten Plots ##################################
############################################################################################
############################################################################################
############################################################################################



##############################################################
## Einlesen der deutschen Karten- und Dateninformationen ##
##############################################################

# Quelle: http://www.r-inla.org/examples/volume-1/code-for-bym-example, Datum: 27.05.2013
# source("http://www.math.ntnu.no/inla/givemeINLA.R")


library(INLA)
library(classInt)


# Lade Germany Datensatz
data.germany <- read.table("Germany.dat", skip=1)

# Übergebe Variablennamen
names(data.germany) <- c("Region", "Expected", "Observed", "Kov.Smoke")

# Zähler der Region um 1 rauf gesetzt, da Start bei 0
data.germany$Region <- data.germany$Region + 1

# Berechne standardisierte Sterblichkeitsrate
data.germany$SMR <- data.germany$Observed/data.germany$Expected 


# Lege logarithmierte Variablen an

# Berechne logarithmierte Variable Observed. 
# Problem Nuller im Datensatz! Deswegen For-Schleife
data.germany$logObserved <- 0

for(i in 1:length(data.germany$Observed)){
  if(data.germany$Observed[i] == 0){
    data.germany$logObserved[i] <- 0
  }else{
    data.germany$logObserved[i] <- log(data.germany$Observed[i])
  }
}

data.germany$logExpected <- log(data.germany$Expected)




# data.germany$Region <- as.factor(data.germany$Region)
graph.germany <- system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))


# eine doppelte Spalte der Variable "Region" anlegen
data.germany <- cbind(data.germany, Region.struct=data.germany$Region)

summary(data.germany)



########################################
## Deskriptive Aufbereitung der Daten ##
########################################

# Lade Funktion zum Plotten der Deutschlandkarte. Übernommen und angepasst aus INLA.
source("germany.map.new.R")

# Plot der Landkreisumrisse
white <- rep.int(0, 544)
germany.map.new(data.germany$SMR, legend=F, axes=F, a=200, b=200, farbe="white")


# Plot der Daten
germany.map.new(data.germany$Observed, legend=T, axes=T)
germany.map.new(data.germany$logObserved, legend=T, axes=T)

germany.map.new(data.germany$Expected, legend=T, axes=T)
germany.map.new(data.germany$logExpected, legend=T, axes=T)

germany.map.new(data.germany$Kov.Smoke, legend=T, axes=T)
germany.map.new(data.germany$SMR, legend=T, axes=T)


# Einfärbung einzelner Distrikte (wichtige Städte)

source("germany.map.new.help.R")

data.germany$Stadt <- rep(0, length(data.germany$Observed))
# Städte: 2=Kiel, 16=Hamburg, 28=Hannover, 64=Bremen, 66=Düsseldorf, 
#         182=Stuttgart, 123=Wiesbaden, 166=Mainz, 227=München, 322=Saarbrücken,
#         328=Westberlin, 329=Ostberlin, 334=Potsdam, 377=Schwerin
#         412=Dresden, 487=Magdeburg, 505=Erfurt
data.germany$Stadt[c(2,16,28,64,66,182,123,166,227,322,328,329,334,377,412,487,505)] <- 100
germany.map.new.help(data.germany$Stadt, legend=F, axes=F, a=230, b=230)

data.germany$Stadt <- rep(0, length(data.germany$Observed))
# Städte: 2=Kiel, 16=Hamburg, 28=Hannover, 64=Bremen, 66=Düsseldorf, 
#         182=Stuttgart, 123=Wiesbaden, 166=Mainz, 227=München, 322=Saarbrücken,
#         328=Westberlin, 329=Ostberlin, 334=Potsdam, 377=Schwerin
#         412=Dresden, 487=Magdeburg, 505=Erfurt
data.germany$Stadt[c(1,2)] <- 100
germany.map.new.help(data.germany$Stadt, legend=F, axes=F, a=230, b=230)




#############################################
## Kartenaufbereitung für Schätzergebnisse ##
#############################################


########## Achtung!!! ########## 
# Beim plotten der Ergebnisse aufpassen! Wieder in ursprüngliche Reihenfolge bringen!!!


# Rücktransformation der Reihenfolge der Ergebnisse

# Reihenfolge der Originalindizes
germanyRedIndizes <- c(229,249,214,227,246,255,228,215,207,206,252,316,309,244,232,239,230,251,250,240,257,259,248,234,
                       222,213,212,211,210,208,256,223,310,320,237,308,241,226,247,245,258,261,254,235,224,218,221,315,
                       219,217,205,209,194,264,454,441,318,317,236,231,238,307,267,243,253,225,216,314,220,202,182,203,
                       197,195,260,268,262,167,170,162,321,461,412,416,444,425,422,417,273,418,447,414,312,311,242,313,
                       265,287,233,183,181,192,184,285,204,185,187,174,196,168,161,164,198,263,266,282,269,271,286,284,
                       291,180,177,163,166,325,175,323,324,322,156,443,434,449,455,411,431,462,456,410,429,430,280,438,
                       463,415,451,278,450,319,293,294,283,186,193,290,188,178,201,199,275,289,277,270,281,272,172,169,
                       173,176,326,160,420,419,446,427,445,423,437,279,433,459,519,505,543,525,533,288,190,189,292,297,
                       200,129,123,274,302,296,276,538,158,171,165,157,148,149,122,458,440,439,436,452,426,453,432,80,
                       528,532,535,509,517,530,539,191,306,119,303,124,295,125,301,305,507,522,521,159,150,179,153,128,
                       131,347,330,368,367,331,448,337,435,421,460,413,85,84,146,87,86,81,82,73,83,531,506,508,485,515,
                       477,523,504,484,304,130,298,121,120,300,299,540,511,526,534,151,154,145,126,135,152,350,343,342,
                       346,345,332,366,334,457,428,351,468,442,424,473,409,373,400,65,77,88,69,91,90,74,72,76,70,541,479,
                       510,516,478,480,127,132,139,518,513,133,155,134,147,338,355,356,348,369,336,344,335,396,466,467,
                       464,471,352,482,406,380,385,383,377,66,71,68,78,112,89,67,107,92,93,111,108,109,536,524,481,512,
                       475,137,140,514,527,136,116,115,353,372,357,333,328,340,370,404,363,469,470,339,476,465,393,379,
                       374,392,387,403,79,75,114,97,118,537,483,142,144,520,138,542,143,113,327,362,360,349,394,498,472,
                       329,341,496,474,499,486,501,398,407,405,397,95,117,96,110,94,529,502,141,20,24,103,106,100,358,
                       359,371,487,365,489,490,494,503,491,56,391,388,399,375,98,101,99,49,21,26,23,17,32,102,104,354,
                       364,500,492,22,16,18,488,493,54,47,52,395,376,402,382,381,378,59,31,25,29,27,105,34,361,495,19,
                       497,60,53,57,62,50,390,401,408,28,30,33,35,38,44,58,51,55,46,48,63,389,384,386,45,42,39,61,40,64,
                       5,2,41,37,36,15,14,7,43,8,12,9,3,1,13,10,4,6,11,0)

# alle Werte +1, da bei Index 0 gestartet wird
germanyRedIndizes <- germanyRedIndizes + 1

# load("ohne Gewichte/Ga.sample.mean.RData")
# load("ohne Gewichte/Ga.sample.mean.zent.RData")
# load("ohne Gewichte/Ga.sample.sign.RData")
# load("ohne Gewichte/alpha.sample.mean.RData")
# load("ohne Gewichte/alpha.sample.mean.zent.RData")
# load("ohne Gewichte/alpha.sample.sign.RData")
# load("ohne Gewichte/eta.sample.mean.RData")


load("mit Gewichte/Ga.sample.mean.RData")
load("mit Gewichte/Ga.sample.mean.zent.RData")
load("mit Gewichte/Ga.sample.sign.RData")
load("mit Gewichte/alpha.sample.mean.RData")
load("mit Gewichte/alpha.sample.mean.zent.RData")
load("mit Gewichte/alpha.sample.sign.RData")
load("mit Gewichte/eta.sample.mean.RData")


# Stelle alte Reihenfolge vor reduzierter Form wieder her
ergebnisse <- data.frame(germanyRedIndizes, Ga.sample.mean, Ga.sample.mean.zent, Ga.sample.sign,
                         alpha.sample.mean, alpha.sample.mean.zent, alpha.sample.sign,
                         eta.sample.mean)

ergebnisse <- ergebnisse[order(ergebnisse$germanyRedIndizes),]

data.germany$Ga.sample.mean <- ergebnisse$Ga.sample.mean
data.germany$Ga.sample.mean.zent <- ergebnisse$Ga.sample.mean.zent
data.germany$Ga.sample.sign <- ergebnisse$Ga.sample.sign
data.germany$alpha.sample.mean <- ergebnisse$alpha.sample.mean
data.germany$alpha.sample.mean.zent <- ergebnisse$alpha.sample.mean.zent
data.germany$alpha.sample.sign <- ergebnisse$alpha.sample.sign
data.germany$eta.sample.mean <- ergebnisse$eta.sample.mean
data.germany$exp.eta.sample.mean <- exp(ergebnisse$eta.sample.mean)
data.germany$sum.Ga.alpha.sample.mean.zent <- data.germany$Ga.sample.mean.zent + data.germany$alpha.sample.mean.zent


summary(data.germany)

source("germany.map.new.ergebnisse.R")
source("germany.map.new.ergebnisse2.R")
source("germany.map.new.ergebnisse3.R")
source("germany.map.new.ergebnisse.sign.R")


# strukturierter Effekt
ogrenze <- max((ceiling(sort(abs(data.germany$Ga.sample.mean), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse(data.germany$Ga.sample.mean, legend=T, axes=T, brks=cuts)

# unstrukturierter Effekt
ogrenze <- max((ceiling(sort(abs(data.germany$alpha.sample.mean), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse(data.germany$alpha.sample.mean, legend=T, axes=T, brks=cuts)


# zentrierter strukturierter Effekt
ogrenze <- max((ceiling(sort(abs(data.germany$Ga.sample.mean.zent), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse(data.germany$Ga.sample.mean.zent, legend=T, axes=T, brks=cuts)

# zentrierter unstrukturierter Effekt
ogrenze <- max((ceiling(sort(abs(data.germany$alpha.sample.mean.zent), decreasing = TRUE)*1000))/1000)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,4)
germany.map.new.ergebnisse3(data.germany$alpha.sample.mean.zent, legend=T, axes=T, brks=cuts)


# summierte zentrierte Effekte
ogrenze <- max((ceiling(sort(abs(data.germany$sum.Ga.alpha.sample.mean.zent), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse(data.germany$sum.Ga.alpha.sample.mean.zent, legend=T, axes=T, brks=cuts)


# Prädiktor
ogrenze <- max((ceiling(sort(abs(data.germany$eta.sample.mean), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse(data.germany$eta.sample.mean, legend=T, axes=T, brks=cuts)

germany.map.new(data.germany$logObserved, legend=T, axes=T)


# exp Prädiktor
ogrenze <- max((ceiling(sort(abs(data.germany$exp.eta.sample.mean), decreasing = TRUE)[5]*10))/10)
cuts <- seq(0,ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse2(data.germany$exp.eta.sample.mean, legend=T, axes=T, brks=cuts)


ogrenze <- max((ceiling(sort(abs(data.germany$Observed), decreasing = TRUE)[3]*10))/10)
cuts <- seq(0,ogrenze, length.out=17)
cuts <- round(cuts,2)
germany.map.new.ergebnisse2(data.germany$Observed, legend=T, axes=T, brks=cuts)



# Signifikanzen strukturierte Effekte
ogrenze <- max((ceiling(sort(abs(data.germany$Ga.sample.sign), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=3)
cuts <- round(cuts,2)
germany.map.new.ergebnisse.sign(data.germany$Ga.sample.sign, legend=T, axes=T, brks=cuts)

# bei Level Problemen Alternativ
Ga.sample.sign.test <- data.germany$Ga.sample.sign
Ga.sample.sign.test <- c(Ga.sample.sign.test, -1)

ogrenze <- max((ceiling(sort(abs(Ga.sample.sign.test), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=3)
cuts <- round(cuts,2)
germany.map.new.ergebnisse.sign(Ga.sample.sign.test, legend=T, axes=T, brks=cuts)


# Signifikanzen unstrukturierte Effekte
ogrenze <- max((ceiling(sort(abs(data.germany$alpha.sample.sign), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=3)
cuts <- round(cuts,2)
germany.map.new.ergebnisse.sign(data.germany$alpha.sample.sign, legend=T, axes=T, brks=cuts)

# bei Level Problemen Alternativ
alpha.sample.sign.test <- data.germany$alpha.sample.sign
alpha.sample.sign.test <- c(alpha.sample.sign.test, -1, 1)

ogrenze <- max((ceiling(sort(abs(alpha.sample.sign.test), decreasing = TRUE)*10))/10)
cuts <- seq(-ogrenze, ogrenze, length.out=3)
cuts <- round(cuts,2)
germany.map.new.ergebnisse.sign(alpha.sample.sign.test, legend=T, axes=T, brks=cuts)



##############################################################
# Plotte auffällige Landkreise bezüglich Gewichtsschätzungen #
##############################################################

hohe.gewichte.redindizes <- c(1,2,4,13,19,23,36,42,53,75,77,78,83,89,107,
                              118,131,153,164,166,184,192,198,216,243,246,
                              252,264,295,304,305,480,484,489,500,543,
                              5,6,14,34,41,43,58,63,73,114,115,122,133,
                              148,157,161,172,188,196,199,221,225,232,254,
                              289,292,300,316,338,348,350,502,503,505,512,544)

germanyRedIndizes[hohe.gewichte.redindizes]



# Einfärbung einzelner Distrikte (Hohe Gewichte)

data.germany$hohe.Gewichte <- rep(0, length(data.germany$Observed))

data.germany$hohe.Gewichte[germanyRedIndizes[hohe.gewichte.redindizes]] <- 100
germany.map.new.help(data.germany$hohe.Gewichte, legend=F, axes=F, a=230, b=230)


