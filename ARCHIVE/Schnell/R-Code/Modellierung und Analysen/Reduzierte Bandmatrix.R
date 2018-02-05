##########################################################################################
### Bestimmung der Nachbarsschaftsmatrix mit reduzierter Bandweite aus der Graph-Datei ###
##########################################################################################

library(Matrix) # wichtig für das erstellen von knappen Matrizen
library(grid) # Wird benötigt für das Zeichnen der Matrizen

source("readandchangegraph.R") # Funktion zum Einlesen einer .graph Datei (Nachbarschaftsmatrix)
source("MatrizenIndizes.R") # Gibt Zeilen- und Spaltenindizes mit Einträgen zurück
source("knappeMatrix.R") # Umformung zu einer sparse Matrix

### Anmerkung: Der Reverse Cuthill Sampler stellt dann einen Vorteil da, wenn gaussian bla bla verwendet wird!!!



###########################################################################
# Bestimme Indexreihenfolge für reduzierte Bandweite für Deutschlanddaten #
###########################################################################

# Laden der deutschen Nachbarschaftsmatrix
germanygraph <- readandchangegraph(file="germany2.graph")
germanyQ <- germanygraph$pmat

# Abspeichern der ursprünglichen Zeilen- und Spaltenindizes
germanyQ.urIndizes <- MatrizenIndizes(germanyQ)
germanyQ.urZeilenIndizes <- germanyQ.urIndizes$Zeilenindizes
germanyQ.urSpaltenIndizes <- germanyQ.urIndizes$Spaltenindizes

# Abspeichern der Diagonaleinträge in der ursprünglichen Reihenfolge
germanyQ.urDiago <- diag(germanyQ)

# Aufbereitung der Nachbarschaftsmatrix für Bandweitenreduzierungsprogramm
germanyQ2 <- germanyQ*(-1) # alle Werte positiv (außerhalb der Diagonale)
diag(germanyQ2) <- 1 # Diagonale mit 1er besetzen
germanyQ2

# Abspeichern der deutschen Nachbarschaftsmatrix
write.table(germanyQ2, "germanyQ.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

# In die txt Datei wird händisch in die erste Zeile die Anzahl der Regionen 544 eingetragen (benötigt das Programm)

# Programm gibt folgende neue Reihenfolge der Indizes aus:
# Reverse Cuthill McKee Algorithmus mit Bandweite 59 (ursprünglich 522):
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


##########################################################################
# Bestimme Indexreihenfolge für reduzierte Bandweite für Schottlanddaten #
##########################################################################

# Laden der schottischen Nachbarschaftsmatrix
scotgraph <- readandchangegraph(file="scotland.graph")
scotQ <- scotgraph$pmat

# Abspeichern der ursprünglichen Zeilen- und Spaltenindizes
scotQ.urIndizes <- MatrizenIndizes(scotQ)
scotQ.urZeilenIndizes <- scotQ.urIndizes$Zeilenindizes
scotQ.urSpaltenIndizes <- scotQ.urIndizes$Spaltenindizes

# Abspeichern der Diagonaleinträge in der ursprünglichen Reihenfolge
scotQ.urDiago <- diag(scotQ)

# Aufbereitung der Nachbarschaftsmatrix für Bandweitenreduzierungsprogramm
scotQ2 <- scotQ*(-1) # alle Werte positiv (außerhalb der Diagonale)
diag(scotQ2) <- 1 # Diagonale mit 1er besetzen
scotQ2

# Abspeichern der schottischen Nachbarschaftsmatrix
write.table(scotQ2, "scotQ.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

# In die txt Datei wird händisch in die erste Zeile die Anzahl der Regionen 56 eingetragen (benötigt das Programm)

# Programm gibt folgende neue Reihenfolge der Indizes aus:
# Reverse Cuthill McKee Algorithmus mit Bandweite 14 (ursprünglich 38):
scotRedIndizes <- c(19,3,31,13,34,26,17,55,32,27,30,36,35,45,23,54,44,46,40,43,52,47,29,37,48,53,41,
                    51,50,39,21,9,1,42,20,25,24,49,14,33,38,15,6,28,22,16,12,8,18,0,4,11,2,10,7,5)   
                    
# alle Werte +1, da bei Index 0 gestartet wird
scotRedIndizes <- scotRedIndizes + 1


########################################################################################
# Funktion für die Erstellung der reduzierten Matrix anhand der neuen Indexreihenfolge #
########################################################################################

new.reduzierteMatrix <- function(urZeilenIndizes, urSpaltenIndizes, urDiago, redIndizes){

  # "urZeilenIndizes": Zeilen Indizes der ursprünglichen Matrix
  # "urSpaltenIndizes": Zeilen Indizes der ursprünglichen Matrix
  # "urDiago": Diagonaleinträge der ursprünglichen Matrix
  # "redIndizes": neu geordnete Indizes für reduzierte Matrix


  # Neue, reduzierte Ordnung der Zeilen- und der Spaltenindizes
  redZeilenIndizes <- vector()
  redSpaltenIndizes <- vector()

  for(j in 1:length(redIndizes)){
    redZeilenIndizes[which(urZeilenIndizes == redIndizes[j])] <- j
    redSpaltenIndizes[which(urSpaltenIndizes == redIndizes[j])] <- j
  }

  
  # Neue, reduzierte ordnung der Diagonaleinträge
  redDiago <- vector()
  
  for(j in 1:length(redIndizes)){
    redDiago[j] <- urDiago[redIndizes[j]]
  }
  
  
  # Bestimme Einträge
  redEintraege <- vector()
  
  for(j in 1:length(redZeilenIndizes)){
    if(redZeilenIndizes[j] != redSpaltenIndizes[j]){
      redEintraege[j] <- -1
    } else if(redZeilenIndizes[j] == redSpaltenIndizes[j]){
      redEintraege[j] <- redDiago[redZeilenIndizes[j]]
    }
  }
  
  
  # Reduzierte Matrix ausgeben
  redQ <- sparseMatrix(i = redZeilenIndizes, j = redSpaltenIndizes, x = redEintraege, dims = c(length(redDiago), length(redDiago)))  
  return(redQ)
}


########################
# Plotten der Matrizen #
########################

# Funktion die Matrix aufbereitet für die Funktion "grid.raster"

hilfe.grid.raster <- function(M){
  M.grid.raster <- matrix(data=NA, nrow=dim(M)[1], ncol=dim(M)[2])
  
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      if(M[i,j] != 0){
        M.grid.raster[i,j] <- 0 
      } else if(M[i,j] == 0){
        M.grid.raster[i,j] <- 1
      }
    }
  }
  
  return(M.grid.raster)
}


# Deutsche Nachbarschaftsmatrix:

# Ursprüngliche Reihenfolge
sparse.Indizes.germanyQ <- knappeMatrix(germanyQ) 

# neue, reduzierte Reihenfolge
red.germanyQ <- new.reduzierteMatrix(urZeilenIndizes=germanyQ.urZeilenIndizes, urSpaltenIndizes=germanyQ.urSpaltenIndizes,
                                     urDiago=germanyQ.urDiago, redIndizes=germanyRedIndizes)

# Umformung von einer knappen Matrix zu einer vollen Matrix
red.germanyQvoll <- as.matrix(red.germanyQ + 1 - 1)

# Abspeichern der reduzierten deutschen Nachbarschaftsmatrix
write.table(red.germanyQvoll, "redgermanyQ.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Anwendung der Hilfe-Funktion
germanyQ.grid.raster <- hilfe.grid.raster(germanyQ)
red.germanyQvoll.grid.raster <- hilfe.grid.raster(red.germanyQvoll)

# Zeichnen der Matrizen
grid.raster(germanyQ.grid.raster, interp=FALSE)
grid.raster(red.germanyQvoll.grid.raster, interp=FALSE)



# Schottische Nachbarschaftsmatrix:

# Ursprüngliche Reihenfolge
sparse.Indizes.scotQ <- knappeMatrix(scotQ) 

# neue, reduzierte Reihenfolge
red.scotQ <- new.reduzierteMatrix(urZeilenIndizes=scotQ.urZeilenIndizes, urSpaltenIndizes=scotQ.urSpaltenIndizes,
                                  urDiago=scotQ.urDiago, redIndizes=scotRedIndizes)

# Umformung von einer knappen Matrix zu einer vollen Matrix
red.scotQvoll <- as.matrix(red.scotQ + 1 - 1)

# Abspeichern der reduzierten deutschen Nachbarschaftsmatrix
write.table(red.scotQvoll, "redscotQ.txt", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Anwendung der Hilfe-Funktion
scotQ.grid.raster <- hilfe.grid.raster(scotQ)
red.scotQvoll.grid.raster <- hilfe.grid.raster(red.scotQvoll)

# Zeichnen der Matrizen
grid.raster(scotQ.grid.raster, interp=FALSE)
grid.raster(red.scotQvoll.grid.raster, interp=FALSE)