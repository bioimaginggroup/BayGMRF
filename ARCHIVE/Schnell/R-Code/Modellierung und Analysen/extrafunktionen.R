###############################################################################
## 1) Funktion zur Bestimmung der Nachbarschaftsmatrix aus einer Graph-Datei ##
###############################################################################

readandchangegraph <- function(file){
  
  # Lese Information aus der Graph-Datei ein
  mapinfo <- scan(file)
  
  # Bestimme die Anzahl der Regionen über das erste Element
  S <- mapinfo[1]  # z.B. 544 Regionen
  mapinfo <- mapinfo[-1] # Entfernen dieser Info aus dem zwischengespeicherten Datensatz
  
  # Lege einen Dummy-Vektor für die Nummerierung der Regionen an
  regions <- rep(0,S)
  
  # Lege eine Dummy-Matrix an 
  pmat <- matrix(0,S,S)
  
  # Schleife über die Regionen für Zuordnung der Nachbarschaftsstruktur
  for(i in 1:S){
    # Extrahiere tatsächliche Nummerierung der i-ten Region
    regions[i] <- mapinfo[1]
    mapinfo <- mapinfo[-1] # Entfernen dieser Info aus dem zwischengespeicherten Datensatz
    
    # Extrahiere Anzahl der Nachbarn der i-ten Region
    nrn <- mapinfo[1]
    if(nrn == 0){ # Ist die Region eine Insel, bzw. hat eine Region keine Nachbarn, bleiben alle Eintrage der zugehörigen Zeile der Matrix bei 0
      mapinfo <- mapinfo[-1] # Entfernen dieser Info aus dem zwischengespeicherten Datensatz
      pmat[i,i] <- nrn # Zuordnung dieser Info zum Diagonalelement der Matrix
      
    } else{ # Ist die Anzahl der Nachbarn einer Region ungleich 0, dann werden die zughörigen Einträge in der Zeile gemacht. 
      mapinfo <- mapinfo[-1] # Entfernen dieser Info aus dem zwischengespeicherten Datensatz
      pmat[i,i] <- nrn # Zuordnung dieser Info zum Diagonalelement der Matrix
      
      # Bestimme Einträge neben der Diagonalen mit -1
      
      # Für alle Nachbarn (1:nrn) mit nrn entpricht der Anzahl an Nachbarn, 
      # d.h. die nächsten nrn Informationen aus dem Datensatz beinhalten die Nummerierung der jeweiligen Nachbarn      
      for(j in 1:nrn){ 
        # An der jeweiligen Stelle der Nummerierung der einzelnen Nachbarn wird eine -1 eingefügt
        pmat[i,mapinfo[j]] <- -1 # mapinfo[1] ist erster Nachbar, mapinfo[2] ist zweiter Nachbar, usw.
      }
      mapinfo <- mapinfo[-(1:nrn)] # Zeile in der die Nachbar standen entfernen
    }
  }
  
  # Ergebnis in einer Liste zurückliefern
  erg <- list(pmat=pmat, regions=regions)
  return(erg)
}



############################################################
# 2) Funktion zur Bestimmung einer knappen (sparse) Matrix ##
############################################################

knappeMatrix <- function(M){
  
  zeilenindizes.l <- c()
  spaltenindizes.k <- c()
  eintrag.x <- c()
  
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      if(M[i,j] != 0){
        zeilenindizes.l <- c(zeilenindizes.l, i)
        spaltenindizes.k <- c(spaltenindizes.k, j)
        eintrag.x <- c(eintrag.x, M[i,j])
      }
    }
  }
  
  
  # l und k sind Vektoren für die Position der Einträge ungleich Null (Zeilen- und Spaltenindizes)
  # x sind die Werte ungleich Null
  # komprimierte Matrix über die Informationen ungleich 0
  M <- sparseMatrix(i = zeilenindizes.l, j = spaltenindizes.k, x = eintrag.x, dims = dim(M))  
  return(M)
}



##########################################################################
## 3) Funktion zur Bestimmung der Indizes und der Einträge einer Matrix ##
##########################################################################

MatrizenIndizes <- function(M){
  
  Indizes <- list()
  
  zeilenindizes.l <- c()
  spaltenindizes.k <- c()
  eintrag.x <- c()
  
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      if(M[i,j] != 0){
        zeilenindizes.l <- c(zeilenindizes.l, i)
        spaltenindizes.k <- c(spaltenindizes.k, j)
        eintrag.x <- c(eintrag.x, M[i,j])
      }
    }
  }
  
  
  # l und k sind Vektoren für die Position der Einträge ungleich Null (Zeilen- und Spaltenindizes)
  # x sind die Werte ungleich Null
  M <- sparseMatrix(i = zeilenindizes.l, j = spaltenindizes.k, x = eintrag.x, dims = dim(M))  # komprimierte Matrix über die Informationen ungleich 0
  Indizes[[1]] <- zeilenindizes.l
  Indizes[[2]] <- spaltenindizes.k
  Indizes[[3]] <- eintrag.x
  names(Indizes)[1] <- "Zeilenindizes"
  names(Indizes)[2] <- "Spaltenindizes"
  names(Indizes)[3] <- "Einträge"
  return(Indizes)
}



#####################################################################
## 4) Funktion zum Ziehen aus einer multivariaten Normalverteilung ##
#####################################################################

rmvnorm_RBA <- function(n, b, P){
  L <- chol(P) # 1. Schritt: Berechne Cholesky-Zerlegung P = LL^T
  u <- solve(t(L), b) # 2. Schritt: Löse Lu = b nach u auf mit u = bL^T
  v <- solve(L, u) # 3. Schritt: Löse L^T u = v nach v auf mit v = uL
  z <- rnorm(n = length(b), mean = 0, sd = 1) # 4. Schritt: Ziehe z als unabhängige standardnormalverteilte Zufallsgröße mit Länge d
  m <- solve(L, z) # 5. Schritt: Löse L^T m = z nach m auf mit m = zL
  Gamma <- v + m # 6. Schritt Berechne die Gamma's als Summe von mu und y
  return(Gamma)
}



##########################################################################
## 5) Funktion für Bestandteile einer "sparseMatrix" in geordneter Form ##
##########################################################################

sortierte.indizes <- function(M){
  
  Sortierte.Indizes <- list()
  
  # Abspeichern der Diagonaleinträge der Ausgangsnachbarschaftsmatrix
  diagonale.M <- diag(M)
  
  # obere Dreiecksmatrix
  oD_M <- matrix(data=0, nrow=dim(M)[1], ncol=dim(M)[2])
  for(r in 1:dim(M)[1]){
    for(c in 1:dim(M)[2]){
      if(r < c){
        oD_M[r,c] <- M[r,c]
      }
    }
  }
  
  
  Indizes.oD_M <- MatrizenIndizes(oD_M)
  
  # Abspeichern der Zeilenindizes der oberen Dreiecksmatrix ohne Diagonale, für die gesampelt werden muss!
  zeilenindizes.oD_M <- Indizes.oD_M[[1]] 
  # Abspeichern der Spaltenindizes der oberen Dreiecksmatrix ohne Diagonale, für die gesampelt werden muss!
  spaltenindizes.oD_M <- Indizes.oD_M[[2]]
  # Abspeichern der Einträge ungleich 0 der oberen Dreiecksmatrix ohne Diagonale, für die gesampelt werden muss!
  eintraege.oD_M <- Indizes.oD_M[[3]]
  
  
  zeilenindizes.M.ordered <- c(zeilenindizes.oD_M, spaltenindizes.oD_M, 1:dim(M)[1])
  spaltenindizes.M.ordered <- c(spaltenindizes.oD_M, zeilenindizes.oD_M, 1:dim(M)[1])
  eintraege.M.ordered <- c(eintraege.oD_M, eintraege.oD_M, diagonale.M)
  
  
  Sortierte.Indizes[[1]] <- zeilenindizes.M.ordered
  Sortierte.Indizes[[2]] <- spaltenindizes.M.ordered
  Sortierte.Indizes[[3]] <- eintraege.M.ordered
  names(Sortierte.Indizes)[1] <- "Zeilenindizes.ordered"
  names(Sortierte.Indizes)[2] <- "Spaltenindizes.ordered"
  names(Sortierte.Indizes)[3] <- "Einträge.ordered"
  return(Sortierte.Indizes)
  
}
