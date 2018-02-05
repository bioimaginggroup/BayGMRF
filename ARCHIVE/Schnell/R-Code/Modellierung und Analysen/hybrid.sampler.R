###########################################################################################
###########################################################################################
###########################################################################################
################################# Hybrid-Sampler Funktion #################################
###########################################################################################
###########################################################################################
###########################################################################################



#################################
## Lade Package und Funktionen ##
#################################

library(Matrix)

source("extrafunktionen.R") # Laden der benötigten Funktionen
# beinhaltet die Funktionen:
# 1) "readandchangegraph" zum Einlesen einer .graph Datei (Nachbarschaftsmatrix)
# 2) "knappeMatrix" zur Umformung einer Matrix zu einer sparse Matrix
# 3) "MatrizenIndizes" zum Ausgeben der besetzten Zeilen- und Spaltenindizes einer Matrix
# 4) "rmvnorm_RBA" zum Ziehen aus einer multivariaten Normalverteilung
# 5) "sortierte.indizes" Funktion für Bestandteile einer "sparseMatrix" in geordneter Form



#############################
## Hybrid-Sampler Funktion ##
#############################

# Funktionseinstellungen:

# size: Iterationsschritte der Sampler
# x: Kovariable
# y: Zielgröße
# e: erwartete Zielgröße
# Q: Nachbarschaftsmatrix
# nu: Parameter der Gewichtsziehung
# burn: Umfang der Burn-In Phase
# thin: Ausdünnungsfaktor
# Qvoll: normale Form der Matrix, keine "sparseMatrix" -> spart Zeit bei den Deutschlanddaten bei Bestimmung der oberen Dreiecksmatrix
# knappeMatrizen: Verwende "sparseMatrix" -> spart Zeit bei Deutschlanddaten
# AdaptiveGewichte: Gewichte werden angepasst
# Intercept: Wird mu mitgeschätzt oder nicht
# IterationPrint: gibt aktuelle Iteration aus
# Hyp.kappa: Wenn TRUE wird kappa normal upgedatet bei jeder Iteration. Bei FALSE ist kappa fest bei 1.
# Hyp.tau: Wenn TRUE wird tau normal upgedatet bei jeder Iteration. Bei FALSE ist tau fest bei 1.
# Koeff: Wenn TRUE wird beta normal upgedatet bei jeder Iteration. Bei FALSE ist beta fest bei 0.
# eta.start.ab: Iteration bis zu der für eta alle Vorschläge akzeptiert werden
# eta.start.logDat: Wenn TRUE, ergibt sich der Startvektor von eta aus den logarithmierten Daten, sonst 1.
# tau.fest.wert: Wird tau nicht geupdated, d.h. Hyp.tau=FALSE, dann kann über diesen Wert der feste Wert bestimmt werden für tau
# kappa.fest.wert: Wird kappa nicht geupdated, d.h. Hyp.kappa=FALSE, dann kann über diesen Wert der feste Wert bestimmt werden für kappa
# GewDiagoVar: Variante der Berechnung der Diagonale in der Gewichts-for-Schleife über "sparseMatrix". Zahlt sich bei Deutschlanddaten aus
# ZS: Zwischenspeicherung alle 5000 Iterationen, wenn ZS=TRUE

hybrid.sampler <- function(size, x, y, e, Q, nu, burn, thin, Qvoll=Q, knappeMatrizen=TRUE, AdaptiveGewichte=TRUE, Intercept=TRUE, IterationPrint=FALSE, Hyp.kappa=TRUE, Hyp.tau=TRUE, Koeff=TRUE, eta.start.ab=0, eta.start.logDat=TRUE,kappa.fest.wert=1, tau.fest.wert=1, GewDiagoVar=FALSE, ZS=FALSE){
  
  # Anzahl an Beobachtungen als Länge für Laufindizes
  d <- length(y)  

  # Bestimme Rang von Q
  rgQ <- d - 1
  
  # Festlegung der Nachbarschaftsmatrix bei AdaptiveGewichte=FALSE. Q_weight.sample.aktuell ist fest.
  if(!AdaptiveGewichte){
    Q_weight.sample.aktuell <- Q
  }

  # Speichere Indizes von der oberen Dreiecksmatrix oD_Q die besetzt sind 
  # für das Sampling der weights jedoch ohne Diagonalelemente:
  # Bilde obere Dreiecksmatrix ohne Diagonale
  if(AdaptiveGewichte){
    oD_Q <- matrix(data=0, nrow=dim(Qvoll)[1], ncol=dim(Qvoll)[2])
    for(r in 1:dim(Qvoll)[1]){
      for(c in 1:dim(Qvoll)[2]){
        if(r < c){
          oD_Q[r,c] <- Qvoll[r,c]
        }
      }
    }
  
  # setze Diagonalelemente 0, da dort die negative Summe der Gewichte der jeweiligen Zeile eingesetzt wird ohne extra Sampling
  diag(oD_Q) <- 0 
  # Abspeichern der Zeilenindizes ohne Diagonale, für die gesampelt werden muss!
  zeilenindizes.oD_Q <- MatrizenIndizes(oD_Q)[[1]] 
  # Abspeichern der Spaltenindizes ohne Diagonale, für die gesampelt werden muss!
  spaltenindizes.oD_Q <- MatrizenIndizes(oD_Q)[[2]]
  
  }
  


  #################################################
  ## Ergebnisvektoren für Samples initialisieren ##
  #################################################

  # Ergebnisvektoren für komplette Samples
  Ga.sample.aktuell <- vector(length = d) # Parametervektor für GMZF (strukturierter Effekt) 
  eta.sample.aktuell <- vector(length = d) # Parametervektor für Prädiktor
  eta.zufall.sample <- rep(x=0, time=d)
  alpha.eta <- 0 # Akzeptanzwahrscheinlichkeit für eta
  alpha.eta.proz.sample <- 0
  alpha.sample.aktuell <- vector(length = d) # Parametervektor für unstrukturierte Effekt
  if(AdaptiveGewichte){
    Q_weight.sample.aktuell <- matrix(data=0, nrow=dim(Q)[1], ncol=dim(Q)[1]) # gewichteter Nachbarschaftsmatrix
    weight_sr.sample.aktuell <- vector(length = length(zeilenindizes.oD_Q)) # Parametervektor der Gewichte
    alpha.weight <- 0  # Akzeptanzwahrscheinlichkeit für weight
    alpha.weight.proz.sample <- 0
  }
  
  
  ###############################
  ## Festlegung der Startwerte ##
  ###############################

  # Priori von mu: mu|. ~ NV(upsilon, sigma2)
  sigma_mu <- 100
  upsilon_mu <- 0
  mu.sample.aktuell <- 1
  if(!Intercept){mu.sample.aktuell <- 0}
  
  # Priori von beta: beta|. ~ NV(mu, sigma2)
  sigma_beta <- 100
  mu_beta <- 0
  beta.sample.aktuell <- 1
  if(!Koeff){beta.sample.aktuell <- 0}
  
  # Priori von kappa: kappa|. ~ G(a, b)
  a_kappa <- 1
  b_kappa <- 1
  kappa.sample.aktuell <- kappa.fest.wert
  
  # Priori von tau: tau|. ~ G(g, h)
  g_tau <- 1
  h_tau <- 0.00001
  tau.sample.aktuell <- tau.fest.wert
  
  # Priori von Gamma: Gamma|. ~ MNV(0, P^(-1))
  Ga.sample.aktuell <- rep(x=1, time=d)  
  
  if(AdaptiveGewichte){
    
    # Priori von weight_sr: weight_sr|. ~ G(e, f)
    e_weight_sr <- nu/2
    f_weight_sr <- nu/2
    # Ziehe negative symmetrische Gewichte für alle Nachbarschaften, also für alle Einträge von Q oberhalb der Diagonalmatrix (verwende oD_Q)
    for(p in 1:length(zeilenindizes.oD_Q)){ 
      Q_weight.sample.aktuell[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]] <- (rgamma(n = 1,  shape = e_weight_sr, rate = f_weight_sr))*(-1) # für Symmetrie
      Q_weight.sample.aktuell[spaltenindizes.oD_Q[p], zeilenindizes.oD_Q[p]] <- Q_weight.sample.aktuell[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]] # für Symmetrie
      weight_sr.sample.aktuell[p] <- (Q_weight.sample.aktuell[zeilenindizes.oD_Q[p], spaltenindizes.oD_Q[p]])*(-1)
    }    
    
    # Berechne Diagonale als negative Summe der jeweiligen Zeile
    diag(Q_weight.sample.aktuell) <- 0 # Setze Diagonale auf 0
    diag(Q_weight.sample.aktuell) <- (rowSums(Q_weight.sample.aktuell))*(-1) # Diagonaleinträge haben Gewicht gleich der negativen Zeilensumme
        
    if(GewDiagoVar){
      # Erstelle sortierte Vektoren einer "sparseMatrix" für spätere Neubestimmung der Einträge. 
      # Sortierung: 1) Indizes der oberen Dreiecksmatrix 2) Indizes der unteren Dreiecksmatrix in gleicher Reihenfolge 3) Indizes der Diagonaleinträge
      Sortierte.Indizes <- sortierte.indizes(Q_weight.sample.aktuell)
      Q_weight.sample.aktuell.sort.ZeilenInd <- Sortierte.Indizes[[1]]
      Q_weight.sample.aktuell.sort.SpaltenInd <- Sortierte.Indizes[[2]]
      Q_weight.sample.aktuell.sort.Eintraege <- Sortierte.Indizes[[3]]
    }    
    
    if(knappeMatrizen){
      Q_weight.sample.aktuell <- knappeMatrix(Q_weight.sample.aktuell) # Erstelle eine sparse Matrix der gewichteten Nachbarschaftsmatrix
    }
    
    # Berechne Startwert der logarithmierten Determinante von Q als ersten aktuellen Zustand
    # Verwende Cholesky-Zerlegung für schnellere Berechnung
    L_Q111 <- chol(Q_weight.sample.aktuell[-dim(Q)[1], -dim(Q)[1]])
    Q_weight.log.determinante.aktuell <- sum(2*log(diag(L_Q111)))  
  }
  
  # Priori von eta_s: eta_s|. ~ NV(mu + x_s*beta + Gamma_s, tau^(-1))
  # Ziehe Startvektor für eta.sample[[1]] aber nicht daraus sondern verwende log(y_s)
  # log(lambda_s) = eta_s -> Verwende Beobachtungen als Erwartungswert, also y_s gleichzusetzen mit lambda_s!
  # Zu beachten ist, dass y auch 0 sein kann, setze dann log(y_s) auf Null!
  if(eta.start.logDat){
    for(u in 1:length(y)){
      if(y[u] == 0){
        eta.sample.aktuell[u] <- 0
      }else{
        eta.sample.aktuell[u] <- log(y[u])
      }
    }
  }else{
    eta.sample.aktuell <- rep(x=1, time=d)
  }

  

  # Berechne den Startwert alpha.sample[1] aus der Reparametrisierung alpha_s = eta_s - mu - x_s*beta - Ga_s
  alpha.sample.aktuell <- eta.sample.aktuell - mu.sample.aktuell - x*beta.sample.aktuell - Ga.sample.aktuell

  
  # Übergebe den aktuellen Zustand an die 1. Position der Samples
  mu.sample <- mu.sample.aktuell
  beta.sample <- beta.sample.aktuell
  kappa.sample <- kappa.sample.aktuell
  tau.sample <- tau.sample.aktuell
  Ga.sample <- Ga.sample.aktuell
  eta.sample <- eta.sample.aktuell
  alpha.sample <- alpha.sample.aktuell
  if(AdaptiveGewichte){
    weight_sr.sample <- weight_sr.sample.aktuell
  }
  

  #####################################################
  ## Umsetzung der Iterations-Schritte ("size"-mal): ##
  #####################################################

    
  for (i in 2:(size+1)){ # gilt sowohl für Gibbs-Sampler als auch für Metropolis-Hastings

    
    ##############################
    ## Gibbs-Sampling-Schritte: ##
    ##############################
    
    if(Intercept){
      # Ziehe mu.sample[i] gegeben tau.sample[i-1], beta.sample[i-1], Ga.sample[,i-1], eta.sample[,i-1], Q_weight.sample[[i-1]] und x:
      sigma2_mu_schlange <- (sigma_mu^(-2) + tau.sample.aktuell*d)^(-1) 
      upsilon_mu_schlange <- (upsilon_mu*(sigma_mu^(-2)) - tau.sample.aktuell*sum(x*beta.sample.aktuell + Ga.sample.aktuell - eta.sample.aktuell))*sigma2_mu_schlange
      mu.sample.aktuell <- rnorm(n = 1, mean = upsilon_mu_schlange, sd = (sigma2_mu_schlange)^(0.5))
    }else{
      mu.sample.aktuell <- 0
    }
    
    if(Koeff){
      # Ziehe beta.sample[i] gegeben mu.sample[i], tau.sample[i-1], Ga.sample[,i-1], eta.sample[,i-1], Q_weight.sample[[i-1]] und x:
      sigma2_beta_schlange <- (sigma_beta^(-2) + tau.sample.aktuell*sum(x^2))^(-1) 
      mu_beta_schlange <- (mu_beta*(sigma_beta^(-2)) - tau.sample.aktuell*sum(x*(Ga.sample.aktuell + mu.sample.aktuell - eta.sample.aktuell)))*sigma2_beta_schlange
      beta.sample.aktuell <- rnorm(n = 1, mean = mu_beta_schlange, sd = (sigma2_beta_schlange)^(0.5))
    }else{
      beta.sample.aktuell <- 0
    }
      
    # Ziehe Ga.sample[,i] gegeben mu.sample[i], beta.sample[i], tau.sample[i-1], kappa.sample[i-1], eta.sample[,i-1], Q_weight.sample[[i-1]] und x:
    # P_weight_schlange1 <- diag(d)*tau.sample[i-1] + kappa.sample[i-1]*Q_weight.sample[[i-1]]
    P_weight_schlange <- kappa.sample.aktuell*Q_weight.sample.aktuell # Zeitersparnis
    diag(P_weight_schlange) <- diag(P_weight_schlange) + tau.sample.aktuell # Zeitersparnis
    b <- tau.sample.aktuell*(eta.sample.aktuell - mu.sample.aktuell - x*beta.sample.aktuell)
    Ga.sample.aktuell <- as.vector(rmvnorm_RBA(n = 1, b = b, P = P_weight_schlange))
  
    
    if(Hyp.kappa){
    # Ziehe kappa.sample[i] gegeben Ga.sample[,i], mu.sample[i], beta.sample[i], tau.sample[i-1], eta.sample[,i-1], Q_weight.sample[[i-1]] und x:
    a_kappa_schlange <- a_kappa + rgQ/2
    b_kappa_schlange <- b_kappa + as.vector(0.5*(t(Ga.sample.aktuell)%*%Q_weight.sample.aktuell%*%Ga.sample.aktuell))
    kappa.sample.aktuell <- rgamma(n = 1,  shape = a_kappa_schlange, rate = b_kappa_schlange)
    }else{
      kappa.sample.aktuell <- kappa.fest.wert
    } 
    
    if(Hyp.tau){
    # Ziehe tau.sample[i] gegeben Ga.sample[,i], kappa.sample[i],  mu.sample[i], beta.sample[i], eta.sample[,i-1], Q_weight.sample[[i-1]] und x:
    g_tau_schlange <- g_tau + d/2
    h_tau_schlange <- h_tau + 0.5*(sum((eta.sample.aktuell - x*beta.sample.aktuell - mu.sample.aktuell - Ga.sample.aktuell)^2))
    tau.sample.aktuell <- rgamma(n = 1,  shape = g_tau_schlange, rate = h_tau_schlange)
    }else{
      tau.sample.aktuell <- tau.fest.wert
    }

  
    
    ###################################
    ## Metropolis-Hastings-Schritte: ##
    ###################################
    
  
    ###################################################################################################################
    #################################################### eta ##########################################################
    ###################################################################################################################
    
   
      #####################################################
      ## 1. Schritt: Wähle aktuellen Zustand eta_s^(t-1) ##
      #####################################################
      
      # Startvektor bereits zu Beginn bestimmt mit eta.sample[,1] = eta_s^(0)
      eta.aktuell <- eta.sample.aktuell 
      

      ######################################################
      ## 2. Schritt: Ziehe eta_s^(*) aus Vorschlagsdichte ##
      ######################################################
      
      # Bestimme den zufälligen neuen Zustand eta_s^(*) gemäß der Vorschlagsdichte 
      # Verwende den aktuellen Zustand eta_s.aktuell als geigneten Wert für eta_s0
      b <- y + (e*exp(eta.aktuell))*(eta.aktuell - 1) + tau.sample.aktuell*(mu.sample.aktuell + x*beta.sample.aktuell + Ga.sample.aktuell)
      c <- e*exp(eta.aktuell) + tau.sample.aktuell # Präzision
      eta.zufall <- rnorm(n = d, mean = b/c, sd = sqrt(c^(-1))) # Normalverteilung über Taylorentwicklung

      
      ######################################################
      ## 3. Schritt: Berechne Akzeptanzwahrscheinlichkeit ##
      ######################################################
      
      
      # Zieldichte p(): vollständig bedingte Dichte p(eta_s | . ) gegeben 
      # tau.sample[i], Ga.sample[[i]], kappa.sample[i], mu.sample[i] und beta.sample[i].
      # Zustände dieser restlichen Parameter sind fest!
      
      # p(eta_s^(t-1) | . ): Zieldichte ausgewertet an der Stelle eta_s^(t-1), also zum aktuellen Zustand.
      log.zieldichte.eta.aktuell <- y*eta.aktuell - e*exp(eta.aktuell) - (tau.sample.aktuell/2)*(eta.aktuell - mu.sample.aktuell - x*beta.sample.aktuell - Ga.sample.aktuell)^2
      
      # p(eta_s^(*) | . ): Zieldichte ausgewertet an der Stelle eta_s^(*), 
      # also zu einem neuen zufällige Zustand gemäß der Vorschlagsdichte.
      log.zieldichte.eta.zufall <- y*eta.zufall - e*exp(eta.zufall) - (tau.sample.aktuell/2)*(eta.zufall - mu.sample.aktuell - x*beta.sample.aktuell - Ga.sample.aktuell)^2
      
      
      # Vorschlagsdichte q(): multivariate Normalverteilung über Taylorentwicklung gegeben 
      # tau.sample[i], Ga.sample[[i]], kappa.sample[i], mu.sample[i] und beta.sample[i].
      # Zustände dieser restlichen Parameter sind fest!
      
      # q(eta_s^(t-1) | eta_s^(*), . ): Vorschlagsdichte ausgewertet an der Stelle eta_s^(t-1), also zum aktuellen Zustand.
      # Wichtig: Taylorentwicklung um eta_s0 = eta_s^(*)
      b.aktuell <- y + (e*exp(eta.zufall))*(eta.zufall - 1) + tau.sample.aktuell*(mu.sample.aktuell + x*beta.sample.aktuell + Ga.sample.aktuell)
      c.aktuell <- e*exp(eta.zufall) + tau.sample.aktuell  # Präzision
      log.vorschlDichte.eta.aktuell <- 0.5*c.aktuell*(eta.aktuell)^2 + b.aktuell*eta.aktuell
      
      
      # q(eta_s^(*) | eta_s^(t-1), . ): Vorschlagsdichte ausgewertet an der Stelle eta_s^(*), 
      # also zu einem neuen zufällige Zustand gemäß der Vorschlagsdichte.
      # Wichtig: Taylorentwicklung um eta_s0 = eta_s^(t-1)
      b.zufall <- y + (e*exp(eta.aktuell))*(eta.aktuell - 1) + tau.sample.aktuell*(mu.sample.aktuell + x*beta.sample.aktuell + Ga.sample.aktuell)
      c.zufall <- e*exp(eta.aktuell) + tau.sample.aktuell  # Präzision
      log.vorschlDichte.eta.zufall <- 0.5*c.zufall*(eta.zufall)^2 + b.zufall*eta.zufall
      
      
      # Explizite Berechnung der Akzeptanzwahrscheinlichkeit
      
      # Berechne Zähler mit p(eta_s.zufall)*q(eta_s.aktuell)
      log.alpha.zaehler.eta <- log.zieldichte.eta.zufall + log.vorschlDichte.eta.aktuell
      
      # Berechne Nenner mit p(eta_s.aktuell)*q(eta_s.zufall)
      log.alpha.nenner.eta <- log.zieldichte.eta.aktuell + log.vorschlDichte.eta.zufall
      
      if(i > eta.start.ab){ 
        for(s in 1:d){
          # Vermeide Subtraktion von -Inf
          if(log.alpha.nenner.eta[s] == -Inf){
            log.alpha.eta <- Inf # d.h. akzeptieren
          }else{
            log.alpha.eta <- log.alpha.zaehler.eta[s] - log.alpha.nenner.eta[s]
          }
          
          # Akzeptiere...
          if(log(runif(1)) < log.alpha.eta){
            eta.sample.aktuell[s] <- eta.zufall[s]
            alpha.eta <- alpha.eta + 1
          }else{ # ...oder verwerfe.
            eta.sample.aktuell[s] <- eta.aktuell[s]
          }
        }
      }else{ # akzeptiere die 1. 10 Vorschläge
        eta.sample.aktuell <- eta.zufall
        alpha.eta <- alpha.eta + d
      }
      
      # Berechne den alpha.sample (unstrukturierter Effekt) aus der Reparametrisierung alpha_s = eta_s - mu - x_s*beta - Ga_s
      alpha.sample.aktuell <- eta.sample.aktuell - mu.sample.aktuell - x*beta.sample.aktuell - Ga.sample.aktuell
      
    
    
    if(AdaptiveGewichte){
    ###################################################################################################################
    ################################################### weight ########################################################
    ###################################################################################################################
    
  
    # Ist für alle Gewichte gleich, deswegen außerhalb der for-Schleife
    e_weight_sr_schlange <- nu/2
    f_weight_sr_schlange <- nu/2 + as.numeric((kappa.sample.aktuell*(t(Ga.sample.aktuell)%*%Ga.sample.aktuell))/2)
    
    ####################################################
    ## 1. Schritt: Wähle aktuellen Zustand w_sr^(t-1) ##
    ####################################################
    
    # Bestimme den aktuellen Zustand von Q_weight.aktuell aus der Matrix Q_weight.aktuell^(t-1) zum aktuellen Zustand
    # weight_sr^(0) wurde bereits zu Beginn bestimmt entspricht Q_weight.sample[[1]]
    Q_weight.aktuell <- Q_weight.sample.aktuell
    
    for(g in 1:length(zeilenindizes.oD_Q)){ # for Schleife um alle Gewichtseinträge oberhalb der Diagonalmatrix zu bestimmen
      
      ##########################################################
      ## 2. Schritt: Ziehe weight_sr^(*) aus Vorschlagsdichte ##
      ##########################################################
      
      # Bestimme den zufälligen neuen Zustand weight_sr^(*) gemäß der Vorschlagsdichte: weight_sr|. ~ G(e_weight_sr_schlange, f_weight_sr_schlange) 
      # Schließe Werte kleiner 10^(-6) aus für neue Gewichte. Setze Werte dann auf Untergrenze.
      weight_sr.zufall <- rgamma(n = 1,  shape = e_weight_sr_schlange, rate = f_weight_sr_schlange)
      if(weight_sr.zufall < 10^(-6)){
        weight_sr.zufall <- 10^(-6)
      }
      
      
      # Erstelle gewichtete Nachbarschaftsmatrizen mit weight_sr^(*)
      # Neue Zufallsmatrix entspricht zunächst dem aktuellen Zustand von Q_w
      Q_weight.zufall <- Q_weight.aktuell 
      
      
      if(GewDiagoVar){
        Q_weight.zufall.sort.Eintraege <- Q_weight.sample.aktuell.sort.Eintraege
        
        # Abspeichern des alten Gewichts an der g-ten Stelle 
        weight_sr.alt <- Q_weight.zufall.sort.Eintraege[g]
        # Änderung des zugehörigen Eintrags
        Q_weight.zufall.sort.Eintraege[g] <- weight_sr.zufall*(-1)
        # für Symmetrie
        Q_weight.zufall.sort.Eintraege[g + length(zeilenindizes.oD_Q)] <- weight_sr.zufall*(-1)
        # Berechnung des zugehörigen Diagonalentrags der Zeile
        Q_weight.zufall.sort.Eintraege[zeilenindizes.oD_Q[g] + 2*length(zeilenindizes.oD_Q)] <- 
          Q_weight.zufall.sort.Eintraege[zeilenindizes.oD_Q[g] + 2*length(zeilenindizes.oD_Q)] + weight_sr.alt + weight_sr.zufall
        # Berechnung des zugehörigen Diagonalentrags der Spalte
        Q_weight.zufall.sort.Eintraege[spaltenindizes.oD_Q[g] + 2*length(zeilenindizes.oD_Q)] <- 
          Q_weight.zufall.sort.Eintraege[spaltenindizes.oD_Q[g] + 2*length(zeilenindizes.oD_Q)] + weight_sr.alt + weight_sr.zufall
        
        Q_weight.zufall <- sparseMatrix(i = Q_weight.sample.aktuell.sort.ZeilenInd, j = Q_weight.sample.aktuell.sort.SpaltenInd, x = Q_weight.zufall.sort.Eintraege, dims = dim(Q_weight.zufall))
        
      }else{        
        # Ordne die beiden (wg. Symmetrie) zufälligen Werte von w_sr in diese Zufallsmatrix ein
        Q_weight.zufall[zeilenindizes.oD_Q[g], spaltenindizes.oD_Q[g]] <- (weight_sr.zufall)*(-1)
        Q_weight.zufall[spaltenindizes.oD_Q[g], zeilenindizes.oD_Q[g]] <- (weight_sr.zufall)*(-1) # Für Symmetrie 
  
        # Berechne Diagonale neu  
        Q_weight.zufall[zeilenindizes.oD_Q[g],zeilenindizes.oD_Q[g]] <- 0 # Setze entsprechenden Diagonaleintrag der Zeile zunächst auf 0
        Q_weight.zufall[spaltenindizes.oD_Q[g],spaltenindizes.oD_Q[g]] <- 0 # Setze entsprechenden Diagonaleintrag der Spalte zunächst auf 0
        Q_weight.zufall[zeilenindizes.oD_Q[g],zeilenindizes.oD_Q[g]] <- (sum(Q_weight.zufall[zeilenindizes.oD_Q[g],]))*(-1) # Ersetze Diagonaleintrag der Zeilen durch negative Zeilensumme
        Q_weight.zufall[spaltenindizes.oD_Q[g],spaltenindizes.oD_Q[g]] <- (sum(Q_weight.zufall[spaltenindizes.oD_Q[g],]))*(-1) # Ersetze Diagonaleintrag der Spalten durch negative Zeilensumme
      }
      
      
      ######################################################
      ## 3. Schritt: Berechne Akzeptanzwahrscheinlichkeit ##
      ######################################################
      
      # Zähler der Akzeptanzwahrscheinlichkeit:

      # Berechne logarithmierte Determinante über Cholesky-Zerlegung, 
      # da möglicherweise Pobleme bei extrem kleinen Diagonaleinträgen entstehen können
      L_Qwz11 <- chol(Q_weight.zufall[-d, -d]) # mit Löschen der letzten Zeile und der letzten Spalte 
      Q_weight.zufall11.logdet <- sum(2*log(diag(L_Qwz11)))
      
      # Logarithmus (für Berechnung von alpha nötig) der Wurzel der neuen exponentierten logarithmierten Determinante als Zähler
      log.alpha.zaehler.Q_weight.zufall11 <- log((exp(Q_weight.zufall11.logdet))^(0.5))
      
      
      # Nenner der Akzeptanzwahrscheinlichkeit:
      # Berechne logarithmierte Determinante nicht neu, sondern verwende die zuvor berechnete für den aktuellen Zustand 
      # Logarithmus (für Berechnung von alpha nötig) der Wurzel der aktuellen exponentierten logarithmierten Determinante als Nenner
      log.alpha.nenner.Q_weight.aktuell11 <-  log((exp(Q_weight.log.determinante.aktuell))^(0.5))
      
      
      # Explizite Berechnung der Akzeptanzwahrscheinlichkeit
      
      # Vermeide Subtraktion von -Inf
      if(log.alpha.nenner.Q_weight.aktuell11 == -Inf){
        log.alpha.Q_weight <- 1 # d.h. akzeptieren
      } else{
        log.alpha.Q_weight <- log.alpha.zaehler.Q_weight.zufall11 - log.alpha.nenner.Q_weight.aktuell11
      }
      
      # Akzeptiere...
      if(log(runif(1)) < log.alpha.Q_weight){
        Q_weight.aktuell <- Q_weight.zufall # dadurch wird immer die aktuell angepasste Gewichtsmatrix verwendet
        Q_weight.log.determinante.aktuell <- Q_weight.zufall11.logdet # Abspeichern der logarithmierten Determinante
        alpha.weight <- alpha.weight + 1
        if(GewDiagoVar){Q_weight.sample.aktuell.sort.Eintraege <- Q_weight.zufall.sort.Eintraege}
      }

      
      # Zusätzliches Abspeichern des Sampling-Pfads der einzelnen Gewichte
      weight_sr.sample.aktuell[g] <- Q_weight.aktuell[zeilenindizes.oD_Q[g], spaltenindizes.oD_Q[g]]*(-1)
      
      
    } # Ende der for-Schleife mit Laufindex g
    
    Q_weight.sample.aktuell <- Q_weight.aktuell
    
    
    } # Ende der Bedingung, ob Gewichte angepasst werden sollen.
    
    #############################################################################
    ## Übergeben bei jeder thin-ten Iteration an den Sample inklusive Burn-In: ##
    #############################################################################
    
    # Berechne Akzeptanzwahrscheinlichkeiten
    alpha.eta.proz <- alpha.eta / (d * (i-1))
    
    if(AdaptiveGewichte){
      alpha.weight.proz <- alpha.weight / (length(zeilenindizes.oD_Q) * (i-1))
    }
    
    
    if ((i%%thin) == 0 && i > burn){ # Thinning=10 (z.B. jede 10. Iteration wird abgespeichert) und Burn-In (z.B. nach der 100. Iteration)
      mu.sample <- c(mu.sample, mu.sample.aktuell)
      beta.sample <- c(beta.sample, beta.sample.aktuell)
      Ga.sample <- cbind(Ga.sample, Ga.sample.aktuell)
      kappa.sample <- c(kappa.sample, kappa.sample.aktuell)
      eta.sample <- cbind(eta.sample, eta.sample.aktuell)
      eta.zufall.sample <- cbind(eta.zufall.sample, eta.zufall)
      alpha.eta.proz.sample <- c(alpha.eta.proz.sample, alpha.eta.proz)
      tau.sample <- c(tau.sample, tau.sample.aktuell)
      alpha.sample <- cbind(alpha.sample, alpha.sample.aktuell)
      if(AdaptiveGewichte){
        weight_sr.sample <- cbind(weight_sr.sample, weight_sr.sample.aktuell)
       alpha.weight.proz.sample <- c(alpha.weight.proz.sample, alpha.weight.proz)
      }
    }

  if(IterationPrint){print(i)}
  
  if(ZS){  
    # Zwischenspeichern
    if((i%%5000) == 0){
      
      zs <- i/5000
      
      save(mu.sample, file=paste("Samples/Zwischenspeicher", zs, "/mu.sample.RData", sep=""))
      save(beta.sample, file=paste("Samples/Zwischenspeicher", zs, "/beta.sample.RData", sep=""))
      save(Ga.sample, file=paste("Samples/Zwischenspeicher", zs, "/Ga.sample.RData", sep=""))
      save(kappa.sample, file=paste("Samples/Zwischenspeicher", zs, "/kappa.sample.RData", sep=""))
      save(eta.sample, file=paste("Samples/Zwischenspeicher", zs, "/eta.sample.RData", sep=""))
      save(eta.zufall.sample, file=paste("Samples/Zwischenspeicher", zs, "/eta.zufall.sample.RData", sep=""))
      save(alpha.eta.proz.sample, file=paste("Samples/Zwischenspeicher", zs, "/alpha.eta.proz.sample.RData", sep=""))
      save(tau.sample, file=paste("Samples/Zwischenspeicher", zs, "/tau.sample.RData", sep=""))
      save(alpha.sample, file=paste("Samples/Zwischenspeicher", zs, "/alpha.sample.RData", sep=""))
      if(AdaptiveGewichte){
        save(weight_sr.sample, file=paste("Samples/Zwischenspeicher", zs, "/weight_sr.sample.RData", sep=""))
        save(Q_weight.sample.aktuell, file=paste("Samples/Zwischenspeicher", zs, "/Q_weight.sample.RData", sep=""))
        save(alpha.weight.proz.sample, file=paste("Samples/Zwischenspeicher", zs, "/alpha.weight.proz.sample.RData", sep=""))
      }
          
    }
  }
 

  } # Ende der for-Schleife mit Laufindex i


  
  ##############################
  ## Abspeichern der Samples: ##
  ##############################
  
  # Anpassen der Form
  colnames(Ga.sample) <- c()
  colnames(eta.sample) <- c()
  colnames(alpha.sample) <- c()
  if(AdaptiveGewichte){colnames(weight_sr.sample) <- c()}
  
  # Funktionsparameter übergeben und abspeichern
  funktionsparameter <- data.frame(data=c(size, nu, burn, thin, kappa.fest.wert, tau.fest.wert, eta.start.ab))
  colnames(funktionsparameter) <- "Funktionsparameter:"
  rownames(funktionsparameter) <- c("Iterationsumfang: size =", "Gewichtungsparameter: nu =", "Burn-In Phase: burn =", "Ausdünnungsfaktor: thin =", "Festgesetzter Wert von kappa =", "Festgesetzter Wert von tau =", "Akzeptanzwahrscheinlichkeit für eta 100% bis")
  write.table(funktionsparameter, "Samples/Funktionsparameter.txt", row.names = TRUE, col.names = TRUE, quote=FALSE)
  
  
  # Vollständige Samples
  save(mu.sample, file="Samples/mu.sample.RData")
  save(beta.sample, file="Samples/beta.sample.RData")
  save(Ga.sample, file="Samples/Ga.sample.RData")
  save(kappa.sample, file="Samples/kappa.sample.RData")
  save(eta.sample, file="Samples/eta.sample.RData")
  save(eta.zufall.sample, file="Samples/eta.zufall.sample.RData")
  save(alpha.eta.proz.sample, file="Samples/alpha.eta.proz.sample.RData")
  save(tau.sample, file="Samples/tau.sample.RData")
  save(alpha.sample, file="Samples/alpha.sample.RData")
  if(AdaptiveGewichte){
    save(weight_sr.sample, file="Samples/weight_sr.sample.RData")
    save(Q_weight.sample.aktuell, file="Samples/Q_weight.sample.aktuell.RData")
    save(alpha.weight.proz.sample, file="Samples/alpha.weight.proz.sample.RData")
  }
  

} # Ende der Hybrid-Sampler Funktion
  



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


###################################
## Anwendung des Hybrid-Samplers ##
###################################

size <- 25000 # Samplingumfang
thin <- 5 # jeder thin-te Wert  
burn <- 0 # erst ab burn-ten Wert



##### Schottland #####

# ohne Gewichtsadaption

# Schottland ohne "knappeMatrix" ohne Gewichtsanpassung ohne mu mit tau=1 fest und kappa=TRUE.
system.time({  # Laufzeit messen
  hybrid.sampler(size=size, x=x.scot, y=y.scot, e=e.scot, Q=Q.scot, nu=1, burn=burn, thin=thin, Qvoll=Q.scot, knappeMatrizen=FALSE, Intercept=FALSE, AdaptiveGewichte=FALSE, IterationPrint=TRUE, Hyp.kappa=TRUE, Hyp.tau=FALSE, eta.start.ab=0, tau.fest.wert=1)
})


# mit Gewichtsadaption

# Schottland ohne "knappeMatrix" mit Gewichtsanpassung ohne mu, nu=1 und tau=1 fest und kappa=0.001 fest ohne Koeff.
system.time({  # Laufzeit messen
  hybrid.sampler(size=size, x=x.scot, y=y.scot, e=e.scot, Q=Q.scot, nu=1, burn=burn, thin=thin, Qvoll=Q.scot, knappeMatrizen=FALSE, Intercept=FALSE, AdaptiveGewichte=TRUE, IterationPrint=TRUE, Hyp.kappa=FALSE, Hyp.tau=FALSE, Koeff=FALSE, eta.start.ab=0, eta.start.logDat=TRUE, kappa.fest.wert=0.001, tau.fest.wert=1)
})



##### Deutschland #####

# ohne Gewichtsadaption

# Deutschland mit "knappeMatrix" ohne Gewichtsanpassung ohne mu mit tau=1 fest und kappa=TRUE.
system.time({  # Laufzeit messen
  hybrid.sampler(size=size, x=x.germany, y=y.germany, e=e.germany, Q=Q.germanyknapp, nu=1, burn=burn, thin=thin, Qvoll=Q.germany, knappeMatrizen=TRUE, Intercept=FALSE, AdaptiveGewichte=FALSE, IterationPrint=TRUE, Hyp.kappa=TRUE, Hyp.tau=FALSE, eta.start.ab=5, eta.start.logDat=TRUE, tau.fest.wert=1)
})


# mit Gewichtsadaption

# Deutschland mit "knappeMatrix" mit Gewichtsanpassung ohne mu mit tau=1 fest und kappa=0.001 fest ohne Koeff, GewDiagoVar=TRUE.
system.time({  # Laufzeit messen
  hybrid.sampler(size=size, x=x.germany, y=y.germany, e=e.germany, Q=Q.germanyknapp, nu=1, burn=burn, thin=thin, Qvoll=Q.germany, knappeMatrizen=TRUE, Intercept=FALSE, AdaptiveGewichte=TRUE, IterationPrint=TRUE, Hyp.kappa=FALSE, Hyp.tau=FALSE, Koeff=FALSE, eta.start.ab=5, eta.start.logDat=TRUE, kappa.fest.wert=0.001, tau.fest.wert=1, GewDiagoVar=TRUE, ZS=TRUE)
})



##### Anwendung auf simulierte Daten #####

load("y.sim10.setting3.RData")
y.sim

size <- 200000 # Samplingumfang
thin <- 200 # jeder thin-te Wert  
burn <- 0 # erst ab burn-ten Wert

system.time({  # Laufzeit messen
  hybrid.sampler(size=size, x=1, y=y.sim, e=1, Q=Q.scot, nu=1, burn=burn, thin=thin, Qvoll=Q.scot, knappeMatrizen=FALSE, Intercept=FALSE, AdaptiveGewichte=TRUE, IterationPrint=TRUE, Hyp.kappa=FALSE, Hyp.tau=FALSE, Koeff=FALSE, eta.start.ab=0, eta.start.logDat=TRUE, kappa.fest.wert=0.001, tau.fest.wert=1)
})

