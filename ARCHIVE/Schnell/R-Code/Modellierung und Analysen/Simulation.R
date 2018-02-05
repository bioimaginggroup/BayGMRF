##################################
## Lade Package und Funktionen ##
##################################

library(Matrix)

source("extrafunktionen.R") # Laden der benötigten Funktionen
# beinhaltet die Funktionen:
# 1) "readandchangegraph" zum Einlesen einer .graph Datei (Nachbarschaftsmatrix)
# 2) "knappeMatrix" zur Umformung einer Matrix zu einer sparse Matrix
# 3) "MatrizenIndizes" zum Ausgeben der besetzten Zeilen- und Spaltenindizes einer Matrix
# 4) "rmvnorm_RBA" zum Ziehen aus einer multivariaten Normalverteilung


# Laden der reduzierten schottischen Nachbarschaftsmatrix der 4 Settings
Q.setting1 <- as.matrix(read.table("redscotQ.stark.setting1.txt"))
colnames(Q.setting1) <- c()
Q.setting1 <- Q.setting1[-(54:56),-(54:56)] # ohne Inseln

Q.setting2 <- as.matrix(read.table("redscotQ.stark.setting2.txt"))
colnames(Q.setting2) <- c()
Q.setting2 <- Q.setting2[-(54:56),-(54:56)] # ohne Inseln

Q.setting3 <- as.matrix(read.table("redscotQ.stark.setting3.txt"))
colnames(Q.setting3) <- c()
Q.setting3 <- Q.setting3[-(54:56),-(54:56)] # ohne Inseln

Q.setting4 <- as.matrix(read.table("redscotQ.stark.setting4.txt"))
colnames(Q.setting4) <- c()
Q.setting4 <- Q.setting4[-(54:56),-(54:56)] # ohne Inseln

Q <- Q.setting4




# Werte der Hyperparameter von kappa und tau sind fest
kappa <- 1
tau <- 100


# Ziehen aus einer multivariaten Normalverteilung
# Priori von Gamma: Gamma|. ~ MNV(0, P^(-1))
P.Gamma <- kappa*Q # Verwende hier ungewichtete Nachbarschaftsmatrix und addiere kleinstmöglichen Wert um Invertierbarkeit zu gewährleisten!!!
diag(P.Gamma) <- diag(P.Gamma) + 10^(-13) 
b_null <- rep(0, nrow(P.Gamma)) # b ist komplett 0, da Erwartungswert 0 ist
Gamma.sim <- as.vector(rmvnorm_RBA(n = 1, b = b_null, P = P.Gamma)) # hier wird Ga komplett gezogen d.h. Vektor der Länge d
#Gamma.sim

# Standardisierung der gezogenen Werte
Gamma.standardisiert <- (Gamma.sim - mean(Gamma.sim)) / sd(Gamma.sim) # Erwartungswert soll 0 sein
#Gamma.standardisiert
#mean(Gamma.standardisiert)
#sd(Gamma.standardisiert)


# Ziehen aus einer multivariaten Normalverteilung
# Priori von alpha: alpha|. ~ MNV(0, T^(-1))

# Diagonalmatrix Tau
Tau <- matrix(data=0, ncol=dim(Q)[1], nrow=dim(Q)[2])
diag(Tau) <- tau

P.alpha <- Tau  # Verwende hier ungewichtete Nachbarschaftsmatrix und addiere kleinstmöglichen Wert um Invertierbarkeit zu gewährleisten!!!
b_null <- rep(0, nrow(P.alpha)) # b ist komplett 0, da Erwartungswert 0 ist
alpha <- as.vector(rmvnorm_RBA(n = 1, b = b_null, P = P.alpha)) # hier wird Ga komplett gezogen d.h. Vektor der Länge d

alpha.standardisiert <- (alpha - mean(alpha)) / sd(alpha)
#mean(alpha.standardisiert)
#sd(alpha.standardisiert)


# Berechne lambda_s
lambda <- exp(Gamma.standardisiert + alpha.standardisiert)

# Ziehe Zufallszahlen aus der Poissonverteilung 

y.sim <- NA

for(i in 1:length(lambda)){
 y.sim[i] <- rpois(1, lambda[i])
}

# y.sim

save(y.sim, file="Simulierte Daten/Setting 4/10. Lauf/y.sim10.Setting4.RData")
save(Gamma.standardisiert, file="Simulierte Daten/Setting 4/10. Lauf/Gamma.sim10.Setting4.RData")
save(alpha.standardisiert, file="Simulierte Daten/Setting 4/10. Lauf/alpha.sim10.Setting4.RData")

y.sim