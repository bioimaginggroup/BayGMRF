# Berechne gleitenden Mittelwert
gleit.Mittelwert <- function(parameter.sample){
  rm <- 0
  for(i in 1:length(parameter.sample[(burn):size])){
    rm[i] <- (sum(parameter.sample[(burn):(burn+i-1)]))/i
  }
  return(rm)
}


for(k in 1:4){
  for(j in 1:10){

    # Laden der Samples
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/mu.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/beta.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/kappa.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/Ga.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/eta.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/tau.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/alpha.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/eta.zufall.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/weight_sr.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/alpha.eta.proz.sample.RData", sep=""))
    load(paste("Simulierte Samples/Setting", k,"/y.sim", j, "/alpha.weight.proz.sample.RData", sep=""))
    
    
    size <- length(mu.sample)
    burn <- 500
    
    
    
    
    #################################
    ## Grafische Veranschaulichung ##
    #################################
    
    
    # Strukturierter Effekt
    
    # apply(X=Ga.sample[,burn:size], MARGIN=1, FUN=mean)
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Samplingpfad - 53 strukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(Ga.sample[,size])){
        plot(Ga.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
         cex.lab=1, cex.axis=1, main=paste("Samplingpfad strukturierter Effekt (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Dichte - 53 strukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(Ga.sample[,size])){
    plot(density(Ga.sample[i,burn:size]), xlab=expression(paste(gamma, " der Region s")), 
         lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte strukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Autokorrelationsfunktion - 53 strukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,3))
    for(i in 1:length(Ga.sample[,size])){
    acf(Ga.sample[i,burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
        cex.lab=1, cex.axis=1 , main=paste("ACF strukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Mittelwertsentwicklung - 53 strukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(Ga.sample[,size])){
    plot(gleit.Mittelwert(Ga.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
         ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. strukturierter Effekt (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    
    
    
    # Unstrukturierter Effekt
    
    #apply(X=alpha.sample[,burn:size], MARGIN=1, FUN=mean)
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Samplingpfad - 53 unstrukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(alpha.sample[,size])){
      plot(alpha.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
           cex.lab=1, cex.axis=1, main=paste("Samplingpfad unstrukturierter Effekt (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Dichte - 53 unstrukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(alpha.sample[,size])){
      plot(density(alpha.sample[i,burn:size]), xlab=expression(paste(alpha, " der Region s")), 
           lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte unstrukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Autokorrelationsfunktion - 53 unstrukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,3))
    for(i in 1:length(alpha.sample[,size])){
      acf(alpha.sample[i,burn:size], lag.max=size, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
          cex.lab=1, cex.axis=1 , main=paste("ACF unstrukturierter Effekt (Region s=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Mittelwertsentwicklung - 53 unstrukturierte Effekte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(alpha.sample[,size])){
      plot(gleit.Mittelwert(alpha.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
           ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. unstrukturierter Effekt (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    
    
    
    
    # Prädiktor eta
    
    #apply(X=eta.sample[,burn:size], MARGIN=1, FUN=mean)
    
    
    #alpha.eta.proz.sample[1] <- 1
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Akzeptanzwahrscheinlichkeit - Prädiktor (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=10, height=8) 
    par(mfrow=c(1,1))
    plot(alpha.eta.proz.sample[-1], type="l", ylim=c(0,1), cex.lab=1.8, cex.axis=1.8, 
         ylab="", xlab="Ausgedünnte Iterationen", lwd=4)
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Samplingpfad - 53 Prädiktorenwerte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(eta.sample[,size])){
      plot(eta.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
           cex.lab=1, cex.axis=1, main=paste("Samplingpfad Prädiktor (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Dichte - 53 Prädiktorenwerte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(eta.sample[,size])){
      plot(density(eta.sample[i,burn:size]), xlab=expression(paste(eta, " der Region s")), 
           lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte Prädiktor (Region s=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Autokorrelationsfunktion - 53 Prädiktorenwerte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,3))
    for(i in 1:length(eta.sample[,size])){
      if(eta.sample[i,burn] != eta.sample[i,size]){
        acf(eta.sample[i,burn:size], lag.max=100, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
            cex.lab=1, cex.axis=1 , main=paste("ACF Prädiktor (Region s=", i, ")", sep=""), ylab="")
      }
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Mittelwertsentwicklung - 53 Prädiktorenwerte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(eta.sample[,size])){
      plot(gleit.Mittelwert(eta.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
           ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. Prädiktor (Region s=", i, ")", sep=""))
    }
    dev.off()
    
    
    
    
    
    # Gewichte
    
    #apply(X=weight_sr.sample[,burn:size], MARGIN=1, FUN=mean)
    
    
    #alpha.weight.proz.sample[1] <- 1
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Akzeptanzwahrscheinlichkeit - Gewichte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=10, height=8) 
    par(mfrow=c(1,1))
    plot(alpha.weight.proz.sample[-1], type="l", ylim=c(0,1), cex.lab=1.8, cex.axis=1.8, 
         ylab="", xlab="Ausgedünnte Iterationen", lwd=4)
    dev.off()
    
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Samplingpfad - 117 Nachbarschaftsgewichte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(weight_sr.sample[,size])){
      plot(weight_sr.sample[i,burn:size], type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", ylab="", 
           cex.lab=1, cex.axis=1, main=paste("Samplingpfad Gewicht (Nachbarschaft g=", i, ")", sep=""))
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Dichte - 117 Nachbarschaftsgewichte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(weight_sr.sample[,size])){
      plot(density(weight_sr.sample[i,burn:size]), xlab="Gewicht der Nachbarschaft g", 
           lwd=4, cex.lab=1, cex.axis=1 , main=paste("Dichte Gewicht (Nachbarschaft g=", i, ")", sep=""), ylab="")
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Autokorrelationsfunktion - 117 Nachbarschaftsgewichte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,3))
    for(i in 1:length(weight_sr.sample[,size])){
      if(weight_sr.sample[i,burn] != weight_sr.sample[i,size]){
        acf(weight_sr.sample[i,burn:size], lag.max=100, xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)", 
            cex.lab=1, cex.axis=1 , main=paste("ACF Gewicht (Nachbarschaft g=", i, ")", sep=""), ylab="")
      }
    }
    dev.off()
    
    
    pdf(paste("Plots/Setting ", k,"/Simulation ", j, "/Mittelwertsentwicklung - 117 Nachbarschaftsgewichte (Setting ", k, " Simulation ", j, ").pdf", sep=""), width=8.27, height=11.69)
    par(mfrow=c(5,2))
    for(i in 1:length(weight_sr.sample[,size])){
      plot(gleit.Mittelwert(weight_sr.sample[i,]), type="l", xlab="Ausgedünnte Iterationen (nach der Burn-In Phase)",
           ylab="", cex.lab=1, cex.axis=1, lwd=2, main=paste("Mittelwertsentw. Gewicht (Nachbarschaft g=", i, ")", sep=""))
    }
    dev.off()

  }
}

