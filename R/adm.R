adm <- function(y, e, coords, nr.it, burnin, thin, nu=2, adaptive=TRUE, method="mh", do.kappa=TRUE,  kappa=1, ab.kappa=c(1,1), tau=1000, print.time=FALSE){
   
  
  d <- length(y)  
  nu.a<-nu.b<-nu/2
  if(length(nu)==2)
  {
  nu.a<-nu[1]
  nu.b<-nu[2]
  }
  p<-max(dim(coords))

  gamma <- vector(length = d) 
  eta <- vector(length = d) 
  if(adaptive){
    weights <- vector(length = p) 
  }
    
  gamma <- rep(x=0, time=d)  
  
 Q<-weightedQ(rep(1,p),coords)

 if (method=="mh2"|method=="gibbs2"|method=="gibbs2beta") Q1<-solve(Q[-d,-d])
 
 # Hier Überprüfung, welche Regionen nur einen Nachbarn haben
 # Es fehlen die Nachbarschaften von Regionen, die auch indirekt keine gemeinsamen Nachbarn
 #  haben, also Regionen trennen (Markoveigenschaft!) 
 if (adaptive){#} & method!="approx"){
  nr.nei<-apply(coords,1,test.for.eins,coords)
  nr.einer <- any(nr.nei==1)
  nr.nei.eins <- which(nr.nei==1)
}

 
 diagM <- array(0,c(p,d))
 for (i in 1:p)
 {
   diagM[i,as.vector(coords[i,])]<-1
 }
 
if(adaptive)if (method=="gibbs.approx.simple")
{
  if(require(parallel))next.n <- mclapply(1:p, find.nn, coords, nr.nei)
  if(!require(parallel))next.n <- lapply(1:p, find.nn, coords, nr.nei)
}
if(adaptive)if (method=="gibbs.approx"|method=="gibbs")
{
  if(require(parallel))Gg <- mclapply(1:p, estimate.g, coords)
  if(!require(parallel))Gg <- lapply(1:p, estimate.g, coords)  
}

  if(adaptive){
   weights <- rep(nu.a/nu.b,p)
   if (method=="gibbs2beta")weights=rep(1,p)
    Q<-weightedQ(weights,coords)
  }

if (method=="gibbs")
{
    #diags<-unlist(mclapply(1:d,function(i,weights,coords){sum(weights[coords[,1]==i])+sum(weights[coords[,2]==i])},weights,coords))

    diags<-(weights%*%diagM)[1,]
}

  for(u in 1:length(y)){
    if(y[u] == 0){
      eta[u] <- 0
      }
    else
      {
        eta[u] <- log(y[u]/e[u])
      }
    }  
  
if (method=="gibbs3"|method=="gibbs3beta")
{
  par.list=par.list2=list()
  counter=0
  for (i in sample(1:p))
  {
    if (!nr.nei[i])
    {
      if (counter==0)
      {
        par.list[[1]]=rep(FALSE,p)
        par.list[[1]][coords[i,]]=TRUE
        par.list2[[1]]=i
        counter=1
      }
      else
      {
        test=FALSE
        j=0
        while((!test)&(j<counter))
        {
          j=j+1
          test=!any(par.list[[j]][coords[i,]])
        }
        if(test)
        {
          par.list[[j]][coords[i,]]=TRUE
          par.list2[[j]]=c(par.list2[[j]],i)
        }
        else
        {
          par.list[[j+1]]=rep(FALSE,p)
          par.list[[j+1]][coords[i,]]=TRUE
          par.list2[[j+1]]=i  
          counter=j+1
        }
        par=par.list
        par2=par.list2
        new.sample<-sample(counter)
        for (k in 1:counter)
        {
        par.list[[new.sample[k]]]<-par[[k]]
        par.list2[[new.sample[k]]]<-par2[[k]]
        }
      }
    }
  }
}

  nr.samples<-floor((nr.it-burnin)/thin)
  gamma.sample <- array(NA,c(nr.samples,d))
  alpha.sample <- array(NA,c(nr.samples,d))
  eta.sample <- array(NA,c(nr.samples,d))
  
 kappa.sample <- rep(NA,nr.samples)
 if(adaptive){
   weights.sample <- array(NA,c(nr.samples,p))
 }
 

  for (iteration in 1:nr.it){
    if(print.time)start.time<-proc.time()
    if(iteration%%10==0){cat(".")}
    # gamma
    #print(kappa)
    #image(as.matrix(Q))
    P <- kappa*Q + diag(tau,d)
    #print(class(P))
    #diag(P) <- diag(P) + tau 
    b <- tau*eta
    gamma <- as.vector(rmvnorm(n = 1, b = b, P = P))
    
    # kappa
    if(do.kappa){
      aa <- ab.kappa[1] + (d-1)/2
      bb <- ab.kappa[2] + as.vector(0.5*(t(gamma)%*%Q%*%gamma))
      kappa <- rgamma(n = 1,  shape = aa, rate = bb)
    }
    
    # eta
    m.temp <- (gamma +log(e))*tau + y
    eta<-sample.eta(eta, m.temp, tau,e)
    
    
  
if(adaptive)
      {

      if (method=="gibbs")
      {
        if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
        {
          nee<-coords[nr.nei.eins,]
          if(is.null(dim(nee)))
          {
            b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
          }
          else
          {
            b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
          }
          a.temp <- nu.a+1/2
          
          
          
          for (i in nr.nei.eins)
            {
            diags[coords[i,1]]<-diags[coords[i,1]]-weights[i]
            diags[coords[i,2]]<-diags[coords[i,2]]-weights[i]
          }
          weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
          for (i in nr.nei.eins)
          {
            diags[coords[i,1]]<-diags[coords[i,1]]+weights[i]
            diags[coords[i,2]]<-diags[coords[i,2]]+weights[i]
          }
        }

#        reorg<-1
        for (i in 1:p)
        {
          #print(system.time({
        if(!(nr.nei[i]==1))
          {
            sortier = c((1:d)[-c(coords[i,1],coords[i,2])],coords[i,1])
            if (coords[i,1]!=1){sortier2 = c(1:(coords[i,1]-1))}
              else{sortier2=c()}
            sortier2=c(sortier2,d-1)
            if ((coords[i,1]+1)!=coords[i,2])sortier2=c(sortier2,coords[i,1]:(coords[i,2]-2))
            if (coords[i,2]!=d)
              sortier2=c(sortier2,NA,(coords[i,2]-1):(d-2))
            coeg<-coords[,2]==coords[i,2]|coords[,1]==coords[i,2]
            coeg<-which(!coeg)
         
        
          
            Q <- sparseMatrix(i=c(1:(d-1),sortier2[coords[coeg,1]],sortier2[coords[coeg,2]]),
                              j=c(1:(d-1),sortier2[coords[coeg,2]],sortier2[coords[coeg,1]]),
                              x=c(diags[sortier],-weights[coeg],-weights[coeg])) 
#            reorg.new<-min(which(Q[coords[i,1],]!=0))
#            reorg<-min(reorg,reorg.new)
      
        #print(system.time({
         # if (reorg==1)
         #   {
              L<-chol(Q)
         #   }
         #   else
         #   {
         #     L[reorg:(d-1),reorg:(d-1)] = chol(Q[reorg:(d-1),reorg:(d-1)])        
          #  }
        #}))
        
      
 #       reorg<-reorg.new
            widot = L[d-1,d-1]^2-weights[i]
            sumgjk = widot+sum(L[1:(d-2),d-1]^2)
        
  #      widot2 = sum(weights[Gg[[i]][["nn"]]])-Gg[[i]][["g0"]]
  #      write(x=c(i,weights[i], widot, widot2, sumgjk), file="widot.txt", append=TRUE)
        
        
            b.temp <- (kappa*(gamma[coords[i,1]]-gamma[coords[i,2]])^2/2 + nu.b)
            a.temp <- nu.a+ifelse(runif(1)<(1/(1+widot*nu.b/nu.a)),1,0)
            diags[coords[i,1]]<-diags[coords[i,1]]-weights[i]
            diags[coords[i,2]]<-diags[coords[i,2]]-weights[i]
            weights[i] <-rgamma(1,a.temp,b.temp) 
            diags[coords[i,1]]<-diags[coords[i,1]]+weights[i]
            diags[coords[i,2]]<-diags[coords[i,2]]+weights[i]
    
      
        }
        #}))
        }
      Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
      }#endif method==gibbs
   
#print(system.time({
if (method=="gibbs2")
{
  if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
  {
    nee<-coords[nr.nei.eins,]
    if(is.null(dim(nee)))
    {
      b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
    }
    else
    {
      b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
    }
    a.temp <- nu.a+1/2
    weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
  }
  
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
  Q1 <- solve(Q[-d,-d])
  Q1 <- as.matrix(Q1)
  
  for (i in 1:p)
  {
    if(!(nr.nei[i]==1))
    {
      xy<-as.vector(coords[i,])
      xy<-xy[xy!=d]
      
      Q0<-Q1[xy,xy]
      
      if (length(xy)==2)
      {
        zaehler = diff(Q0)
        zaehler = weights[i]*t(zaehler)%*%zaehler
        qtilde0=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
        Q0=Q0+zaehler/(1-weights[i]*qtilde0)
        qtilde=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
      }
      else
      { 
        zaehler = weights[i]*t(Q0)%*%Q0
        qtilde0=Q0
        qtilde=Q0+zaehler/(1-weights[i]*Q0)
      }
      
      b.temp <- (kappa*(gamma[xy[1]]-gamma[xy[2]])^2/2 + nu.b)
      a.temp <- 1.5
      #print(-1/qtilde)
      weights.proposed <-rstgamma(n=1,shape=a.temp,rate=b.temp,shift=-1/qtilde,a=0,b=1) 
      
      #update Q0  
      qtemp=diff(Q1[xy,])
      qtemp1 = (weights.proposed-weights[i])
      qtemp = qtemp1*t(qtemp)%*%qtemp
      qtemp2 = 1-qtemp1*qtilde0
      Q1 <- Q1 + qtemp/qtemp2
      weights[i]=weights.proposed
    }
  }#endfor
  
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
}#endif method==gibbs2
#}))
#write("\n",file=paste0(path,"/a1.txt"),append=TRUE)

#print(system.time({
if (method=="gibbs2beta")
{
  if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
  {
    nee<-coords[nr.nei.eins,]
    if(is.null(dim(nee)))
    {
      b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2)
    }
    else
    {
      b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2)
    }
    a.temp <- 1
    weights[nr.nei.eins] <-rtrunc(length(nr.nei.eins),"gamma",shape=a.temp,rate=b.temp,a=0,b=1)
  }
  
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
  Q1 <- solve(Q[-d,-d])
  Q1 <- as.matrix(Q1)
  
  for (i in 1:p)
  {
    if(!(nr.nei[i]))
      {
      xy<-as.vector(coords[i,])
      xy<-xy[xy!=d]
      
      Q0<-Q1[xy,xy]
      
      if (length(xy)==2)
      {
        zaehler = diff(Q0)
        zaehler = weights[i]*t(zaehler)%*%zaehler
        qtilde0=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
        Q0=Q0+zaehler/(1-weights[i]*qtilde0)
        qtilde=Q0[1,1]+Q0[2,2]-2*Q0[1,2]
        b.temp <- (kappa*(gamma[xy[1]]-gamma[xy[2]])^2/2)
        qtemp=diff(Q1[xy,])
        
      }
      else
      { 
        zaehler = weights[i]*t(Q0)%*%Q0
        qtilde0=Q0
        qtilde=Q0+zaehler/(1-weights[i]*Q0)
        qtilde=qtilde[1,1]
        b.temp <- (kappa*(gamma[xy]-gamma[d])^2/2)
        qtemp=t(Q1[xy,])
        
      }
        a.temp <- 1.5
        #print(c(b.temp,-1/qtilde,xy,gamma[xy]))
        weights.proposed <-rstgamma(n=1,shape=a.temp,rate=b.temp,shift=-1/qtilde,a=0,b=1) 
      
      #update Q0  
        qtemp1 = (weights.proposed-weights[i])
        qtemp = qtemp1*t(qtemp)%*%qtemp
        qtemp2 = 1-qtemp1*qtilde0
        Q1 <- Q1 + qtemp/qtemp2
        weights[i]=weights.proposed
    }
    }#endfor

diags<-(weights%*%diagM)[1,]
Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
}#endif method==gibbs2
#}))
#write("\n",file=paste0(path,"/a1.txt"),append=TRUE)

if (method=="gibbs3")
{
  if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
  {
    nee<-coords[nr.nei.eins,]
    if(is.null(dim(nee)))
    {
      b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
    }
    else
    {
      b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
    }
    a.temp <- nu.a+1/2
    weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
  }
  
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
  Q1 <- solve(Q[-d,-d])
  Q1 <- as.matrix(Q1)
  
  for (i in 1:length(par.list2))
  {
    #weights[par.list2[[i]]]<-unlist(parallelLapply(par.list2[[i]],single.gibbs.update,coords,d,Q1,weights,kappa,gamma,nu.b))
    weights[par.list2[[i]]]<-unlist(mclapply(par.list2[[i]],single.gibbs.update,coords,d,Q1,weights,kappa,gamma,nu.b))
    diags<-(weights%*%diagM)[1,]
    Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
    if (i!=length(par.list2))
      {
      Q1 <- solve(Q[-d,-d])
      Q1 <- as.matrix(Q1)
    }
  }
}#endif method==gibbs3
#}))
#write("\n",file=paste0(path,"/a1.txt"),append=TRUE)

if (method=="gibbs3beta")
{
  #print(summary(weights))
  if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
  {
    nee<-coords[nr.nei.eins,]
    if(is.null(dim(nee)))
    {
      b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2)
    }
    else
    {
      b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2)
    }
    a.temp <- 1
    weights[nr.nei.eins] <-rtrunc(length(nr.nei.eins),"gamma",shape=a.temp,rate=b.temp,a=0,b=1)
  }
  
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
  Q1 <- solve(Q[-d,-d])
  Q1 <- as.matrix(Q1)
  
  for (i in 1:length(par.list2))
  {
    #weights[par.list2[[i]]]<-unlist(parallelLapply(par.list2[[i]],single.gibbs.update,coords,d,Q1,weights,kappa,gamma,nu.b))
    weights[par.list2[[i]]]<-unlist(mclapply(par.list2[[i]],single.gibbs.update,coords,d,Q1,weights,kappa,gamma,nu.b))
    diags<-(weights%*%diagM)[1,]
    Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
    if (i!=length(par.list2))
    {
      Q1 <- solve(Q[-d,-d])
      Q1 <- as.matrix(Q1)
    }
  }
}#endif method==gibbs3beta

if (method=="gibbs.approx.simple")
      {
        if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
        {
          nee<-coords[nr.nei.eins,]
          if(is.null(dim(nee)))
          {
            b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
          }
          else
          {
            b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
          }
          a.temp <- nu.a+1/2
          weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
        }
        for (i in 1:p)
        {
          if(!(nr.nei[i]==1))
          {
            widot = sum(weights[next.n[[i]]])
            widot = 0.46 + 0.46*widot - 0.01114*widot*widot
            b.temp <- (kappa*(gamma[coords[i,1]]-gamma[coords[i,2]])^2/2 + nu.b)
            a.temp <- nu.a+ifelse(runif(1)<(1/(1+widot*nu.b/nu.a)),1,0)
            weights[i] <-rgamma(1,a.temp,b.temp) 
          }
        }
      }#endif method==gibbs.approx.simple

      if (method=="gibbs.approx")
       {
         if (nr.einer)# nr.nei=1 -> use simple Gibbs sampler
        {
          nee<-coords[nr.nei.eins,]
          if(is.null(dim(nee)))
          {
            b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
          }
          else
          {
            b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
          }
          a.temp <- nu.a+1/2
          weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
        }
        
       for (i in 1:p)
         {
         if(!(nr.nei[i]==1))
          {
             widot = sum((1+(1/Gg[[i]][["g.nn"]]))*weights[Gg[[i]][["nn"]]])+Gg[[i]][["g0"]]
             b.temp <- (kappa*(gamma[coords[i,1]]-gamma[coords[i,2]])^2/2 + nu.b)
             a.temp <- nu.a+ifelse(runif(1)<(1/(1+widot*nu.b/nu.a)),1,0)
             weights[i] <-rgamma(1,a.temp,b.temp) 
           }
         }
       }#endif method==gibbs.approx


      if(method=="approx")
      {
        b.temp <- (kappa*(gamma[coords[,1]]-gamma[coords[,2]])^2/2 + nu.b)
        a.temp <- nu.a+1/2
        weights <-rgamma(p,a.temp,b.temp)
        #diags<-unlist(mclapply(1:d,function(i,weights,coords){sum(weights[coords[,1]==i])+sum(weights[coords[,2]==i])},weights,coords))
        diags<-(weights%*%diagM)[1,]
        Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
      }

if (method=="mh")
{
 # print(system.time({
#   if (nr.einer)# nr.nei=1 -> use Gibbs sampler
#   {
#     nee<-coords[nr.nei.eins,]
#     if(is.null(dim(nee)))
#     {
#       b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
#     }
#     else
#     {
#       b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
#     }
#     a.temp <- nu.a+1/2
#     weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
#   }
  #diags<-unlist(mclapply(1:d,function(i,weights,coords){sum(weights[coords[,1]==i])+sum(weights[coords[,2]==i])},weights,coords))
  diags<-(weights%*%diagM)[1,]
  Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
  #print(system.time({
    ev<-log(diag(chol(Q[-d,-d])))
  #}))
  for (i in p:1)
  {
    #if(!(nr.nei[i]==1))
    {
      b.temp <- (kappa*(gamma[coords[i,1]]-gamma[coords[i,2]])^2/2 + nu.b)
      a.temp <- nu.a+1/2 
      weights.proposed <- weights
      weights.proposed[i] <-rgamma(1,a.temp,b.temp) 
      #Q.proposed <- weightedQ(weights.proposed,coords)  #besser: Nur den Eintrag ändern!!        
      #          Q.proposed <- weightedQ.change(coords[i,],i,diags,weights,weights.proposed)  
      diags.proposed<-diags
      diags.proposed[coords[i,1]]<-diags.proposed[coords[i,1]]-weights[i]+weights.proposed[i]
      diags.proposed[coords[i,2]]<-diags.proposed[coords[i,2]]-weights[i]+weights.proposed[i]
      Q.proposed<-sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags.proposed,-weights.proposed),symmetric=TRUE)
      
      start<-min(coords[i,])
      new.ev<-log(diag(chol(Q.proposed[start:(d-1),start:(d-1)])))
      a<-sum(ev[start:(d-1)]-new.ev)
      a<-exp(a)*sqrt(weights.proposed[i]/weights[i])
      if(runif(1)<a)
      {
        diags<-diags.proposed
        #diags[coords[i,1]]<-diags[coords[i,1]]-weights[i]+weights.proposed[i]
        #diags[coords[i,2]]<-diags[coords[i,2]]-weights[i]+weights.proposed[i]
        
        weights<-weights.proposed
        Q<-Q.proposed
        ev[start:(d-1)]<-new.ev
        
      }
    }
  }
#}))
}#endif mh
 
#print(system.time({
  if (method=="mh2")
      {
        if (nr.einer)# nr.nei=1 -> use Gibbs sampler
        {
          nee<-coords[nr.nei.eins,]
          if(is.null(dim(nee)))
          {
            b.temp <- (kappa*(gamma[nee[1]]-gamma[nee[2]])^2/2 + nu.b)
          }
          else
          {
            b.temp <- (kappa*(gamma[nee[,1]]-gamma[nee[,2]])^2/2 + nu.b)
          }
          a.temp <- nu.a+1/2
          weights[nr.nei.eins] <-rgamma(length(nr.nei.eins),a.temp,b.temp)
        }
        diags<-(weights%*%diagM)[1,]
        Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 
        Q1 <- solve(Q[-d,-d])
        Q0 <- as.matrix(Q1)
        
        for (i in p:1)
        {
            if(!(nr.nei[i]==1))
          {
              b.temp <- (kappa*(gamma[coords[i,1]]-gamma[coords[i,2]])^2/2 + nu.b)
              a.temp <- nu.a
              weights.proposed <-rgamma(1,a.temp,b.temp) 
              
              xy<-as.vector(coords[i,])
              xy<-xy[xy!=d]

              if (length(xy)==2)
                {
                #Q11 = diff(Q1[xy,]) #  v%*%Q1
                #Q111 = diff(Q11[1,xy])
                #a <- 1-(weights.proposed-weights[i])*Q111
                #a <- 1-(weights.proposed-weights[i])*(Q1[xy[1],xy[1]]+Q1[xy[2],xy[2]])
                a <- 1-(weights.proposed-weights[i])*(Q0[xy[1],xy[1]]+Q0[xy[2],xy[2]])
              }
              else
              {
                #a <- 1+(weights.proposed-weights[i])*Q1[xy,xy]
                a <- 1+(weights.proposed-weights[i])*Q0[xy,xy]
              }
              #print(a)
              if(a>1){updated=TRUE}
              else
                if(runif(1)<a){updated=TRUE}

              if (updated==TRUE)
                {
  if (length(xy)==2)
    {
    Q11 = diff(Q0[xy,]) #  v%*%Q1
    Q0 <- Q0 + (weights.proposed-weights[i])*(t(Q11)%*%Q11)/a  
  }
  else
  {
    Q11 = Q0[xy,]
    Q0 <- Q0 - (weights.proposed-weights[i])*(Q11%*%t(Q11))/a
  }
  # update Q1
  weights[i]<-weights.proposed  
}
          }
        }
        
        diags<-(weights%*%diagM)[1,]
        Q <- sparseMatrix(i=c(1:d,coords[,1]),j=c(1:d,coords[,2]),x=c(diags,-weights),symmetric=TRUE) 

}
#}))#endif mh2

}#endif adaptive

    sample<-(iteration-burnin)/thin
    if ((iteration%%thin) == 0 && sample > 0 && sample <= nr.samples){ 
      gamma.sample[sample,] <- gamma
      kappa.sample[sample] <- kappa
      eta.sample[sample,] <- eta
      alpha.sample[sample,] <- eta-gamma
      if(adaptive){
        weights.sample[sample,] <- weights
      }
    }    
     
if(print.time)print((proc.time()-start.time)[3])
}#end iteration
  
returnlist<-list("gamma"=gamma.sample,"kappa"=kappa.sample,"eta"=eta.sample,
                 "alpha"=alpha.sample)
if (adaptive)returnlist[["w"]]=weights.sample
return(returnlist)
  
} 