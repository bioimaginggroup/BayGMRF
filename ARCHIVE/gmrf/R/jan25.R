library(Matrix)
library(BayesLogit)
library(GIGrvg)
source("R/rmvnorm.R"); source("R/rmvnorm-modified.R"); source("R/qtrunc.R"); source("R/rtrunc.R");

N=40000; mcmc=1000; sim=30; nstates<-16;
r<-rep(0.7,nstates); rnot<-rep(0.1,nstates); a<-rep(1,nstates); b<-rep(1,nstates); c<-rep(0.95,nstates); 
s<-0.5; snot<-0.005; as<-1; bs<-0.001; cs=0.95; 

myfunc <- function(nstates,s,as,bs,cs,N,mcmc){
  loc="data/germany.adjacency.arranged";  n <- as.numeric(readLines(loc, n = 1)); nnodes <- id <- numeric(n); adj <- list();
  for (i in 1:n) {
    tmp <- as.numeric(scan(loc, skip=i, nlines = 1, quiet = T, what = list(rep("", 13)))[[1]])
    id[i] <- tmp[1]; nnodes[i] <- tmp[2]; adj[[i]] <- tmp[-c(1:2)]
  }
  adj1 <- as.integer(unlist(adj)); col <- c(1:n,  adj1); row <- c(1:n,  rep(1:n,nnodes));
  entries <- c(nnodes, rep(-1,length(adj1)))
  remove <- col<row; entries <- entries[!remove]
  col <- col[!remove]
  row <- row[!remove] 
  lrq = sparseMatrix(i=row,j=col,x=entries,symmetric=TRUE)
  frq <- lrq+0.5*diag(n);
  
  states <- list(); states[[1]] <- 1:15; states[[2]] <- 16; states[[3]] <- 17:63; states[[4]] <- 64:65;
  states[[5]] <- 66:119; states[[6]] <- 120:145; states[[7]] <- 146:181; states[[8]] <- 182:225;
  states[[9]] <- 226:321; states[[10]] <- 322:327; states[[11]] <- 328:329; states[[12]] <- 330:373;
  states[[13]] <- 374:410; states[[14]] <- 411:464; states[[15]] <- 465:504; states[[16]] <- 505:544;

  qur <- sparseMatrix(i=row,j=col,x=0,symmetric=TRUE);
  for(i in 1:nstates){
    qur[states[[i]],states[[i]]] <- r[i]*lrq[states[[i]],states[[i]]]
  }
  qu <- qur+s*lrq+0.5*diag(n); 

  qur2<-vector(mode="list",length=16)
  for(i in 1:nstates){
    qur2[[i]] <- sparseMatrix(i=row,j=col,x=0,symmetric=TRUE);
    qur2[[i]][states[[i]],states[[i]]] <- lrq[states[[i]],states[[i]]]
    for (j in states[[i]])
      qur2[[i]][j,j]<-qur2[[i]][j,j]-sum(qur2[[i]][j,])+1
  }
  qur_complete <- qur2[[1]]
  for (p in 2:nstates)qur_complete<-qur_complete+qur2[[p]]
  
  x <- xmeans <- mser <- matrix(data=0, nrow=sim, ncol=n); mset <- cases80 <- cases90 <- cases95 <- numeric(sim);
# TODO: specify different values for b, e.g.:
  be <- rep(0,n); be[328:544]<-1
  x <- rmvnorm(1, b=be, P=qu);  p <- 1/(1+exp(-x));  y <- rbinom(n,N,p);
  si <- rks(n); lambdaprior <- 1/(2*si)^2; wprior <- rnorm(n,x,1/lambdaprior);
  ynull <- numeric(n); ynull[wprior>0] <- y[wprior>0]; y <- ynull; wpriorvec <- matrix(wprior,ncol=1,byrow=TRUE);
  mu <- diag(lambdaprior)%*%wpriorvec; prec <- qu+diag(lambdaprior);
  
### x, w and lambda from posterior distribution
  xpost <- rmvnormmodi(n, b=mu, P=prec); wpostnull <- numeric(n);
  wpostpve <- rtrunc(1, spec="logis", a=0, b=Inf, location=xpost, scale=1);
  wpostnve <- rtrunc(1, spec="logis", a=-Inf, b=0, location=xpost, scale=1);
  wpostnull[y>0] <- wpostpve[y>0]; wpostnull[y<1] <- wpostnve[y<1]; wpost <- wpostnull;
  lambdapost <- rgig(n, lambda=0.5, chi=1, psi=(wpost-xpost)^2);

### generating r proposal, acception/rejection based on alpha    
  rprop<-numeric(nstates); qu_r_not<-qu_r_prop<-sparseMatrix(i=row,j=col,x=0,symmetric=TRUE);
  for(i in 1:nstates) { 
    rprop[i] <- rgamma(1,a[i]+c[i],as.numeric(b[i]+0.5*t(as.vector(xpost[states[[i]]]))%*%lrq[states[[i]],states[[i]]]%*%as.vector(xpost[states[[i]]])));
    qu_r_not[states[[i]],states[[i]]] <- rnot[i]*lrq[states[[i]],states[[i]]]; qu_r_prop[states[[i]],states[[i]]] <- rprop[i]*lrq[states[[i]],states[[i]]];
  }
  qur_not<-qu_r_not+s*lrq+0.5*diag(n); qur_prop<-qu_r_prop+s*lrq+0.5*diag(n);
  ###prod(c) because all cis are same
  log_apr <- prod(c)*log(prod(rnot))-prod(c)*log(prod(rprop))+0.5*determinant(qur_prop, logarithm=TRUE)[[1]][1]-0.5*determinant(qur_not, logarithm=TRUE)[[1]][1]; apr<-exp(log_apr);
  ur <- runif(1); rpost<-matrix(0,1,nstates);    if(ur < apr){rpost <- rprop} else {rpost <- rnot}
  
### putting r values in precision matrix partly  
  qur_post<-sparseMatrix(i=row,j=col,x=0,symmetric=TRUE); 
  for(i in 1:nstates) {
    qur_post[states[[i]],states[[i]]] <- rpost[i]*lrq[states[[i]],states[[i]]];
  }
  
### generating s proposal, acception/rejection based on alpha      
  sprop <- rgamma(1,as+cs,as.numeric(bs+0.5*t(as.vector(xpost))%*%lrq%*%as.vector(xpost)));
  qus_prop<-qur_post+sprop*lrq+0.5*diag(n); qus_not<-qur_post+snot*lrq+0.5*diag(n);
  log_aps <- cs*log(snot)-cs*log(sprop)+0.5*determinant(qus_prop, logarithm=TRUE)[[1]][1]-0.5*determinant(qus_not, logarithm=TRUE)[[1]][1]; aps<-exp(log_aps);
  us <- runif(1); if(us < aps){spost <- sprop} else {spost <- snot}
  
  xvec <- wvec <- lambdavec <- matrix(data=0, nrow=10, ncol=n);
  svec <- matrix(data=0, nrow=10, ncol=1); rvec <- matrix(data=0, nrow=10, ncol=nstates);
  xvec[1,] <- xpost; wvec[1,] <- wpost; lambdavec[1,] <- lambdapost; svec[1,] <- spost; rvec[1,] <- rpost;
  counts<-0; reqno<-mcmc/10; xmcmc <- wmcmc <- lambdamcmc <- matrix(data=0, nrow=mcmc/10, ncol=n);
  smcmc <- matrix(data=0, nrow=mcmc/10, ncol=1); rmcmc <- matrix(data=0, nrow=mcmc/10, ncol=nstates);
  
  qur_vec<-rvec[1,1]*qur2[[1]]
  for (p in 2:nstates)
  {
    qur_r_vec<-qur_vec+rvec[1,p]*qur2[[p]]
  }
  
  for(counts in 1:reqno) {   
  for(j in 1:9) {
    mu <- diag(lambdavec[j,])%*%wvec[j,]
    
    # not needed her, qu_r_vec is still valid
    #qu_r_vec <- sparseMatrix(i=row,j=col,x=0,symmetric=TRUE);
    #for(p in 1:nstates) {
    #  qu_r_vec[states[[p]],states[[p]]] <- rvec[j,p]*lrq[states[[p]],states[[p]]];
    #}

    qu_vec<-qu_r_vec+svec[j,]*lrq+0.5*diag(n)
    
    prec <- qu_vec+diag(lambdavec[j,]);
    xvec[j+1,] <- rmvnormmodi(n, b=mu, P=prec); 
    # wpostnull <- numeric(n);
    # update w
    wpostpve <- rtrunc(1, spec="logis", a=0, b=Inf, location=xvec[j+1,y>0], scale=1)
    wpostnve <- rtrunc(1, spec="logis", a=-Inf, b=0, location=xvec[j+1,y<1], scale=1);
    wpostnull[y>0] <- wpostpve
    wpostnull[y<1] <- wpostnve
    wvec[j+1,] <- wpostnull; 
    # update lambda
    lambdavec[j+1,] <- rgig(n, lambda=0.5, chi=1, psi=(wvec[j+1,]-xvec[j+1,])^2);
    
    # update r
    rprop<-r_prev<-rvec[j,]
    qu_r_prop<-qu_vec
    
    for(i in 1:nstates) { 
      rprop[i] <- rgamma(1,a[i]+c[i],as.numeric(b[i]+0.5*t(as.vector(xvec[j+1,states[[i]]]))%*%qur2[[i]][states[[i]],states[[i]]]%*%as.vector(xvec[j+1,states[[i]]])))
      qu_r_prop<-qu_r_prop+(rprop[i]-r_prev[i])*qur2[[i]]
      log_apr <- c[i]*log(r_prev[i])-c[i]*log(rprop[i])+0.5*determinant(qu_r_prop, logarithm=TRUE)[[1]][1]-0.5*determinant(qur_prev, logarithm=TRUE)[[1]][1]
      apr<-exp(log_apr)
      ur <- runif(1)
      if(ur < apr){
        r_prev <- rprop
        qu_r_prev <- qu_r_prop
      } else {
          rprop <- r_prev
          qu_r_prop <- qu_r_prev
      }
    }
    rvec[j+1,] <- r_prev
    
    #qur_vec<-sparseMatrix(i=row,j=col,x=0,symmetric=TRUE);  for(i in 1:nstates) {
    #  qur_vec[states[[i]],states[[i]]] <- rvec[j+1,i]*lrq[states[[i]],states[[i]]];
    #}
    # faster: 
    qur_vec<-rvec[j,1]*qur2[[1]]
    for (p in 2:nstates)
    {
      qur_vec<-qur_vec+rvec[j,p]*qur2[[p]]
    }
    
    s_prev<-svec[j,]
    sprop <- rgamma(1,as+cs,as.numeric(bs+0.5*t(as.vector(xvec[j+1,]))%*%lrq%*%as.vector(xvec[j+1,])))
    qus_prop<-qur_vec+sprop*lrq+0.5*diag(n)
    qus_prev<-qur_vec+s_prev*lrq+0.5*diag(n)
    log_aps <- cs*log(s_prev)-cs*log(sprop)+0.5*determinant(qus_prop, logarithm=TRUE)[[1]][1]-0.5*determinant(qus_prev, logarithm=TRUE)[[1]][1]
    aps<-exp(log_aps)
    us <- runif(1)
    if(us < aps){svec[j+1,] <- sprop} else {svec[j+1,] <- s_prev}
    
    if(counts==reqno) {xmcmc[counts,] <- xvec[10,]; wmcmc[counts,] <- wvec[1,] <- wvec[10,];
    lambdamcmc[counts,] <- lambdavec[1,] <- lambdavec[10,]; smcmc[counts,] <- svec[1,] <- svec[10,]; rmcmc[counts,] <- rvec[1,] <- rvec[10,]
    } else 
      if(j==9) {xmcmc[counts,] <- xvec[1,] <- xvec[10,]; wmcmc[counts,] <- wvec[1,] <- wvec[10,];
      lambdamcmc[counts,] <- lambdavec[1,] <- lambdavec[10,]; smcmc[counts,] <- svec[1,] <- svec[10,]; rmcmc[counts,] <- rvec[1,] <- rvec[10,]
      }
  }                      
  }
  
  xmeans<-apply(xmcmc[101:1000,], 2, mean);  mset<-mean((xmeans-x)^2);
  mser<-apply(xmcmc[101:1000,], 2, function(x)  {sumx<-sum(x); lengthx<-length(x); meanx<-sumx/lengthx; msex<-sum((x-meanx)^2)/lengthx})
  CIs <- apply(xmcmc[101:1000,], 2, quantile, probs=c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975)); cases95<-length(which(CIs[1,] < x & x < CIs[6,]));
  cases90<-length(which(CIs[2,] < x & x < CIs[5,]));   cases80<-length(which(CIs[3,] < x & x < CIs[4,]));
  return(list(x,xmcmc,wmcmc,lambdamcmc,smcmc,rmcmc,xmeans,mser,mset,cases80,cases90,cases95,r,rnot))}

filenames<-c("1.RData","2.RData","3.RData","4.RData","5.RData","6.RData","7.RData","8.RData","9.RData","10.RData",
             "11.RData","12.RData","13.RData","14.RData","15.RData","16.RData","17.RData","18.RData","19.RData","20.RData",
             "21.RData","22.RData","23.RData","24.RData","25.RData","26.RData","27.RData","28.RData","29.RData","30.RData");

for(i in 1:sim){ emp<-myfunc(nstates,s,as,bs,cs,N,mcmc); etc <- warnings(); out <- capture.output(etc);
cat("Errors", out, file="etc.txt", sep="\n ", append=TRUE); save(emp, file=filenames[i], compress=FALSE)
}

