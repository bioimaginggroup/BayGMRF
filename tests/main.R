#main<-function(data,x)
#{
  require(Matrix)
  dims<-dim(data)
  Q1<-sparseMatrix(i = c(1:dims[1],2:dims[1]),j=c(1:dims[1],1:(dims[1]-1)),x = c(1,rep(2,dims[1]-2),1,rep(-1,dims[1]-1)),symmetric = TRUE)
  Q2<-sparseMatrix(i = c(1:dims[2],2:dims[2]),j=c(1:dims[2],1:(dims[2]-1)),x = c(1,rep(2,dims[2]-2),1,rep(-1,dims[2]-1)),symmetric = TRUE)
  Q<-kronecker(Q2,diag(nrow=dims[1]))+kronecker(diag(nrow=dims[2]),Q1)

  tau<-1



  for (iter in 1:iterations)
  {

  }
#}


  Q.from.c<-function(c)
  {
    dims<-dim(c)
    Q1<-sparseMatrix(i = c(1:dims[1],2:dims[1]),j=c(1:dims[1],1:(dims[1]-1)),x = c(1,rep(2,dims[1]-2),1,rep(-1,dims[1]-1)),symmetric = TRUE)
    Q2<-sparseMatrix(i = c(1:dims[2],2:dims[2]),j=c(1:dims[2],1:(dims[2]-1)),x = c(1,rep(2,dims[2]-2),1,rep(-1,dims[2]-1)),symmetric = TRUE)
    Q<-kronecker(Q2,diag(nrow=dims[1]))+kronecker(diag(nrow=dims[2]),Q1)
    for (i in 1:dims[1])
      for (j in 1:dims[2])
        if (c[i,j]==0)
        {
          if (i!=1){
            Q[index(i,j,dims[1]),index(i-1,j,dims[1])]=Q[index(i,j,dims[1]),index(i-1,j,dims[1])]+1
            Q[index(i-1,j,dims[1]),index(i,j,dims[1])]=Q[index(i-1,j,dims[1]),index(i,j,dims[1])]+1
            Q[index(i-1,j,dims[1]),index(i-1,j,dims[1])]=Q[index(i-1,j,dims[1]),index(i-1,j,dims[1])]-1
            Q[index(i,j,dims[1]),index(i,j,dims[1])]=Q[index(i,j,dims[1]),index(i,j,dims[1])]-1
          }
          # if (i!=dims[1]){
          #   Q[index(i,j,dims[1]),index(i+1,j,dims[1])]=0
          #   Q[index(i+1,j,dims[1]),index(i,j,dims[1])]=0
          #   Q[index(i+1,j,dims[1]),index(i+1,j,dims[1])]=Q[index(i+1,j,dims[1]),index(i+1,j,dims[1])]-1
          #   Q[index(i,j,dims[1]),index(i,j,dims[1])]=Q[index(i,j,dims[1]),index(i,j,dims[1])]-1
          # }
          if (j!=1){
            Q[index(i,j,dims[1]),index(i,j-1,dims[1])]=Q[index(i,j,dims[1]),index(i,j-1,dims[1])+1
            Q[index(i,j-1,dims[1]),index(i,j,dims[1])]=Q[index(i,j-1,dims[1]),index(i,j,dims[1])]+1
            Q[index(i,j-1,dims[1]),index(i,j-1,dims[1])]=Q[index(i,j-1,dims[1]),index(i,j-1,dims[1])]-1
            Q[index(i,j,dims[1]),index(i,j,dims[1])]=Q[index(i,j,dims[1]),index(i,j,dims[1])]-1
          }
          # if (j!=dims[2]){
          #   Q[index(i,j,dims[1]),index(i,j+1,dims[1])]=0
          #   Q[index(i,j+1,dims[1]),index(i,j,dims[1])]=0
          #   Q[index(i,j+1,dims[1]),index(i,j+1,dims[1])]=Q[index(i,j+1,dims[1]),index(i,j+1,dims[1])]-1
          #   Q[index(i,j,dims[1]),index(i,j,dims[1])]=Q[index(i,j,dims[1]),index(i,j,dims[1])]-1
          # }
        }#end loop
    return(Q)
  }


        }

        }

    index<-function(i,j,dim1){return((j-1)*dim1+i)}
