

Xb<-function(x, b){
  sum(x*b)
}


simulateData <- function(seedNum){
  
  set.seed(seedNum)
  
  ##generate simulated data
  
  # number of records
  gamma0<-0
  Ntr<-50
  Nts<-100
  

  nTF<-50

  
  # number of target genes
  nTarget<-10
  
  # generate the netwok information
  netwk<-matrix(0,nTF*nTarget,2)
  for(i in 1:nTF) {
    for(j in 1:nTarget) {
      netwk[(i-1)*nTarget+j,1] = (i-1)*(nTarget+1)+1
      netwk[(i-1)*nTarget+j,2] = (i-1)*(nTarget+1)+1+j
    }
  }
  
  corR =  seq(from=0.95,to=0.05,by=-0.1) 
  corR_edges = rep(corR, nTF)
  netwk = cbind(netwk,corR_edges)
  
  b<-c(9, rep(9/sqrt(nTarget),nTarget),
       -7, rep(-7/sqrt(nTarget),nTarget), 
       5, rep(5/sqrt(nTarget),nTarget), 
       -3, rep(-3/sqrt(nTarget),nTarget),
       1, rep(1/sqrt(nTarget),nTarget),
       -9, rep(-9/sqrt(nTarget),nTarget),
       7, rep(7/sqrt(nTarget),nTarget),
       -5, rep(-5/sqrt(nTarget),nTarget),
       3, rep(3/sqrt(nTarget),nTarget),
       -1, rep(-1/sqrt(nTarget),nTarget),
       rep(0, (nTarget+1)*(nTF-10))  )
  
  Ve<-sum(b*b)/2
  
  # training data
  Ytr<-rep(0, Ntr)
  Xtr<-matrix(0, nrow=nTF*(nTarget+1), ncol=Ntr)
  
  j<-1
  for(i in 1:nTF){
    Xtr[j,]<-rnorm(Ntr, 0, 1)
    j<-j+1
    for(k in 1:nTarget){
      Xtr[j,]<-rnorm(Ntr, corR[k]*Xtr[(i-1)*(nTarget+1)+1,], sqrt( 1-corR[k]^2 ) )
      j<-j+1
    }
  }
  Ytr<-rnorm(Ntr, 0, sqrt(Ve)) + apply(Xtr, 2, Xb, b)
  Xtr = t(Xtr)
  
  # testing data
  Yts<-rep(0, Nts)
  Xts<-matrix(0, nrow=nTF*(nTarget+1), ncol=Nts)
  
  j<-1
  for(i in 1:nTF){
    Xts[j,]<-rnorm(Nts, 0, 1)
    j<-j+1
    for(k in 1:nTarget){
      Xts[j,]<-rnorm(Nts, corR[k]*Xts[(i-1)*(nTarget+1)+1,], sqrt( 1-corR[k]^2 ) )
      # Xts[j,]<-rnorm(Nts, corR[k]*Xts[(i-1)*(nTarget+1)+1,] )
      j<-j+1
    }
  }
  Yts<-rnorm(Nts, 0, sqrt(Ve)) + apply(Xts, 2, Xb, b)
  Xts <- t(Xts)
  
  simData = list( Xtr=Xtr, Ytr=Ytr, Xts=Xts, Yts=Yts , netwk = netwk , betas=b)
  return(simData)
}



