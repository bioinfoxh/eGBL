


predictGBLasso <- function(Yts, Xts, result ){
  if( is.null(result) ){
    preResult = NA
  }else{
    beta = result$BetaSelection$betas
    intercept = result$BetaSelection$intercept
    
    Ypre = Xts %*% as.matrix(beta) + intercept
    deltaY = Ypre - Yts
    mse = sqrt(  sum( deltaY * deltaY  ) / length(deltaY)  )
    
    preResult = list(Ypredict=Ypre, MSEpredict = mse)
  } 
  return(preResult)
}



get_result_statistic <- function(b, result, Yts, Xts){
  GSP = which(b!=0)
  GSN = which(b==0)
  
  if( is.null(result$BetaSelection) ){
    res_stat = c(length(GSP),length(GSN),rep(NA,13))
  }else{
    preBeta = result$BetaSelection$betas
    pmse = result$BetaSelection$pmse
    lambda = result$BetaSelection$lambda
    
    preP = which( preBeta!=0 )
    preN = which( preBeta==0 )
    
    TP = intersect(preP, GSP)
    FP = intersect(preP, GSN)
    TN = intersect(preN, GSN)
    FN = intersect(preN, GSP)

    if(length(TP)==0){
        precision = 0
        recall = 0
        corPCC = 0
        corSPM = 0
    }else{
        precision = length(TP)/length(preP)
        recall = length(TP)/length(GSP)
        corPCC = cor( b, preBeta, method = "pearson")
        corSPM = cor( b, preBeta, method = "spearman")
    }
    
    MSEpredict = predictGBLasso(Yts, Xts, result )$MSEpredict
    
    res_stat = c( length(GSP),length(GSN),length(TP),length(FP),length(preP),length(TN),length(FN),length(preN),
                  precision,recall,corPCC,corSPM,lambda,pmse,MSEpredict )
  }
  return(res_stat)
}





GBLasso_simulation <- function(gamma = 2, epsilon=0.01, ksi=0.01, simTimes=100 ){
  
  for( nSim in 1:simTimes){
    
    print( paste( "running simulation",nSim,sep=" " )  )
    
    simData_file = paste("/media/hx239/Disk1/projects/GBLasso/setup2_110feature/simulateData/simData_",nSim,".RData",sep="")
    load(simData_file)
    
 #   simData = simulateData(simSeed,0)
    Ytr = simData$Ytr
    Xtr = simData$Xtr
    Yts = simData$Yts
    Xts = simData$Xts
    netwk = simData$netwk
    b = simData$betas
    
    temp = rbind( netwk, netwk[, c(2,1,3)]  )
    
    temp1 = data.frame( table(temp[,1]))
    temp1[,1] = as.numeric(temp1[,1])
    temp1 = temp1[ order(temp1[,1])  ,]
    degree = temp1[,2]
    
    temp1 = aggregate( data.frame(abs(temp[,3]),stringsAsFactors = F), list(factor( temp[,1],levels = 1:ncol(Xtr) ) ), sum )
    temp2 = data.frame( node=as.character(temp1[,1]), w=temp1[,2], stringsAsFactors = F )
    weightSum = temp2[,2]
    
    # gamma = 2
    # epsilon=0.001
    # ksi=0.001
    
    
    source("/media/hx239/Disk1/projects/GBLasso/GBLasso_code/NetworkRegression_GBLasso.R")
    res_D=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # res_D_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # res_D_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    
    # res_WS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # res_WS_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # res_WS_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    
    
    source("/media/hx239/Disk1/projects/GBLasso/GBLasso_code/NetworkRegression_GBLasso_EB.R")
    # resEB_D=NetworkRegression_GBLasso_EB(Ytr, Xtr, netwk, gamma, degree, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # resEB_D_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # resEB_D_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    
    resEB_WS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # resEB_WS_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=5000)
    # resEB_WS_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    
    
    # source("/Users/Hui/Dropbox/my/Rscripts/GBLasso/NetworkRegression_GBLasso_EBsqrt.R")
    # resEBs_D=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # resEBs_D_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # resEBs_D_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, degree^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # 
    # resEBs_WS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # resEBs_WS_G=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^gamma, epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    # resEBs_WS_GS=NetworkRegression_GBLasso(Ytr, Xtr, netwk, gamma, weightSum^((gamma+1)/2), epsilon, ksi, MAXITER=1000, MAXLOOP=10000)
    
    
    list_res = ls(pattern = "^res_|^resEB_|^resEBs_")
    
    temp_result_stat = c()
    for( i in 1:length(list_res)  ){
      
      temp_stat = get_result_statistic(b, get(list_res[i]), Yts, Xts)
      temp_result_stat = rbind( temp_result_stat, temp_stat  )
    }
    result_stat = data.frame( list_res, temp_result_stat, stringsAsFactors = F  )
    colnames(result_stat)=c("model","GSP","GSN","TP","FP","preP","TN","FN","preN","precision","recall","corPCC","corSPM","lambda","PMSE","MSEpredict")
    
    RDataFile = paste( "GBLasso_sim_",nSim,".RData",sep=""  )
    result_stat_file = paste("resStat_GBLasso_sim_",nSim,".txt",sep="")
    save.image(RDataFile)
    write.table( result_stat, result_stat_file, sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE)
  }
}













