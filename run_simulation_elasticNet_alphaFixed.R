

run_simulation_elasticNet <- function(alpha,simTimes){ 
  
  library(glmnet,quietly = T)
  
  sim100_res_elasticNet=c()
  
  for(i in 1:simTimes){

  simData_file = paste("/Users/Hui/myFiles/cambridge/projects/GBLasso/setup2-110feature/simulateData/simData_",i,".RData",sep="")
  load(simData_file)
  

  Xtr = simData$Xtr
  Ytr = simData$Ytr
  b = simData$betas
  
  Xts = simData$Xts
  Yts = simData$Yts



  glmnetCV = cv.glmnet(Xtr,Ytr,alpha = alpha)
  lambdaCV = glmnetCV$lambda.min
  glmnetFit = glmnet(Xtr,Ytr,alpha = alpha, lambda=lambdaCV)
  
  betaPredict = as.matrix( predict(glmnetFit, type = "coefficient") )[,1]
  
  preIntercept = betaPredict[1]
  preBeta = betaPredict[2:length(betaPredict) ]
  
  preYtr = predict( glmnetFit, type="response", newx = Xtr )
  preYts = predict( glmnetFit, type="response", newx = Xts )
  
  
  delta_Ytr = preYtr - Ytr
  pmse_Ytr = sqrt(  sum( delta_Ytr * delta_Ytr  ) / length(delta_Ytr)  )
  
  delta_Yts = preYts - Yts
  pmse_Yts = sqrt(  sum( delta_Yts * delta_Yts  ) / length(delta_Yts)  )
  

  GSP = which(b!=0)
  GSN = which(b==0)
  
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
  
  res_stat = c( length(GSP),length(GSN),length(TP),length(FP),length(preP),length(TN),length(FN),length(preN),
                precision,recall,corPCC,corSPM,lambdaCV,pmse_Ytr,pmse_Yts )

  sim100_res_elasticNet = rbind(sim100_res_elasticNet, res_stat )
  
  } 
  
  colnames(sim100_res_elasticNet)=c("GSP","GSN","TP","FP","preP","TN","FN","preN","precision","recall","corPCC","corSPM","lambda","PMSE","MSEpredict")
  
  return(sim100_res_elasticNet)
}
