
library(igraph, quietly = T)
library(glmnet, quietly = T)

source("/Users/Hui/Dropbox/my/Rscripts/grace/cvGrace_1.R")
source("/Users/Hui/Dropbox/my/Rscripts/grace/grace.R")
source("/Users/Hui/Dropbox/my/Rscripts/grace/predictGrace.R")


run_simulation_Grace <- function(lambda.1,lambda.2,lambda.L,simTimes){ 

  sim100_res_Grace=c()
  
  for(i in 1:simTimes){
    
    print( paste("run simulation data ",i,sep="") )
  
    simData_file = paste("/Users/Hui/myFiles/cambridge/projects/GBLasso/setup2-110feature/simulateData/simData_",i,".RData",sep="")
    load(simData_file)

    Xtr = simData$Xtr
    Ytr = simData$Ytr
    b = simData$betas
    netwk = simData$netwk
    
    Xts = simData$Xts
    Yts = simData$Yts
    
    netwk_g = graph_from_data_frame(as.data.frame(netwk[,1:2]), directed = F, vertices = 1:ncol(Xtr))
    edge_attr(netwk_g,'weight') <- abs(netwk[,3])
    netwk_adj = as_adjacency_matrix(netwk_g, type="both", names=T, sparse = F, attr = "weight")
    

    graceFit = grace(Ytr, Xtr, netwk_adj, lambda.L, lambda.1, lambda.2, normalize.L=FALSE, K=10  )
    
    best_lambda.1 = graceFit$ParameterEstimation$parameterMin[1]
    best_lambda.L = graceFit$ParameterEstimation$parameterMin[2]
    best_lambda.2 = graceFit$ParameterEstimation$parameterMin[3]
    
    preBeta = graceFit$Beta$beta
    
    preYtr = predictGrace( Xtr, graceFit )
    preYts = predictGrace( Xts, graceFit )
    
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
                  precision,recall,corPCC,corSPM,best_lambda.1,best_lambda.2,best_lambda.L,pmse_Ytr,pmse_Yts )
    
    sim100_res_Grace = rbind(sim100_res_Grace, res_stat )
    
    RDataFile = paste("Grace_sim_",i,".RData",sep="")
    save(graceFit,file=RDataFile)
    
  }

  colnames(sim100_res_Grace)=c("GSP","GSN","TP","FP","preP","TN","FN","preN","precision","recall","corPCC","corSPM","lambda.1","lambda.2","lambda.L","PMSE","MSEpredict")
  
  return(sim100_res_Grace)
  
}








