

GBLasso_correlationFiltering_result_statistic<-function( result_RData ){

  load( result_RData)
  
  score_list = c("EXP","ICM","TE","T5","CC2","S2","CC3")


result_statistic_all = c()
for (i in 1:length(score_list)){
  
  temp_result = get( paste("result",score_list[i],sep="_")   )
  
  temp_result_statistic=c()
  for (j in 1:length( temp_result )){
#    print( paste( "i=",i,", j=",j,sep=""    )  )
    
    temp_result_j = temp_result[[j]]
    num_corGene = length( temp_result_j$corGenes )
    num_netEdge = nrow( temp_result_j$network )
    
    if( is.null( temp_result_j$GBLassoFit$BetaSelection  )   ){
      lambda = NA
      num_beta = NA
      pmse = NA
    }else{
      lambda = temp_result_j$GBLassoFit$BetaSelection$lambda
      num_beta = sum( temp_result_j$GBLassoFit$BetaSelection$betas!=0 )
      pmse = temp_result_j$GBLassoFit$BetaSelection$pmse
    }
 #   print( c(num_corGene,num_netEdge,lambda,num_beta,pmse ) )
    temp_result_statistic = rbind( temp_result_statistic, c(num_corGene,num_netEdge,lambda,num_beta,pmse )  )
  }
  param = names(temp_result)
  
  score_idx = rep( score_list[i], length(temp_result)  )
  
  temp_result_statistic = data.frame( score_idx, param, temp_result_statistic, stringsAsFactors = F )
  
  result_statistic_all = rbind( result_statistic_all, temp_result_statistic )
  
}

colnames(result_statistic_all)=c( "score","GBLassoParam","num_corGenes","num_netEdges","lambda","num_betas","PMSE")

return( result_statistic_all )

}


