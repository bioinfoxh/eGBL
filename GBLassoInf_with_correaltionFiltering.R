

get_corGenesIndex <- function(Y,X,corMethod){
  temp_p = apply( X, 2 , function(x) { cor.test( x, Y,method=corMethod,exact=F)$p.value } )
  temp_r = apply( X, 2 , function(x) { cor.test( x, Y,method=corMethod,exact=F)$estimate } )
  cor_Y = cbind( temp_r, temp_p )
  rownames(cor_Y) = colnames(X)
  corGene_Y = cor_Y[ cor_Y[,2]<=0.05,]
#  colnames(corGene_Y) = c("corR","p_value")
#  X_corGene = X[, colnames(X) %in% rownames(corGene_Y) ]
#  X_corGene_vid = which(  colnames(X) %in% colnames(X_corGene)   )
  X_corGene_vid = which( colnames(X) %in% rownames(corGene_Y) )
  return( X_corGene_vid )
}



GBLasso_corGenes_Wdegree <- function(Y,X,corMethod,Graph_X,gamma,epsilon,ksi,MAXITER,MAXLOOP){
  
  X_corGene_vid = get_corGenesIndex(Y,X,corMethod)
  
  print( paste( length(X_corGene_vid), corMethod, "correlated genes identified", sep=" " )  )
  
  network_g_Y = induced.subgraph(Graph_X, X_corGene_vid )
  netwk_Y = get.edgelist(network_g_Y, names=FALSE)
  temp = get.data.frame(network_g_Y)
  netwk_Y = cbind( netwk_Y, temp[,3]  )
  
  ppi = get.data.frame( Graph_X )
  
  temp1 = as.data.frame( table(c(ppi[,1],ppi[,2])), stringsAsFactors = F)
  temp2 = temp1[ match( V(network_g_Y)$name, temp1[,1] ) , ]
  degree_node = temp2[,2]
  degree_node[is.na(degree_node)]=0
  names(degree_node) = V(network_g_Y)$name
  
  Y0=Y
  X0=X[ , X_corGene_vid  ]
  
  res_GBLasso = NetworkRegression_GBLasso(Y0, X0, netwk_Y, gamma, degree_node, epsilon, ksi, MAXITER, MAXLOOP )
  
  corGenes = colnames(X)[X_corGene_vid]
  network = netwk_Y
  weight_node = degree_node
  
  result = list( corGenes=corGenes,
                 network = netwk_Y,
                 weight_node = degree_node,
                 GBLassoFit = res_GBLasso )
  
}



GBLasso_corGenes_Wcoexp <- function(Y,X,corMethod,Graph_X,gamma,epsilon,ksi,MAXITER,MAXLOOP){
  
  X_corGene_vid = get_corGenesIndex(Y,X,corMethod)
  
  print( paste( length(X_corGene_vid), corMethod, "correlated genes identified", sep=" " )  )
  
  network_g_Y = induced.subgraph(Graph_X, X_corGene_vid )
  netwk_Y = get.edgelist(network_g_Y, names=FALSE)
  temp = get.data.frame(network_g_Y)
  netwk_Y = cbind( netwk_Y, temp[,3]  )
  
  ppi = get.data.frame( Graph_X )
  gene_exp = colnames(X)
  
  temp =  ppi[, c(2,1,3)] 
  colnames(temp) = colnames(ppi)
  temp = rbind(ppi,temp)
  
  temp1 = aggregate( data.frame(abs(temp[,3]),stringsAsFactors = F), list(factor( temp[,1],levels = gene_exp)) , sum  )
  temp2 = data.frame( gene=as.character(temp1[,1]), w=temp1[,2], stringsAsFactors = F )
  temp3 = temp2[ match( V(network_g_Y)$name, temp2[,1] ) , ]
  weight_coexp = temp3[,2]
  weight_coexp[is.na(weight_coexp)]=0
  names(weight_coexp) = V(network_g_Y)$name
  
  Y0=Y
  X0=X[ , X_corGene_vid  ]
  
  res_GBLasso = NetworkRegression_GBLasso(Y0, X0, netwk_Y, gamma, weight_coexp, epsilon, ksi, MAXITER, MAXLOOP )
  
  corGenes = colnames(X)[X_corGene_vid]
  network = netwk_Y
  weight_node = weight_coexp
  
  result = list( corGenes=corGenes,
                 network = netwk_Y,
                 weight_node = weight_coexp,
                 GBLassoFit = res_GBLasso )
}






