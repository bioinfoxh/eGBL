

source("/Users/Hui/Dropbox/my/Rscripts/GBLasso/NetworkRegression_GBLasso_EB.R")
source("GBLassoInf_with_correaltionFiltering.R")
source("GBLassoInf_batch_correlationFiltering.R")

library(igraph,quietly = TRUE)

network_file_sig="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_network_analysis/coexp_network_significant.txt"
network_file_all="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_network_analysis/coexp_network_all.txt"
exp_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/parameter_analysis/gene_expression_data/gene_expression_RMA_call_PM1.txt"
parameter_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/parameter_analysis/gene_parameter_association/arrayID_mapping_score.txt"  
  

ppi = read.table(network_file_sig,header = TRUE, stringsAsFactors = F ) 
ppi[,1] = as.character(ppi[,1])
ppi[,2] = as.character(ppi[,2])
ppi_sig = ppi

ppi = read.table(network_file_all,header = TRUE, stringsAsFactors = F ) 
ppi[,1] = as.character(ppi[,1])
ppi[,2] = as.character(ppi[,2])
ppi_all = ppi

rm(ppi)

expRMA_gene = read.table(exp_file, header = TRUE , stringsAsFactors=F )

gene2symbol = expRMA_gene[,1:2]
gene2symbol[,1]=as.character(gene2symbol[,1])

expRMA_gene = expRMA_gene[,3:ncol(expRMA_gene)]
rownames(expRMA_gene) = gene2symbol[,1]

##### filtering expression and ppi, only keep the overlapped genes

ppi_sig = ppi_sig[ ppi_sig[,1]%in%rownames(expRMA_gene) & ppi_sig[,2]%in%rownames(expRMA_gene) ,]
keepGene = union( ppi_sig[,1],ppi_sig[,2]  ) 

ppi_all = ppi_all[ ppi_all[,1]%in%keepGene & ppi_all[,2]%in%keepGene ,]


expRMA_gene = expRMA_gene[ rownames(expRMA_gene) %in% keepGene ,]

gene_exp =  rownames(expRMA_gene)

##################
 
expData = expRMA_gene
ppi = ppi_all

ppiData = ppi[ (ppi[,1]%in%gene_exp)+(ppi[,2]%in%gene_exp) == 2 , ]
ppiData_g = graph_from_data_frame(as.data.frame(ppiData[,1:2]), directed = F, vertices = gene_exp)
edge_attr(ppiData_g,'weight') <- abs(ppiData[,3])
edge_attr(ppiData_g,'direction') <- ppiData[,5]

#### get scores

score = read.table(parameter_file,header=TRUE, stringsAsFactors=F )
score[,1] = paste("X",score[,1],sep="")
score = score[  match( colnames(expData), score[,1]  ) ,  ]

x10 = as.matrix(t(expData))
scoreEXP = score[,3]
scoreICM = score[,4]
scoreTE = score[,5]

del_sample="X0613_Suek15_Ep_H_SMH_E6005_10.CEL"
x9 = as.matrix(t( expData[, colnames(expData)!=  del_sample ]   ))
scoreT5 = score[  score[,1]!=del_sample  , colnames(score)=="t5"  ]
scoreCC2 = score[ score[,1]!=del_sample , colnames(score)=="cc2"  ]
scoreS2 = score[ score[,1]!=del_sample  , colnames(score)=="s2"  ]
scoreCC3 = score[ score[,1]!=del_sample , colnames(score)=="cc3"  ]

############
result_EXP = GBLasso_batch_correlationFiltering(scoreEXP,x10,ppiData_g)
result_ICM = GBLasso_batch_correlationFiltering(scoreICM,x10,ppiData_g)
result_TE = GBLasso_batch_correlationFiltering(scoreTE,x10,ppiData_g)
result_T5 = GBLasso_batch_correlationFiltering(scoreT5,x9,ppiData_g)
result_CC2 = GBLasso_batch_correlationFiltering(scoreCC2,x9,ppiData_g)
result_S2 = GBLasso_batch_correlationFiltering(scoreS2,x9,ppiData_g)
result_CC3 = GBLasso_batch_correlationFiltering(scoreCC3,x9,ppiData_g)


save.image("result_GBLassoEBinf_call_PM1_geneOvlp_maxloop10000.RData")







source("GBLasso_correlationfiltering_result_statistic.R")
resultStatistic = GBLasso_correlationFiltering_result_statistic("result_GBLassoEBinf_call_PM1_geneOvlp_maxloop10000.RData")
write.table(resultStatistic, "result_GBLassoEBinf_call_PM1_geneOvlp_statistic_maxloop10000.txt", quote=F, sep="\t", row.names=FALSE, col.names=TRUE)

###########################

dir.create("featureGenes")

setwd("featureGenes")

for (i in 1:nrow(resultStatistic)){

  temp_score = resultStatistic[i,1]
  temp_param = resultStatistic[i,2]
  temp_numBeta = resultStatistic[i,6]
  
  if( !is.na(temp_numBeta) ){
    
    result_name = paste("result_",temp_score,"$",temp_param ,sep="" )
    temp_result = eval( parse(text = result_name) )
    # temp_result = result_CC2$pcc_W_g0e3k3
    temp_corGene = temp_result$corGenes
    temp_beta = temp_result$GBLassoFit$BetaSelection$betas
    temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
    temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
    featureGene = temp_gene[ temp_gene[,3]!=0,  ]
    outputFile = paste( "featureGene_",temp_score,"_",temp_param,sep="")
    
    write.table( featureGene, outputFile, sep="\t", quote=F, row.names=F, col.names=F )
    
  }
}

setwd("..")

# save.image("result_GBLassoEBinf_call_PM1_geneOvlp.RData")

############################

setwd("featureGenes_selected")

featureGene_CC2 = read.table("featureGene_CC2_spm_W_g0e3k1", header = F, stringsAsFactors = F)
featureGene_CC3 = read.table("featureGene_CC3_spm_W_g0e3k1", header = F,stringsAsFactors = F)
featureGene_EXP = read.table("featureGene_EXP_spm_W_g0e3k1", header = F,stringsAsFactors = F)
featureGene_ICM = read.table("featureGene_ICM_spm_W_g0e3k1", header = F,stringsAsFactors = F)
featureGene_S2 = read.table("featureGene_S2_spm_W_g0e3k1", header = F,stringsAsFactors = F)
featureGene_T5 = read.table("featureGene_T5_spm_W_g0e3k1", header = F,stringsAsFactors = F)
featureGene_TE = read.table("featureGene_TE_spm_W_g0e3k1", header = F,stringsAsFactors = F)


featureGene_list = list( CC2=featureGene_CC2[,1],
                         CC3=featureGene_CC3[,1],
                         EXP=featureGene_EXP[,1],
                         ICM=featureGene_ICM[,1],
                         S2=featureGene_S2[,1],
                         T5=featureGene_T5[,1],
                         TE=featureGene_TE[,1])

source("/Users/Hui/Dropbox/my/Rscripts/enrichmentGO/FDRoptional/enrichmentGO_human_EntrezGeneID.R")

dir.create("featureGene_enrichment")
enrichmentGO( featureGene_list,  category="BP", minDEG=3, cutoff=0.05, FDRcorrection=FALSE, DelRedundency=FALSE, "featureGene_enrichment" )

dir.create("featureGene_enrichment_childNode")
enrichmentGO( featureGene_list,  category="BP", minDEG=3, cutoff=0.05, FDRcorrection=FALSE, DelRedundency=TRUE, "featureGene_enrichment_childNode" )





source("../getEnrichmentBarplot_ggplot2.R")

dir.create("featureGene_enrichment_figure")

temp_files = list.files( "featureGene_enrichment" )

for(i in 1:length(temp_files)){
  temp_inputFile = paste("featureGene_enrichment/", temp_files[i], sep="" )
  temp_outputFile = paste( "featureGene_enrichment_figure/" ,sub("txt","pdf",temp_files[i] ) ,sep="")
  getEnrichmentBarplot_ggplot2( temp_inputFile, 10, temp_outputFile)
}





#################### plot feature genes in network

library(igraph)

dir.create("featureGene_display_figure")

pdf("featureGene_display_figure/TE.pdf")

temp_feature = featureGene_TE
temp_gene = as.character(temp_feature[,1])
temp_vid = match( temp_gene, gene_exp )
temp_graph = induced.subgraph(ppiData_g, temp_vid   )

temp_v_sign = temp_feature[,3]
temp_v_sign[temp_v_sign<=0]=0
temp_v_sign[temp_v_sign>0]=1

vertex_attr(temp_graph,"symbol") = temp_feature[,2]
vertex_attr(temp_graph,"beta") = abs( temp_feature[,3] )
vertex_attr(temp_graph,"direction") = temp_v_sign

v_size = V(temp_graph)$beta
v_size = (( v_size/max(v_size) )^0.25) * 20

v_color = as.character( V(temp_graph)$direction )
v_color[v_color=="0"] = "steel blue"
v_color[v_color=="1"] = "firebrick3"

e_color = as.character( E(temp_graph)$direction )
e_color[e_color=="0"] = "steel blue"
e_color[e_color=="1"] = "firebrick3"

e_width = E(temp_graph)$weight
e_width = sqrt(e_width) * 5

plot(temp_graph,
     layout = layout.fruchterman.reingold,
     vertex.shape="circle",
     vertex.size=v_size,
     vertex.color=v_color,
     vertex.frame.color=v_color,
     vertex.label=V(temp_graph)$symbol,
     vertex.label.dist=1,
     vertex.label.font=3,
     edge.color=e_color,
     edge.width=e_width   )

dev.off()









setwd("..")
save.image("result_GBLassoEBinf_call_PM1_geneOvlp_maxloop10000_featureGenes.RData")



# temp_result = result_CC3$pcc_W_g0e3k3
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_CC3 = temp_gene[ temp_gene[,3]!=0,  ]
# 
# temp_result = result_EXP$spm_W_g0e3k3
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_EXP = temp_gene[ temp_gene[,3]!=0,  ]
# 
# temp_result = result_ICM$spm_W_g0e3k4
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_ICM = temp_gene[ temp_gene[,3]!=0,  ]
# 
# temp_result = result_S2$pcc_W_g0e3k3
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_S2 = temp_gene[ temp_gene[,3]!=0,  ]
# 
# temp_result = result_T5$pcc_W_g0e3k2
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_T5 = temp_gene[ temp_gene[,3]!=0,  ]
# 
# temp_result = result_TE$spm_W_g0e3k2
# temp_corGene = temp_result$corGenes
# temp_beta = temp_result$GBLassoFit$BetaSelection$betas
# temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
# temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
# featureGene_TE = temp_gene[ temp_gene[,3]!=0,  ]




# ###############
# 
# source("GBLasso_batch_correaltionFiltering_g8.R")
# 
# result_EXP = GBLasso_batch_correlationFiltering_g8(scoreEXP,x10,ppiData_g,ppi)
# result_ICM = GBLasso_batch_correlationFiltering_g8(scoreICM,x10,ppiData_g,ppi)
# result_TE = GBLasso_batch_correlationFiltering_g8(scoreTE,x10,ppiData_g,ppi)
# result_T5 = GBLasso_batch_correlationFiltering_g8(scoreT5,x9,ppiData_g,ppi)
# result_CC2 = GBLasso_batch_correlationFiltering_g8(scoreCC2,x9,ppiData_g,ppi)
# result_S2 = GBLasso_batch_correlationFiltering_g8(scoreS2,x9,ppiData_g,ppi)
# result_CC3 = GBLasso_batch_correlationFiltering_g8(scoreCC3,x9,ppiData_g,ppi)
# 
# 
# save.image("result_GBLassoEB_call_0_geneOvlp_g8.RData")
# 
# source("GBLasso_correlationfiltering_result_statistic.R")
# resultStatistic = GBLasso_correlationFiltering_result_statistic("result_GBLassoEB_call_0_geneOvlp_g8.RData")
# write.table(resultStatistic, "result_GBLassoEB_call_0_geneOvlp_statistic_g8.txt", quote=F, sep="\t", row.names=F, col.names=T)
# 
# 
# 
# 
# ###############
# 
# source("GBLasso_batch_correaltionFiltering_e3k1.R")
# 
# result_EXP = GBLasso_batch_correlationFiltering_e3k1(scoreEXP,x10,ppiData_g,ppi)
# result_ICM = GBLasso_batch_correlationFiltering_e3k1(scoreICM,x10,ppiData_g,ppi)
# result_TE = GBLasso_batch_correlationFiltering_e3k1(scoreTE,x10,ppiData_g,ppi)
# result_T5 = GBLasso_batch_correlationFiltering_e3k1(scoreT5,x9,ppiData_g,ppi)
# result_CC2 = GBLasso_batch_correlationFiltering_e3k1(scoreCC2,x9,ppiData_g,ppi)
# result_S2 = GBLasso_batch_correlationFiltering_e3k1(scoreS2,x9,ppiData_g,ppi)
# result_CC3 = GBLasso_batch_correlationFiltering_e3k1(scoreCC3,x9,ppiData_g,ppi)
# 
# 
# save.image("result_GBLassoEB_call_0_geneOvlp_e3k1.RData")
# 
# source("GBLasso_correlationfiltering_result_statistic.R")
# resultStatistic = GBLasso_correlationFiltering_result_statistic("result_GBLassoEB_call_0_geneOvlp_e3k1.RData")
# write.table(resultStatistic, "result_GBLassoEB_call_0_geneOvlp_statistic_e3k1.txt", quote=F, sep="\t", row.names=F, col.names=T)

