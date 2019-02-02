

### get feature gene list from the results of call_1 maxloop10000

setwd("/Users/Hui/myFiles/cambridge/projects/epihealthnet/ManchesterSecondment/EmbryoScope/parameter_analysis/gene_expression_association/GBLasso_inf/")
load("result_GBLassoEB_call_1_geneOvlp_maxloop10000.RData")

dir.create("featureGenes_call_1_maxloop10000")
setwd("featureGenes_call_1_maxloop10000")

gene2symbol = expRMA_gene[,1:2]
gene2symbol[,1] = as.character(gene2symbol[,1])

temp_result = result_CC2$pcc_W_g0e3k3
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_CC2 = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_CC3$pcc_W_g0e3k3
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_CC3 = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_EXP$spm_W_g0e3k3
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_EXP = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_ICM$spm_W_g0e3k4
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_ICM = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_S2$pcc_W_g0e3k3
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_S2 = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_T5$pcc_W_g0e3k2
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_T5 = temp_gene[ temp_gene[,3]!=0,  ]

temp_result = result_TE$spm_W_g0e3k2
temp_corGene = temp_result$corGenes
temp_beta = temp_result$GBLassoFit$BetaSelection$betas
temp_symbol = gene2symbol[match(temp_corGene, gene2symbol[,1]),2]
temp_gene = data.frame( Gene= temp_corGene,Symbol=temp_symbol, beta=temp_beta ,stringsAsFactors=FALSE)
featureGene_TE = temp_gene[ temp_gene[,3]!=0,  ]




# library(org.Hs.eg.db)
# columns(org.Hs.eg.db)
# 
# #uniKeys = keys(org.Hs.eg.db, keytype="ENTREZID")
# uniKeys = as.character(exp_gene[,1])
# cols = c("SYMBOL")
# gene2symbol = select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENTREZID")
# 
# for (i in 1:nrow(gene2symbol)){
#   if( is.na(gene2symbol[i,2]) ){
#     gene2symbol[i,2] = gene2symbol[i,1]
#   }
# }









write.table(featureGene_CC2,"featureGene_CC2.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_CC3,"featureGene_CC3.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_EXP,"featureGene_EXP.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_ICM,"featureGene_ICM.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_S2,"featureGene_S2.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_T5,"featureGene_T5.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.table(featureGene_TE,"featureGene_TE.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)






###############

pivot = read.table("/Users/Hui/myFiles/cambridge/projects/epihealthnet/art_exp_development_pathwayCommon/test/linkcomm_network_sig/linkcomm_cluster/analysis_max200_one2one/hub2mod_pivot.txt",header=TRUE,stringsAsFactors = FALSE)

pivot[,1] = as.character(pivot[,1])







