GBLasso_batch_correlationFiltering <- function(Y,X,ppiData_g){

# print("parameter g0e2k2")
# #pcc_D_g0e2k2 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.01,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# pcc_W_g0e2k2 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.01,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# #spm_D_g0e2k2 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.01,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# spm_W_g0e2k2 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.01,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# 
# print("parameter g0e2k3")
# #pcc_D_g0e2k3 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.01,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# pcc_W_g0e2k3 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.01,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# #spm_D_g0e2k3 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.01,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# spm_W_g0e2k3 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.01,ksi=0.001,MAXITER=1000,MAXLOOP=10000)

print("parameter g0e3k1")
#pcc_D_g0e3k1 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.1,MAXITER=1000,MAXLOOP=10000)
pcc_W_g0e3k1 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.1,MAXITER=10000,MAXLOOP=10000)
#spm_D_g0e3k1 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.1,MAXITER=1000,MAXLOOP=10000)
spm_W_g0e3k1 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.1,MAXITER=10000,MAXLOOP=10000)

# print("parameter g0e3k2")
# #pcc_D_g0e3k2 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# pcc_W_g0e3k2 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.01,MAXITER=10000,MAXLOOP=10000)
# #spm_D_g0e3k2 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.01,MAXITER=1000,MAXLOOP=10000)
# spm_W_g0e3k2 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.01,MAXITER=10000,MAXLOOP=10000)

# print("parameter g0e3k3")
# #pcc_D_g0e3k3 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# pcc_W_g0e3k3 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# #spm_D_g0e3k3 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# spm_W_g0e3k3 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.001,MAXITER=1000,MAXLOOP=10000)
# 
# print("parameter g0e3k4")
# #pcc_D_g0e3k4 <- GBLasso_corGenes_Wdegree(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.0001,MAXITER=1000,MAXLOOP=10000)
# pcc_W_g0e3k4 <- GBLasso_corGenes_Wcoexp(Y,X,"pearson",ppiData_g,gamma=0,epsilon=0.001,ksi=0.0001,MAXITER=1000,MAXLOOP=10000)
# #spm_D_g0e3k4 <- GBLasso_corGenes_Wdegree(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.0001,MAXITER=1000,MAXLOOP=10000)
# spm_W_g0e3k4 <- GBLasso_corGenes_Wcoexp(Y,X,"spearman",ppiData_g,gamma=0,epsilon=0.001,ksi=0.0001,MAXITER=1000,MAXLOOP=10000)



result = list( pcc_W_g0e3k1 = pcc_W_g0e3k1,
               spm_W_g0e3k1 = spm_W_g0e3k1
               )

return(result)



}