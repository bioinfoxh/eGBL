
setwd("/Users/Hui/myFiles/cambridge/projects/GBLasso/setup2-110feature/")

source("simulation_data_550gene110feature.R")

dir.create("simulateData")
setwd("simulateData")

nDataset = 120
for (i in 1:nDataset){
  simData = simulateData(i)
  filename = paste("simData_",i,".RData",sep="")
  save( simData, file=filename )
}

setwd("..")


