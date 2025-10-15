



# this function will read corregulatory network module saved as a R.Data file and perform WGCNA based hierarchical clustering 
# to obtain gene modules. The module membership file will be returned as the output. 

getGeneModule = function(readPath, readKey, writePath, writeKey){
  library(WGCNA)
  library(parallel)  
  load(paste0(readPath,"/",readKey,'_coregulationNet.RData'))
  outfile = paste0(writePath,'/', writeKey,'_coregulationModule_membership.csv')
  
  allgenes = colnames(mat)
  TOM <- TOMsimilarity(as.matrix(mat))
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average") # modify as see fit
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 100) # modify as see fit

  df = data.frame(Gene = row.names(mat),'module' = dynamicMods)
  write.csv(df,paste0(outfile))
}
