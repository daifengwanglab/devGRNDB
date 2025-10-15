

# this function calculates AUCell values for the subnetworks inferred at the previous step.

# inputs 
# adata = anndata (h5ad object with gene expression)
# make sure the adata$var_names have same gene name format as your geneList
# geneList = list of subnetwork genes for AUCell calculations


calcAUC = function(adata,geneList, savePath, saveKey){
  library(anndata)
  library(Matrix)
  # reading expression matrix from anndata object
  exprMat = as.matrix(t(adata$X))
  
  cells_AUC <- AUCell_run(exprMat, geneSets = geneList)
  AUCs = data.frame(t(cells_AUC@assays@data$AUC))
  

  AUCs = merge(AUCs, meta, by.x = 'row.names',by.y = 'row.names')
  row.names(AUCs) = AUCs$Row.names
  AUCs = subset(AUCs, select = -"Row.names")
  
  write.csv(AUCs,paste0(savePath, saveKey,'.csv'))
}