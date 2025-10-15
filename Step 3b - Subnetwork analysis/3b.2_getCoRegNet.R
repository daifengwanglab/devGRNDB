


# this file inputs a GRN and save the co-regulation module on user specified location.
# net = GRN dataframe with columns, TF, TG, ctx_score
# savePath = path for the folder where the network needs to be saved
# saveKey = file name

getCoRegNet = function(net,savePath, saveKey){
  
  library(WGCNA)
  library(parallel)
  
  jaccard_similarity <- function(a, b) {
    length(intersect(a, b)) / length(union(a, b))
  }
  
  allgenes = unique(c(net$TF,net$TG))
  
  mat = matrix(0,nrow = length(allgenes), ncol = length(allgenes))
  colnames(mat) = allgenes
  row.names(mat) = allgenes
  

  # Create list of TFs for each TG
  gene_to_tfs <- setNames(lapply(allgenes, function(g) net$TF[net$TG == g]), allgenes)
  
  
  mat_list <- mclapply(seq_along(allgenes), function(i) {
    gene_i <- allgenes[i]
    tf_set_i <- gene_to_tfs[[gene_i]]
    sapply(gene_to_tfs, function(tf_set_j) jaccard_similarity(tf_set_i, tf_set_j))
  }, mc.cores = 20)
  
  mat <- do.call(rbind, mat_list)
  
  rownames(mat) <- colnames(mat) <- allgenes
  
  save(mat,allgenes, mat_list,file = paste0(savePath,"/",saveKey,'_coregulationNet.RData'))
  return(mat_list)
}