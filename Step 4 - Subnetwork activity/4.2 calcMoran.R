
# This function calculates Moran's I scores for given the adata file containing pseudotime information in adata$obs. 
# if the dataset contains more than 10000 cells, we will divide the dataset into multiple groups and calculate the Moran's I. 
# the cleanMoran function will combine the scores and adjust for multiple testing.

cleanMoran = function(savePath,saveKey){
  
  files = list.files(paste0(savePath))
  files = files[grep(pattern = paste0('MoransI_',saveKey,'_subsample_'), files)]
  print(files)
  
  
  meta = c()
  for (file in files){
    temp = read.csv(paste0(savePath,file))
    temp$padj = p.adjust(temp$p, method ='BH')
    meta = rbind(meta, temp)
    
  }
  print(dim(meta))
  
  library(dplyr)
  
  library(dplyr)
  metaF <- meta %>%
    filter(padj < 0.05) %>%
    group_by(X) %>%
    summarise(
      mean_observed = mean(observed),
      mean_sd = mean(sd),
      combined_padj = p.adjust(min(padj), method = "BH"),  # Conservative: min padj per group, adjusted across groups
      .groups = 'drop'
    )
  
  metaF = as.data.frame(metaF)
  colnames(metaF) = c('Subnetwork','DynamicScore','SD','padj')
  write.csv(metaF,paste0(savePath,'MoransI_',saveKey,'_cleaned.csv'))
}


calcMoran  = function(adata, AUCs, pseudotimeSpaceColumns = c('pr_pseudotime', 'subsample_0', 'subsample_1', 'subsample_2', 'subsample_3'), savePath, saveKey){
  require(anndata)
  require(ape)
  require(doParallel)
  
  obs = adata$obs
  obs = obs[,pseudotimeSpaceColumns]
  
  AUCOrig = AUCs
  AUCOrig = merge(AUCs, obs, by.x = 'row.names',by.y = 'row.names', all.x = T)
  
  print(paste0('Number of cells',dim(AUCOrig)[1]))
  if(dim(AUCOrig)[1] < 10000){
    n = 1
  }else{
    n = ceiling(dim(AUCOrig)[1]/10000)
    print(paste0('AUCOrig file is too long to handle, changing n value to ',n))
  }
  
  sampleOrder  = sample(dim(AUCOrig)[1],dim(AUCOrig)[1],replace = F)
  numSamples =c(seq(from = 1, to = dim(AUCOrig)[1], by = ceiling(dim(AUCOrig)[1]/n)),dim(AUCOrig)[1])
  
  for (iter in 1:n){
    
    AUCs = AUCOrig[sampleOrder[numSamples[iter]:numSamples[iter+1]],]
    print(dim(AUCs))
    # calculating Moran's I scores for combined pseudotimes
    dists <- as.matrix(dist(cbind(AUCs[,pseudotimeSpaceColumns])))
    
    dists.inv <- 1/dists
    diag(dists.inv) <- 0
    dists.inv[is.infinite(dists.inv)] <- 0
    
    gc()
    print('obtaining dists.inv complete')
    
    runMoran = function(AUCs, i, dists.inv, alternative = 'greater',scaled = T, na.rm = T){
      # require(spdep)
      require(ape)
      gc()
      res = Moran.I(AUCs[,i],dists.inv,alternative = alternative,scaled = scaled,na.rm = na.rm)
      gc()
      return(c(res$observed, res$sd, res$p.value))
      
    }
    library(doParallel)
    
    N = dim(AUCs)[2]-length(pseudotimeSpaceColumns)
    print('Parallel run initiating')
    cl = makeCluster(30,type = 'PSOCK')
    registerDoParallel(cl)
    
    
    outcome = foreach (i = c(1:N), .combine = rbind) %dopar% { 
      runMoran(AUCs, i, dists.inv = dists.inv)
    }
    
    stopCluster(cl)
    print('parallel run ending')
    outcome = as.data.frame(outcome)
    row.names(outcome) = colnames(AUCs)[1:N]
    colnames(outcome) = c('observed','sd','p.value')
    
    write.csv(outcome,paste0(savePath,'/MoransI_',saveKey,'_subsample_',iter,'.csv'))
  }
  
  cleanMoran(savePath = savePath , saveKey = saveKey)
}

