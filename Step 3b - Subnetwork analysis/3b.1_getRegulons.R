

# GRN = dataframe with columns TF, TG and Value
# minTarget = 5 (default) - the number of targets required to define a regulon
# savePath - path to save the regulon list
getRegulons = function(GRN ,minTarget = 5,savePath){
    GRN = GRN[,c('TF','TG','Value')]
    GRN = GRN[order(GRN$Value,decreasing = T),]
    
    require(dplyr)
    regulons = unique(GRN$TF)
    
    dex = list(unique(GRN$TG[GRN$TF == regulons[1]]))
    for (i in 2:length(regulons)){
      dex[[i]] = unique(GRN$TG[GRN$TF == regulons[i]])
    }
    names(dex) = regulons
    
    # removing regulons that has less than minTarget
    toRemove = c()
    for(i in 1:length(dex)){
      if(length(dex[[i]]) < minTarget){
        toRemove = c(toRemove,i)
      }
    }
    for (i in rev(toRemove)){
      dex[[i]] = NULL
    }
    print('Number of regulons:')
    print(length(dex))
    save(dex,file = paste0(savePath,'.RData'))
    return(dex)
}
  
  
