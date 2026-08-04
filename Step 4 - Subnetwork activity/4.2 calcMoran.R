

## Moran's I calculations --------------------------------------------------

# The following functions calculates the Moran's I score based dynamic scores. 

#' Internal function of calculating Moran's I using "spdep" R package
#'
#' @param i iterating index
#' @param AUC a dataframe with AUC scores (cells by AUC). each column represents a regulon/gene module
#' @param case_names_df A data.frame containing metadata matching each file
#' 
.moran_worker <- function(i, AUC, lw, nperm, returnPerm = FALSE) {
  
  message("Working on the gene ", colnames(AUC)[i])
  
  x <- AUC[, i]
  
  res <- spdep::moran.mc(
    x,
    lw,
    nsim = nperm,
    alternative = "greater"
  )
  
  I_obs  <- res$statistic
  I_perm <- res$res[1:nperm]
  
  # mt_rand <- moran.test(x, lw, randomisation = TRUE) 
  
  summarydf = data.frame(
    gene = colnames(AUC)[i],
    I_obs = as.numeric(I_obs),
    meanPerm = mean(I_perm),
    sdPerm = sd(I_perm),
    Z_I = (I_obs - mean(I_perm)) / sd(I_perm),
    NES = (I_obs - mean(I_perm)) /
      max(abs(I_perm - mean(I_perm))),
    Pval = res$p.value,
    maxAUC = max(x, na.rm = T),
    meanAUC = mean(x),
    sdAUC = sd(x, na.rm = T),
    maxAUCZ = (max(x,na.rm = T)-mean(x))/sd(x, na.rm = T),
    stringsAsFactors = FALSE
  )
  
  
  if (returnPerm) {
    return(list(summary = summarydf, perm = I_perm))
  } else {
    return(list(summary = summarydf))
  }
  
  
}


#' Executable function to estimate Moran's I 
#' @param AUCInput AUCell score input. It could be provided as a dataframe
#' (cell by regulon/cell by gene module) or the path to a csv file with 1st column with cell IDs.
#' @param metaInput cell metadata. if NULL, the meta data are included as the last columns of the AUCInput.
#' otherwise it could be provided as a dataframe or a csv file with 1st column with cell IDs. 
#' @param savePath folder path to store results, folder will be created if not available.
#' @param saveKey a string variable to distinguish saved result files
#' @param k number of nearest neighbors. if k = NULL, value of 100 will be assigned by default.
#' @param numOfMetaColumns number of metadata columns included in AUCInput. 
#' @pseudotimeColumns column names of pseudotime values in metadata
#' @numProcessorsTobeRemoved number of processors removed from the parallel computing. 
#' @param padjustMethod name of the p.value adjustment method
#' @param nperm number of permutations to be performed for p-value calculations
#' @param savePerm whether to save the permuted Moran's I values
#' @param permFile the folder where the permutations are saved. 
#' 
getMoran <- function(AUCInput,
                     metaInput = NULL,
                     savePath,
                     saveKey,
                     k = NULL,
                     numOfMetaColumns = 8,
                     pseudotimeColumns = c(
                       "pr_pseudotime",
                       "subsample_0",
                       "subsample_1",
                       "subsample_2",
                       "subsample_3"
                     ),
                     numProcessorsTobeRemoved = 6,
                     padjustMethod = "BH",
                     nperm = 999,
                     savePerm = FALSE,
                     permFile = NULL
                     
) {
  
  library(FNN)          # Fast nearest neighbors
  library(spdep)        # Spatial dependence
  library(future.apply) 
  
  if (!dir.exists(savePath)) dir.create(savePath, recursive = TRUE)
  
  # helpers
  
  
  .is_file <- function(x) {
    is.character(x) && length(x) == 1L && file.exists(x)
  }
  
  
  # reading AUCInput
  if (is.data.frame(AUCInput)) {
    AUC = AUCInput
    message('Input is provided as a data.frame')
  } else if (.is_file(AUCInput)) {
    message('Input is provided as a path to the csv file')
    AUC = read.csv(AUCInput,row.names = 1)
  } else{
    stop('Please provide a valid path or a data.frame to AUCInput')
  }
  
  # Reading metadata input
  if(!is.null(metaInput)){
    if(is.data.frame(metaInput)){
      meta = metaInput
      message('Metadata are provided as a data.frame')
    }else if(.is_file(metaInput)){
      meta = read.csv(metaInput, row.names = 1)
      message('meta data is provided as a csv file')
    }else{
      stop('Please provide a valid path or a data.frame to metaInput')
    }
  }else{
    message('meta data are included with the AUCInput.')
    meta = AUC[, tail(names(AUC),numOfMetaColumns)]
    AUC = AUC[,1:(ncol(AUC)-numOfMetaColumns)]
  }
  
  # Ensure metadata cellIDs matches with the AUCInput cell IDs.
  if (!setequal(rownames(AUC), rownames(meta))) {
    stop("Error: Row names in AUC and meta do not match.")
  }
  # Reorder meta to match AUC
  meta <- meta[rownames(AUC), , drop = FALSE]
  # Optional sanity check to ensure the cells are ordered correctly.
  if (!identical(rownames(AUC), rownames(meta))) {
    stop("Error: Failed to align row names.")
  }
  
  if(!colnames(meta) %in% pseudotimeColumns){
    stop('the provided pseudotime columns cannot be found in meta data')
  }
  embed <- meta[, pseudotimeColumns]
  
  embed <- as.matrix(embed)
  if(is.null(k)){
    k = 100
    message('k is not provided, using a default value of ', k)
  }
  knn_res <- get.knn(embed, k = k)
  
  nb_list <- lapply(seq_len(nrow(knn_res$nn.index)), function(i) {
    knn_res$nn.index[i, ]
  })
  
  nb <- structure(
    nb_list,
    class = "nb",
    region.id = rownames(embed),
    call = match.call(),
    type = "knn"
  )
  
  
  
  glist <- lapply(seq_len(nrow(knn_res$nn.dist)), function(i) {
    d <- knn_res$nn.dist[i, ]
    w <- 1/(d^2 + 1e-9)
    # w <- 1/(d + 1e-9)
    w / sum(w)
  })
  
  

  lw <- nb2listw(
    nb,
    glist = glist,
    style = "W"
  )
  lwU <- listw2U(lw)  # recommended for kNN symmetry assumptions
  
  #   #Constants
  n <- nrow(embed)

  #Parallel setup
  plan(
    multisession,
    workers = max(1, parallel::detectCores() - numProcessorsTobeRemoved)
  )
  set.seed(123)
  
  gene_idx <- seq_len(ncol(AUC))
  
  #Parallel Moran's I
  results_list <- future_lapply(
    gene_idx,
    .moran_worker,
    AUC = AUC,
    lw = lwU,
    nperm = nperm,
    returnPerm = savePerm,
    future.seed = T
  )
  
  results <- do.call(rbind, lapply(results_list, `[[`, "summary"))
  
  # Multiple testing correction
  results$padj <- p.adjust(results$Pval, method = padjustMethod)
  
  # Save summary output
  summary_file <- file.path(savePath, paste0("MoransI_", saveKey, ".csv"))
  write.csv(results, summary_file, row.names = FALSE)
  
  
  # Save permutation outputs if requested
  if (savePerm) {
    if (is.null(permFile)) {
      permFile <- file.path(savePath, paste0("Perm_", saveKey, ".csv"))
    }
    
    # Build permutation matrix (genes x nperm)
    perm_mat <- do.call(rbind, lapply(results_list, function(x) x$perm))
    genes <- results$gene
    
    colnames(perm_mat) <- paste0("perm_", seq_len(nperm))
    perm_df <- data.frame(gene = genes, perm_mat, check.names = FALSE)
    write.csv(perm_df, permFile, row.names = FALSE)
    
    
    message("Permutation values saved to: ", permFile)
  }
  
  
  invisible(results)
}



# This function will preprocess the raw results obtained in getMoran function

#' @param filePath folder path to the results section
#' @param fileName name of the file without the prefix "moranI_" and suffix ".csv"
#' @param getColumns specific columns need to be kept from the original results file
#' @param suffix description
#' @param filterRegulonCutoff optional cutoff value
#' @param filterRegulonsby optional cutoff variable
#' @param addKeyColumn adding a new key column to the results section
#' @param pvalCutoff cutoff for the adjusted p-value

readAndProcessMoran <- function(
    filePath,
    fileName,
    getColumns = NULL,
    suffix = NULL,
    filterRegulonsCutoff = NULL,
    filterRegulonsby = "Percentile",
    addKeyColumn = FALSE,
    pvalCutoff = 0.05,

    
){
  
  inputfile <- file.path(filePath, paste0("moranI_", fileName, ".csv"))
  
  if (!is.character(inputfile) || !file.exists(inputfile)) {
    stop("Input file not available: ", inputfile)
  }
  
  m <- read.csv(inputfile, check.names = FALSE)
  
  valid <- c(
    "Regulon",
    "Moran's I",
    "Mean - permutations",
    "S.D - permutations",
    "Z-score",
    "NES",
    "P-value",
    "Adjusted P-value",
    "Normalized Z score",
    "Percentile"
  )
  
  missing_cols <- setdiff(valid, colnames(m))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.null(getColumns)) {
    bad <- setdiff(getColumns, valid)
    if (length(bad) > 0) {
      stop("Invalid getColumns: ", paste(bad, collapse = ", "))
    }
  }
  
  # keep required columns
  m <- m[, valid, drop = FALSE]
  
  # remove NA only in relevant columns
  m <- m[complete.cases(m), , drop = FALSE]
  
  # significance filter
  m <- m[m$`Adjusted P-value` < pvalCutoff, , drop = FALSE]
  
  # optional filtering
  if (!is.null(filterRegulonsCutoff)) {
    
    if (!filterRegulonsby %in% colnames(m)) {
      stop("filterRegulonsby must be one of: ", paste(colnames(m), collapse = ", "))
    }
    
    if (!is.numeric(m[[filterRegulonsby]])) {
      stop(filterRegulonsby, " must be numeric")
    }
    
    m <- m[m[[filterRegulonsby]] >= filterRegulonsCutoff, , drop = FALSE]
  }
  
  # sort after filtering
  m <- m[order(m$`Moran's I`, decreasing = TRUE), , drop = FALSE]
  
  # optional column selection
  if (!is.null(getColumns)) {
    
    m <- m[, c("Regulon", getColumns), drop = FALSE]
    colnames(m)[1] <- "gene"
    
    if (!is.null(suffix) & !isTRUE(addKeyColumn)) {
      colnames(m)[-1] <- paste0(colnames(m)[-1], "_", suffix)
    }
  }
  
  # optional key column
  if (isTRUE(addKeyColumn)) {
    if (is.null(suffix)) {
      stop("suffix is required when addKeyColumn = TRUE")
    }
    m$key <- suffix
  }
  
  return(m)
}

