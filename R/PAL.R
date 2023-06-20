PAL = function(data, info, grouplabels, adjust=NULL, mainfeature=NULL, userandom=NULL, adjustmentformula=NULL, pathwayformula=NULL, 
               adjustmentmodel="lmer", pathwaymodel="lmer",  pathwayadress=NULL, useKEGG=TRUE, score="activity", nodemin=5, seed=1234){
  
  oldwd = getwd()
  
  # Check and initialise PAL specific input parameters
  mainfeature = CheckInput_PAL(data, info, grouplabels, adjust, mainfeature, userandom, adjustmentformula, pathwayformula, adjustmentmodel, pathwaymodel)
  
  # Extract arguments from info
  grouplabelindex = match(grouplabels, colnames(info), 0)
  grouplabels = info[,grouplabelindex]
  if(ncol(info) < 2){
    info = NULL
  } else info = info[,-grouplabelindex,drop=FALSE]
  
  # Check and initialize input parameters related to general pathway analysis
  originalgenedata = PASI::CheckInput(data, grouplabels, pathwayadress, useKEGG, score, nodemin)
  
  # Read in and preprocess pathways from local file (if available) and KEGG API
  cat("\n") + cat("Accessing and preprocessing pathways, this may take a while.") + cat("\n")
  pathways = PASI:::PreProcessPathways(pathwayadress, useKEGG, data, nodemin, score)
  
  # Extract genes appearing in at least one pathway
  pathwaygenes = unlist(lapply(pathways, function(x){
    entrez = unlist(strsplit(x$nodeinfo$Entrez, split="_")) 
    return(entrez)
  }))
  
  # Adjust for the effect of wanted coefficients in 'info' from pathway genes
  if(!(is.null(neutralise))){
    cat("\n") + cat("Adjusting for coefficients (takes time)...") + cat("\n")
    originalgenedata = NeutraliseCoefficients(originalgenedata, info, adjust, userandom, grouplabels, pathwaygenes, adjustmentformula, adjustmentmodel)
  } 
  
  # Normalize the measurements
  cat("\n") + cat("Scaling data...") + cat("\n")
  scaleddata = PASI::ScaleData(originalgenedata, grouplabels, score)
  
  # Drop genes that don't appear in any pathway
  genedata = scaleddata[rownames(scaleddata) %in% pathwaygenes, ]
  
  cat("\n") + cat("Starting to process the scaled measurements...") + cat("\n")
  
  # Calculate initial node values (NA for missing values)
  nodevalues = PASI::MeasurementsToNodes(pathways, genedata)
  
  # Drop pathways with less than nodemin measured nodes for any sample
  droppathways = unlist(lapply(nodevalues, function(x) all(colSums(!is.na(x)) < nodemin)))
  dropindex = which(droppathways)
  if(length(dropindex)!=0){
    nodevalues = nodevalues[-dropindex]
    cat("Following pathways are left out from analysis due to too few measured nodes:")+cat("\n")
    cat(names(pathways)[dropindex], sep="\n")
    pathways = pathways[-dropindex]
  } 
  
  # Detect structural info from pathways
  pathwaystatistics = PASI::ExtractPathwayStatistics(pathways)
  topologyvalues = PASI::CalculateTopologyFactors(pathways, pathwaystatistics)
  
  # Process the node values with feedback and calculate pathway scores
  if(score=="deregulation"){
    nodevalues = PASI::ConsidereFeedBack(nodevalues, pathways, pathwaystatistics$occurrences)
    results = PASI::CalculatePathwayValues_Deregulation(nodevalues, topologyvalues)
  }
  
  # Calculate relation values and pathway scores
  if(score=="activity"){
    relationvalues = PASI::MeasurementsToRelations(pathways, genedata)
    results = PASI::CalculatePathwayValues_Activity(nodevalues, relationvalues, topologyvalues)
  }
  
  # Add row and column names to pathway results
  rownames(results) = paste(names(pathways), unlist(lapply(pathways, function(p) p$pathwayname)), sep=": ")
  colnames(results) = colnames(genedata)
  
  cat("\n") + cat("Pathway scores calculated, extracting pathway information for reporting") + cat("\n")
  
  # Extract pathway info
  if(score=="activity"){
    pathwayinfo = ExtractPathwayInfo(pathways, results, nodevalues, relationvalues)
  } else pathwayinfo = ExtractPathwayInfo(pathways, results, nodevalues)
  
  # Write the pathway values into a text file
  outputfile = paste("PAL_", ".txt", sep=as.character(Sys.Date()))
  setwd(oldwd)
  write.table(results, file=outputfile, quote=FALSE, sep="\t")
  toreturn = list(pathwayscores=results, pathwayinfo=pathwayinfo)
  
  # Calculate significance level for the feature of interest
  if(!is.null(mainfeature)){
    cat("\n") + cat("Calculating significance levels (takes time)...") + cat("\n")
    significance = switch(pathwaymodel, 
                          rlm = PickImportant_rlm(results, mainfeature, userandom, info, pathwayformula, seed, n=1000),
                          rlmer = PickImportant_rlmer(results, mainfeature, userandom, info, pathwayformula),
                          lmer = PickImportant_lmer(results, mainfeature, userandom, info, pathwayformula, seed, n=1000))
    toreturn = list(pathwayscores=results, significance=significance, pathwayinfo=pathwayinfo)
  }
  
  cat("\n") + cat("Done!") + cat("\n")
  return(toreturn)
}
