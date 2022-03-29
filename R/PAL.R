# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector indicating samples groups (control: 0),
#        "pathwayadress" is path to pathways,
#        "datatype" is should be either "microarray" or "rnaseq",
#        "noisedefault" is the default cutoff (numeric value) for real signal,
#        "score" defines whether the returned values reflect activity (default) or deregulation,
#        "nodemin" is the smallest number of measured nodes accepted for analysis,
#        "info" includes variables (cols) to be used in neutralisation step (rows are samples),
#        "neutralize" is a logical vector indicating the variables from info to be actually neutralised,
#        "mainfeature" is a vector (corresp. to samples) including the feature whose significance is calculated from the ready pathway scores
# OUTPUT: Returns a list of 1) pathway scores, 2) pathways' significance levels, and 3) pathway info. 

PAL = function(data, grouplabels, pathwayadress=NULL, datatype="rnaseq",
               noisedefault="automatic", score="activity", nodemin=5, 
               info=NA, neutralize=NA, mainfeature=NA){
  
  oldwd = getwd()
  
  # Ensure that functions from packages MASS and PASI are ready to be used
  if(!("MASS" %in% (.packages()))) library(MASS)
  if(!("PASI" %in% (.packages()))) library(PASI)
  
  # Check and initialize input parameters
  parameters = CheckInput(data, grouplabels, pathwayadress, datatype, noisedefault, score, 
                          nodemin, info, mainfeature)
  originalgenedata = parameters[[1]]
  noisedefault = parameters[[2]]
  mainfeature = CheckInput_PAL(data, grouplabels, info, neutralize, mainfeature)
  
  # Detect noise
  if(is.na(noisedefault) | is.numeric(noisedefault)){
    noiselevel = noisedefault
  } else noiselevel = DetectNoise(originalgenedata, noisedefault)
  
  # Read in and preprocess pathways from local file (if available) and KEGG API
  pathways = PreProcessPathways(pathwayadress, data, nodemin, score)
  
  # Extract genes appearing in at least one pathway
  pathwaygenes = unlist(lapply(pathways, function(x){
    entrez = unlist(strsplit(x$nodeinfo$Entrez, split="_")) 
    return(entrez)
  }))
  
  # Neutralize effect of coefficient in 'info' from pathway genes
  if(!is.null(nrow(info))){
    cat("\n") + cat("Neutralizing coefficients...") + cat("\n")
    originalgenedata = NeutralizeCoefficients(originalgenedata, info, neutralize, grouplabels, pathwaygenes, mainfeature)
  } 
  
  # Normalize the measurements
  scaleddata = ScaleData(originalgenedata, grouplabels, noiselevel, score)
  
  # Drop genes that don't appear in any pathway
  genedata = scaleddata[rownames(scaleddata) %in% pathwaygenes, ]
  
  # Calculate initial node values (NA for missing values)
  nodevalues = MeasurementsToNodes(pathways, genedata)
  
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
  pathwaystatistics = ExtractPathwayStatistics(pathways)
  topologyvalues = CalculateTopologyFactors(pathways, pathwaystatistics)
  
  # Process the node values with feedback and calculate pathway scores
  if(score=="deregulation"){
    nodevalues = ConsidereFeedBack(nodevalues, pathways, pathwaystatistics$occurrences)
    results = CalculatePathwayValues_Deregulation(nodevalues, topologyvalues)
  }
  
  # Calculate relation values and pathway scores
  if(score=="activity"){
    relationvalues = MeasurementsToRelations(pathways, genedata)
    results = CalculatePathwayValues_Activity(nodevalues, relationvalues, topologyvalues)
  }
  
  # Add row and column names to pathway results
  rownames(results) = paste(names(pathways), unlist(lapply(pathways, function(p) p$pathwayname)), sep=": ")
  colnames(results) = colnames(genedata)
  
  # Extract pathway info
  if(score=="activity"){
    pathwayinfo = ExtractPathwayInfo(pathways, results, nodevalues, relationvalues)
  } else pathwayinfo = ExtractPathwayInfo(pathways, results, nodevalues)
  
  # Write the pathway values into a text file
  outputfile = paste("PAL_", ".txt", sep=as.character(Sys.Date()))
  setwd(oldwd)
  write.table(results, file=outputfile, quote=F, sep="\t")
  toreturn = list(pathwayscores=results, pathwayinfo=pathwayinfo)
  
  # Calculate significance level for the feature of interest
  if(!all(is.na(mainfeature))){
    if(length(unique(mainfeature[!is.na(mainfeature)])) == 1){
      mainfeature[is.na(mainfeature)] = "Control"
    } 
    significance = PickImportant(results, mainfeature)
    toreturn = list(pathwayscores=results, significance=significance, pathwayinfo=pathwayinfo)
  }
  
  return(toreturn)
}