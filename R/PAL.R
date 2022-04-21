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

PAL = function(data, info, grouplabels, pathwayadress=NULL, useKEGG=TRUE, score="activity", nodemin=5, 
               neutralize=NA, mainfeature=NA, seed=1234){
  
  oldwd = getwd()
  
  # Check that data frame 'info' is properly given
  if(!("data.frame" %in% class(info))) stop("Argument 'info' should be a data frame.")
  message = "Arguments 'grouplabels', 'neutralize' and 'mainfeature' should be either NA or names of columns in argument 'info'."
  if(!all(c(grouplabels,neutralize,mainfeature) %in% c(NA, colnames(info)))) stop(message)
  if(length(mainfeature) > 1) stop("Only one 'mainfeatrure' can be used.")
  
  # Extract arguments from info
  grouplabelindex = match(grouplabels, colnames(info), 0)
  grouplabels = info[,grouplabelindex]
  if(ncol(info) < 2){
    info = NA
  } else info = info[,-grouplabelindex,drop=FALSE]
  if(!is.na(mainfeature)){
    mainfeatureindex = match(mainfeature, colnames(info), 0)
    mainfeature = info[,mainfeatureindex]
    if(ncol(info) < 2){
      info = NA
    } else info = info[,-mainfeatureindex,drop=FALSE]
  }
  if(!is.na(neutralize)){
    neutralize = c(colnames(info) %in% neutralize)
  }
  
  # Ensure that functions from packages MASS and PASI are ready to be used
  if(!("MASS" %in% (.packages()))) library(MASS)
  if(!("PASI" %in% (.packages()))) library(PASI)
  
  # Check and initialize input parameters
  originalgenedata = CheckInput(data, grouplabels, pathwayadress, useKEGG, score, nodemin)
  mainfeature = CheckInput_PAL(data, grouplabels, info, neutralize, mainfeature)
  
  # Read in and preprocess pathways from local file (if available) and KEGG API
  cat("\n") + cat("Accessing and preprocessing pathways, this may take a while.") + cat("\n")
  pathways = PreProcessPathways(pathwayadress, useKEGG, data, nodemin, score)
  
  # Extract genes appearing in at least one pathway
  pathwaygenes = unlist(lapply(pathways, function(x){
    entrez = unlist(strsplit(x$nodeinfo$Entrez, split="_")) 
    return(entrez)
  }))
  
  # Neutralize effect of coefficient in 'info' from pathway genes
  if(TRUE %in% neutralize){
    cat("\n") + cat("Neutralizing coefficients...") + cat("\n")
    originalgenedata = NeutralizeCoefficients(originalgenedata, info, neutralize, grouplabels, pathwaygenes, mainfeature)
  } 
  
  # Normalize the measurements
  cat("\n") + cat("Scaling data...") + cat("\n")
  scaleddata = ScaleData(originalgenedata, grouplabels, score)
  
  # Drop genes that don't appear in any pathway
  genedata = scaleddata[rownames(scaleddata) %in% pathwaygenes, ]
  
  cat("\n") + cat("Starting to process the scaled measurements...") + cat("\n")
  
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
  if(!all(is.na(mainfeature))){
    
	cat("\n") + cat("Calculating significance levels...") + cat("\n")
	
    if(length(unique(mainfeature[!is.na(mainfeature)])) == 1){
      mainfeature[is.na(mainfeature)] = "Control"
    } 
    
    if(!is.na(neutralize[1])){
      
      index = which(neutralize == FALSE)
      
      # For categorical main variable, only variables that overlap between the categories can be used
      if(all(is.na(as.numeric(mainfeature))) & (length(index)>0)){
        keepvariables = apply(info[,index,drop=FALSE], 2, function(f){
          groups = split(f, mainfeature)
          keep = TRUE
          for(i in 1:length(groups)){
            if(length(intersect(groups[[i]], unlist(groups[-i]))) == 0){
              keep = FALSE
              break
            }
          }
          return(keep)
        })
        index = index[keepvariables]
      }
      
      if(length(index) > 0){
        significance = PickImportant(results, mainfeature, seed, info[,index,drop=FALSE])
      } else significance = PickImportant(results, mainfeature, seed, NA)
    } else{
      significance = PickImportant(results, mainfeature, seed, NA)
    }
    toreturn = list(pathwayscores=results, significance=significance, pathwayinfo=pathwayinfo)
  }
  
  cat("\n") + cat("Done!") + cat("\n")
  return(toreturn)
}