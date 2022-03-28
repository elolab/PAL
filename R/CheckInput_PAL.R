# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector (starting from 0) indicating which samples belong to the same group,
#        "pathwayadress" is path to .xml pathways downloaded from KEGG, 
#        "datatype" is a character like "rnaseq_tpm" describing the type of genomic input data in file "data",
#        "noisedefault" is the default cutoff (numeric value) for real signal, and
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# OUTPUT: Checks that all arguments are valid (gives an error if not) and sets some default values. 

CheckInput_PAL = function(data, grouplabels, info, neutralise, mainfeature){
  
  # Check that 'mainfeature' is either NA or a vector
  if(length(mainfeature) == length(grouplabels)){
    if(length(unique(grouplabels)) > 1) mainfeature[grouplabels == 0] = NA
  } else{
    if(length(mainfeature) == 1){
      if(!is.na(mainfeature)) stop("Argument 'mainfeature' should be either NA or a vector of length ncol(data).")
    } else stop("Argument 'mainfeature' should be either NA or a vector of length ncol(data).")
  }
  
  if(is.data.frame(info)){
    if(nrow(info) != ncol(data)){
      stop("The number of rows in 'info' should equal to the number of columns in 'data'.")
    }
    
    # Check that info is ok
    if(any(is.na(info))){
      stop("'info' should not contain missing values (NA) if it is a data frame.")
    }
  } else{
    if(!is.na(info)) stop("'info' should be either NA or a data frame.")
  }
  
  # Check that neutralise is NA or a logical vector and its length equals to the number of columns in 'info'
  if(!is.na(neutralise[1])){
    if(!is.logical(neutralise)) stop("Argument 'neutralize' should be either NA or a logical vector.")
    if(length(neutralise) != ncol(info)) stop("Length of argument 'neutralize' (if not NA) should equal to the number of columns in 'info'.")
  }
  
  return(mainfeature)
}