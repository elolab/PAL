# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector (starting from 0) indicating which samples belong to the same group,
#        "pathwayadress" is path to .xml pathways downloaded from KEGG, 
#        "datatype" is a character like "rnaseq_tpm" describing the type of genomic input data in file "data",
#        "noisedefault" is the default cutoff (numeric value) for real signal, and
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# OUTPUT: Checks that all arguments are valid (gives an error if not) and sets some default values. 

CheckInput_PAL = function(data, grouplabels, info, neutralise, mainfeature){

  # Check that 'info' is a data frame
  if(!("data.frame" %in% class(info))) stop("Argument 'info' should be a data frame.")
  
  # Check that the rows of 'info' correspond to the columns of 'data'
  if(nrow(info) != ncol(data)) stop("The number of rows in 'info' should equal to the number of columns in 'data'.")
  
  # Check that 'neutralize' and 'mainfeature' are either NA or columns of 'info'
  message = "Arguments 'neutralize' and 'mainfeature' should be either NA or names of columns in argument 'info'."
  if(!all(c(neutralise,mainfeature) %in% c(NA, colnames(info)))) stop(message)
  if((NA %in% neutralise) & any(colnames(info) %in% neutralise)) stop("Argument 'neutralize' can be either NA or name(s) of column(s) in 'info', but not mix both.")
  
  # Check that there is maximum of one main feature
  if(length(mainfeature) > 1) stop("Only one 'mainfeature' can be used.")
  
  # Check that sample groups are defined
  if(!(grouplabels %in% colnames(info))) stop("Argument 'grouplabels' should be a name of a column in argument 'info'.")
  
  # Check that no other column than possibly 'mainfeature' in 'info' has missing values
  keepcols = setdiff(colnames(info), mainfeature)
  if(any(is.na(info[,keepcols,drop=FALSE]))) stop("data frame 'info' can contain missing values (NA) only in the column identified by 'mainfeature'.")
 
  
  # Check that 'mainfeature' is either NA or a vector (done already in PAL)
  #if(length(mainfeature) == length(grouplabels)){
  #  if(length(unique(grouplabels)) > 1) mainfeature[grouplabels == 0] = NA
  #} 
  
  return(mainfeature)
}