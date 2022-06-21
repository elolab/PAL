# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector (starting from 0) indicating which samples belong to the same group,
#        "pathwayadress" is path to .xml pathways downloaded from KEGG, 
#        "datatype" is a character like "rnaseq_tpm" describing the type of genomic input data in file "data",
#        "noisedefault" is the default cutoff (numeric value) for real signal, and
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# OUTPUT: Checks that all arguments are valid (gives an error if not) and sets some default values. 

CheckInput_PAL = function(data, info, grouplabels, neutralise, mainfeature, userandom, neutralisationformula, pathwayformula, neutralisationmodel, pathwaymodel){

  # Check that 'info' is a data frame
  if(!("data.frame" %in% class(info))) stop("Argument 'info' should be a data frame.")
  
  # Check that column names in info do not include spaces or notation used in formulas
  includesspace = grep(" ", colnames(info), value=T)
  if(length(includesspace) > 0){
    message = paste("Please remove spaces from column names of 'info'", paste(includesspace, collapse=", "), sep=": ")
    stop(message)
  }
  includessyntax = grep("\\+|\\-|\\*|\\(|\\)|\\:|\\||\\%", colnames(info), value=T)
  if(length(includessyntax) > 0){
    message = paste("Please remove formula syntax from following column names of 'info'", paste(includessyntax, collapse=", "), sep=": ")
    stop(message)
  }
  
  # Check that the rows of 'info' correspond to the columns of 'data'
  if(nrow(info) != ncol(data)) stop("The number of rows in 'info' should equal to the number of columns in 'data'.")
  
  # Check that sample groups are defined
  if(!(grouplabels %in% colnames(info))) stop("Argument 'grouplabels' should be a name of a column in argument 'info'.")
  if(all(info[,grouplabels] != 0)) stop("In 'info', the column 'grouplabels' should include 0's to indicate control samples.")
  
  # Check that variables used in neutralisation step have some variation among control samples
  if(is.null(neutralisationformula)){
    usevariables = c(neutralise, userandom)
    if(length(usevariables) > 0){
      uniquevals = apply(info[info[,grouplabels]==0,usevariables,drop=FALSE], 2, function(v){return(length(unique(v)))})
      if(any(uniquevals < 2)){
        message = paste("Following variables can not be used in neutralisation step as they do not differ within control samples", 
                        paste(usevariables[uniquevals < 2], collapse=", "), sep=": ")
        stop(message)
      }
    }
  }
  
  # Check argument 'neutralise'
  if(!is.null(neutralise)){
    if(!all(neutralise %in% colnames(info))) stop("Argument 'neutralise' should be either NULL or a vector of name(s) of a column(s) in argument 'info' (check spelling).")
    # Are missing values ok?
    if(any(info[,grouplabels] != 0)){ # Check that cases do not include categorical variables to be neutralised that are not represented among controls
      categorical = neutralise[apply(info[,neutralise,drop=FALSE],2,function(x){all(is.na(as.numeric(x)))})]
      if(length(categorical) > 0){
        controlcategories = unique(unlist(info[info[,grouplabels]==0,categorical]))
        othercategories = unique(unlist(info[info[,grouplabels]!=0,categorical]))
        notrepresented = setdiff(unique(othercategories), unique(controlcategories))
        if(length(notrepresented) > 0){
          message = paste("Following categories are not presented in control samples, but should be neutralised for (not allowed)", paste(notrepresented,collapse=", "), sep=": ")
          stop(message)
        }
      }
    }
  }
  
  # Check argument 'userandom'
  if(!is.null(userandom)){
    if(!all(userandom %in% colnames(info))) stop("Argument 'userandom' should be either NULL or a vector of name(s) of a column(s) in argument 'info' (check spelling).")
    # Are missing values ok?
  }
  
  # Check argument 'mainfeature'
  if(!is.null(mainfeature)){
    #if(length(mainfeature) > 1) stop("Only one 'mainfeature' (name of the column in 'info') can be used in this version of PAL.")
    mainpieces = unlist(strsplit(mainfeature, split=":"))
    missingpieces = unique(setdiff(mainpieces, colnames(info)))
    if(length(missingpieces) > 0){
      message = paste("Argument 'mainfeature' contains following variables not present in 'info'", paste(missingpieces,collapse=", "), sep=": ")
      stop(message)
    } 
  }
  
  # Check 'neutralisationformula' (note: neutralise can include only fixed effects)
  if(!is.null(neutralisationformula)){
    if(!any(class(neutralisationformula) %in% c("formula","character"))) stop("Argument 'neutralisationformula' should be either NULL, a character, or a formula.")
    if(is.null(neutralise)) stop("If argument 'neutralisationformula' is not NULL, also argument 'neutralise' should be defined.")
    sides = as.character(as.formula(neutralisationformula))
    if(sides[2] != "Expression") stop("If defined, argument 'neutralisationformula' should have left hand side of Expression~")
    fixed = unlist(strsplit(gsub("\\(.*\\)","",gsub(" ","",sides[3])), split="\\+"))
    if(!all(neutralise %in% fixed)) stop("Variables in argument 'neutralise' should be fixed effects in argument 'neutralisationformula'.")
    #random = unlist(strsplit(gsub("\\).*\\(","",gsub(" ","",sides[3])), split="\\+")) # needs fixing
  }
 
  # Check 'pathwayformula'
  if(!is.null(pathwayformula)){
    if(!any(class(pathwayformula) %in% c("formula","character"))) stop("Argument 'pathwayformula' should be either NULL, a character, or a formula.")
    if(is.null(mainfeature)) stop("If argument 'pathwayformula' is not NULL, also argument 'mainfeature' should be defined.")
    sides = as.character(as.formula(pathwayformula))
    if(sides[2] != "Score") stop("If defined, argument 'pathwayformula' should have left hand side of Score~")
    if(length(grep(mainfeature, sides[3])) < 1) stop("Argument 'mainfeature' should be present in 'pathwayformula', if both are given.")
  }
  
  # Check 'neutralisationmodel'
  if(!(neutralisationmodel %in% c("rlm","lmer","rlmer"))) stop("Argument 'neutralisationmodel' should be one of the following characters: 'rlm', 'lmer', 'rlmer'.")
  
  # Check 'pathwaymodel'
  if(!(pathwaymodel %in% c("rlm","lmer","rlmer"))) stop("Argument 'pathwaymodel' should be one of the following characters: 'rlm', 'lmer', 'rlmer'.")

  
  return(mainfeature)
}