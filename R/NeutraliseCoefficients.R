# Input: 'data' is a matrix of values (rows: genes/pathways/whatever, cols: samples)
#        'times' is a named (col names in data) numeric vector indicating ages/times at samples
#        'controls' is an integer vector indicating which samples in 'data' are controls
#        'donor' is a character vector telling a donor for each sample in 'data'. Can be NA
# Output: Returns a version of 'data' so that all pathway genes have been neutralised from the 
#         effect of given variables (e.g. age).
NeutraliseCoefficients = function(data, info, neutralise, userandom, labels, pathwaygenes, neutralisationformula=NULL, neutralisationmodel="rlm"){
  
  # Preprocess inputs
  pathwaygenes = intersect(pathwaygenes, rownames(data))
  
  # Extract control data for pathway genes and control info
  controlindex = which(labels == 0)
  controldata = data[pathwaygenes,controlindex]
  controlinfo = info[colnames(controldata),,drop=FALSE]
  
  # initialize coefficient data frame (only for variables to be neutralised)
  categoricalindex = which(apply(info[,neutralise,drop=FALSE],2,function(x){all(is.na(as.numeric(x)))}))
  numericindex = setdiff(1:length(neutralise), categoricalindex)
  cols = NULL
  if(length(numericindex) > 0) cols = c(cols, neutralise[numericindex])
  if(length(categoricalindex) > 0){
    categorical = neutralise[categoricalindex]
    for(cn in neutralise[categoricalindex]) cols = c(cols, paste(cn,setdiff(unique(info[,cn]),NA),sep=""))
  }
  coeftemplate = matrix(0, nrow=length(pathwaygenes), ncol=length(cols))
  colnames(coeftemplate) = cols
  rownames(coeftemplate) = pathwaygenes
  
  # Construct formula (or use the one given in argument 'neutralisationformula')
  if(is.null(neutralisationformula)){
    if(!is.null(userandom)){
      if(neutralisationmodel == "rlm"){
        prefixes = rep(neutralise, times=length(userandom))
        suffixes = rep(userandom, each=length(neutralise))
        righthand = paste(prefixes, suffixes, sep="*")
        useformula = paste("Expression", righthand, sep="~")
      } else{
        useformula = paste("Expression~", paste(neutralise, collapse="+"), sep="")
        for(r in userandom){
          addrandom = paste(paste(neutralise, r, sep="|"),collapse=")+(")
          useformula = paste(c(useformula,"+(",addrandom,")"), collapse="")
        }
      }
    } else{
      useformula = paste("Expression~", paste(neutralise, collapse="+"), sep="")
      if(neutralisationmodel != "rlm"){
        print("\n")
        print("NOTE: as no random effect variables were given, neutralisation model is changed to 'rlm'.")
        print("\n")
        neutralisationmodel = "rlm"
      } 
    }
    useformula = as.formula(useformula)
  } else useformula = as.formula(neutralisationformula)
  
  # Estimate coefficients to be neutralised
  coefficients_cntrl = switch(neutralisationmodel, 
                              rlm = NeutraliseCoefficients_rlm(controldata, controlinfo, useformula, coeftemplate),
                              rlmer = NeutraliseCoefficients_rlmer(controldata, controlinfo, useformula, coeftemplate),
                              lmer = NeutraliseCoefficients_lmer(controldata, controlinfo, useformula, coeftemplate))
  
  # Remove the effect of given coefficients to be neutralized (NOTE: intercept is not to be neutralised)
  residuals = data[pathwaygenes,]
  if(length(numericindex) > 0){          # neutralise numerical coefficients
    for(n in neutralise[numericindex]) residuals = residuals - coefficients_cntrl[,n,drop=FALSE] %*% t(info[,n,drop=FALSE])
  }
  if(length(categoricalindex) > 0){ # neutralise categorical coefficients
    for(m in categorical){
      residuals = residuals - coefficients_cntrl[,paste(m,info[,m],sep="")]
    }
  }
  data[pathwaygenes,] = residuals
  
  return(data)
}