# Input: 'data' is a matrix of values (rows: genes/pathways/whatever, cols: samples)
#        'times' is a named (col names in data) numeric vector indicating ages/times at samples
#        'controls' is an integer vector indicating which samples in 'data' are controls
#        'donor' is a character vector telling a donor for each sample in 'data'. Can be NA
# Output: Returns a scaled version of 'data' so that all values are between -1 and 1. The values
#         have been neutralized from effect of donor (if available) and age.
NeutralizeCoefficients = function(data, info, neutralize, labels, pathwaygenes, mainfeature){
  
  # Preprocess inputs
  pathwaygenes = intersect(pathwaygenes, rownames(data))
  neutralizedcoefnames = colnames(info)[neutralize]
  
  # Extract info for numeric and categorical features
  info_numeric = info[,apply(info, 2, function(f){any(!is.na(as.numeric(f)))}),drop=F]
  info_categorical = info[,setdiff(colnames(info),colnames(info_numeric)),drop=F]
  
  # Extract control data for pathway genes
  controlindex = which(labels == 0)
  controldata = data[pathwaygenes,controlindex]
  
  # initialize coefficient data frame
  cathegorynames = unlist(lapply(intersect(colnames(info_categorical), neutralizedcoefnames), function(g){paste(g,unique(info[,g]),sep="")}))
  cols = c("Intercept_control", intersect(colnames(info_numeric), neutralizedcoefnames), cathegorynames)
  coefficients = matrix(0, nrow=length(pathwaygenes), ncol=length(cols))
  colnames(coefficients) = cols
  rownames(coefficients) = pathwaygenes
  
  controlinfo = cbind(info_numeric[controlindex,,drop=F],info_categorical[controlindex,,drop=F])
  colnames(controlinfo) = c(colnames(info_numeric),colnames(info_categorical))
  controlinfo = controlinfo[,which(apply(controlinfo,2,function(x){length(unique(x))})>1),drop=F]
  
  # Record coefficients 
  for(p in pathwaygenes){
    
    # Construct data frame (controls) for the model fitting
    expressionandinfo = as.data.frame(cbind(as.numeric(t(controldata[p,,drop=F])), controlinfo))
    colnames(expressionandinfo) = c("Expression",colnames(controlinfo))
    
    # Calculate coefficients based on control samples
    coefs_cntrl = MASS::rlm(Expression~., data=expressionandinfo)$coefficients
    
    # Record those coefficients
    shared = intersect(names(coefs_cntrl), colnames(coefficients))
    coefficients[p,"Intercept_control"] = coefs_cntrl["(Intercept)"]
    coefficients[p,shared] = coefs_cntrl[shared] 
  }
  
  # Remove the effect of given coefficients to be neutralized (NOTE: intercept is not to be neutralized)
  residuals = data[pathwaygenes,]
  if(length(intersect(colnames(info_numeric), neutralizedcoefnames)) > 0){          # neutralize numerical coefficients
    for(n in intersect(colnames(info_numeric), neutralizedcoefnames)) residuals = residuals - coefficients[,n,drop=F] %*% t(info_numeric[,n,drop=F])
  }
  if(length(intersect(colnames(info_numeric), neutralizedcoefnames)) < sum(neutralize)){ # neutralize categorical coefficients
    for(m in setdiff(neutralizedcoefnames,colnames(info_numeric))){
      residuals = residuals - coefficients[,paste(m,info[,m],sep="")]
    }
  }
  data[pathwaygenes,] = residuals
  
  return(data)
}