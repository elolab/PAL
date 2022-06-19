# Input: 'data' is a matrix of values (rows: genes/pathways/whatever, cols: samples)
#        'times' is a named (col names in data) numeric vector indicating ages/times at samples
#        'controls' is an integer vector indicating which samples in 'data' are controls
#        'donor' is a character vector telling a donor for each sample in 'data'. Can be NA
# Output: Returns a version of 'data' so that all pathway genes have been neutralised from the 
#         effect of given variables (e.g. age).
NeutraliseCoefficients_rlmer = function(usedata, useinfo, useformula, coefficients){
  
  # Record coefficients 
  for(p in rownames(usedata)){
    
    # Construct data frame (controls) for the model fitting
    expressionandinfo = as.data.frame(cbind(as.numeric(t(usedata[p,,drop=FALSE])), useinfo))
    colnames(expressionandinfo) = c("Expression",colnames(useinfo))
    
    # Calculate coefficients based on control samples
    modelres = robustlmm::rlmer(useformula, data=expressionandinfo)
    fixedcoefs_cntrl = coef(summary((modelres)))[,"Estimate"]
    #randomcoefs_cntrl = lme4::ranef(modelres) # Does it make sense to allow random effects to be neutralised?
    
    # Record those coefficients
    coefnames = intersect(colnames(coefficients), names(fixedcoefs_cntrl))
    coefficients[p,coefnames] = fixedcoefs_cntrl[coefnames]
  }
  
  return(coefficients)
}