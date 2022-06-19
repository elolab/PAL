# Input: 'data' is a numeric matrix (cols: samples, rows: pathways) with wanted effects neutralized
#
#
#
PickImportant_rlm = function(data, mainvariable, userandom, info, pathwayformula, seed, n=1000){ 
  
  # Pick samples with the main feature defined
  keepindex = which(apply(info[,mainvariable,drop=FALSE],1,function(s){any(!is.na(s))}))
  usedata = data[,keepindex]
  info = info[keepindex,,drop=FALSE]
  
  # Construct formula
  if(is.null(pathwayformula)){
    if(!is.null(userandom)){
      prefixes = rep(mainvariable, times=length(userandom))
      suffixes = rep(userandom, each=length(mainvariable))
      righthand = paste(prefixes, suffixes, sep="*")
      useformula = paste("Expression", righthand, sep="~")
    } else useformula = paste("Expression~", paste(mainvariable, collapse="+"), sep="")
    useformula = as.formula(useformula)
  } else useformula = as.formula(neutralisationformula)
  
  # Extract main variables (otherwise as is, but in case of categorical ones, different classes get their own coefficients and p-values)
  categoricalindex = which(apply(info[,mainvariable,drop=FALSE],2,function(m){all(is.na(as.numeric(m)))}))
  numericindex = setdiff(1:length(mainvariable), categoricalindex)
  maincols = NULL
  if(length(numericindex) > 0) maincols = c(maincols, mainvariable[numericindex])
  if(length(categoricalindex) > 0){
    for(cv in mainvariable[categoricalindex]) maincols = c(maincols, paste(cv,setdiff(unique(info[,cv]),NA),sep=""))
  }
  
  # Initialise coefficient and p-value tables
  coefficients = matrix(NA, nrow=nrow(usedata), ncol=length(maincols))
  colnames(coefficients) = maincols
  rownames(coefficients) = rownames(usedata)
  pvalues = coefficients
  
  # Initialise coefficient matrix (empty template) for p-value calculation
  coefstemplate = matrix(NA, ncol=length(maincols), nrow=n)
  colnames(coefstemplate) = maincols
  
  # Analyse one pathway at time
  for(i in 1:nrow(data)){
    
    # Construct data frame for the model fitting
    scoreandinfo = as.data.frame(cbind(as.numeric(usedata[i,]), info))
    colnames(scoreandinfo) = c("Score",colnames(info))
    
    # Fit model
    fittedmodel = MASS::rlm(useformula, data=scoreandinfo)
    coefficients_row = fittedmodel$coefficients
    pickcols = intersect(maincols, names(coefficients_row))
    coefficients[i, pickcols] = coefficients_row[pickcols]
    
    # Generate random sampling indices for the 'mainvariable'
    set.seed(seed+i)
    randomindices = t(data.frame(lapply(1:n, function(r){sample(1:ncol(usedata), size=ncol(usedata), replace=FALSE)})))
    
    # Initialise coefficient matrix for p-value calculation
    coefs = coefstemplate
    
    # collect coefficients from different random samplings
    sampleddata = scoreandinfo
    for(j in 1:n){
      sampleddata[,mainvariable] = scoreandinfo[randomindices[j,],mainvariable]
      coefs[j,pickcols] = MASS::rlm(useformula, data=sampleddata)$coefficients[pickcols]
    }
    
    # Calculate and store p-values
    for(m in pickcols){
      realgreater = sum(coefficients_row[m] > coefs[,m])
      pvalues[i,m] = 2*pmin(realgreater/n, (n-realgreater)/n)
    } 
    
  }
  
  # Calculate p-values and convert them into fdr
  fdr = apply(pvalues, 2, p.adjust, method="fdr")
  
  # Combine p-values, FDR and coefficients into one table
  significance = cbind(pvalues, cbind(fdr, coefficients))
  colnames(significance) = paste(rep(c("Pval","FDR","coef"),each=length(maincols)), rep(maincols, times=3), sep="_")
  rownames(significance) = rownames(usedata)
  
  return(significance)
}