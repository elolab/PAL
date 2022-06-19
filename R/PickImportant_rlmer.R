# Input: 'data' is a numeric matrix (cols: samples, rows: pathways) with wanted effects neutralized
#
#
#
PickImportant_rlmer = function(data, mainvariable, userandom, info, pathwayformula){ 
  
  # Pick samples with the main feature defined
  keepindex = which(apply(info[,mainvariable,drop=FALSE],1,function(s){any(!is.na(s))}))
  usedata = data[,keepindex]
  info = info[keepindex,,drop=FALSE]
  
  # Construct formula
  if(is.null(pathwayformula)){
    useformula = paste("Score~", paste(mainvariable, collapse="+"), sep="")
    if(!is.null(userandom)){
      for(r in userandom){
        addrandom = paste(paste(mainvariable, r, sep="|"),collapse=")+(")
        useformula = paste(c(useformula,"+(",addrandom,")"), collapse="")
      }
    } 
    useformula = as.formula(useformula)
  } else useformula = as.formula(pathwayformula)
  
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
  
  # Analyse one pathway at time
  for(i in 1:nrow(usedata)){
    
    # Construct data frame for the model fitting
    scoreandinfo = as.data.frame(cbind(as.numeric(usedata[i,]), info))
    colnames(scoreandinfo) = c("Score",colnames(info))
    
    # Fit model and extract coefficients
    fittedmodel = robustlmm::rlmer(useformula, data=scoreandinfo)
    coefs = coef(summary(fittedmodel))
    pickcols = intersect(maincols, rownames(coefs))
    coefficients[i,pickcols] = coefs[pickcols,"Estimate"] 
    
    # fit a non-robust model and extract degree of freedom for p-value calculation purposes
    nonrobustmodel = lmerTest::lmer(useformula, data=scoreandinfo)
    dfs = coef(summary(nonrobustmodel))[pickcols,"df"]
    
    # calculate p-values based on robust t-values and non-robust approximated degree of freedoms
    pvalues[i,pickcols] = 2*pt(abs(coefs[pickcols,"t value"]), dfs, lower=FALSE)
    # NOTE: Package sjPlot has function tab_model, which can also calculate p-values based on Satterthwaite-approximated dfs.
    # However, as this implementation was simple, I decided to not use sjPlot. 
  } 
  
  # Calculate FDR
  fdr = apply(pvalues, 2, p.adjust, method="fdr")
  
  # Combine p-values, FDR and coefficients into one table
  significance = cbind(pvalues, cbind(fdr, coefficients))
  colnames(significance) = paste(rep(c("Pval","FDR","coef"),each=length(maincols)), rep(maincols, times=3), sep="_")
  rownames(significance) = rownames(usedata)
  
  return(significance)
}