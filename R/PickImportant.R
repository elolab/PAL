# Input: 'data' is a numeric matrix (cols: samples, rows: pathways) with wanted effects neutralized
#
#
#
PickImportant = function(data, mainfeature, seed, info){ 
  
  # Pick samples with the main feature defined
  index = which(!is.na(mainfeature))
  feature = mainfeature[index]
  usedata = data[,index]
  if(!is.na(info[1])) info = info[index,,drop=F]
  
  # Initialize coefficient matrix
  if(is.numeric(feature)){
    cols = 1
    coefnames = "feature"
  } else{
    cols = length(unique(feature))
    coefnames = paste("feature", unique(feature), sep="")
  }
  coefficient = matrix(0, nrow=nrow(data), ncol=cols)
  rownames(coefficient) = rownames(data)
  colnames(coefficient) = coefnames
  
  # Initialize coefficient matrix for randomised main feature
  coefs = array(0, dim=c(1000, nrow(usedata), ncol(coefficient)))
  featurenames = paste("temp",colnames(coefficient),sep="")
  dimnames(coefs) = list(paste("R",1:1000,sep=""),rownames(usedata),featurenames)
  
  # Also other variables (e.g. donor) than the main one  
  if(!is.na(info[1])){
    
    # Fit model and collect the coefficient
    for(i in 1:nrow(data)){
      
      # Construct data frame for the model fitting
      expressionandinfo = as.data.frame(cbind(cbind(as.numeric(usedata[i,]), feature), info))
      colnames(expressionandinfo) = c("Expression","feature",colnames(info))
      
      modelres = MASS::rlm(Expression~., data=expressionandinfo)
      coefficient[i,coefnames] = modelres$coefficients[coefnames]
    }
    
    # Calculate coefficients with randomly sampled main variable (takes forever)
    for(i in 1:nrow(coefs)){
      set.seed(seed+i)
      tempfeature = sample(feature, size=length(feature), replace=FALSE)
      for(j in 1:nrow(usedata)){
        tempdata = as.data.frame(cbind(cbind(as.numeric(usedata[j,]), tempfeature), info))
        colnames(tempdata) = c("Expression","tempfeature",colnames(info))
        coefs[i,j,featurenames] = MASS::rlm(Expression~., data=tempdata)$coefficients[featurenames]
      }
    }
    
  # No other variables than the main one  
  } else{
    
    # Fit model and collect the coefficient
    for(i in 1:nrow(data)){
      
      modelres = MASS::rlm(as.numeric(usedata[i,])~feature)
      coefficient[i,coefnames] = modelres$coefficients[coefnames]
    }
    
    # Calculate coefficients with randomly sampled main variable (takes forever)
    for(i in 1:nrow(coefs)){
      set.seed(seed+i)
      tempfeature = sample(feature, size=length(feature), replace=FALSE)
      for(j in 1:nrow(usedata)){
        coefs[i,j,featurenames] = MASS::rlm(as.numeric(usedata[j,])~tempfeature)$coefficients[featurenames]
      }
    }
  }
  
  # Initialize a matrix for p-value calculation
  realgreater = matrix(0, nrow=nrow(usedata), ncol=ncol(coefficient))
  rownames(realgreater) = rownames(usedata)
  colnames(realgreater) = gsub("feature","",colnames(coefficient))
  
  # Calculate p-values and convert them into fdr
  for(n in 1:dim(coefs)[1]) realgreater = realgreater + (coefficient >= as.matrix(coefs[n,,]))
  realsmaller = nrow(coefs) - realgreater
  pvals = 2*pmin(realgreater/nrow(coefs), realsmaller/nrow(coefs))
  fdr = apply(pvals, 2, p.adjust, method="fdr")
  
  # Combine p-values and fdr into one significance matrix
  significance = cbind(coefficient, cbind(pvals, fdr))
  rownames(significance) = rownames(pvals)
  colnames(significance) = c(gsub("feature","coefficient",colnames(coefficient)), paste("pval",colnames(pvals),sep="_"), paste("fdr",colnames(pvals),sep="_"))

  return(significance)
}