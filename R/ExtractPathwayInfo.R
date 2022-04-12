#
#
#
#
ExtractPathwayInfo = function(pathways, results, nodevalues, relationvalues=NA){
  
  pathways = pathways[gsub(":.*", "", rownames(results))]
  
  # Extract the number of measured nodes (median over samples) from each pathway and the proportion of group nodes
  measurednodes = unlist(lapply(nodevalues, function(x){median(colSums(!is.na(x)))}))
  groupnodeproportion = mapply(function(x,y){
    ids = rownames(y)
    nodes = x$nodeinfo[x$nodeinfo$Id %in% ids,]
    return(sum(nodes$Type != "gene")/nrow(nodes))
  }, x=pathways, y=nodevalues)
  
  # If available, extract the number of analysed relations and proportion of inhibiting relations
  if(is.list(relationvalues)){
    measuredrelations = unlist(lapply(relationvalues, function(x){
      if(any(!is.na(x))){
        return(median(colSums(!is.na(x))))
      } else return(0)
      }))
    inhibitingrelations = mapply(function(x,y){
      if(any(!is.na(y))){
        indices = which(rowSums(!is.na(y)) > 0)
        relations = x$relationinfo[indices,]
        return(sum(relations$Direction == "inhibition")/nrow(relations))
      } else return(NA)
    }, x=pathways, y=relationvalues)
  }
  
  # Collect extracted info into a matrix
  if(!is.list(relationvalues)){
    info = cbind(measurednodes, groupnodeproportion)
    colnames(info) = c("MeasuredNodes", "GroupNodeProportion")
  } else{
    info = cbind(cbind(measurednodes, measuredrelations), cbind(groupnodeproportion, inhibitingrelations))
    colnames(info) = c("MeasuredNodes", "MeasuredRelations", "GroupNodeProportion", "InhibitingRelationProportion")
  }
  rownames(info) = rownames(results)
  
  return(info)
}