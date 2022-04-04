#
#
#
#
ExtractPathwayInfo = function(pathways, results, nodevalues, relationvalues=NA){
  
  pathways = pathways[gsub(":.*", "", rownames(results))]
  
  # Extract the number of measured nodes (median over samples) from each pathway and the proportion of inhibitory nodes
  measurednodes = unlist(lapply(nodevalues, function(x){median(colSums(!is.na(x)))}))
  inhibitingnodes = mapply(function(x,y){
    ids = rownames(y)
    nodes = x$nodeinfo[x$nodeinfo$Id %in% ids,]
    return(sum(nodes$Role == -1)/nrow(nodes))
  }, x=pathways, y=nodevalues)
  
  # If available, extract the number of analysed relations and proportion of relations that should be inactive when the pathway is active
  if(!is.na(relationvalues)){
    measuredrelations = unlist(lapply(relationvalues, function(x){median(colSums(!is.na(x)))}))
    inhibitingrelations = mapply(function(x,y){
      indices = which(rowSums(!is.na(y)) > 0)
      relations = x$relationinfo[indices,]
      return(sum(relations$Role == -1)/nrow(relations))
    }, x=pathways, y=relationvalues)
  }
  
  # Collect extracted info into a matrix
  if(is.na(relationvalues)){
    info = cbind(measurednodes, inhibitingnodes)
    colnames(info) = c("MeasuredNodes", "InhibitingNodeProportion")
  } else{
    info = cbind(cbind(measurednodes, measuredrelations), cbind(inhibitingnodes, inhibitingrelations))
    colnames(info) = c("MeasuredNodes", "MeasuredRelations", "InhibitingNodeProportion", "InhibitingRelationProportion")
  }
  rownames(info) = rownames(results)
  
  return(info)
}