#'@title Tracing nodes along paths
#'@description Given a vector of nodes, the function collates nodes along
#'  individual lineages from the youngest (i.e. furthest from the tree root) to
#'  the oldest.
#'@usage node.paths(tree, vec)
#'@param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'  fully dichotomous.
#'@param vec a vector of node numbers
#'@export
#'@return A list of node paths, each starting from the youngest node (i.e.
#'  furthest from the tree root) and ending to the oldest along the path.
#'@author Silvia Castiglione, Pasquale Raia
#' @examples
#' require(ape)
#'
#' rtree(100)->tree
#' sample(seq(Ntip(tree)+1,Ntip(tree)+Nnode(tree)),20)->nods
#' plot(tree,show.tip.label=FALSE)
#' nodelabels(node=nods,frame="n",col="red")
#' node.paths(tree=tree, vec=nods)

node.paths<-function(tree,vec){
  # require(ape)
  # require(geiger)

  if(!is.numeric(vec)) as.numeric(vec)->vec
  # vec[order(dist.nodes(tree)[vec,(Ntip(tree)+1)],decreasing=TRUE)]->vec
  vec[order(sapply(vec,function(x) length(tips(tree,x))))]->vec

  lapply(vec,function(x) c(x,vec[which(vec%in%getMommy(tree,x))]))->nodP

  k=1
  while(k<=length(nodP)){
    nodP[[k]]->nodPk
    any(sapply(nodP[-k], function(w) all(nodPk%in%w)))->rem
    if(rem) nodP[-k]->nodP else k=k+1
  }

  return(nodP)
}
