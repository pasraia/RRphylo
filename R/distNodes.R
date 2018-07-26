#' @title Finding distance between nodes
#' @description The function returns the distance between pairs of nodes. The distance is meant as both patristic distance and the number of nodes intervening between the pair.
#' @usage distNodes(tree,n)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and fully dichotomous.
#' @param n either a single node or a pair of nodes.
#' @export
#' @importFrom data.tree Distance FindNode as.Node
#' @return The function returns a data frame with distances between the focal node and the other nodes on the tree (or for the selected pair only).
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' distNodes(tree=Tstage,n=64)
#' distNodes(tree=Tstage,n=c(64,48))

distNodes<-function(tree,n){
  #require(phytools)
  #require(data.tree)
  #require(ape)

  getDescendants(tree,(Ntip(tree)+1))->sis
  c((Ntip(tree)+1),sis)->sis
  sis[-which(sis<=Ntip(tree))]->sis
  if(length(which(sis==n[1]))>0) sis[-which(sis==n[1])]->sis else sis->sis

  as.Node(tree)->treeN
  distN<-array()
  if(length(n)>1) {
    Distance(FindNode(treeN,as.character(n[1])),FindNode(treeN,as.character(n[2])))->distN
    dist.nodes(tree)[which(as.numeric(colnames(dist.nodes(tree)))==n[1]),which(as.numeric(colnames(dist.nodes(tree)))==n[2])]->distL
    data.frame(n.nodes=distN,length=distL)->dfN
  }else{
    for(i in 1:length(sis)){
      Distance(FindNode(treeN,as.character(n)),FindNode(treeN,as.character(sis[i])))->distN[i]
    }
    data.frame(node=sis,n.nodes=distN)->dfN
    dist.nodes(tree)[which(as.numeric(colnames(dist.nodes(tree)))>Ntip(tree)),which(as.numeric(colnames(dist.nodes(tree)))==n[1])]->distL
    if(length(which(names(distL)==n))>0) distL[-which(names(distL)==n)]->distL else distL->distL
    data.frame(dfN,length=distL[match(dfN[,1],names(distL))])->dfN
  }
  return(dfN)
}


