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

distNodes<-function(tree,n)
{
  #require(phytools)
  #require(data.tree)
  #require(ape)

  tree$node.label<-NULL
  sis <- getDescendants(tree, (Ntip(tree) + 1))
  sis <- c((Ntip(tree) + 1), sis)
  sis <- sis[-which(sis <= Ntip(tree))]
  if (length(which(sis == n[1])) > 0)
    sis <- sis[-which(sis == n[1])] else sis <- sis
  treeN <- as.Node(tree)
  distN <- array()
  if (length(n) > 1) {
    distN <- Distance(FindNode(treeN, as.character(n[1])),
                      FindNode(treeN, as.character(n[2])))
    distL <- dist.nodes(tree)[which(as.numeric(colnames(dist.nodes(tree))) ==
                                      n[1]), which(as.numeric(colnames(dist.nodes(tree))) ==
                                                     n[2])]
    dfN <- data.frame(n.nodes = distN, length = distL)
  }else {
    for (i in 1:length(sis)) {
      distN[i] <- Distance(FindNode(treeN, as.character(n)),
                           FindNode(treeN, as.character(sis[i])))
    }
    dfN <- data.frame(node = sis, n.nodes = distN)
    distL <- dist.nodes(tree)[which(as.numeric(colnames(dist.nodes(tree))) >
                                      Ntip(tree)), which(as.numeric(colnames(dist.nodes(tree))) ==
                                                           n[1])]
    if (length(which(names(distL) == n)) > 0)
      distL <- distL[-which(names(distL) == n)]else distL <- distL
    dfN <- data.frame(dfN, length = distL[match(dfN[, 1],
                                                names(distL))])
  }
  return(dfN)
}


