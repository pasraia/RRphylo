#' @title Upward tip or node to root path
#' @description This function is a wrapper around \pkg{phytools}
#'   \code{getDescendants} (\cite{Revell 2012}). It returns the node path from a
#'   given node or species to the root of the phylogeny.
#' @usage getMommy(tree,N)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param N the number of node or tip to perform the function on. The function also works with tip labels.
#' @export
#' @return The function produces a vector of node numbers as integers, collated
#'   from a node or a tip towards the tree root.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Revell, L. J. (2012). phytools: An R package for phylogenetic
#' comparative biology (and other things). \emph{Methods in Ecology and
#' Evolution}, 3: 217-223.doi:10.1111/j.2041-210X.2011.00169.x
#' @examples
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' getMommy(tree=Tstage,N=12)

getMommy<-function(tree,N){
  if (is.character(N)) if(N%in%tree$tip.label) N <- which(tree$tip.label == N) else N<-as.numeric(N)

  N->node
  curr<-vector()
  daughters<-tree$edge[which(tree$edge[,2]==node),1]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) curr<-c(curr,getMommy(tree,daughters[w[i]]))
  return(curr)
}


