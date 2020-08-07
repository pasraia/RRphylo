#' @title Upward tip or node to root path
#' @description This function is a wrapper around \pkg{phytools}
#'   \code{getDescendants} (\cite{Revell 2012}). It returns the node path from a
#'   given node or species to the root of the phylogeny.
#' @usage getMommy(tree,N,curr=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param N the number of node or tip to perform the function on. Notice the
#'   function only works with number, not tip labels.
#' @param curr has not to be provided by the user.
#' @export
#' @details The object \code{'curr'} is created inside the function in order to
#'   produce an array of nodes on the path.
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

getMommy<-function(tree,
                   N,
                   curr=NULL){

  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }

  N->node
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,2]==node),1]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getMommy(tree,daughters[w[i]],curr)
  return(curr)
}
