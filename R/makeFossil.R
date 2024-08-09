#' @title Make fossil species on a phylogeny
#' @description This function takes an object of class \code{'phylo'} and
#'   randomly changes the lengths of the leaves.
#' @usage makeFossil(tree,p=0.5,ex=0.5)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param p the proportion of tips involved. By default it is half of the number
#'   of tips.
#' @param ex the multiplying parameter to change the leaf lengths. It is set at
#'   0.5 by default.
#' @export
#' @return The function produces a phylogeny having the same backbone of the
#'   original one.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' makeFossil(tree=Tstage)


makeFossil<-function(tree,
                     p=0.5,
                     ex=0.5){
  #require(ape)

  which(tree$edge[,2]<(Ntip(tree)+1))->leaves
  Ntip(tree)*p->prop
  sample(leaves,prop)->ext.leaves
  tree$edge.length[ext.leaves]*ex->tree$edge.length[ext.leaves]
  return(tree)
}
