#' @title Get sister clade
#' @description The function identifies and returns the sister clade of a given
#'   node/tip.
#' @usage getSis(tree,n,printZoom=TRUE)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param n number of focal node or name of focal tip.
#' @param printZoom if \code{TRUE} the function plots the tree section of
#'   interest.
#' @return The sister node number or sister tip name. In case of polytomies, the
#'   function returns a vector.
#' @export
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @importFrom ape zoom
#' @examples
#' data(DataOrnithodirans)
#' DataOrnithodirans$treedino->treedino
#' getSis(tree=treedino,n=677,printZoom=FALSE)
#' getSis(tree=treedino,n="Shenzhoupterus_chaoyangensis",printZoom=FALSE)



getSis<-function(tree,n,printZoom=TRUE)
{
  #require(ape)
  #require(geiger)
  if(class(n)=="character"){
    which(tree$tip.label==n)->n}
  tree$edge[which(tree$edge[,2]==n),1]->mom
  tree$edge[which(tree$edge[,1]==mom),2]->daug
  daug[which(daug!=n)]->sis
  if(length(which(sis<=Ntip(tree)))>0){
    tree$tip.label[sis[which(sis<=Ntip(tree))]]->sis[which(sis<=Ntip(tree))]
  }
  if(printZoom==TRUE) zoom(tree,tips(tree,mom))
  return(sis)
}
