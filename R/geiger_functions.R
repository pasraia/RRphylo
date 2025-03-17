#' @title Cross-reference tree and data
#' @description The function matches data names with tree tips. If either there
#'   is no data for a tip or it is not present on the tree, the function removes
#'   the entry from both.
#' @usage treedataMatch(tree,y)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param y named variable. It can be a vector or a multivariate dataset or a 3D
#'   array. Alternatively, \code{y} can also be a vector of species names.
#' @return The function returns a \code{list} object. If no mismatch between
#'   \code{tree} and \code{y} is detected, the list only includes the matrix of
#'   \code{y} ordered according to the order of tips on the tree (\code{$y}). If
#'   some tips on the \code{tree} are missing from \code{y}, they are removed
#'   from the phylogeny. Thus, the list also includes the pruned tree
#'   (\code{$tree}) and the vector of dropped tips (\code{$removed.from.tree}).
#'   Similarly, if some entries in \code{y} are missing from the \code{tree},
#'   the list also includes the vector of mismatching entry names
#'   (\code{$removed.from.y}). In this latter case, the first element of the
#'   list (\code{$y}) does not include the entries \code{$removed.from.y}, so
#'   that it perfectly matches the phylogeny.
#' @export
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @examples
#' data(DataCetaceans)
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#'
#' treedataMatch(tree=treecet,y=masscet)
#' treedataMatch(tree=treecet,y=brainmasscet)
#' treedataMatch(tree=treecet,y=names(brainmasscet))


treedataMatch<-function (tree, y){
  #if(!inherits(y,"matrix")&!inherits(y,"data.frame")) as.matrix(y)->y
  if(length(dim(y))<3){
    if(is.null(ncol(y))||is.na(ncol(y))) as.matrix(y)->y
    rownames(y)->ynams
  }else dimnames(y)[[3]]->ynams
  if(is.null(ynams)){
    if(any(y%in%tree$tip.label)) rownames(y)<-ynams<-y else stop("y needs to be named")
  }
  if(all(!tree$tip.label%in%ynams)) stop("There is no match between tree tip labels and y names")

  if(!all(tree$tip.label%in%ynams)){
    tree$tip.label[which(!tree$tip.label%in%ynams)]->rem.tree
    drop.tip(tree,which(!tree$tip.label%in%ynams))->tree
  }else rem.tree<-NULL

  if(!all(ynams%in%tree$tip.label)){
    ynams[which(!ynams%in%tree$tip.label)]->rem.y
    #y[which(ynams%in%tree$tip.label),,drop=FALSE]->y
  }else rem.y<-NULL

  if(length(dim(y))<3) y[match(tree$tip.label,rownames(y)),,drop=FALSE]->y else
    y[,,match(tree$tip.label,dimnames(y)[[3]])]->y

  if(all(is.character(y))&all(y%in%tree$tip.label)) rownames(y)<-NULL
  list(y=y)->res
  if(!is.null(rem.tree)) c(res,list(tree=tree,removed.from.tree=rem.tree))->res
  if(!is.null(rem.y)) c(res,removed.from.y=list(rem.y))->res
  return(res)
}


#' @title Brownian Motion rate computation
#' @description The function computes rate of phenotypic evolution along a phylogeny assuming Brownian Motion model of evolution.
#' @usage sig2BM(tree,y)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param y either a single vector variable or a multivariate dataset. In any
#'   case, \code{y} must be named.
#' @return The Brownian Motion rate of phenotypic evolution for each variable in \code{y}.
#' @export
#' @importFrom ape pic
#' @author Pasquale Raia, Silvia Castiglione
#' @examples
#'
#' ### Univariate data ###
#' data(DataCetaceans)
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' sig2BM(tree=treecet,y=masscet)
#'
#' ### Multivariate data ###
#' data(DataUng)
#' DataUng$treeung->treeung
#' DataUng$PCscoresung->PCung
#' sig2BM(tree=treeung,y=PCung)

sig2BM<-function(tree,y){
  if (is.binary(tree))
    tree <- tree else tree <- multi2di(tree, random = FALSE)
  as.matrix(y)->y
  apply(y,2,function(e){
    pic(e,tree)->picc
    sum(picc^2)/(length(picc))
  })->sigg
  if(!is.null(colnames(y))) colnames(y)->names(sigg) else
    sapply(1:ncol(y),function(j) paste("y",j,sep=""))->names(sigg)
  return(sigg)
}



#' @title Get descending tips
#' @description The function returns the numbers or labels of tips descending from a given node.
#' @usage tips(tree,node,labels=TRUE)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param node the number of focal node
#' @param labels if \code{TRUE} (default) the function returns the labels of descending tips.
#' @return The tips, either labels or numbers depending on the argument \code{labels}, descending from the \code{node}.
#' @export
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @examples
#' data(DataOrnithodirans)
#' DataOrnithodirans$treedino->treedino
#' tips(tree=treedino,node=677,labels=FALSE)
#' tips(tree=treedino,node=677,labels=TRUE)

tips<-function(tree,node,labels=TRUE){
  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }
  getDescendants(tree,node)->des
  if(isTRUE(labels)) tree$tip.label[sort(des[which(des<=Ntip(tree))])] else
    sort(des[which(des<=Ntip(tree))])
}
