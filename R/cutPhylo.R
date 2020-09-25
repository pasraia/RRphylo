#' @title Cut the phylogeny at a given age or node
#' @description The function cuts all the branches of the phylogeny which are
#'   younger than a specific age or node (i.e. the age of the node).
#' @usage cutPhylo(tree,age=NULL,node=NULL,keep.lineage=TRUE)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param age the age (in terms of time distance from the recent) at which the
#'   tree must be cut
#' @param node the node whose age must be used as cutting limit.
#' @param keep.lineage logical specifying whether lineages with no descendant tip must be retained (see example below). Default is \code{TRUE}.
#' @export
#' @seealso \href{../doc/Tree-Manipulation.html#cutPhylo}{\code{cutPhylo} vignette}
#' @importFrom phytools drop.clade
#' @importFrom ape axisPhylo
#' @details When an entire lineage is cut (i.e. one or more nodes along a path) and \code{keep.lineages = TRUE},
#'   the leaves left are labeled as "l" followed by a number.
#' @return The function returns the cut phylogeny and plots it into the graphic
#'   device. The time axis keeps the root age of the original tree.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' \dontrun{
#' library(ape)
#'
#' set.seed(22)
#' rtree(100)->tree
#' 3->age
#'
#' cutPhylo(tree,age=age)->t1
#' cutPhylo(tree,age=age,keep.lineage=FALSE)->t1
#' cutPhylo(tree,node=151)->t2
#' cutPhylo(tree,node=151,keep.lineage=FALSE)->t2
#' }

cutPhylo<-function(tree,age=NULL,node=NULL,keep.lineage=TRUE){
  # require(ape)
  # require(geiger)
  # require(phytools)
  # require(picante)

  if (!requireNamespace("picante", quietly = TRUE)) {
    stop("Package \"picante\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }

  distNodes(tree,(Ntip(tree)+1))->dN
  if(is.null(node)) max(nodeHeights(tree))-age->cutT else dN[match(node,rownames(dN)),2]->cutT
  dN[,2]->dd
  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter

  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
  #}
  }else{
  ### Tips and nodes ###
  #if(any(suppressWarnings(as.numeric(cutter))>Ntip(tree))){
    as.numeric(cutter[which(suppressWarnings(as.numeric(cutter))>Ntip(tree))])->cutn
    i=1
    while(i<=length(cutn)){
      if(any(cutn%in%getDescendants(tree,cutn[i]))) cutn[-which(cutn%in%getDescendants(tree,cutn[i]))]->cutn
      i=i+1
    }

    tree->tt
    i=1
    while(i<=length(cutn)){
      getMRCA(tt,tips(tree,cutn[i]))->nn
      drop.clade(tt,tips(tt,nn))->tt
      tt$tip.label[which(tt$tip.label=="NA")]<-paste("l",i,sep="")
      i=i+1
    }

    diag(vcv(tt))->times
    if(any(times>cutT)){
      times[which(times>cutT)]->times
      tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]<-
        tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]-(times-cutT)
    }

  }

  if(isFALSE(keep.lineage)){
    drop.tip(tt,which(!tt$tip.label%in%tree$tip.label))->tt
    if(Ntip(tt)<=100) plot(tt,cex=.8) else plot(tt,cex=.5 )
  }else{
    tip.col<-rep("black",Ntip(tt))
    tip.col[which(!tt$tip.label%in%tree$tip.label)]<-"red"
    if(Ntip(tt)<=100) plot(tt,cex=.8,tip.color=tip.col) else plot(tt,cex=.5,tip.color=tip.col)
  }

  axisPhylo(root.time = max(nodeHeights(tree)))

  return(tt)
}
