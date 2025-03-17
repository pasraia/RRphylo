#' @title Taxonomic inspection of the tree at the genus level
#' @description The function returns the most recent common ancestor and the
#'   number of species belonging to each or some user-specified genera within
#'   the phylogenetic tree.
#' @usage getGenus(tree,genera=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous. Generic name and specific epithet must be separated by
#'   '_'.
#' @param genera a character vector including one or more genera to focus on.
#'   Please notice the function is case sensitive.
#' @export
#' @return The function returns a data-frame including the number of species,
#'   the most recent common ancestor of each genera (if a genus includes one
#'   species this is the species tip number), and
#'   whether the genera form monophyletic or paraphyletic clades.
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @examples
#' DataCetaceans$treecet->treecet
#'
#' getGenus(treecet)->gg1
#' getGenus(treecet,c("Mesoplodon","Balaenoptera"))->gg2

getGenus<-function(tree,genera=NULL){
  # require(ape)
  # require(RRphylo)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  if(any(!grepl("_",tree$tip.label))) stop("Generic name and specific epithet must be separated by _ .")

  table(sapply(strsplit(tree$tip.label,"_"),"[[",1))->gens
  sapply(1:length(gens),function(x) paste(names(gens)[x],"_",sep=""))->names(gens)
  if(!is.null(genera)) {
    sapply(genera,function(x) paste(x,"_",sep=""))->genera
    gens[match(genera,names(gens))]->gens

    if(all(is.na(gens))){
      stop("Required genera not on the phylogenetic tree")
    }else if(any(is.na(gens))){
      warning(paste(genera[which(!genera%in%names(gens))]," not on the phylogenetic tree"),immediate. = TRUE)
      gens[which(!is.na(gens))]->gens
    }
  }
  sapply(1:length(gens),function(x){
    if(gens[x]>1) getMRCA(tree,tree$tip.label[grep(names(gens)[x],tree$tip.label)]) else
      grep(names(gens)[x],tree$tip.label)
  })->nns
  cond<-mapply(a=nns,b=gens,function(a,b) ifelse(a<=Ntip(tree),"mono",ifelse(length(tips(tree,a))>b,"para","mono")))
  gsub("_","",names(gens))->names(gens)
  data.frame(as.table(gens),nns,cond)->df
  colnames(df)<-c("genus","nspecies","mrca","condition")
  return(df)
}
