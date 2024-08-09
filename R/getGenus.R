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
#' @return The function returns a data-frame including the number of species and
#'   the most recent common ancestor of each genera.
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @examples
#' DataCetaceans$treecet->tree
#'
#' getGenus(tree)
#' getGenus(tree,c("Mesoplodon","Balaenoptera"))

getGenus<-function(tree,genera=NULL){
  # require(ape)
  # require(RRphylo)
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
      getMommy(tree,grep(names(gens)[x],tree$tip.label))[1]
  })->nns
  gsub("_","",names(gens))->names(gens)
  data.frame(as.table(gens),nns)->df
  colnames(df)<-c("genus","nspecies","mrca")
  return(df)
}
