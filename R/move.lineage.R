#' @title Move tips or clades
#' @description Move a single tip or an entire clade to a different position
#'   within the tree.
#' @usage move.lineage(tree,focal,sister,poly=FALSE,rescale=TRUE,rootage=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param focal the lineage to be moved. It can be either a tip name/number or a
#'   node number. If \code{tree$node.label} is not \code{NULL}, a \code{focal}
#'   clade can be indicated as "Clade NAMEOFTHECLADE" when appropriate.
#'   Similarly, an entire genus can be indicated as "Genus NAMEOFTHEGENUS" (see
#'   examples below).
#' @param sister the sister tip/node where the \code{focal} must be attached. It
#'   can be tip name/number or node number. If \code{tree$node.label} is not
#'   \code{NULL}, a \code{focal} clade can be indicated as "Clade
#'   NAMEOFTHECLADE" when appropriate. Similarly, an entire genus can be
#'   indicated as "Genus NAMEOFTHEGENUS" (see examples below).
#' @param poly logical indicating whether the \code{focal} and the \code{sister}
#'   should form a polytomous clade.
#' @param rescale logical. If the most recent common ancestor of the
#'   \code{focal} clade is older than its new ancestor (i.e. the node right
#'   above \code{sister}), the user can choose whether the height of the
#'   \code{focal} clade must be rescaled on the height of the new ancestor
#'   (\code{rescale=TRUE}), or the topology of the tree must be modified to
#'   accommodate the height of \code{focal} as it is (rescale=FALSE, in this
#'   case \code{\link{scaleTree}} is applied). This is ignored under \code{poly
#'   = TRUE}.
#' @param rootage the age of the tree root to be supplied if \code{focal} must
#'   be attached to it (and \code{poly=FALSE}). If \code{rootage=NULL} the total
#'   height of the tree increases by 10\%.
#' @return The phylogenetic tree with required topological changes.
#' @export
#' @author  Silvia Castiglione, Pasquale Raia
#' @examples
#' require(phytools)
#' DataCetaceans$tree->treecet
#'
#' ### Case 1. Moving a single tip
#' # sister to a tip
#' move.lineage(treecet,focal="Orcinus_orca",sister="Balaenoptera_musculus")->mol1
#' # sister to a clade
#' move.lineage(treecet,focal="Orcinus_orca",sister=131)->mol2
#' # sister to a clade by using treecet$node.label
#' move.lineage(treecet,focal="Balaenoptera_musculus",sister="Clade Delphinida")->mol3
#' # sister to a specific genus
#' move.lineage(treecet,focal="Orcinus_orca",sister="Genus Balaenoptera")->mol4
#' # sister to the tree root with and without rootage
#' move.lineage(treecet,focal="Balaenoptera_musculus",sister=117)->mol5
#' move.lineage(treecet,focal="Balaenoptera_musculus",sister=117,rootage=max(diag(vcv(treecet))))->mol6
#'
#' ### Case 2. Moving a clade
#' # sister to a tip
#' move.lineage(treecet,focal="Genus Mesoplodon",sister="Balaenoptera_musculus")->mol7
#' move.lineage(treecet,focal="Clade Delphinida",sister="Balaenoptera_musculus")->mol8
#' move.lineage(treecet,focal=159,sister="Balaenoptera_musculus")->mol9
#' # sister to a clade
#' move.lineage(treecet,focal="Genus Mesoplodon",sister=131)->mol10
#' move.lineage(treecet,focal="Clade Delphinida",sister=131)->mol11
#' move.lineage(treecet,focal=159,sister=131)->mol12
#' # sister to a clade by using treecet$node.label
#' move.lineage(treecet,focal="Genus Mesoplodon",sister="Clade Plicogulae")->mol13
#' move.lineage(treecet,focal="Clade Delphinida",sister="Clade Plicogulae")->mol14
#' move.lineage(treecet,focal=159,sister="Clade Plicogulae")->mol15
#' # sister to a specific genus
#' move.lineage(treecet,focal="Genus Mesoplodon",sister="Genus Balaenoptera")->mol16
#' move.lineage(treecet,focal="Clade Delphinida",sister="Genus Balaenoptera")->mol17
#' move.lineage(treecet,focal=159,sister="Genus Balaenoptera")->mol18
#' # sister to the tree root with and without rootage
#' move.lineage(treecet,focal="Genus Mesoplodon",sister=117)->mol19
#' move.lineage(treecet,focal="Clade Delphinida",sister=117)->mol20
#' move.lineage(treecet,focal=159,sister=117)->mol21
#' move.lineage(treecet,focal="Genus Mesoplodon",
#'              sister=117,rootage=max(diag(vcv(treecet))))->mol22
#' move.lineage(treecet,focal="Clade Delphinida",
#'              sister=117,rootage=max(diag(vcv(treecet))))->mol23
#' move.lineage(treecet,focal=159,sister=117,rootage=max(diag(vcv(treecet))))->mol24

move.lineage<-function(tree,focal,sister,poly=FALSE,rescale=TRUE,rootage=NULL){

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  rootname<-Ntip(tree)+1
  tipages<-max(diag(vcv(tree)))-diag(vcv(tree))

  if(grepl("Genus",focal)){
    trimws(gsub("Genus ","",focal))->genref
    getGenus(tree,genref)->getgenref
    if(getgenref[2]>1) focal<-getgenref[,3] else
      focal<-grep(genref,tree$tip.label,value=TRUE)
  }else if(grepl("Clade",focal)){
    if(!gsub("Clade ","",focal)%in%tree$node.label) stop("Required node.label not indicated on the tree")
    focal<-Ntip(tree)+which(tree$node.label==gsub("Clade ","",focal))
  }

  if(is.character(focal)|focal<=Ntip(tree)){
    if(is.character(focal)) ifelse(focal%in%tree$tip.label,focal<-which(tree$tip.label==focal),stop("The focal tip is missing from the tree"))
    if(is.numeric(sister)) tips(tree,sister)->sistips
    bind.tip(tree,tip.label = "sis_abcde",edge.length = 0.0000001,where=focal,position=0.0000001)->tree
    focal<-getMRCA(tree,c(tree$tip.label[focal],"sis_abcde"))
    if(is.numeric(sister)) ifelse(length(sistips)>1,sister<-getMRCA(tree,sistips),sister<-which(tree$tip.label==sistips))
  }

  if(grepl("Genus",sister)){
    trimws(gsub("Genus ","",sister))->genref
    getGenus(tree,genref)->getgenref
    if(getgenref[2]>1) sister<-getgenref[,3] else
      sister<-grep(genref,tree$tip.label,value=TRUE)
  }else if(grepl("Clade",sister)){
    if(!gsub("Clade ","",sister)%in%tree$node.label) stop("Required node.label not indicated on the tree")
    sister<-Ntip(tree)+which(tree$node.label==gsub("Clade ","",sister))
  }

  if(is.character(sister)){
    if(!sister%in%tree$tip.label) stop("The sister tip is missing from the tree")
    sister<-which(tree$tip.label==sister)
  }

  if(sister==(Ntip(tree)+1)){
    if(!poly){
      if(is.null(rootage)){
        rootageSc<-rootage<-max(diag(vcv(tree)))+max(diag(vcv(tree)))/10
        warning("Argument 'rootage' is missing, the tree root will be arbitrarily moved back in time",immediate. = TRUE)
      }else rootageSc<-rootage

      if(rootage==max(diag(vcv(tree)))) rootage<-rootage+0.0000001
      tree$root.edge<-rootage-max(diag(vcv(tree)))
      tree<-bind.tip(tree,tip.label="tip_abcde",edge.length=max(diag(vcv(tree)))+tree$root.edge,where=Ntip(tree)+1,
                     position=tree$root.edge)
      sister<-which(tree$tip.label=="tip_abcde")
      focal<-focal+ifelse(poly,1,2)
    }else rootageSc<-max(diag(vcv(tree)))
  }else{
    rootdesc<-getDescendants(tree,Ntip(tree)+1)[1:2]
    if(!focal%in%rootdesc) rootageSc<-max(diag(vcv(tree))) else
      rootageSc<-max(diag(vcv(tree)))-dist.nodes(tree)[(Ntip(tree)+1),rootdesc[which(rootdesc!=focal)]]
  }
  names(rootageSc)<-rootname

  momsist<-getMommy(tree,sister)[1]
  rootdist<-dist.nodes(tree)[Ntip(tree)+1,c(focal,sister,momsist)]
  claHei<-max(diag(vcv(tree)))-rootdist
  sis.el<-tree$edge.length[which(tree$edge[,2]==sister)]

  extract.clade(tree,focal)->bind.clade
  extages<-tipages[bind.clade$tip.label]
  tree$tip.label[which(tree$tip.label%in%bind.clade$tip.label)]<-
    paste(tree$tip.label[which(tree$tip.label%in%bind.clade$tip.label)],"abcde",sep="_")

  if(poly) {
    pos<-0
    bind.clade$root.edge<-NULL
    if((claHei[1]+rootdist[2])>max(diag(vcv(tree)))) bind.clade<-rescaleRR(bind.clade,height=claHei[2])
  }else{
    if(claHei[1]<=claHei[2]){
      pos<-sis.el/2
      bind.clade$root.edge<-pos+claHei[2]-claHei[1]
    }else if(claHei[1]<claHei[3]){
      pos<-sis.el-(claHei[3]-claHei[1])/2
      bind.clade$root.edge<-(claHei[3]-claHei[1])/2
    }else{
      if(any(extages>claHei[3])|!isTRUE(rescale)){
        bind.clade$root.edge<-min(tree$edge.length)
        pos<-sis.el-min(tree$edge.length)
      }else{
        rescaleRR(bind.clade,height=claHei[3]-min(tree$edge.length)*2)->bind.clade
        pos<-sis.el-min(tree$edge.length)
        bind.clade$root.edge<-min(tree$edge.length)
      }
    }
  }
  bind.tree(tree,bind.clade,where=sister,position=pos)->tt
  drop.tip(tt,grep("_abcde",tt$tip.label))->newtree

  scaleTree(newtree,tip.ages = tipages,node.ages = rootageSc)->newtree

  return(newtree)
}
