#' @title Resolving polytomies to non-zero length branches
#' @description The function either collapses clades under a polytomy or
#'   resolves polytomous clades to non-zero length branches, dichotomous clades.
#' @usage fix.poly(tree,type=c("collapse","resolve"),node=NULL)
#' @param tree a phylogenetic tree.
#' @param type either 'collapse' to create or 'resolve' to resolve (fix) a
#'   polytomy to a specific node indicated by the argument \code{node}.
#' @param node the node in the tree where a polytomy should be resolved of
#'   fixed, either. If \code{type='resolve'} and \code{node} is left unspecified
#'   all the polytomies present in the tree are resolved.
#' @return A phylogenetic tree with randomly fixed (i.e. \code{type='resolve'})
#'   polytomies or created polytomies (i.e. \code{type='collapse'}).
#' @importFrom ape di2multi
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio
#' @details Under \code{type='resolve'} polytomous clades are resolved adding
#'   non-zero length branches to each new node. The evolutionary time attached
#'   to the new nodes is partitioned equally below the dichotomized clade.
#' @export
#' @seealso \href{../doc/Tree-Manipulation.html#fix.poly}{\code{fix.poly} vignette};
#' @references Castiglione, S., Serio, C., Piccolo, M., Mondanaro, A.,
#'   Melchionna, M., Di Febbraro, M., Sansalone, G., Wroe, S., & Raia, P.
#'   (2020). The influence of domestication, insularity and sociality on the
#'   tempo and mode of brain size evolution in mammals. \emph{Biological Journal
#'   of the Linnean Society},in press. doi:10.1093/biolinnean/blaa186
#' @examples
#' \dontrun{
#'  require(ape)
#'
#'  data("DataCetaceans")
#'  DataCetaceans$treecet->treecet
#'
#'  # Resolve all the polytomies within Cetaceans phylogeny
#'  fix.poly(treecet,type="resolve")->treecet.fixed
#'  par(mfrow=c(1,2))
#'  plot(treecet,no.margin=TRUE,show.tip.label=FALSE)
#'  plot(treecet.fixed,no.margin=TRUE,show.tip.label=FALSE)
#'
#'  # Resolve the polytomies pertaining the genus Kentriodon
#'  fix.poly(treecet,type="resolve",node=221)->treecet.fixed2
#'  par(mfrow=c(1,2))
#'  plot(treecet,no.margin=TRUE,show.tip.label=FALSE)
#'  plot(treecet.fixed2,no.margin=TRUE,show.tip.label=FALSE)
#'
#'  # Collapse Delphinidae into a polytomous clade
#'  fix.poly(treecet,type="collapse",node=179)->treecet.collapsed
#'  par(mfrow=c(1,2))
#'  plot(treecet,no.margin=TRUE,show.tip.label=FALSE)
#'  plot(treecet.collapsed,no.margin=TRUE,show.tip.label=FALSE)
#' }

fix.poly<-function(tree,type=c("collapse","resolve"),node=NULL){
  # require(ape)
  # require(phytools)
  # require(geiger)

  tree->treeO

  if(is.null(node)==FALSE & length(node)>1){
    sapply(node,function(x) tips(tree,x),simplify=FALSE)->xx
    sapply(node,function(x) length(tips(tree,x)))->lls
    names(xx)<-node
    combn(names(xx),2)->cb
    for(y in 1:ncol(cb)){
      if(lls[match(cb[1,y],names(lls))]>=lls[match(cb[2,y],names(lls))]) cb[,y]->cb[,y] else rev(cb[,y])->cb[,y]
    }

    outs<-array()
    for(y in 1:ncol(cb)){
      if(length(which(xx[match(cb[2,y],names(xx))][[1]]%in%xx[match(cb[1,y],names(xx))][[1]]))==length(xx[match(cb[2,y],names(xx))][[1]])){
        warning(paste("node", names(xx[match(cb[2,y],names(xx))]), "is nested within", names(xx[match(cb[1,y],names(xx))]), "and will be removed"))
        names(xx[match(cb[2,y],names(xx))])->outs[y]
      } else {
        outs[y]<-NA
      }
    }
    outs[!is.na(outs)]->outs
    if(length(outs)==0) node->node else node[-which(node==unique(outs))]->node
  }

  if(type=="collapse"){
    if(is.null(node)) stop("node must be supplied for type 'collapse' ")
    diag(vcv(tree))->las
    mean(tree$edge.length)/10e6->xme
    for(x in 1:length(node)) tree$edge.length[match(getDescendants(tree,node[x]),tree$edge[,2])]<-xme
    di2multi(tree,xme*1.000001)->treeN
    unique(unlist(sapply(node, function(x) tips(treeO,x))))->focal
    max(diag(vcv(treeN)))-las[focal]->fixer
    suppressWarnings(scaleTree(treeN,tip.ages = fixer)->treeN)
    return(treeN)



  }else{
    if(is.binary.tree(tree))
    {
      if(min(tree$edge.length)==0) di2multi(tree,tol=1e-06)->tree else stop("binary tree provided, no polytomies to resolve")
    }

    tree->treeN
    max(diag(vcv(treeN)))-diag(vcv(treeN))->f2
    if(is.null(node)){
      table(treeN$edge[,1])->tt
      if(any(tt>2)) names(which(tt>2))->nn
    } else {node->nn}


    which(sapply(nn, function(x) length(which(getDescendants(treeN,x)>Ntip(treeN))))>=2)->checker
    if(length(which(checker==TRUE)>0)) nn[-checker]->nn
    if(length(nn)==0) stop("non polytomous clade selected to resolve")

    sapply(nn, function(y) any(getMommy(treeN,as.numeric(y))%in%nn))->nest
    if(length(which(nest==TRUE)>0)) nn[-which(nest==TRUE)]->nn

    for(i in 1:length(nn)){
      if(i==1) treeN->xtree else xtree->xtree
      extract.clade(treeN,as.numeric(nn[i]))->tar
      multi2di(tar)->trx
      diag(vcv(trx))->dtar
      max(dtar)-dtar->fixer

      names(which(table(tar$edge[,1])==2))->fixnode
      max(dtar)-nodeHeights(tar)[match(fixnode,tar$edge[,2]),2]->nodage
      names(nodage)<-fixnode

      trx$edge.length+.01->trx$edge.length
      makeL(trx)->L->Lx
      Lx[,1]<-L[,1]<-0
      xtime<-array(); for(e in 1:Ntip(trx)) (sum(Lx[e,])-length(getMommy(trx,e))*.01)/(length(getMommy(trx,e)))->xtime[e]

      for(o in 1:Ntip(trx)){
        Lx[o,][which(Lx[o,]!=0)]<-xtime[o]
      }
      suppressWarnings(apply(Lx,2,function(x) min(x[-x!=0]))->xtar)
      xtar[-1]->xtar


      trx$edge.length->egg
      names(egg)<-trx$edge[,2]
      names(egg)[which(as.numeric(names(egg))<=Ntip(trx))]<-trx$tip.label[trx$edge[which(as.numeric(names(egg))<=Ntip(trx)),2]]
      xtar[match(names(egg),names(xtar))]->trx$edge.length
      geiger::rescale(trx,"depth",max(nodeHeights(tar)))->trx

      if(length(nodage)==0)
      {
        suppressWarnings(scaleTree(trx,tip.ages=fixer)->aa)
      }else{
      newn<-array(); for(e in 1:length(nodage)) getMRCA(trx,tips(tar,as.numeric(names(nodage)[e])))->newn[e]
      names(nodage)<-newn
      suppressWarnings(scaleTree(trx,tip.ages=fixer,node.ages = nodage)->aa)
      }

      for(e in 1:Nnode(aa)) aa$node.label[e]<-paste("new.node",i,sep = "")
      aa$node.label[1]<-nn[i]
      match(aa$tip.label,tar$tip.label)->m
      n.lab<-array(); for(o in 1:Ntip(aa)) paste("tip",o,sep = "_")->n.lab[o]
      n.lab->aa$tip.label

      if(i==1) bind.tree(xtree,aa,as.numeric(nn[i]))->temp else bind.tree(xtree,aa,getMRCA(xtree,tips(treeN,as.numeric(nn[i]))))->temp
      tips(treeN,as.numeric(nn[i]))->nams
      drop.tip(temp,nams)->tree1
      tree1$tip.label[match(aa$tip.label,tree1$tip.label)]<-nams[m]
      geiger::rescale(tree1,"depth",max(nodeHeights(treeN)))->tree1
      suppressWarnings(scaleTree(tree1,tip.ages=f2)->tree1)
      if(is.binary.tree(tree1)==FALSE) multi2di(tree1)->xtree1 else tree1->xtree1
      xtree1->xtree
    }
    return(xtree)
  }
}
