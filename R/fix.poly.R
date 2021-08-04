#' @title Resolving polytomies to non-zero length branches
#' @description The function either collapses clades under a polytomy or
#'   resolves polytomous clades to non-zero length branches, dichotomous clades.
#' @usage
#' fix.poly(tree,type=c("collapse","resolve"),node=NULL,tol=1e-10,random=TRUE)
#' @param tree a phylogenetic tree.
#' @param type either 'collapse' to create polytomies to one or more specific
#'   nodes or 'resolve' to resolve (fix) all the polytomies within the tree or
#'   to one or more specific nodes.
#' @param node the node in the tree where a polytomy should be created or fixed,
#'   either. If \code{type='resolve'} and \code{node=NULL} all the polytomies
#'   present in the tree are resolved.
#' @param tol the tolerance to consider a branch length significantly greater
#'   than zero, set at 1e-10 by default. If \code{type='resolve'}, all the
#'   branch lengths smaller than \code{tol} are treated as polytomies.
#' @param random a logical value specifying whether to resolve the polytomies
#'   randomly (the default) or in the order they appear in the tree (if
#'   \code{random = FALSE}).
#' @return A phylogenetic tree with randomly fixed (i.e. \code{type='resolve'})
#'   polytomies or created polytomies (i.e. \code{type='collapse'}).
#' @importFrom ape di2multi
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @details Under \code{type='resolve'} polytomous clades are resolved adding
#'   non-zero length branches to each new node. The evolutionary time attached
#'   to the new nodes is partitioned equally below the dichotomized clade.
#' @export
#' @seealso \href{../doc/Tree-Manipulation.html#fix.poly}{\code{fix.poly}
#'   vignette};
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
fix.poly<-function(tree,type=c("collapse","resolve"),node=NULL,tol=1e-10,random=TRUE){
  # require(ape)
  # require(phytools)
  # require(geiger)

  tree->treeN

  if(type=="collapse"){
    if(is.null(node)) stop("node must be supplied for type 'collapse' ")
    diag(vcv(treeN))->las
    min(treeN$edge.length)/10->xme
    for(x in 1:length(node)) treeN$edge.length[match(getDescendants(treeN,node[x]),treeN$edge[,2])]<-xme
    di2multi(treeN,xme*1.000001)->treeN
    unique(unlist(sapply(node, function(x) tips(tree,x))))->focal
    max(diag(vcv(treeN)))-las[focal]->fixer
    suppressWarnings(scaleTree(treeN,tip.ages = fixer)->xtree)
  }else{
    max(diag(vcv(treeN)))-diag(vcv(treeN))->f2

    if(is.null(node)){
      if(any(treeN$edge.length<=tol)) di2multi(treeN,tol=tol)->treeN
      if(is.binary(treeN)) stop("binary tree provided, no polytomies to resolve")
      table(treeN$edge[,1])->tt
      if(any(tt>2)) names(which(tt>2))->nn
      max(diag(vcv(treeN)))-dist.nodes(treeN)[(Ntip(treeN)+1),(Ntip(treeN)+1):(Ntip(treeN)+Nnode(treeN))]->nodages
    } else {
      unlist(lapply(node,function(x){
        extract.clade(treeN,as.numeric(x))->cla
        if(any(cla$edge.length<=tol)) di2multi(cla,tol=tol)->cla
        table(cla$edge[,1])->tt
        if(any(tt>2)){
          names(which(tt>2))->nncla
          unname(sapply(nncla,function(w) getMRCA(treeN,tips(cla,as.numeric(w)))))
        } else NULL
      }))->nn

      if(is.null(nn)) stop("no polytomous clade provided") else{
        if(any(!node%in%nn)) sapply(node[which(!node%in%nn)],function(q){
          getDescendants(treeN,q)->desnode
          if(!any(desnode%in%nn)) warning("node ",q," does not subtend to a polytomous clade")
        })
      }

      max(diag(vcv(treeN)))-dist.nodes(treeN)[(Ntip(treeN)+1),(Ntip(treeN)+1):(Ntip(treeN)+Nnode(treeN))]->nodages
      unlist(unname(sapply(node,function(w) getDescendants(treeN,w)[which(getDescendants(treeN,w)>Ntip(treeN))])))->nod.des
      nodages[which(!nodages%in%nod.des[which(!nod.des%in%nn)])]->nodages
    }

    sapply(nn,function(k) length(tips(treeN,k)))->ntips
    nn[order(ntips)]->nn

    for(i in 1:length(nn)){
      if(i==1) treeN->xtree
      extract.clade(treeN,as.numeric(nn[i]))->tar
      extract.clade(xtree,getMRCA(xtree,tips(treeN,as.numeric(nn[i]))))->trx1
      if(!is.null(node)&any(trx1$edge.length<=tol)) di2multi(trx1,tol=tol)->trx1
      if(all(table(trx1$edge[,1])<3)) next else{
        multi2di(trx1,random=random)->trx

        trx$edge.length+min(trx$edge.length[which(trx$edge.length!=0)])->trx$edge.length
        makeL(trx)->Lx
        Lx[,1]<-0

        for(e in 1:Ntip(trx))
          (sum(Lx[e,])-length(getMommy(trx,e))*min(trx$edge.length[which(trx$edge.length!=0)]))/(length(getMommy(trx,e)))->Lx[e,][which(Lx[e,]!=0)]

        suppressWarnings(apply(Lx,2,function(x) min(x[-x!=0]))->xtar)
        xtar[-1]->xtar

        egg<-trx$edge[,2]
        egg[which(egg<=Ntip(trx))]<-trx$tip.label[egg[which(egg<=Ntip(trx))]]
        xtar[match(egg,names(xtar))]->trx$edge.length
        root.trx<-max(diag(vcv(tar)))
        names(root.trx)<-Ntip(trx)+1

        for(e in 1:Nnode(trx)) trx$node.label[e]<-paste("new.node",i,sep = "")
        trx$node.label[1]<-nn[i]
        unname(sapply(trx$tip.label, function(o) paste("RRphylofixed",o,sep = "_")))->trx$tip.label

        bind.tree(xtree,trx,getMRCA(xtree,tips(treeN,as.numeric(nn[i]))))->temp
        drop.tip(temp,tips(treeN,as.numeric(nn[i])))->tree1
        gsub("RRphylofixed_","",tree1$tip.label)->tree1$tip.label
        tree1->xtree
      }
    }

    unname(sapply(names(nodages), function(e) getMRCA(xtree,tips(treeN,as.numeric(e)))))->newn
    names(nodages)<-newn

    nodages[order(sapply(names(nodages),function(w) length(tips(xtree,as.numeric(w)))))]->nodages

    fix.val<-1e-19
    while(max(nodages)==max(nodages)+fix.val) fix.val<-fix.val*10

    for(i in 1:length(nodages)){
      getMommy(xtree,names(nodages[i]))->moms
      moms[which(moms%in%names(nodages))][1]->mom
      if(mom%in%names(nodages)&&nodages[i]>=nodages[which(names(nodages)==mom)])
        nodages[which(names(nodages)==mom)]<-nodages[i]+fix.val
    }

    suppressWarnings(scaleTree(xtree,tip.ages=f2,node.ages=nodages,min.branch=min(xtree$edge.length)/max(apply(makeL(xtree),1,function(x) length(which(x!=0)))))->xtree)
  }
  return(xtree)
}
