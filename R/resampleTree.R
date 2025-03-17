#' @title Altering phylogenetic trees
#' @description The function alters the topology and randomly removes a
#'   user-specified proportion of species from a phylogenetic tree.
#' @usage resampleTree(tree,s=0.25,sdata=NULL,nodes=NULL,categories=NULL,
#'              swap.si=0.1,swap.si2=0.1,swap.node=NULL,nsim=1)
#' @param tree	a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous.
#' @param s the percentage of tips to be cut off. It is set at 25\% by default.
#' @param sdata to be supplied to condition the species sampling. It can be
#'   either a named vector or a data.frame/matrix having the species names as
#'   first column. In case of stratified random sampling, \code{sdata} should
#'   contain the strata. Otherwise, the user can provide a sampling probability
#'   (meant as the probability to be removed from the tree) for each species.
#' @param nodes the clades to be preserved. In this case the function maintains
#'   no less than 5 species at least in each of them.
#' @param categories the categories to be preserved. In this case the function
#'   maintains no less than 5 species at least in each of them.
#' @param nsim number of phylogenies to return. It is set at 1 by default.
#' @param swap.si,swap.si2,swap.node arguments \code{si, si2, node} as passed to
#'   \code{\link{swapONE}}. The default for both \code{si} and \code{si2} is
#'   0.1.
#' @return The function returns \code{phylo} or \code{multiPhylo} object.
#' The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Giorgia Girardi
#' @export
#' @seealso \href{../doc/search.conv.html}{\code{search.conv} vignette};
#' \href{../doc/overfitRR.html}{\code{overfitRR} vignette};
#' \href{../doc/Alternative-trees.html}{\code{Alternative-trees} vignette}
#' @examples
#' \dontrun{
#' DataCetaceans$treecet->treecet
#' plot(treecet,show.tip.label = FALSE,no.margin = TRUE)
#' nodelabels(frame="n",col="red")
#'
#' # Select two clades for stratified random sampling
#' clanods=c("crown_Odo"=150,"crown_Mysti"=131)
#' sdata1<-do.call(rbind,lapply(1:length(clanods),function(w)
#'   data.frame(species=tips(treecet,clanods[w]),group=names(clanods)[w])))
#'
#' # generate a vector of probabilities based on body mass
#' prdata<-max(DataCetaceans$masscet)-DataCetaceans$masscet
#'
#' # select two nodes to be preserved
#' nn=c(180,159)
#'
#' # generate two fictional categorical vectors to be preserved
#' cat1<-sample(rep(c("a","b","c"),each=39),Ntip(treecet))
#' names(cat1)<-treecet$tip.label
#' cat2<-rep(c("d","e"),each=100)
#' names(cat2)<-sample(treecet$tip.label,100)
#'
#' # 1. Random sampling
#' resampleTree(treecet,s=0.25,swap.si=0.3)->treecet1
#'
#' # 1.1 Random sampling preserving clades
#' resampleTree(treecet,s=0.25,nodes=nn)->treecet2
#'
#' # 2. Stratified random sampling
#' resampleTree(treecet,sdata = sdata1,s=0.25)->treecet3
#'
#' # 2.1 Stratified random sampling preserving clades and categories
#' resampleTree(treecet,sdata = sdata1,s=0.25,nodes=nn,categories = list(cat1,cat2))->treecet4
#'
#' # 3. Sampling conditioned on probability
#' resampleTree(treecet,sdata = prdata,s=0.25,nsim=5)->treecet5
#' }




resampleTree<-function(tree,s=0.25,sdata=NULL,nodes=NULL,categories=NULL,
                       swap.si=0.1,swap.si2=0.1,swap.node=NULL,nsim=1){
  # require(ape)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  funcall <- match.call()
  if(!is.null(sdata)){
    if(is.null(ncol(sdata))){
      sdata<-data.frame(species=names(sdata),sdata)
      ifelse(is.numeric(sdata[,2]),"probs","group")->colnames(sdata)[2]
    }else{
      if(!all(colnames(sdata)%in%c("species","group","probs"))){
        warning("Colnames not matching: columns assumed to be ordered as 'species','group/probs' \n",immediate. = TRUE)
        colnames(sdata)[1]<-"species"
        ifelse(is.numeric(sdata[,2]),"probs","group")->colnames(sdata)[2]
      }
    }
    # if(is.null(ncol(sdata))) sdata<-data.frame(species=names(sdata),group=sdata,probs=NA) else{
    #   if(ncol(sdata)==2) sdata<-data.frame(sdata,probs=NA)
    # }
    # if(!all(colnames(sdata)%in%c("species","group","probs"))) {
    #   warning("Colnames not matching: columns assumed to be ordered as 'species','group','probs\n'",immediate. = TRUE)
    #   colnames(sdata)[1:2]<-c("species","group")
    #   if(ncol(sdata)==3) colnames(sdata)[3]<-"probs"
    # }
    if(any(duplicated(sdata$species))) stop("Species",paste(sdata$species[duplicated(sdata$species)],sep=" "),"in sdata are duplicated")
    if(!all(tree$tip.label%in%sdata$species)){
      if(!is.null(sdata$probs))
        stop("Please provide sampling probabilities for all the species on the tree")
      sdata<-rbind(sdata,data.frame(species=tree$tip.label[which(!tree$tip.label%in%sdata$species)],group="others"))
    }

    # if(all(is.na(sdata$probs))){
      # table(sdata$group)/Ntip(tree)->grprob
      # # if(any(sdata$group=="others")) grprob["others"]<-0.5
      # sdata$probs<-1-grprob[match(sdata$group,names(grprob))]
      # if(any(grprob<=0.05)){
      #   sapply(names(grprob)[which(grprob<=0.05)],function(j) sdata$probs[which(sdata$group%in%j)][sample(1:sum(sdata$group%in%j),2)]<<-0)
      # }
    # }
  }else data.frame(species=tree$tip.label,probs=1)->sdata
  tree$node.label<-(Ntip(tree)+1:Nnode(tree))

  tree.list<-list()
  for(k in 1:nsim){
    if(!is.null(nodes)){
      unlist(lapply(nodes,function(x){
        length(tips(tree,x))->lenx
        if(lenx<=5) tips(tree,x) else sample(tips(tree,x),5)
      }))->out.nodes
    } else out.nodes<-NULL

    if(!is.null(categories)){
      # if(is.list(categories)) do.call(cbind,categories)->categories
      # categories<-as.matrix(categories)
      # categories <- treedataMatch(tree, categories)[[1]]
      #
      # out.cat<-unlist(lapply(1:ncol(categories),function(w){
      #   table(categories[,w])->tab.ss
      #   unlist(lapply(1:length(tab.ss),function(x){
      #     if(tab.ss[x]<=5) names(which(categories[,w]==names(tab.ss)[x])) else
      #       sample(names(which(categories[,w]==names(tab.ss)[x])),5)
      #   }))
      # }))

      if(!is.list(categories)){
        categories<-as.matrix(categories)
        asplit(categories,2)->categories
      }
      categories<-lapply(categories,function(q) treedataMatch(tree, q)[[1]])

      out.cat<-unlist(lapply(categories,function(w){
        table(w)->tab.ss
        unlist(lapply(1:length(tab.ss),function(x){
          if(tab.ss[x]<=5) rownames(w)[which(w==names(tab.ss)[x])] else
            sample(rownames(w)[which(w==names(tab.ss)[x])],5)
        }))
      }))
    } else out.cat<-NULL

    unique(c(out.nodes,out.cat))->outs
    sx<-s
    repeat(if((Ntip(tree)-length(outs))>Ntip(tree)*sx) break else s*0.9->sx)

    if(is.null(sdata$probs)){
      grsam<-round(table(sdata$group)*sx)
      if(any(grsam==table(sdata$group))) grsam[which(grsam==table(sdata$group))]<-grsam[which(grsam==table(sdata$group))]-1

      datagr<-split(sdata,sdata$group)
      offs<-unlist(mapply(a=datagr,b=grsam,function(a,b){
        samrow<-seq(1,nrow(a))[which(!a$species%in%outs)]
        if(length(samrow)>1) a[sample(samrow,b),1] else a[samrow,1]
      },SIMPLIFY = FALSE))
    }else{
      samdata<-sdata
      if(length(outs)>0) samdata$probs[which(samdata$species%in%outs)]<-0
      samtips<-sum(samdata$probs>0)
      sample(samdata$species,round(Ntip(tree)*sx),prob=samdata$probs)->offs
    }

    if(is.null(swap.node)) swap.node<-NULL

    suppressWarnings(swapONE(tree,si=swap.si,si2=swap.si2,node=swap.node,plot.swap=FALSE)[[1]])->tree.swap

    tree.swap$edge[tree.swap$edge[,1]==(Ntip(tree.swap)+1),2]->rootdesc
    if(any(rootdesc<=Ntip(tree.swap))){
      tree.swap$tip.label[rootdesc[which(rootdesc<=Ntip(tree.swap))]]->saver
      offs[which(!offs%in%saver)]->offs
    }
    drop.tip(tree.swap,offs)->tree.list[[k]]
  }

  if(length(tree.list)>1) class(tree.list)<-"multiPhylo" else tree.list[[1]]->tree.list
  attr(tree.list,"Call")<-funcall

  return(tree.list)
}
