#' @title Visualize the difference between phylogenetic trees
#' @description The function scans a pair of phylogenetic trees to find
#'   topological differences.
#' @param tree,tree1 a phylogenetic tree. The tree needs not to be ultrametric
#'   and fully dichotomous. Generic name and specific epithet must be separated
#'   by '_'.
#' @param plot if \code{TRUE}, the function produces an interactive plotting
#'   device to check differences in species placement at the genus level.
#' @export
#' @return The function returns a data-frame indicating for each un-matching
#'   species its sister species/clades on both trees.
#' @author Silvia Castiglione, Carmela Serio, Antonella Esposito
#' @examples
#' \dontrun{
#' DataSimians$tree->tree
#' set.seed(22)
#' swapONE(tree,si=0.5)[[1]]->tree1
#'
#' treeCompare(tree,tree1)
#'
#' }

treeCompare<-function(tree,tree1,plot=TRUE){
  if(isTRUE(plot)){
    mars <- par("mar")
    on.exit(par(mar = mars))
  }

  if(!identical(tree1$edge[tree1$edge[,2]<=Ntip(tree1),2],seq(1,Ntip(tree1)))){
    tree1$tip.label<-tree1$tip.label[tree1$edge[tree1$edge[,2]<=Ntip(tree1),2]]
    tree1$edge[tree1$edge[,2]<=Ntip(tree1),2]<-seq(1,Ntip(tree1))
  }

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  tree->treeO
  tree1->tree1O

  if(any(!tree$tip.label%in%tree1$tip.label)) drop.tip(tree,which(!tree$tip.label%in%tree1$tip.label))->tree
  if(any(!tree1$tip.label%in%tree$tip.label)) drop.tip(tree1,which(!tree1$tip.label%in%tree$tip.label))->tree1

  data.frame(row.names=tree$tip.label)->sispair
  sapply(tree$tip.label,function(j) getSis(tree,j,printZoom = FALSE))->sispair$sist
  sapply(sispair$sist,function(ss) {
    if(!all(ss%in%tree$tip.label)){
      c(ss[which(ss%in%tree$tip.label)],unlist(sapply(ss[which(!ss%in%tree$tip.label)],function(k) tips(tree,k))))
    }else ss
  },simplify = FALSE)->sispair$sis.list
  sapply(sispair$sist,length)->sispair$poly

  sapply(tree1$tip.label,function(j) getSis(tree1,j,printZoom = FALSE))->sist1
  sist1[match(rownames(sispair),names(sist1))]->sispair$sist1
  sapply(sispair$sist1,function(ss) {
    if(!all(ss%in%tree1$tip.label)){
      c(ss[which(ss%in%tree1$tip.label)],unlist(sapply(ss[which(!ss%in%tree1$tip.label)],function(k) tips(tree1,k))))
    }else ss
  },simplify = FALSE)->sispair$sis.list1
  sapply(sispair$sist1,length)->sispair$poly1

  sispair[which((sispair$poly>1&sispair$poly1==1)|(sispair$poly==1&sispair$poly1>1)),]->nomatch.sp
  sispair[which(!rownames(sispair)%in%rownames(nomatch.sp)),]->sispair
  rbind(nomatch.sp,sispair[which(!apply(sispair,1,function(x) all(x[[2]]%in%x[[5]])&all(x[[5]]%in%x[[2]]))),])->nomatch.sp

  nomatch.sp[,c(1,4)]->nomatch.sp

  if(isTRUE(plot)){
    sapply(unique(sapply(strsplit(rownames(nomatch.sp),"_"),"[[",1)),function(j) getGenus(treeO,j)[,3])->gentree
    sapply(unique(sapply(strsplit(rownames(nomatch.sp),"_"),"[[",1)),function(j) getGenus(tree1O,j)[,3])->gentree1

    gentree1[order(names(gentree))]->gentree1
    gentree[order(names(gentree))]->gentree

    font.tree<-rep(3,length(tree$tip.label))
    font.tree[which(tree$tip.label%in%rownames(nomatch.sp))]<-4
    names(font.tree)<-tree$tip.label

    plotClades<-function(cla.list,cla1.list,colo,colo1,font.cla,font.cla1){
      par(mfrow=c(1,2))
      plot(cla.list,tip.col = colo,no.margin=TRUE,font=font.cla)
      plot(cla1.list,tip.col = colo1,no.margin = TRUE,direction="leftward",font=font.cla1)
    }

    cla.list<-cla1.list<-colo<-colo1<-font.cla<-font.cla1<-list()
    for(w in 1:length(gentree)){
      extract.clade(treeO,gentree[w])->cla->cla.list[[w]]
      extract.clade(tree1O,gentree1[w])->cla1->cla1.list[[w]]
      match(rownames(nomatch.sp)[grep(paste(names(gentree)[w],"_",sep=""),rownames(nomatch.sp))],cla$tip.label)->spcla
      match(rownames(nomatch.sp)[grep(paste(names(gentree1)[w],"_",sep=""),rownames(nomatch.sp))],cla1$tip.label)->spcla1

      colo[[w]]<-rep("black",Ntip(cla))
      colo1[[w]]<-rep("black",Ntip(cla1))
      colo[[w]][spcla]<-colo1[[w]][spcla1]<-scales::hue_pal()(length(spcla))
      font.cla[[w]]<-font.tree[match(cla$tip.label,names(font.tree))]
      font.cla1[[w]]<-font.tree[match(cla1$tip.label,names(font.tree))]
    }
    names(cla.list)<-names(cla1.list)<-names(colo)<-names(colo1)<-names(font.cla)<-names(font.cla1)<-names(gentree)

    taxon<-manipulate::picker(as.list(names(gentree)))
    manipulate::manipulate(do.call(plotClades,list(cla.list=cla.list[[which(names(cla.list)==taxon)]],
                                                   cla1.list=cla1.list[[which(names(cla1.list)==taxon)]],
                                                   colo=colo[[which(names(colo)==taxon)]],
                                                   colo1=colo1[[which(names(colo1)==taxon)]],
                                                   font.cla=font.cla[[which(names(font.cla)==taxon)]],
                                                   font.cla1=font.cla1[[which(names(font.cla1)==taxon)]])),taxon=taxon)
  }

  colnames(nomatch.sp)<-c("tree","tree1")
  nomatch.sp[order(rownames(nomatch.sp)),]->nomatch.sp
  return(nomatch.sp)
}
