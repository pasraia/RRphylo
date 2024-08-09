#' @title Create alternative phylogenies from a given tree
#' @usage swapONE(tree,node=NULL,si=0.5,si2=0.5,plot.swap=FALSE)
#' @description The function produces an alternative phylogeny with altered
#'   topology and branch length, and computes the Kuhner-Felsenstein (Kuhner &
#'   Felsenstein 1994) distance between original and 'swapped' tree.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous.
#' @param node if specified, the clades subtended by such \code{node(s)} are
#'   imposed to be monophyletic. In this case, the function can still swap tips
#'   \emph{within} the clade.
#' @param si the proportion of tips whose topologic arrangement will be swapped.
#' @param si2 the proportion of nodes whose age will be changed.
#' @param plot.swap if \code{TRUE}, the function plots the swapped tree. Swapped
#'   positions appear in red. Nodes with altered ages appear in green.
#' @details \code{swapONE} changes the tree topology and branch lengths. Up to
#'   half of the tips, and half of the branch lengths can be changed randomly.
#'   Each randomly selected node is allowed to move up to 2 nodes apart from its
#'   original position.
#' @export
#' @importFrom stats runif
#' @importFrom ape cophenetic.phylo
#' @return The function returns a list containing the 'swapped' version of the
#'   original tree, and the Kuhner-Felsenstein distance between the trees. Note,
#'   tip labels are ordered according to their position in the tree.
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Kuhner, M. K. & Felsenstein, J. (1994). A simulation comparison
#'   of phylogeny algorithms under equal and unequal evolutionary rates,
#'   \emph{Molecular Biology and Evolution}, 11: 459-468.
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#'
#' ## Case 1. change the topology and the branch lengths for the entire tree
#' swapONE(tree=treedino,si=0.5,si2=0.5,plot.swap=FALSE)
#'
#' ## Case 2. change the topology and the branch lengths of the
#' ##         tree by keeping the monophyly of a specific clade
#' swapONE(tree=treedino,node=422,si=0.5,si2=0.5,plot.swap=FALSE)
#' }


swapONE<-function(tree,
                  node=NULL,
                  si=0.5,
                  si2=0.5,
                  plot.swap=FALSE){

  #require(phangorn)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  tree1<-tree
  maxN <- function(x, N=2){
    len <- length(x)
    if(N>len){
      warning("N greater than length(x).  Setting N=length(x)")
      N <- length(x)
    }
    sort(x,partial=len-N+1)[len-N+1]
  }

  ### swap tips ###
  if(si>0){
    if(Ntip(tree)*si==1) 2/Ntip(tree)->si
    apply(vcv(tree),1,function(x) which(x==maxN(x,N=3)))->shifts
    if(!is.null(node)){
      for(i in 1:length(node)){
        tips(tree,node[i])->nod.tip
        lapply(shifts[which(names(shifts)%in%nod.tip)],function(x) which(!names(x)%in%nod.tip))->rem.tip
        which(sapply(rem.tip,length)>0)->rr
        if(any(rr)) for(j in 1:length(rr)){
          shifts[which(names(shifts)%in%nod.tip)][[rr[j]]][-rem.tip[[rr[j]]]]->shifts[which(names(shifts)%in%nod.tip)][[rr[j]]]
        }
        lapply(shifts[which(!names(shifts)%in%nod.tip)],function(x) which(names(x)%in%nod.tip))->rem.tip1
        which(sapply(rem.tip1,length)>0)->rr1
        if(any(rr1)) for(k in 1:length(rr1)){
          shifts[which(!names(shifts)%in%nod.tip)][[rr1[k]]][-rem.tip1[[rr1[k]]]]->shifts[which(!names(shifts)%in%nod.tip)][[rr1[k]]]
        }
      }
    }

    cophenetic.phylo(tree)->cop
    if(any(sapply(shifts,length)==0)) {
      cop[-which(sapply(shifts,length)==0),]->cop
      shifts[-which(sapply(shifts,length)==0)]->shifts
    }

    lapply(1:length(shifts),function(x) {
      cop[x,which(names(cop[x,])%in%names(shifts[[x]]))]->shift.dist
      names(shift.dist)<-colnames(cop)[which(names(cop[x,])%in%names(shifts[[x]]))]
      return(shift.dist)
    })->patr.dist
    names(patr.dist)<-names(shifts)
    sd(sapply(patr.dist,mean))*2->lim

    for(x in 1:length(patr.dist)){
      dN<-c()
      getSis(tree,names(shifts)[x],printZoom = FALSE)->sis
      suppressWarnings(as.numeric(sis)->nsis)
      if(any(is.na(nsis))) which(tree$tip.label%in%sis[which(is.na(nsis))])->sis[which(is.na(nsis))]
      as.numeric(sis)->sis
      if(length(sis)<2){
        if(sis<=(Ntip(tree))) c(sis,dN)->dN else {
          tree$edge[tree$edge[,1]==sis,2]->sis2
          if(any(sis2<=Ntip(tree))) c(dN,sis2[which(sis2<=Ntip(tree))])->dN
        }

      }else{
        for(y in sis){
          if(y<=(Ntip(tree))) c(y,dN)->dN else {
            tree$edge[tree$edge[,1]==y,2]->sis2
            if(any(sis2<=Ntip(tree))) c(dN,sis2[which(sis2<=Ntip(tree))])->dN
          }
        }
      }

      getMommy(tree,which(tree$tip.label==names(shifts)[x]))[1]->mom
      getSis(tree,mom,printZoom = FALSE)->sismom
      suppressWarnings(as.numeric(sismom)->nsismom)
      if(any(is.na(nsismom))) c(dN,which(tree$tip.label%in%sismom[which(is.na(nsismom))]))->dN
      tree$tip.label[dN]->dN

      shifts[[x]][unique(c(match(dN,names(shifts[[x]]),nomatch=0),which(patr.dist[[x]]<lim)))]->shifts[[x]]
    }
    names(shifts)<-names(patr.dist)

    if(any(sapply(shifts,length)==0)) shifts[-which(sapply(shifts,length)==0)]->shifts

    if((Ntip(tree)*si)>length(shifts)) shifts->t.shifts else sample(shifts,Ntip(tree)*si)->t.shifts
    lapply(t.shifts,function(x) if(length(x)==1) x<-x else sample(x,1))->t.change
    sapply(1:length(t.change), function(k) paste(names(t.change)[k],"><",names(t.change[[k]]),sep=""))->names(t.change)
    unlist(lapply(t.change,unname))->t.change

    diag(vcv(tree))->ages
    data.frame(tree$edge[,2],tree$edge.length)->DF
    DF[which(DF[,1]<=Ntip(tree)),]->DF
    data.frame(tree$tip.label,DF,ages-DF[,2],ages)->DF
    colnames(DF)<-c("tip","Ntip","leaf","age.node","age")

    check<-array()
    for(i in 1:length(t.change)){
      if(DF[DF[,1]==strsplit(names(t.change),"><")[[i]][1],5]<DF[DF[,1]==strsplit(names(t.change),"><")[[i]][2],4]|
         DF[DF[,1]==strsplit(names(t.change),"><")[[i]][2],5]<DF[DF[,1]==strsplit(names(t.change),"><")[[i]][1],4]) check[i]<-"bar" else check[i]<-"good"
    }
    if(length(which(check=="bar"))>0) t.change[-which(check=="bar")]->t.change

    if(length(t.change)<1) {
      warning("no swap implemented, fuzzy node aging only will be performed")
      tree1->tree1
    } else {
      repeat({
        tree->tree1
        sw.tips<-c()
        for(i in 1:length(t.change)){
          c(sw.tips,c(which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][1]),
                      which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][2])))->sw.tips
          tree1$tip.label[replace(seq(1:Ntip(tree1)),
                                  c(which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][1]),
                                    which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][2])),
                                  c(which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][2]),
                                    which(tree1$tip.label==strsplit(names(t.change),split="><")[[i]][1])))]->tree1$tip.label
        }


        data.frame(tree1$edge[,2],tree1$edge.length)->DF1
        DF1[which(DF1[,1]<=Ntip(tree1)),]->DF1
        data.frame(tree1$tip.label[DF1[,1]],DF1,diag(vcv(tree1))[DF1[,1]],ages[match(tree1$tip.label[DF1[,1]],names(ages))])->DF1
        colnames(DF1)<-c("tip","Ntip","leaf","age","age.real")
        DF1$age-DF1$age.real->DF1$corr
        as.character(DF1[which((DF1$leaf-DF1$corr)<0),1])->no.change

        if(length(no.change)>0) t.change[-unlist(lapply(no.change, function(k) grep(k,names(t.change))))]->t.change else break
      })

      data.frame(tree1$edge[,2],tree1$edge.length)->DF1
      DF1[which(DF1[,1]<=Ntip(tree1)),]->DF1
      data.frame(tree1$tip.label[DF1[,1]],DF1,diag(vcv(tree1))[DF1[,1]],ages[match(tree1$tip.label[DF1[,1]],names(ages))])->DF1
      colnames(DF1)<-c("tip","Ntip","leaf","age","age.real")
      DF1$age-DF1$age.real->DF1$corr
      DF1$leaf-DF1$corr->DF1$new.leaf

      tree1$edge.length[match(DF1[,2],tree1$edge[,2])]<-DF1$new.leaf
      # tree1$edge.length[as.numeric(rownames(DF1))]<-(DF1[,2]-corr)
      # scaleTree(tree1,max(ages)-ages[match(tree1$tip.label,names(ages))])->tree1
    }
  }else sw.tips<-NULL
  ### change node ages ######
  if(si2>0){
    if(Ntip(tree)*si2==1) 2/Ntip(tree)->si2
    data.frame(tree1$edge,nodeHeights(tree1),tree1$edge.length)->nodedge
    sample(nodedge[nodedge[,2]>Ntip(tree1)+1,2],(Nnode(tree1)-1)*si2)->N

    for(i in 1:length(N)){
      runif(1,nodedge[nodedge[,2]==N[i],3],min(nodedge[nodedge[,1]==N[i],4]))->new.pos
      nodedge[nodedge[,2]==N[i],4]-new.pos->Xcorr
      nodedge[nodedge[,2]==N[i],4]<-new.pos
      nodedge[nodedge[,2]==N[i],5]-Xcorr-> nodedge[nodedge[,2]==N[i],5]
      nodedge[nodedge[,1]==N[i],5]+Xcorr->nodedge[nodedge[,1]==N[i],5]
      nodedge[nodedge[,1]==N[i],3]-Xcorr->nodedge[nodedge[,1]==N[i],3]
    }
    nodedge[,5]->tree1$edge.length
  }else N<-NULL
  phangorn::KF.dist(tree,tree1)->KF

  if(isTRUE(plot.swap)){
    colo<-rep("black",nrow(tree$edge))
    if(length(sw.tips)>0) colo[which(tree$edge[,2]%in%unique(sw.tips))]<-"red"
    plot(tree1,edge.color = colo,cex=.5,no.margin=TRUE)
    if(length(N)>0) nodelabels(bg="w",frame="n",col="green",node=N)
  }

  return(list("modified tree"=tree1,"Kuhner-Felsenstein distance"=KF))
}
