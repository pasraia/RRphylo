#' @title Create alternative phylogenies from a given tree
#' @usage swapONE(tree,si=0.5,si2=0.5)
#' @description The function produces an alternative phylogeny with altered topology and branch length, and computes the Kuhner-Felsenstein (Kuhner & Felsenstein 1994) distance between original and 'swapped' tree.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous.
#' @param si the proportion of tips whose topologic arrangement will be swapped.
#' @param si2 the proportion of nodes whose age will be changed.
#' @details \code{swap.one} changes the tree topology and branch lengths. Up to half of the tips, and half of the branch lengths can be changed randomly. Each randomly selected node is allowed to move up to 2 nodes apart from its original position.
#' @export
#' @importFrom stats runif
#' @importFrom phangorn KF.dist
#' @return The function returns a list containing the 'swapped' version of the original tree, and the Kuhner-Felsenstein distance between the trees.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references Kuhner, M. K. & Felsenstein, J. (1994). A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates, \emph{Molecular Biology and Evolution}, \strong{11}, 459-–468
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#'
#' swapONE(tree=treedino)


swapONE<-function(tree,
                  si=0.5,
                  si2=0.5){

  tree->tree1

  maxN <- function(x, N=2){
    len <- length(x)
    if(N>len){
      warning("N greater than length(x).  Setting N=length(x)")
      N <- length(x)
    }
    sort(x,partial=len-N+1)[len-N+1]
  }


  ### swap tips ####
  apply(vcv(tree),1,function(x) which(x==maxN(x,N=3)))->shifts
  sample(shifts,length(shifts)*si)->t.shifts
  sapply(t.shifts,function(x) if(length(x)==1) x<-x else sample(x,1))->t.change

  diag(vcv(tree))->ages
  data.frame(tree$edge[,2],tree$edge.length)->DF
  DF[which(DF[,1]<=Ntip(tree)),]->DF
  data.frame(tree$tip.label,DF,ages-DF[,2],ages)->DF
  colnames(DF)<-c("tip","Ntip","leaf","age.node","age")


  check<-array()
  for(i in 1:length(t.change)){
    if(DF[DF[,1]==strsplit(names(t.change),"\\.")[[i]][1],5]<DF[DF[,1]==strsplit(names(t.change),"\\.")[[i]][2],4] | DF[DF[,1]==strsplit(names(t.change),"\\.")[[i]][2],5]<DF[DF[,1]==strsplit(names(t.change),"\\.")[[i]][1],4]) check[i]<-"bar" else check[i]<-"good"
  }
  t.change[which(check=="bar")]<-"NULL"

  if(length(which(check=="bar"))>0) t.change[-which(t.change=="NULL")]->t.change else t.change->t.change
  if(length(t.change)<1)
  {
    warning("no swap implemented, fuzzy node aging only will be performed")
    tree1->tree1
  } else {
    for(i in 1:length(t.change)){
      tree1$tip.label[replace(seq(1:Ntip(tree1)),c(which(tree1$tip.label==strsplit(names(t.change),split="\\.")[[i]][1]),which(tree1$tip.label==strsplit(names(t.change),split="\\.")[[i]][2])),c(which(tree1$tip.label==strsplit(names(t.change),split="\\.")[[i]][2]),which(tree1$tip.label==strsplit(names(t.change),split="\\.")[[i]][1])))]->tree1$tip.label
    }


    diag(vcv(tree))-ages[match(tree$tip.label,names(ages))]->corr
    data.frame(tree$edge[,2],tree$edge.length)->DF1
    DF1[which(DF1[,1]<=Ntip(tree)),]->DF1

    tree1$edge.length[as.numeric(rownames(DF1))]<-(DF1[,2]-corr)
  }

  ### change node ages ######
  data.frame(tree1$edge,nodeHeights(tree1),tree1$edge.length)->nodedge
  sample(nodedge[nodedge[,2]>Ntip(tree1)+1,2],Nnode(tree1)*si2)->N

  for(i in 1:length(N)){
    runif(1,nodedge[nodedge[,2]==N[i],3],min(nodedge[nodedge[,1]==N[i],4]))->new.pos
    nodedge[nodedge[,2]==N[i],4]-new.pos->Xcorr
    nodedge[nodedge[,2]==N[i],4]<-new.pos
    nodedge[nodedge[,2]==N[i],5]-Xcorr-> nodedge[nodedge[,2]==N[i],5]
    nodedge[nodedge[,1]==N[i],5]+Xcorr->nodedge[nodedge[,1]==N[i],5]
    nodedge[nodedge[,1]==N[i],3]-Xcorr->nodedge[nodedge[,1]==N[i],3]
  }
  nodedge[,5]->tree1$edge.length
  KF.dist(tree,tree1)->KF
  return(list("modified tree"=tree1,"Kuhner-Felsenstein distance"=KF))
}
