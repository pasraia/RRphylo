#' @title Phylogenetic tree calibration
#' @description The function is a wrapper around the functions "scalePhylo", "assign.ages", and "assign.brlen" written by Gene Hunt (http://paleobiology.si.edu/staff/individuals/hunt.cfm) . It rescales tree branch lengths according to given calibration dates.
#' @usage scaleTree(tree, tip.ages, node.ages=NULL, min.branch=0.1)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and fully dichotomous.
#' @param tip.ages a named vector of tip ages in terms of distance from the youngest tip.
#' @param node.ages a named vector of node ages in terms of distance from the youngest tip. If no calibration date for nodes is available except for the tree root, the vector must contain as only entry the tree height with the root number as name.
#' @param min.branch the minimum branch length that will be imposed when no age is specified.
#' @export
#' @return Rescaled phylogentic tree.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#' \donttest{
#' library(ape)
#' library(phytools)
#' library(geiger)
#' set.seed(76)
#' rtree(100)->tree
#'
#' max(nodeHeights(tree))->H
#'
#' #### Example 1 ####
#' H-diag(vcv(tree))->tipAges
#' tipAges[tips(tree,160)]<-min(tipAges)
#' H->nodeAges
#' names(nodeAges)<-Ntip(tree)+1
#'
#' scaleTree(tree,tipAges,nodeAges)->treeS1
#'
#' edge.col<-rep("black",nrow(tree$edge))
#' edge.col[which(treeS1$edge[,2]%in%getDescendants(treeS1,160))]<-"red"
#'
#' layout(2:1)
#' plot(tree,edge.color = edge.col,show.tip.label=F)
#' plot(treeS1,edge.color = edge.col,show.tip.label=F)
#'
#' #### Example 2 ####
#' nodeAges<-c(H,4.5,6,6)
#' names(nodeAges)<-c(101,143,171,121)
#' tipAges<-H-diag(vcv(tree))
#' scaleTree(tree,tipAges,nodeAges)->treeS2
#'
#' edge.col<-rep("black",nrow(tree$edge))
#' edge.col[which(treeS1$edge[,2]%in%c(getDescendants(treeS1,143),
#'                                     getDescendants(treeS1,171),
#'                                     getDescendants(treeS1,121)))]<-"red"
#'
#' layout(2:1)
#' plot(tree,edge.color = edge.col,show.tip.label=F)
#' plot(treeS2,edge.color = edge.col,show.tip.label=F)
#' }

scaleTree<- function(tree, tip.ages, node.ages=NULL, min.branch=0.1)
{
  # require(ape)
  # require(phytools)

  tree->tr
  abs(min.branch)->min.diff

  node.mins<-array(dim=Nnode(tr))
  names(node.mins)<-seq((Ntip(tr)+1),(Ntip(tr)+Nnode(tr)))
  if(!is.null(node.ages)) node.mins[match(names(node.ages),names(node.mins))]<-node.ages

  if(any(is.na(charmatch(tr$tip.label, names(tip.ages))))) stop("names of tip.ages do not match the tree tip labels")
  tip.ages[charmatch(tr$tip.label, names(tip.ages))]->tip.ages

  c(tip.ages,node.mins)->aa
  names(aa)[1:Ntip(tree)]<-seq(1,Ntip(tree))
  names(is.na(node.mins))->ii
  tn<- tr$edge

  # go through internal nodes, assign ages to them
  while (sum(is.na(aa))>=1){ # loop through as long as ages for some inodes not yet known
    for (i in ii){
      ci<- i
      dec<- tn[which(tn[,1]==i),2]# direct descendants of i
      aad<- aa[as.character(dec)] # ages of these direct descendants

      if (sum(is.na(aad))==0){ # if all ages are known.
        if(!is.na(node.mins[ci])) aa[ci]<- node.mins[ci] else aa[ci]<- max(aad) + min.diff
      }
    }
  }

  aa[-which(as.numeric(names(aa))==(Ntip(tree)+1))]->aa.nods
  for(i in 1:length(aa.nods)){
    aa.nods[i]->x
    aa[as.character(getMommy(tr,names(x))[1])]->agemom
    if(as.numeric(names(x))>Ntip(tree)) aa[as.character(getDescendants(tr,names(x))[1:2])]->agedes else agedes<-0
    if(x>agemom|any(x<agedes)) stop("min.branch is too big")
  }

  ne<- nrow(tr$edge)
  bl<- array(dim=ne)

  for (i in 1:ne){
    anc<- tr$edge[i,1]
    dec<- tr$edge[i,2]
    bl[i]<- aa[anc] - aa[dec]
  }

  t2<- tr
  t2$edge.length<- bl
  return(t2)
}

