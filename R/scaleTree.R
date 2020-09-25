#'@title Phylogenetic tree calibration
#'@description The function is a wrapper around the functions "scalePhylo",
#'  "assign.ages", and "assign.brlen" written by Gene Hunt
#'  (http://paleobiology.si.edu/staff/individuals/hunt.cfm). It rescales tree
#'  branch lengths according to given calibration dates.
#'@usage scaleTree(tree, tip.ages, node.ages=NULL, min.branch=0.1)
#'@param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'  fully dichotomous.
#'@param tip.ages a named vector including the ages (i.e. distance from the
#'  youngest tip within the tree) of the tips to be changed. If unspecified, the
#'  function assumes all the tips are correctly placed with respect to the root.
#'@param node.ages a named vector including the ages (i.e. distance from the
#'  youngest tip within the tree) of the nodes to be changed. If no calibration
#'  date for nodes is supplied, the function shifts node position only where
#'  needed to fit tip ages.
#'@param min.branch the minimum branch length that will be imposed for shifted
#'  nodes.
#'@export
#'@seealso \href{../doc/Tree-Manipulation.html#scaleTree}{\code{scaleTree} vignette}
#'@importFrom geiger rescale
#'@return Rescaled phylogentic tree.
#'@author Silvia Castiglione, Pasquale Raia, Carmela Serio, Alessandro
#'  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'  Carotenuto
#' @examples
#' \donttest{
#' library(ape)
#' library(phytools)
#' library(geiger)
#'
#' data("DataFelids")
#' DataFelids$treefel->tree
#'
#' max(nodeHeights(tree))->H
#'
#' #### Example 1 ####
#' rep(0,4)->tipAges
#' names(tipAges)<-tips(tree,146)
#' scaleTree(tree,tipAges)->treeS1
#'
#' edge.col<-rep("black",nrow(tree$edge))
#' edge.col[which(treeS1$edge[,2]%in%getDescendants(treeS1,146))]<-"red"
#'
#' layout(2:1)
#' plot(tree,edge.color = edge.col,show.tip.label=FALSE)
#' plot(treeS1,edge.color = edge.col,show.tip.label=FALSE)
#'
#' #### Example 2 ####
#' nodeAges<-c(23.5,15.6)
#' names(nodeAges)<-c(85,139)
#' scaleTree(tree,node.ages=nodeAges)->treeS2
#'
#' edge.col<-rep("black",nrow(tree$edge))
#' edge.col[which(treeS1$edge[,2]%in%c(getDescendants(treeS1,85),
#'                                     getDescendants(treeS1,139)))]<-"red"
#'
#' layout(2:1)
#' plot(tree,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=c(85,139),col="green")
#' plot(treeS2,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=c(85,139),col="green")
#'
#' #### Example 3 ####
#' 16->nodeAges
#' names(nodeAges)<-"145"
#' tipAges<-19
#' names(tipAges)<-tree$tip.label[1]
#' scaleTree(tree,tip.ages = tipAges,node.ages=nodeAges)->treeS3
#'
#' edge.col<-rep("black",nrow(tree$edge))
#' edge.col[which(treeS3$edge[,2]%in%c(1,getMommy(tree,1),
#'                                     getDescendants(treeS3,145)))]<-"red"
#'
#' layout(2:1)
#' plot(tree,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=145,col="green")
#' plot(treeS3,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=145,col="green")
#'}


scaleTree<- function(tree, tip.ages=NULL, node.ages=NULL, min.branch=0.1)
{
  # require(ape)
  # require(phytools)

  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }

  abs(min.branch)->min.diff

  if(is.null(tip.ages)) max(diag(vcv(tree)))-diag(vcv(tree))->tip.ages else{
    max(diag(vcv(tree)))-diag(vcv(tree))->ta
    ta[match(names(tip.ages),names(ta))]<-tip.ages
    ta->tip.ages
    if(any(is.na(charmatch(tree$tip.label, names(tip.ages))))) stop("names of tip.ages do not match the tree tip labels")
    tip.ages[charmatch(tree$tip.label, names(tip.ages))]->tip.ages
  }

  if(!is.null(node.ages)&(Ntip(tree)+1)%in%names(node.ages)) rescale(tree,"depth",node.ages[as.character(Ntip(tree)+1)])->tree

  c(max(nodeHeights(tree)),max(nodeHeights(tree))-nodeHeights(tree)[which(tree$edge[,2]>Ntip(tree)),2])->node.mins
  names(node.mins)<-c(Ntip(tree)+1,tree$edge[which(tree$edge[,2]>Ntip(tree)),2])

  if(!is.null(node.ages)) node.mins[match(names(node.ages),names(node.mins))]<-node.ages

  sapply((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)),function(x) length(tips(tree,x)))->node.tips
  names(node.tips)<-(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
  sort(node.tips)->node.tips
  node.mins[match(names(node.tips),names(node.mins))]->nn
  names(tip.ages)<-which(tree$tip.label%in%names(tip.ages))
  c(tip.ages,nn)->age.vec

  i=1
  while(i<=length(age.vec)){
    age.vec[as.character(getMommy(tree,names(age.vec[i])))]->moms
    if(any(moms<age.vec[i])){
      if(as.numeric(names(moms[which(moms<age.vec[i])][1]))==(Ntip(tree)+1)) warning("The tree root has been moved, tree height has changed")
      if(!is.null(node.ages)&names(moms[which(moms<age.vec[i])][1])%in%names(node.ages)){
        names(age.vec[i])->nam
        if(as.numeric(nam)<=Ntip(tree)|nam%in%names(node.ages)) stop(paste("The age for",nam,"cannot be older than the age indicated for",names(moms[which(moms<age.vec[i])][1])))
        tree$edge[which(tree$edge[,1]==nam),2]->des
        age.vec[i]<-mean(c(max(age.vec[as.character(des)]),moms[which(moms<age.vec[i])][1]))
      }else age.vec[names(moms[which(moms<age.vec[i])][1])]<- age.vec[i] + min.diff
    }
    i=i+1
  }

  ne<- nrow(tree$edge)
  bl<- array(dim=ne)

  for (i in 1:ne){
    anc<- tree$edge[i,1]
    dec<- tree$edge[i,2]
    bl[i]<- age.vec[as.character(anc)] - age.vec[as.character(dec)]
  }

  t2<- tree
  t2$edge.length<- bl
  return(t2)
}
