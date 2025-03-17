#'@title Phylogenetic tree calibration
#'@description The function is a wrapper around the functions "scalePhylo",
#'  "assign.ages", and "assign.brlen" written by Gene Hunt
#'  (http://paleobiology.si.edu/staff/individuals/hunt.cfm). It rescales tree
#'  branch lengths according to given calibration dates.
#'@usage scaleTree(tree, tip.ages=NULL, node.ages=NULL)
#'@param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'  fully dichotomous.
#'@param tip.ages a named vector including the ages (i.e. distance from the
#'  youngest tip within the tree) of the tips to be changed. If unspecified, the
#'  function assumes all the tips are correctly placed with respect to the root.
#'  Names can be either tip labels or numbers.
#'@param node.ages a named vector including the ages (i.e. distance from the
#'  youngest tip within the tree) of the nodes to be changed. If no calibration
#'  date for nodes is supplied, the tree root is fixed and the function shifts
#'  node position only where needed to fit tip ages.
#'@export
#'@seealso \href{../doc/Tree-Manipulation.html#scaleTree}{\code{scaleTree}
#'  vignette}
#'@return Rescaled phylogentic tree with tip labels ordered according to their
#'  position in the tree.
#'@author Silvia Castiglione, Pasquale Raia, Carmela Serio, Alessandro
#'  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'  Carotenuto
#' @examples
#' \donttest{
#' library(ape)
#' library(phytools)
#'
#' data("DataFelids")
#' DataFelids$treefel->treefel
#'
#' max(nodeHeights(treefel))->H
#'
#' #### Example 1 ####
#' rep(0,4)->tipAges
#' names(tipAges)<-tips(treefel,146)
#' scaleTree(treefel,tipAges)->treeS1
#'
#' edge.col<-rep("black",nrow(treefel$edge))
#' edge.col[which(treeS1$edge[,2]%in%getDescendants(treeS1,146))]<-"red"
#'
#' layout(2:1)
#' plot(treefel,edge.color = edge.col,show.tip.label=FALSE)
#' plot(treeS1,edge.color = edge.col,show.tip.label=FALSE)
#'
#' #### Example 2 ####
#' nodeAges<-c(23.5,15.6)
#' names(nodeAges)<-c(85,139)
#' scaleTree(treefel,node.ages=nodeAges)->treeS2
#'
#' edge.col<-rep("black",nrow(treefel$edge))
#' edge.col[which(treeS1$edge[,2]%in%c(getDescendants(treeS1,85),
#'                                     getDescendants(treeS1,139)))]<-"red"
#'
#' layout(2:1)
#' plot(treefel,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=c(85,139),col="green")
#' plot(treeS2,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=c(85,139),col="green")
#'
#' #### Example 3 ####
#' 16->nodeAges
#' names(nodeAges)<-"145"
#' tipAges<-19
#' names(tipAges)<-treefel$tip.label[1]
#' scaleTree(treefel,tip.ages = tipAges,node.ages=nodeAges)->treeS3
#'
#' edge.col<-rep("black",nrow(treefel$edge))
#' edge.col[which(treeS3$edge[,2]%in%c(1,getMommy(treefel,1),
#'                                     getDescendants(treeS3,145)))]<-"red"
#'
#' layout(2:1)
#' plot(treefel,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=145,col="green")
#' plot(treeS3,edge.color = edge.col,show.tip.label=FALSE)
#' nodelabels(bg="w",frame="n",node=145,col="green")
#'}


scaleTree<- function(tree, tip.ages=NULL, node.ages=NULL)
{
  # require(ape)
  # require(phytools)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  if(!is.null(ncol(tip.ages))){
    tag<-tip.ages[,1]
    names(tag)<-rownames(tip.ages)
    tip.ages<-tag
  } # stop("tip.ages must be a vector")
  if(!is.null(ncol(node.ages))) {
    nag<-node.ages[,1]
    names(nag)<-rownames(node.ages)
    node.ages<-nag
  } #stop("node.ages must be a vector")

  if(is.null(tip.ages)) max(diag(vcv(tree)))-diag(vcv(tree))->tip.ages else{

    if(all(!is.na(suppressWarnings(as.numeric(names(tip.ages))))))
      names(tip.ages)<-tree$tip.label[as.numeric(names(tip.ages))]

    if(anyDuplicated(names(tip.ages))){
      sapply(names(which(table(names(tip.ages))>1)),function(k)
        if(length(unique(tip.ages[which(names(tip.ages)%in%k)]))>1)  NA else k)->dup.tip
      if(any(is.na(dup.tip)))
        stop(paste("More than one age value supplied for",paste(names(dup.tip)[which(is.na(dup.tip))],collapse=", ")))
      tip.ages[which(!duplicated(names(tip.ages)))]->tip.ages
    }

    max(diag(vcv(tree)))-diag(vcv(tree))->ta
    ta[match(names(tip.ages),names(ta))]<-tip.ages
    ta->tip.ages
    if(any(is.na(charmatch(tree$tip.label, names(tip.ages))))) stop("names of tip.ages do not match the tree tip labels")
    tip.ages[charmatch(tree$tip.label, names(tip.ages))]->tip.ages
  }

  if(anyDuplicated(names(node.ages))){
    sapply(names(which(table(names(node.ages))>1)),function(k)
      if(length(unique(node.ages[which(names(node.ages)%in%k)]))>1) NA else k)->dups
    if(any(is.na(dups)))
      stop(paste("More than one age value supplied for nodes",paste(names(dups)[which(is.na(dups))],collapse=", ")))
    node.ages[-which(duplicated(node.ages))]->node.ages
  }


  if(!(Ntip(tree)+1)%in%names(node.ages)&!any(node.ages>=max(diag(vcv(tree))))){
    # If the tree root is not fixed and no other node is older than the current root age
    c(max(diag(vcv(tree))),node.ages)->node.ages
    names(node.ages)[1]<-(Ntip(tree)+1)
  }else{
    if((Ntip(tree)+1)%in%names(node.ages)&any(node.ages[-which(names(node.ages)==as.character(Ntip(tree)+1))]>=node.ages[as.character(Ntip(tree)+1)]))
      # If the tree root is fixed and some other nodes are older than the fixed root age
      stop("The age for ",paste(names(node.ages)[which(node.ages[-which(names(node.ages)==as.character(Ntip(tree)+1))]>=node.ages[as.character(Ntip(tree)+1)])],collapse=", "),
           " cannot be older the the age indicated for the tree root")

    if(!(Ntip(tree)+1)%in%names(node.ages)&any(node.ages>=max(diag(vcv(tree))))){
      # If the tree root is not fixed and any other nodes are older than the current root age
      node.ages[which(node.ages>=max(diag(vcv(tree))))]->max.ages
      warning("The age for ",paste(names(max.ages),collapse=", "),
              " is older than the tree root. The tree root has been moved, tree height has changed.")
      c(max(max.ages)+min(tree$edge.length),node.ages)->node.ages
      names(node.ages)[1]<-(Ntip(tree)+1)
    }
    # rescaleRR(tree,height=node.ages[as.character(Ntip(tree)+1)])->tree
    tree$edge.length[which(tree$edge==(Ntip(tree)+1))]<-tree$edge.length[which(tree$edge==(Ntip(tree)+1))]+(node.ages[as.character(Ntip(tree)+1)]-max(diag(vcv(tree))))
  }

  c(max(nodeHeights(tree)),max(nodeHeights(tree))-nodeHeights(tree)[which(tree$edge[,2]>Ntip(tree)),2])->node.mins
  names(node.mins)<-c(Ntip(tree)+1,tree$edge[which(tree$edge[,2]>Ntip(tree)),2])

  if(!is.null(node.ages)) node.mins[match(names(node.ages),names(node.mins))]<-node.ages

  sapply((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)),function(x){
    c(length(tips(tree,x)),node.mins[as.character(x)],length(getMommy(tree,x)))
  })->node.tips
  colnames(node.tips)<-(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
  node.tips[,order(node.tips[1,],node.tips[2,],node.tips[3,],decreasing = c(FALSE,TRUE,TRUE))]->node.tips
  node.mins[match(colnames(node.tips),names(node.mins))]->nn
  names(tip.ages)<-which(tree$tip.label%in%names(tip.ages))
  c(sort(tip.ages,decreasing = TRUE),nn)->age.vec
  names(age.vec)->nam.vec

  data.frame(ed=tree$edge,bl=NA)->eds.dat

  lenpb<-length(nam.vec)
  pb = txtProgressBar(min = 0, max = lenpb, style = 3)
  on.exit(close(pb))
  while(length(nam.vec)>0){
    setTxtProgressBar(pb,lenpb-length(nam.vec))
    nam.vec[1]->nnfoc
    age.vec[as.character(nnfoc)]->agefoc
    age.vec[as.character(getMommy(tree,nnfoc))]->moms
    if(any(moms<=agefoc)){

      which(names(moms)%in%names(node.ages))[1]->mompos
      node.ages[which(names(node.ages)==names(moms[mompos]))]->momage

      if(nnfoc%in%names(node.ages)|as.numeric(nnfoc)<=Ntip(tree)){
        if(agefoc>=momage) stop(paste("The age for",nnfoc,"cannot be older than the age indicated for",names(momage)))
        if(mompos==1){
          nam.vec[-which(nam.vec==nnfoc)]->nam.vec
          next
        }

        (momage-agefoc)/mompos->bl1
        rep(bl1,mompos)->bl2
        names(bl2)<-c(nnfoc,names(moms)[(1:mompos)-1])

        eds.dat[match(names(bl2),eds.dat$ed.2),]->bldat
        if(any(!is.na(bldat$bl))){
          bldat[which(!is.na(bldat$bl)),]->bldat
          if(any(bldat$bl<bl2[match(bldat$ed.2,names(bl2))])){
            bldat[which(bldat$bl<bl2[match(bldat$ed.2,names(bl2))]),]->bldat.cut
            bl2[match(bldat.cut$ed.2,names(bl2))]-bldat.cut$bl->bldiff
            bldat.cut$bl->bl2[match(bldat.cut$ed.2,names(bl2))]
            sapply(1:length(bldiff),function(k)
              bl2[which(names(bl2)%in%getDescendants(tree,names(bldiff)[k])[1:2])]<<-
                bl2[which(names(bl2)%in%getDescendants(tree,names(bldiff)[k])[1:2])]+bldiff[k])
          }
        }

        agefoc+cumsum(bl2[-length(bl2)])->new.ages
        names(new.ages)<-names(bl2)[2:length(bl2)]
        eds.dat$bl[match(names(bl2),eds.dat$ed.2)]<-bl2
        age.vec[match(names(bl2)[2:length(bl2)],names(age.vec))]<-new.ages
      }else {
        getDescendants(tree,nnfoc)[1:sum(tree$edge[,1]==as.numeric(nnfoc))]->des

        max(age.vec[as.character(des)])->maxdes
        min(age.vec[as.character(des)])->mindes
        names(maxdes)<-names(which.max(age.vec[as.character(des)]))
        names(mindes)<-names(which.min(age.vec[as.character(des)]))

        (momage-maxdes)/(mompos+1)->bl1
        c((maxdes+bl1-mindes),rep(bl1,(mompos+1)))->bl2
        names(bl2)<-c(names(mindes),names(maxdes),nnfoc,names(moms)[(1:mompos)-1])
        eds.dat[match(names(bl2),eds.dat$ed.2),]->bldat
        if(any(!is.na(bldat$bl))){
          bldat[which(!is.na(bldat$bl)),]->bldat
          if(any(bldat$bl<bl2[match(bldat$ed.2,names(bl2))])){
            bldat[which(bldat$bl<bl2[match(bldat$ed.2,names(bl2))]),]->bldat.cut
            bl2[match(bldat.cut$ed.2,names(bl2))]-bldat.cut$bl->bldiff
            bldat.cut$bl->bl2[match(bldat.cut$ed.2,names(bl2))]
            sapply(1:length(bldiff),function(k)
              bl2[which(names(bl2)%in%getDescendants(tree,names(bldiff)[k])[1:2])]<<-
                bl2[which(names(bl2)%in%getDescendants(tree,names(bldiff)[k])[1:2])]+bldiff[k])
          }
        }
        maxdes+cumsum(bl2[2:(length(bl2)-1)])->new.ages
        names(new.ages)<-names(bl2)[-c(1,2)]
        eds.dat$bl[match(names(bl2),eds.dat$ed.2)]<-bl2
        age.vec[match(names(bl2)[-c(1,2)],names(age.vec))]<-new.ages
      }

      unlist(lapply(1:length(new.ages),function(k) {
        age.vec[as.character(getDescendants(tree,names(new.ages)[k])[1:2])]->agedes
        if(any(agedes>new.ages[k])) names(agedes)[which(agedes>new.ages[k])]
      }))->namadd
      if(any(!namadd%in%nam.vec)) c(namadd[which(!namadd%in%nam.vec)],nam.vec)->nam.vec
      nam.vec[-which(nam.vec==nnfoc)]->nam.vec

    }else{
      eds.dat$bl[match(nnfoc,eds.dat$ed.2)]<-moms[1]-agefoc
      nam.vec[-which(nam.vec==nnfoc)]->nam.vec
    }

  }

  t2<- tree
  apply(eds.dat,1,function(j){
    age.vec[as.character(j[1])]-age.vec[as.character(j[2])]
  })->t2$edge.length

  return(t2)
}
