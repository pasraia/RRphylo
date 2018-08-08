#' @title Searching morphological convergence
#' @description The function scans a phylogenetic tree looking for morphological convergence between pairs of distant nodes. (For further details see \code{\link{evo.dir}})
#' @usage search.conv(RR,y,nodes=NULL,min.dim=NULL,max.dim=NULL,min.dist=NULL,nsim=1000,clus=.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a multivariate dataset of phenotypic variables.
#' @param nodes node pair to be tested. If unspecified, the function automatically searches for converging subtrees.
#' @param min.dim the minimum size of the subtrees to be compared. When \code{nodes} is indicated, it coincides with the size of the smallest clades in \code{nodes}, otherwise it is set at one tenth of the tree size.
#' @param max.dim the maximum size of the subtrees to be compared. When \code{nodes} is indicated, it is \code{min.dim}*2 if the largest clade in \code{nodes} is smaller than this value, otherwise it corresponds to the size of the largest clade. Whitout \code{nodes} it is set at one third of the tree size.
#' @param min.dist the minimum distance, in terms of number of nodes, between the subtrees to be compared. When \code{nodes} is indicated, it is the distance between the pair, otherwise it is set at 10 nodes.
#' @param nsim number of simulations to be performed. It is set at 1000 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @export
#' @return A data frame including for each pair of nodes:
#' @return \strong{real.diff} the difference between the mean evolutionary direction of the tips of the respective nodes divided by squared \code{nod.dist}
#' @return \strong{diff} the difference between the mean evolutionary direction of the tips of the respective nodes plus the angle between ancestral shapes (\code{ang}), divided by squared \code{nod.dist}
#' @return \strong{ang} the angle between ancestral shapes
#' @return \strong{nod.dist} the distance in number of nodes between the focal nodes pair
#' @return \strong{p.real.diff} the p-value for \code{real.diff}
#' @return \strong{p.diff} the p-value for \code{diff}
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references
#' Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method for testing evolutionary rate variation and shifts in phenotypic evolution. \emph{Methods in Ecology and Evolution}, in press.doi:10.1111/2041-210X.12954
#' @examples
#'  \donttest{
#' ### Tree and data preparation ###
#' library(RRphylo)
#' library(ape)
#' library(geiger)
#' library(phytools)
#'
#' repeat({
#'   rtree(100)->tree
#'   sizedsubtree(tree,11)->n
#'   if(length(tips(tree,n))<15) break
#' })
#'
#'
#' setBM(tree,type="brown",nY=4)->y
#' y[match(tree$tip.label,rownames(y)),]->y
#' y[match(tips(tree,n),rownames(y)),]->a
#'
#'
#' apply(y,2,range)->m.a
#' m.a[2,]*.7->m.a
#'
#' a1<-matrix(ncol=dim(y)[2],nrow=dim(a)[1])
#' for(m in 1:dim(a)[1])
#' {
#'   v<-array()
#'   for(i in 1:length(m.a)) jitter(m.a[i],amount=(sd(a[,i])*1))->v[i]
#'   v->a1[m,]
#' }
#'
#' rownames(a1)<-rownames(a)
#' y[match(tips(tree,n),rownames(y)),]<-a1
#'
#'
#' extract.clade(tree,n)->t1
#' swapONE(t1)[[1]]->t1
#' drop.tip(t1,t1$tip.label[c(1,length(t1$tip.label))])->t1
#' t1$root.edge<-data.frame(tree$edge,tree$edge.length)[which(
#'   data.frame(tree$edge,tree$edge.length)[,2]==n),3]
#'
#' apply(y[match(t1$tip.label,rownames(y)),],2,jitter)->y.t1
#'
#'
#' tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
#' rownames(y)<-tree$tip.label
#' t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
#' rownames(y.t1)<-t1$tip.label
#'
#'
#' distNodes(tree,n)->dfN
#' dfN[which(dfN[,2]<=Ntip(tree)*0.1),1]->bar
#' dfN[-which(dfN[,2]<=Ntip(tree)*0.1),1]->nodN
#' sample(nodN,1)->tar
#'
#' if(length(which(getDescendants(tree,tar)%in%bar))>0)
#' {
#'   getDescendants(tree,tar)->ta
#'   ta[which(ta%in%bar)]->ex
#'   exa<-list()
#'   for(w in 1:length(ex)) tips(tree,ex[i])->exa[[w]]
#'   unlist(exa)->exil
#'   tips(tree,tar)->tipper
#'   tipper[-which(exil%in%tipper)]->tipper
#'   sample(tipper,1)->at.tip
#' }else{
#'   sample(tips(tree,tar),1)->at.tip
#' }
#'
#' at<-which(tree$tip.label==at.tip)
#' data.frame(tree$edge,tree$edge.length)[which(data.frame(tree$edge,
#'   tree$edge.length)[,2]==at),3]->pos
#' bind.tree(tree,t1,where=at,position=pos/2)->tree1
#' c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
#'   getMRCA(tree1,t1$tip.label))->nod.par
#'
#' ### nod.par is the node pair set to converge
#'
#' rbind(y,y.t1)->y
#' RRphylo(tree1,y)->RR
#' length((tips(tree1,nod.par[2])))->min.dim
#'
#' search.conv(RR,y,min.dim = min.dim,min.dist = NULL)
#'
#' search.conv(RR,y,nodes=nod.par)
#'
#'     }


search.conv<-function(RR,y,nodes=NULL,min.dim=NULL,max.dim=NULL,min.dist=NULL,nsim=1000,clus=.5)
{
  #require(ape)
  #require(geiger)
  #require(phytools)
  #require(foreach)
  #require(doParallel)
  #require(parallel)

  unitV <- function(x) {
    sum(x^2)^0.5
  }

  deg2rad <- function(deg) {
    (deg * pi)/(180)
  }
  rad2deg <- function(rad) {
    (rad * 180)/(pi)
  }

  RR$tree->tree1
  RR$aces->RRaces
  RR$tip.path->L

  if(is.null(nodes)){
    if(is.null(min.dim)) min.dim<-Ntip(tree1)*0.1 else min.dim<-min.dim
    if(is.null(min.dist)) min.dist<-10 else min.dist<-min.dist
    if(is.null(max.dim))  max.dim<-Ntip(tree1)/3 else max.dim<-max.dim
  }else{
    distNodes(tree1,nodes)->matDist
    if(is.null(min.dim)) min.dim<-min(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2])))) else min.dim<-min.dim
    if(is.null(min.dist)) min.dist<-matDist[,1] else min.dist<-min.dist
    if(is.null(max.dim)){
      if (max(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2]))))>min.dim*2) max.dim<-max(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2])))) else max.dim<-min.dim*2
    } else max.dim<-max.dim
  }

  subtrees(tree1)->subT
  if(length(unlist(lapply(subT,Ntip))>max.dim)>0){
    subT[-c(which(unlist(lapply(subT,Ntip))<min.dim),which(unlist(lapply(subT,Ntip))>max.dim))]->subT

  }else{
    subT[-which(unlist(lapply(subT,Ntip))<min.dim)]->subT
  }
  unlist(lapply(subT,function(x) getMRCA(tree1,x$tip.label)))->names(subT)
  as.numeric(names(subT))->nod
  c(nod,Ntip(tree1)+1)->nod
  subT[[length(subT)+1]]<-tree1
  names(subT)[length(subT)]<-Ntip(tree1)+1

  RRaces[match(nod,rownames(RR$aces)),]->aces
  rbind(RRaces,y)->phen

  ang.tot<-list()
  for(i in 1:length(nod)){
    nod[i]->sel
    RRaces[match(sel,rownames(RRaces)),]->ace.sel
    Lsub <- L[match(tips(tree1, sel), rownames(L)),]
    subN<-sel
    savenode <- getDescendants(tree1, subN)
    savenode <- c(subN, savenode)
    cutnode <- which(apply(Lsub, 2, sum) == 0)
    if (length(which(names(cutnode) %in% savenode)) == 0) cutnode <- cutnode else cutnode <- cutnode[-which(names(cutnode) %in% savenode)]
    Lsub <- Lsub[, which(!colnames(Lsub) %in% names(cutnode))]
    Lsub <- Lsub[, seq(match(subN, colnames(Lsub)),dim(Lsub)[2])]
    Htree <- subT[[i]]
    Htree$node.label <- colnames(Lsub)[seq(1:(Ntip(Htree) - 1))]
    subetas <- phen[match(colnames(Lsub), rownames(phen)),]
    vec.len <- apply(subetas, 1, function(x) sum(x^2)^0.5)
    vec.len <- as.matrix(vec.len)
    tips(tree1,sel)->tippa

    matrix(ncol=2,nrow=length(tippa))->ang2sel

    for (p in 1:length(tippa)){
      a <- Lsub[match(tippa[p], rownames(Lsub)), ]
      node.a <- savenode[which(savenode %in% getMommy(tree1,which(tree1$tip.label == tippa[p])))]
      a <- names(a[c(which(names(a) %in% node.a), unname(which(a !=0)))])
      a <- a[!duplicated(a)]

      a.betas <- phen[match(a, rownames(phen)), ]
      a.betas <- as.matrix(a.betas)
      a.betas<- rbind(rep(1,ncol(a.betas)),a.betas)
      if (dim(a.betas)[1] == 2) a.resultant <- apply(a.betas, 2, diff) else a.resultant <- apply(a.betas, 2, function(x) x[1]-sum(x[2:length(x)]))

      angA <- rad2deg(acos((a.betas[1,] %*% a.resultant)/(unitV(a.betas[1,]) *unitV(a.resultant))))
      a.size <- unitV(a.resultant)
      angA->ang2sel[p,1]
      a.size->ang2sel[p,2]
    }
    rownames(ang2sel)<-tippa
    colnames(ang2sel)<-c("angle","size")
    data.frame(node=rep(sel,dim(ang2sel)[1]),species=rownames(ang2sel),ang2sel)->ang2nod

    ang2nod->ang.tot[[i]]
  }

  do.call(rbind,ang.tot)->ang
  subset(ang,ang[,1]==Ntip(tree1)+1)->angR
  subset(ang,ang[,1]!=Ntip(tree1)+1)->ang

  apply(L,1,function(x) length(which(x!=0))-1)->par
  par[match(angR[,2],names(par))]->par
  data.frame(angR,dist=unname(par))->angR
  apply(angR,1,function(x) as.numeric(x[3])/as.numeric(x[5])^2)->pang
  data.frame(angR,pang)->angR
  angR[,c(1,2,6,3,4,5)]->angR

  tapply(as.numeric(as.character(ang[,3])),ang[,1],mean)->mean.ang

  ang.nRoot<-array()
  for (i in 1:length(mean.ang)){
    aces[which(rownames(aces)%in%c(names(mean.ang[i]),(Ntip(tree1)+1))),]->acRoot
    rad2deg(acos((acRoot[1,] %*% acRoot[2,])/(unitV(acRoot[1,]) *unitV(acRoot[2,]))))->ang.nRoot[i]
  }
  names(ang.nRoot)<-names(mean.ang)
  mean(ang.nRoot)->mean.aR

  if(is.null(nodes)){
    res<-list()
    cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res <- foreach(i = 1:length(mean.ang),
                   .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
                   {
                     gc()

                     mean.ang[i]->sel1
                     subset(ang,ang[,1]==names(sel1))->ang.sel


                     if(length(which(names(mean.ang)%in%getDescendants(tree1,names(sel1))))>0)
                       mean.ang[-which(names(mean.ang)%in%getDescendants(tree1,names(sel1)))]->mean.sel else
                         mean.ang->mean.sel

                     if(length(which(names(mean.sel)%in%getMommy(tree1,names(sel1))))>0)
                       mean.sel[-which(names(mean.sel)%in%getMommy(tree1,names(sel1)))]->mean.sel else
                         mean.sel->mean.sel
                     mean.sel[-which(names(mean.sel)==names(sel1))]->mean.sel

                     distNodes(tree1,names(sel1))->matDist
                     matDist->matNod
                     if (length(which(names(mean.sel) %in% matDist[which(matDist[,2]<min.dist),1]))>0) mean.sel[-which(names(mean.sel) %in% matDist[which(matDist[,2]<min.dist),1])]->mean.sel else mean.sel->mean.sel



                     if(length(mean.sel)==0) {
                       c(dir.diff=NULL,diff=NULL,ang=NULL)->diff.p
                       nDD<-NULL
                     }else{
                       nDD<-array()
                       diff.p<-matrix(ncol=4,nrow=length(mean.sel))
                       for(k in 1:length(mean.sel)){
                         matDist[match(as.numeric(as.character(names(mean.sel[k]))),matDist[,1]),2]->nD
                         nD->nDD[k]
                         aces[which(rownames(aces)%in%c(names(sel1),names(mean.sel[k]))),]->ac
                         rad2deg(acos((ac[1,] %*% ac[2,])/(unitV(ac[1,]) *unitV(ac[2,]))))->ang.ac
                         c(dir.diff=abs(sel1-mean.sel[k])/(nD^2),diff=(abs(sel1-mean.sel[k])/(nD^2))+(ang.ac)/(nD^2),ang=ang.ac,nD=nD)->diff.p[k,]
                       }
                       rownames(diff.p)<-names(mean.sel)
                       colnames(diff.p)<-c("real.diff","diff","ang","nod.dist")

                     }

                     diff.p->mean.diff

                     res[[i]]<-list(matNod,mean.diff,mean.sel,nDD)
                   }

    stopCluster(cl)



    lapply(res,"[[",1)->matNod
    lapply(res,"[[",2)->mean.diff
    lapply(res,"[[",3)->mean.sel
    lapply(res,"[[",4)->nD.sel

    names(mean.ang)->names(mean.diff)->names(matNod)->names(mean.sel)->names(nD.sel)


    sapply(mean.diff,is.null)->nulls
    if(length(which(nulls==TRUE)>0)) mean.diff[-which(nulls==TRUE)]->mean.diff
    if(length(which(nulls==TRUE)>0)) mean.sel[-which(nulls==TRUE)]->mean.sel
    #if(length(which(nulls==TRUE)>0)) ang.nRoot[-which(nulls==TRUE)]->ang.nRoot
    if(length(which(nulls==TRUE)>0)) matNod[-which(nulls==TRUE)]->matNod
    if(length(which(nulls==TRUE)>0)) nD.sel[-which(nulls==TRUE)]->nD.sel



    res.ran <- list()
    cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res.ran <- foreach(i = 1:nsim,
                       .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
                       {
                         gc()

                         mean.diffR<-list()
                         for(u in 1:length(mean.diff)){
                           names(mean.diff)[[u]]->sel1
                           mean.sel[[u]]->msel
                           nD.sel[[u]]->ndsel
                           diff.pR<-list()
                           for(j in 1:length(msel)){
                             names(msel)[j]->sel2
                             if(length(tips(tree1,sel1))>dim(angR)[1]){
                               mean(angR[,3])->mean.angR1
                             }else{
                               sample(seq(1,dim(angR)[1]),length(tips(tree1,sel1)))->sam1
                               mean(angR[sam1,3])->mean.angR1
                             }
                             if(length(tips(tree1,sel2))>dim(angR)[1]){
                               mean(angR[,3])->mean.angR2
                             } else {
                               if(length(sam1)>(dim(angR)[1]-length(sam1))){
                                 sample(seq(1,dim(angR)[1]),length(tips(tree1,sel2)),replace=TRUE)->sam2
                               } else {
                                 sample(seq(1,dim(angR)[1])[-sam1],length(tips(tree1,sel2)),replace=TRUE)->sam2
                                 mean(angR[sam2,3])->mean.angR2
                               }
                             }



                             data.frame(dir.diff=abs(mean.angR1-mean.angR2),diff=(abs(mean.angR1-mean.angR2)+(mean.aR/ndsel[j]^2)))->diff.pR[[j]]
                           }
                           do.call(rbind,diff.pR)->mean.diffR[[u]]
                           rownames(mean.diffR[[u]])<-names(msel)
                         }
                         names(mean.diffR)<-names(mean.diff)
                         mean.diffR->res.ran[[i]]

                       }
    stopCluster(cl)



    diff.rank<-list()
    for(i in 1:length(mean.diff)){
      rnk<-matrix(ncol=2,nrow=dim(mean.diff[[i]])[1])
      for(k in 1:dim(mean.diff[[i]])[1]){
        apply(rbind(mean.diff[[i]][k,c(1,2)],do.call(rbind,lapply(lapply(res.ran,"[[",i),function(x) x[k,c(1,2)]))[1:(nsim-1),]),2, function(x) rank(x)[1]/nsim )->rnk[k,]

      }
      data.frame(mean.diff[[i]],rnk)->diff.rank[[i]]
      colnames(diff.rank[[i]])[c(5,6)]<-c("p.real.diff","p.diff")

    }
    names(mean.diff)->names(diff.rank)

    lapply(diff.rank,function(x) abs(x[order(abs(x[,1])),]))->diff.rank

    as.data.frame(do.call(rbind,lapply(diff.rank,function(x) cbind(rownames(x)[1],x[1,]))))->df
    colnames(df)[1]<-"node"
    df[order(df[,7],df[,6]),]->df

  }else{

    mean.ang[which(names(mean.ang)%in%nodes)]->mean.angN

    matDist[,1]->nD
    aces[which(rownames(aces)%in%nodes),]->ac
    rad2deg(acos((ac[1,] %*% ac[2,])/(unitV(ac[1,]) *unitV(ac[2,]))))->ang.ac
    data.frame(dir.diff=abs(diff(mean.angN))/(nD^2),diff=(abs(diff(mean.angN))/(nD^2))+(ang.ac)/(nD^2),ang=ang.ac,nD=nD)->mean.diff

    rownames(mean.diff)<-paste(nodes[1],nodes[2],sep="-")
    colnames(mean.diff)<-c("real.diff","diff","ang","nod.dist")


    res.ran <- list()
    cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res.ran <- foreach(i = 1:nsim,
                       .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
                       {
                         gc()

                         sample(seq(1,dim(angR)[1]),length(tips(tree1,nodes[1])))->sam1
                         mean(angR[sam1,3])->mean.angR1

                         sample(seq(1,dim(angR)[1])[-sam1],length(tips(tree1,nodes[2])),replace=TRUE)->sam2
                         mean(angR[sam2,3])->mean.angR2

                         data.frame(dir.diff=abs(mean.angR1-mean.angR2),diff=(abs(mean.angR1-mean.angR2)+(mean.aR/nD^2)))->diff.pR
                         colnames(diff.pR)<-c("real.diff","diff")
                         diff.pR->res.ran[[i]]

                       }
    stopCluster(cl)

    apply(rbind(mean.diff[,c(1,2)],do.call(rbind,res.ran)[1:(nsim-1),]),2,function(x) rank(x)[1]/nsim)->rnk
    data.frame(mean.diff,t(as.data.frame(rnk)))->df
    colnames(df)[c(5,6)]<-c("p.real.diff","p.diff")

  }
  return(df)

}
