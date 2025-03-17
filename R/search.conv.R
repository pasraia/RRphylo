#' @title Searching for morphological convergence among species and clades
#' @description The function scans a phylogenetic tree looking for morphological
#'   convergence between entire clades or species evolving under specific
#'   states.
#' @usage search.conv(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
#'   min.dim=NULL,max.dim=NULL,min.dist=NULL,declust=FALSE,nsim=1000,rsim=1000,
#'   clus=0.5)
#' @param RR an object produced by \code{\link{RRphylo}}. This is not indicated
#'   if convergence among states is tested.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous. This is not indicated if convergence among clades is
#'   tested.
#' @param y a multivariate phenotype. The object \code{y} should be either a
#'   matrix or dataframe with species names as rownames.
#' @param nodes node pair to be tested. It can be either a vector of two, or a
#'   two-columns matrix/data.frame of node pairs to be tested. Notice the node
#'   number must refer to the dichotomic version of the original tree, as
#'   produced by \code{\link{RRphylo}}. If unspecified, the function automatically
#'   searches for convergence among clades.
#' @param state the named vector of tip states. The function tests for
#'   convergence within a single state or among different states (this latter
#'   case is especially meant to test for iterative evolution as for example the
#'   appearance of repeated morphotypes into different clades). In both cases,
#'   the state for non-focal species (i.e. not belonging to any convergent
#'   group) must be indicated as "nostate".
#' @param aceV phenotypic values at internal nodes. The object \code{aceV}
#'   should be either a matrix or dataframe with nodes (referred to the
#'   dichotomic version of the original tree, as produced by \code{\link{RRphylo}}) as
#'   rownames. If \code{aceV} are not indicated, ancestral phenotypes are
#'   estimated via \code{\link{RRphylo}}.
#' @param min.dim the minimum size of the clades to be compared. When
#'   \code{nodes} is indicated, it is the minimum size of the smallest clades in
#'   \code{nodes}, otherwise it is set at one tenth of the tree size.
#' @param max.dim the maximum size of the clades to be compared. When
#'   \code{nodes} is indicated, it is \code{min.dim}*2 if the largest clade in
#'   \code{nodes} is smaller than this value, otherwise it corresponds to the
#'   size of the largest clade. Without \code{nodes} it is set at one third of
#'   the tree size.
#' @param min.dist the minimum distance between the clades to be compared. When
#'   \code{nodes} is indicated, it is the distance between the pair. Under the
#'   automatic mode, the user can choose whether time distance or node distance
#'   (i.e. the number of nodes intervening between the pair) should be used. If
#'   time distance has to be considered, \code{min.dist} should be a character
#'   argument containing the word "time" and then the actual time distance to be
#'   used. The same is true for node distance, but the word "node" must precede
#'   the node distance to be used. For example, if the user want to test only
#'   clades more distant than 10 time units, the argument should be "time10". If
#'   clades separated by more than 8 nodes have to be tested, the argument
#'   \code{min.dist} should be "node8". If left unspecified, it automatically
#'   searches for convergence between clades separated by a number of nodes
#'   larger than one tenth of the tree size.
#' @param declust if species under a given state (or a pair of states) to be
#'   tested for convergence are phylogenetically closer than expected by chance,
#'   trait similarity might depend on proximity rather than true convergence. In
#'   this case, by setting \code{declust = TRUE}, tips under the focal state (or
#'   states) are removed randomly until clustering disappears. A minimum of 3
#'   species per state is enforced to remain anyway.
#' @param nsim number of simulations to perform sampling within the theta random
#'   distribution. It is set at 1000 by default.
#' @param rsim number of simulations to be performed to produce the random
#'   distribution of theta values. It is set at 1000 by default.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{search.conv} set \code{clus} = 0.
#' @export
#' @seealso \href{../doc/search.conv.html}{\code{search.conv} vignette}
#' @seealso \code{\link{overfitSC}}; \href{../doc/overfit.html#overfitSC}{\code{overfitSC} vignette}
#' @seealso \code{\link{plotConv}}; \href{../doc/Plotting-tools.html#plotConv}{\code{plotConv} vignette}
#' @importFrom grDevices chull rgb
#' @importFrom graphics axis layout lines segments
#' @importFrom stats TukeyHSD aov princomp
#' @importFrom ape vcv cophenetic.phylo
#' @return If convergence between clades is tested, the function returns a list
#'   including:
#' @return \itemize{\item\strong{$node pairs}: a dataframe containing for each
#'   pair of nodes: \itemize{\item ang.bydist.tip: the mean theta angle between
#'   clades divided by the time distance. \item ang.conv: the mean theta angle
#'   between clades plus the angle between aces, divided by the time distance.
#'   \item ang.ace: the angle between aces. \item ang.tip: the mean theta angle
#'   between clades. \item nod.dist: the distance intervening between clades in
#'   terms of number of nodes. \item time.dist: the time distance intervening
#'   between the clades. \item p.ang.bydist: the p-value computed for
#'   ang.bydist.tip. \item p.ang.conv: the p-value computed for ang.conv. \item
#'   clade.size: the size of clades. } \item\strong{$node pairs comparison}:
#'   pairwise comparison between significantly convergent pairs (all pairs if no
#'   instance of significance was found) performed on the distance from group
#'   centroids (the mean phenotype per clade). \item\strong{$average distance
#'   from group centroids}: smaller average distances mean less variable
#'   phenotypes within the pair. }
#' @return If convergence between (or within a single state) states is tested,
#'   the function returns a list including: \itemize{\item\strong{state.res} a
#'   dataframe including for each pair of states (or single state): \itemize{
#'   \item ang.state: the mean theta angle between species belonging to
#'   different states (or within a single state). \item ang.state.time: the mean
#'   of theta angle between species belonging to different states (or within a
#'   single state) divided by time distance. \item p.ang.state: the p-value
#'   computed for ang.state. \item p.ang.state.time: the p-value computed for
#'   ang.state.time.}} \itemize{\item\strong{plotData} a dataframe including
#'   data to plot the results via \code{\link{plotConv}}}.
#' @return The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto, Paolo Piras, Davide Tamagnini
#' @references Castiglione, S., Serio, C., Tamagnini, D., Melchionna, M.,
#'   Mondanaro, A., Di Febbraro, M., Profico, A., Piras, P.,Barattolo, F., &
#'   Raia, P. (2019). A new, fast method to search for morphological convergence
#'   with shape data. \emph{PLoS ONE}, 14, e0226949.
#'   https://doi.org/10.1371/journal.pone.0226949
#' @examples
#' \dontrun{
#' data("DataFelids")
#' DataFelids$PCscoresfel->PCscoresfel
#' DataFelids$treefel->treefel
#' DataFelids$statefel->statefel
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(treefel,PCscoresfel,clus=cc)->RRfel
#'
#'
#' ## Case 1. searching convergence between clades
#' # by setting min.dist as node distance
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",clus=cc)->sc.clade
#' # by setting min.dist as time distance
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="time38",clus=cc)->sc.clade.time
#' # by setting node pairs to be tested
#' nodpairs<-rbind(c(85,145),c(85,155))
#' search.conv(RR=RRfel, y=PCscoresfel, nodes=nodpairs ,clus=cc)->sc.node
#'
#' ## Case 2. searching convergence within a single state
#' search.conv(tree=treefel, y=PCscoresfel, state=statefel,declust=TRUE,clus=cc)->sc.state
#'   }

search.conv<-function(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
                      min.dim=NULL,max.dim=NULL,min.dist=NULL,
                      declust=FALSE,nsim=1000,rsim=1000,clus=0.5)
{
  # require(ape)
  # require(phytools)
  # require(foreach)
  # require(doParallel)
  # require(parallel)
  # require(vegan)
  # require(cluster)
  # require(plotrix)
  # require(RColorBrewer)
  # require(geomorph)

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package \"vegan\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  funcall <- match.call()
  cldef<-({
    'cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
     clusterEvalQ(cl, {library(ape)\nlibrary(phytools)})
    '
  })

  if(is.null(state)){
    if(is.null(RR)) stop("missing RRphylo object")

    RR$tree->tree
    RR$aces->RRaces
    if(!is.null(aceV)){
      if (inherits(aceV,"data.frame")) aceV <- as.matrix(aceV)
      RRaces[match(rownames(aceV),rownames(RRaces)),]<-aceV
    }

    #if (inherits(y,"data.frame"))
    # y <- treedata(tree, y, sort = TRUE)[[2]]
    y <- as.matrix(treedataMatch(tree, y)[[1]])

    # RR$tip.path->L
    RR$rates->betas
    RR$node.path->L1
    rbind(RRaces,y)->phen

    if(is.null(nodes)){
      if(is.null(min.dim)) min.dim<-Ntip(tree)*0.1
      if(is.null(min.dist)){
        min.dist<-Ntip(tree)*0.1
        dist.type<-"node"
      }else{
        if(grepl("time",min.dist)){
          min.dist<-as.numeric(gsub("time","",min.dist))
          dist.type<-"time"
        }else{
          min.dist<-as.numeric(gsub("node","",min.dist))
          dist.type<-"node"
        }
      }
      if(is.null(max.dim))  max.dim<-Ntip(tree)/5

      subtrees(tree)->subT->subTN
      # if(any(sapply(subT,Ntip)>max.dim))
      subT[c(which(sapply(subT,Ntip)>=min.dim&sapply(subT,Ntip)<=max.dim))]->subT # else
      #    subT[which(unlist(lapply(subT,Ntip))>=min.dim)]->subT

      # c(as.numeric(names(subT)),Ntip(tree)+1)->nod
      # subT[[length(subT)+1]]<-tree
      sapply(subT,function(x) getMRCA(tree,x$tip.label))->nod
      # names(subT)[length(subT)]<-Ntip(tree)+1

      t(combn(nod,2))->nodpairs
      split(nodpairs[,2],nodpairs[,1])->nodpairs
      output.sel<-lapply(1:length(nodpairs),function(w){
        nodpairs[[w]][which(!nodpairs[[w]]%in%getDescendants(tree,as.numeric(names(nodpairs)[w])))]->mean.sel
        mean.sel[which(!mean.sel%in%getMommy(tree,as.numeric(names(nodpairs)[w])))]->mean.sel
        if(length(mean.sel)>0){
          distNodes(tree,as.numeric(names(nodpairs)[w]),clus=clus)[1:Nnode(tree),]->distN
          ifelse(dist.type=="time",distN[,2]->matNod,distN[,1]->matNod)
          mean.sel[which(!mean.sel%in%names(matNod[which(matNod<min.dist)]))]
          list(mean.sel,distN)
        } else NULL
      })
      names(output.sel)<-names(nodpairs)

      if(all(sapply(output.sel,is.null)))
        stop("There are no nodes distant enough to search convergence, consider using a smaller min.dist")

      output.sel<-output.sel[which(!sapply(output.sel,is.null))]
      ifelse(dist.type=="time",sapply(output.sel,function(j)j[[2]][,2])->matNod,
             sapply(output.sel,function(j)j[[2]][,1])->matNod)

      RRaces[match(nod,rownames(RR$aces)),]->aces

    }else{
      if(is.null(ncol(nodes))) matrix(nodes,ncol=2)->nodes
      if(is.data.frame(nodes)) as.matrix(nodes)->nodes
      matDist<-do.call(rbind,apply(nodes,1,function(k) distNodes(tree,k),simplify = FALSE))

      if(any(matDist[,1]<Ntip(tree)*0.15)){
        ddpcr::quiet(apply(nodes[which(matDist[,1]<Ntip(tree)*0.15),,drop=FALSE],1, function(k)
          warning(paste("Clades",paste(k,collapse=" and "),"are close to each other. Similarity might not represent convergence"),immediate.=TRUE)),all=FALSE)
      }

      nod<-lapply(unique(c(nodes)),function(x){
        as.numeric(names(L1[,which(colnames(L1)==x)][which(L1[,which(colnames(L1)==x)]!=0)]))[-1]->des
        ldes<-sapply(1:length(des),function(i) length(which(des%in%getMommy(tree,des[i]))))
        des[which(ldes<=1)]->des1
        c(getMommy(tree,x)[1:2],des1,x)->nod
        nod[which(!is.na(nod))]->nod
        nod[which(sapply(nod,function(x) length(tips(tree,x)))>=(length(tips(tree,x))/2))]->nod
        nod[which(sapply(nod,function(x) length(tips(tree,x)))<=(length(tips(tree,x))*1.5))]
      })
      names(nod)<-unique(c(nodes))

      nodpairs<-do.call(rbind,apply(nodes,1,function(w) expand.grid(nod[[match(w[1],names(nod))]],nod[[match(w[2],names(nod))]]),simplify = FALSE))
      split(nodpairs[,2],nodpairs[,1])->nodpairs
      output.sel<-lapply(1:length(nodpairs),function(w){
        nodpairs[[w]][which(!nodpairs[[w]]%in%getDescendants(tree,as.numeric(names(nodpairs)[w])))]->mean.sel
        mean.sel[which(!mean.sel%in%getMommy(tree,as.numeric(names(nodpairs)[w])))]->mean.sel
        distNodes(tree,as.numeric(names(nodpairs)[w]),clus=clus)[1:Nnode(tree),]->distN
        list(mean.sel,distN)
      })
      names(output.sel)<-names(nodpairs)
      RRaces[match(unlist(nod),rownames(RR$aces)),]->aces
    }
{
  # if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else
    #   cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
    # registerDoParallel(cl)
    # res <- foreach(i = 1:length(output.sel),
    #                .packages = c("ape", "phytools", "doParallel","RRphylo")) %dopar%
    #   {
    #     # for(i in 1:(length(nod)-1)){
    #     # print(i)
    #     gc()
    #     # nod[i]->sel1
    #
    #     output.sel[[i]][[1]]->mean.sel
    #     output.sel[[i]][[2]]->distN
    #
    #     as.numeric(names(output.sel)[i])->sel1
    #     tips(tree,sel1)->tt1
    #
    #     if(length(mean.sel)>0){
    #       nDD<-nTT<-array()
    #       mean.diff<-matrix(ncol=6,nrow=length(mean.sel))
    #       for(k in 1:length(mean.sel)){
    #         # print(k)
    #         distN[match(as.numeric(as.character(mean.sel[k])),rownames(distN)),1]->nD->nDD[k]
    #         distN[match(as.numeric(as.character(mean.sel[k])),rownames(distN)),2]->nT->nTT[k]
    #
    #         tips(tree,mean.sel[k])->TT
    #         expand.grid(tt1,TT)->ctt
    #         ang.tip<-mean(sapply(1:nrow(ctt),function(g){
    #           as.matrix(phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),])->ppTT
    #           angle.vecs(ppTT[1,],ppTT[2,])
    #         }))
    #
    #         angle.vecs(aces[which(rownames(aces)%in%sel1),],
    #                    aces[which(rownames(aces)%in%mean.sel[k]),])->ang.ac
    #
    #         c(dir.diff=ang.tip/nT,diff=(ang.tip+ang.ac)/nT,ang=ang.ac,ang.tip=ang.tip,nD=nD,nT=nT)->mean.diff[k,]
    #       }
    #       rownames(mean.diff)<-mean.sel
    #       colnames(mean.diff)<-c("ang.bydist.tip","ang.conv","ang.ace","ang.tip","nod.dist","time.dist")
    #
    #     }else{
    #       c(dir.diff=NULL,diff=NULL,ang=NULL)->mean.diff
    #       nDD<-nTT<-NULL
    #     }
    #     # diff.p->mean.diff
    #
    #     list(mean.diff,mean.sel,nDD,nTT)
    #   }
    #
    # stopCluster(cl)
  }

    core.chunk1<-  expression({

      output.sel[[i]][[1]]->mean.sel
      output.sel[[i]][[2]]->distN

      as.numeric(names(output.sel)[i])->sel1
      tips(tree,sel1)->tt1

      if(length(mean.sel)>0){
        nDD<-nTT<-array()
        mean.diff<-matrix(ncol=6,nrow=length(mean.sel))
        for(k in 1:length(mean.sel)){
          # print(k)
          distN[match(as.numeric(as.character(mean.sel[k])),rownames(distN)),1]->nD->nDD[k]
          distN[match(as.numeric(as.character(mean.sel[k])),rownames(distN)),2]->nT->nTT[k]

          tips(tree,mean.sel[k])->TT
          expand.grid(tt1,TT)->ctt
          ang.tip<-mean(sapply(1:nrow(ctt),function(g){
            as.matrix(phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),])->ppTT
            angle.vecs(ppTT[1,],ppTT[2,])
          }))

          angle.vecs(aces[which(rownames(aces)%in%sel1),],
                     aces[which(rownames(aces)%in%mean.sel[k]),])->ang.ac

          c(dir.diff=ang.tip/nT,diff=(ang.tip+ang.ac)/nT,ang=ang.ac,ang.tip=ang.tip,nD=nD,nT=nT)->mean.diff[k,]
        }
        rownames(mean.diff)<-mean.sel
        colnames(mean.diff)<-c("ang.bydist.tip","ang.conv","ang.ace","ang.tip","nod.dist","time.dist")

      }else{
        c(dir.diff=NULL,diff=NULL,ang=NULL)->mean.diff
        nDD<-nTT<-NULL
      }
      # diff.p->mean.diff

      list(mean.diff,mean.sel,nDD,nTT)
    })

    i=NULL
    res=NULL
    if(round((detectCores() * clus), 0)>1)
      eval(parse(text=paste0(cldef,'\nres<-parLapply(cl=cl,1:length(output.sel),function(i)',
                             core.chunk1,")\n stopCluster(cl)"))) else
                               eval(parse(text=paste0('res<-lapply(1:length(output.sel),function(i)',core.chunk1,")")))

    lapply(res,"[[",1)->mean.diff
    lapply(res,"[[",2)->mean.sel
    lapply(res,"[[",3)->nD.sel
    lapply(res,"[[",4)->nT.sel
    names(mean.diff)<-names(nD.sel)<-names(nT.sel)<-names(output.sel)

    mean.aRd<-array()
    AD<-matrix(ncol=5,nrow=rsim)
    suppressWarnings(
      for(h in 1:rsim){
        # print(h)
        tree$tip.label[sample(seq(1:length(tree$tip.label)),2)]->tsam

        angle.vecs(phen[match(as.character(tsam[1]),rownames(phen)),],
                   phen[match(as.character(tsam[2]),rownames(phen)),])->tt

        getMommy(tree,tsam[1])[1]->n1
        getMommy(tree,tsam[2])[1]->n2

        angle.vecs(phen[match(n1,rownames(phen)),],phen[match(n2,rownames(phen)),])->aa
        dist.nodes(tree)[n1,n2]->dt

        paste(n1,n2,sep="/")->cb
        c(diff=(aa+tt)/dt,dist=dt,n1=n1,n2=n2,combo=cb)->AD[h,]
        tt/dt->mean.aRd[h]
      })

    as.data.frame(AD)->AD
    apply(AD[,1:2],2,function(k) as.numeric(as.character(k)))->AD[,1:2]
    # as.numeric(as.character(AD[,1]))->AD[,1]
    # as.numeric(as.character(AD[,2]))->AD[,2]
    if(any(is.na(AD[,1]))){
      AD[which(!is.na(AD[,1])),]->AD
      mean.aRd[which(!is.na(AD[,1]))]->mean.aRd
    }

    if(any(AD[,1]=="Inf")>0){
      AD[which(!AD[,1]=="Inf"),]->AD
      mean.aRd[which(!AD[,1]=="Inf")]->mean.aRd
    }
    colnames(AD)<-c("ang.by.dist","dist","n1","n2","combo")

    diff.rank<-list()
    for(u in 1:length(mean.diff)){
      names(mean.diff)[[u]]->sel1
      mean.sel[[u]]->msel
      nT.sel[[u]]->ntsel
      diff.pR<-list()
      for(j in 1:length(msel)){
        msel[j]->sel2

        as.numeric(names(L1[,which(colnames(L1)==sel1)][which(L1[,which(colnames(L1)==sel1)]!=0)]))[-1]->des
        ldes<-sapply(1:length(des),function(i) length(which(des%in%getMommy(tree,des[i]))))
        des[which(ldes<=1)]->des1

        as.numeric(names(L1[,which(colnames(L1)==sel2)][which(L1[,which(colnames(L1)==sel2)]!=0)]))[-1]->des
        ldes<-sapply(1:length(des),function(i) length(which(des%in%getMommy(tree,des[i]))))
        des[which(ldes<=1)]->des2

        expand.grid(c(getMommy(tree,sel1)[1:2],des1,sel1),c(getMommy(tree,sel2)[1:2],des2,sel2))->ee
        apply(ee,2,function(i) as.numeric(as.character(i)))->ee

        ee[,c(2,1)]->e2
        colnames(e2)<-colnames(ee)
        rbind(ee,e2)->ex
        apply(ex,1, function(x) paste(x[1],x[2],sep="/"))->exx
        as.data.frame(exx)->exx
        if (any(AD$combo%in%exx[,1])) AD[which(!AD$combo%in%exx[,1]),]->ADD else AD->ADD
        replicate(nsim-1,sample(seq(1,nrow(ADD)),1))->s

        c(rank(c(mean.diff[[u]][j,1],mean.aRd[s]))[1]/nsim,
          rank(c(mean.diff[[u]][j,2],ADD[s,1]))[1]/nsim)->diff.pR[[j]]
      }

      data.frame(mean.diff[[u]],do.call(rbind,diff.pR))->diff.rank[[u]]
      colnames(diff.rank[[u]])[c(7,8)]<-c("p.ang.bydist","p.ang.conv")
      abs(diff.rank[[u]][order(abs(diff.rank[[u]][,1])),])->diff.rank[[u]]
    }

    names(mean.diff)->names(diff.rank)

    diff.rank<-lapply(1:length(diff.rank),function(x){
      dd<-cbind(as.numeric(rownames(diff.rank[[x]])),as.matrix(diff.rank[[x]]))
      rownames(dd)<-rep(names(diff.rank)[x],nrow(dd))
      dd
    })


    diff.rank<-lapply(diff.rank,function(dd){
      if(!is.null(nodes)){
        which(sapply(1:nrow(dd),function(k) paste(rownames(dd)[k],dd[k,1],sep="-"))%in%apply(nodes,1,paste,collapse="-"))->ddpos
        if(length(ddpos)>0){
          dd[ddpos,,drop=FALSE]->ddnpairs
          dd[-ddpos,,drop=FALSE]->dd
        } else ddnpairs<-NULL
      }else ddnpairs<-NULL
      if(any(dd[,8]<=0.05 | dd[,9]<=0.05)) dd[which(dd[,8]<=0.05 | dd[,9]<=0.05),,drop=FALSE]->dd
      do.call(rbind,lapply(RRphylo::node.paths(tree,dd[,1]),function(k)
        dd[match(k,dd[,1]),,drop=FALSE][which.min(dd[match(k,dd[,1]),2]),,drop=FALSE]))->dd
      if(!is.null(ddnpairs)) rbind(ddnpairs,dd) else dd
    })

    do.call(rbind,diff.rank)->sc
    # if(!is.null(nodes))
    #    else
    #     do.call(rbind,lapply(diff.rank,function(k) k[1,,drop=FALSE]))->sc
    colnames(sc)[1]<-"node"
    sc[order(sc[,8],sc[,7]),]->sc
    # sc<<-sc
    ######################################################################
    if(is.null(nodes)) sc[which(sc[,8]<=0.05 | sc[,9]<=0.05),,drop=FALSE]->sc.sel else{
      scind<-unique(c(which(sapply(1:nrow(sc),function(k) paste(rownames(sc)[k],sc[k,1],sep="-"))%in%apply(nodes,1,paste,collapse="-")),
                      which(sc[,8]<=0.05 | sc[,9]<=0.05)))
      sc.sel<-sc[scind,,drop=FALSE]
    }

    if(any(sc.sel[,9]<=0.05)|any(sc.sel[,8]<=0.05)){
      if(any(sc.sel[,9]<=0.05)) print("parallel and convergent trajectories") else {
        if(any(sc.sel[,4]/sc.sel[,5]>1.1)) print("convergent trajectories") else
          print("parallel and convergent trajectories")
      }
    }else print("no candidate node pairs selected, no convergence found")

    ### phenotypic vectors descending from node pairs resulting from search.conv ###
    phen2<-list()
    for(i in 1:nrow(sc.sel)){
      c(rownames(sc.sel)[i],getDescendants(tree,rownames(sc.sel)[i]),
        as.numeric(as.character(sc.sel[i,1])),getDescendants(tree,as.numeric(as.character(sc.sel[i,1]))))->d2
      d2[which(as.numeric(d2)<(Ntip(tree)+1))]<-tree$tip.label[as.numeric(d2[which(as.numeric(d2)<(Ntip(tree)+1))])]
      c(phen2,list(phen[which(rownames(phen)%in%d2),]))->phen2
    }

    sapply(phen2,nrow)->len
    unlist(lapply(1:length(len),function(i) rep(paste(rownames(sc.sel)[i],sc.sel[i,1],sep="/"),len[i])))->gr

    ### perform betadisper and TukeyHSD ###
    dist(do.call(rbind,phen2))->dis
    vegan::betadisper(dis,gr)->bd
    data.frame(bd$group,dist=bd$distance)->M

    data.frame(gr=paste(rownames(sc.sel),sc.sel[,1],sep="/"),dtime=sc.sel[,7])->xdist
    M[,2]/xdist[match(M[,1],xdist[,1]),2]->M[,3]
    tapply(M[,3],M[,1],mean)->mean.dist

    if(nrow(sc.sel)>1&any(sc.sel[,8]<=0.05|sc.sel[,9]<=0.05)){
      MM<-M[which(M[,1]%in%sapply(which(sc.sel[,8]<=0.05|sc.sel[,9]<=0.05),function(x) paste(rownames(sc.sel)[x],sc.sel[x,1],sep="/"))),,drop=FALSE]
      TukeyHSD(aov(MM[,3]~MM[,1]))[[1]]->BD
    }else NULL->BD

    cbind(sc.sel,
          clade.size.n1=sapply(as.numeric(rownames(sc.sel)),function(x) length(tips(tree,x))),
          clade.size.n2=sapply(as.numeric(as.character(sc.sel[,1])),function(x) length(tips(tree,x))))->sc.sel
    list(sc.sel,BD,mean.dist)->res.tot
    names(res.tot)<-c("node pairs","node pairs comparison","average distance from group centroids")

  }else{

    if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
      tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
      tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
    }

    y <- as.matrix(treedataMatch(tree, y)[[1]])

    state[match(rownames(y),names(state))]->state
    if("nostate"%in%state) state[which(state!="nostate")]->state.real else state->state.real
    if(max(table(state.real))>Ntip(tree)*0.5)
      warning("One or more states apply to a large portion of the tree,
              this might be inappropriate for testing convergence",.immediate=TRUE)


    if(length(unique(state.real))>1)
      combn(unique(state.real),2)->stcomb1 else
        combn(rep(unique(state.real),2),2)->stcomb1

    state2decl<-lapply(1:ncol(stcomb1),function(k){
      if(any(table(state[which(state%in%stcomb1[,k])])<4)) NULL else{
        replace(state,which(state%in%stcomb1[,k]),"A")->f
        replace(f,which(!state%in%stcomb1[,k]),"B")->f
        f
      }
    })


    if(declust){
      if(ncol(stcomb1)==1&!("nostate"%in%state))
        decl<-rep(list(NULL),ncol(stcomb1)) else {
          decl<-list()
          for(k in 1:ncol(stcomb1)){
            if(!is.null(state2decl[[k]])){
              setTimeLimit(elapsed=60,transient=TRUE)
              try(repeat({
                phyloclust(tree,state2decl[[k]],"A")$declusterized$removed.tips->dec
                if(!is.null(dec)){
                  if(all(table(state[which(!names(state)%in%dec)])>3)&
                     length(table(state))==length(table(state[which(!names(state)%in%dec)]))) break
                }else break
              }),silent=TRUE)->rep.try
              if(inherits(rep.try,"try-error")) decl[[k]]<-list(NULL) else decl[[k]]<-dec
              setTimeLimit(elapsed=Inf)
            }else decl[k]<-list(NULL)
          }
        }
    }else{
      sapply(state2decl,function(x) ifelse(!is.null(x),phylo.run.test(tree,x,"A")$p,1))->prt
      if(any(prt<=0.05)){
        if(length(unique(state.real))>1){
          prnames<-sapply(which(prt<=0.05),function(q) paste(stcomb1[,q],collapse=" and "))
        }else prnames<-unique(state.real)

        for(q in prnames)
          print(paste("Species in state",prnames,
                      "are phylogenetically clustered; consider running search.conv with declust=TRUE"))
      }
      decl<-rep(list(NULL),ncol(stcomb1))
    }

    cophenetic.phylo(tree)->cop
    ang.by.state<-matrix(ncol=2,nrow=ncol(stcomb1))
    vs<-list()
    for(i in 1:ncol(stcomb1)){
      y[which(!rownames(y)%in%decl[[i]]),]->yd
      state[which(!names(state)%in%decl[[i]])]->stated

      tt1<-lapply(stcomb1[,i], function(w) yd[which(stated==w),])
      vs[[i]]<-sapply(tt1,function(w) mean(apply(w,1,unitV)))

      if(length(unique(stcomb1[,i]))==1) t(combn(rownames(tt1[[1]]),2))->ctt else
        expand.grid(rownames(tt1[[1]]),rownames(tt1[[2]]), stringsAsFactors = FALSE)->ctt
      ang.tip<-sapply(1:nrow(ctt),function(g){
        as.matrix(yd[match(as.character(ctt[g,1:2]),rownames(yd)),])->ppTT
        angle.vecs(ppTT[1,],ppTT[2,])
      })
      dt<-apply(ctt,1,function(w) cop[match(w[1],rownames(cop)),match(w[2],rownames(cop))])
      c(mean(ang.tip),mean(ang.tip/dt))->ang.by.state[i,]
    }
    data.frame(state1=stcomb1[1,],state2=stcomb1[2,],
               ang.state=ang.by.state[,1],ang.state.time=ang.by.state[,2])->ang2state

    {
      # if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else
      #   cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
      # registerDoParallel(cl)
      # ang2stateR <- foreach(j = 1:(nsim-1)) %dopar%
      #   {
      #     gc()
      #     ang.by.stateR<-matrix(ncol=2,nrow=ncol(stcomb1))
      #     for(i in 1:ncol(stcomb1)){
      #       y[which(!rownames(y)%in%decl[[i]]),]->yR
      #       state[which(!names(state)%in%decl[[i]])]->stateR
      #       sample(stateR)->stateR
      #       names(stateR)<-names(yR)
      #
      #       tt1R<-lapply(stcomb1[,i], function(w) yR[which(stateR==w),])
      #       if(length(unique(stcomb1[,i]))==1) t(combn(rownames(tt1R[[1]]),2))->cttR else
      #         expand.grid(rownames(tt1R[[1]]),rownames(tt1R[[2]]), stringsAsFactors = FALSE)->cttR
      #       ang.tipR<-sapply(1:nrow(cttR),function(g){
      #         as.matrix(yR[match(as.character(cttR[g,1:2]),rownames(yR)),])->ppTT
      #         angle.vecs(ppTT[1,],ppTT[2,])
      #       })
      #       dtR<-apply(cttR,1,function(w) cop[match(w[1],rownames(cop)),match(w[2],rownames(cop))])
      #       c(mean(ang.tipR),mean(ang.tipR/dtR))->ang.by.stateR[i,]
      #     }
      #     data.frame(state1=stcomb1[1,],state2=stcomb1[2,],ang.state=ang.by.stateR[,1],
      #                ang.state.time=ang.by.stateR[,2])
      #   }
      # stopCluster(cl)
}

    core.chunk2<-expression({
      ang.by.stateR<-matrix(ncol=2,nrow=ncol(stcomb1))
      for(i in 1:ncol(stcomb1)){
        y[which(!rownames(y)%in%decl[[i]]),]->yR
        state[which(!names(state)%in%decl[[i]])]->stateR
        sample(stateR)->stateR
        names(stateR)<-names(yR)

        tt1R<-lapply(stcomb1[,i], function(w) yR[which(stateR==w),])
        if(length(unique(stcomb1[,i]))==1) t(combn(rownames(tt1R[[1]]),2))->cttR else
          expand.grid(rownames(tt1R[[1]]),rownames(tt1R[[2]]), stringsAsFactors = FALSE)->cttR
        ang.tipR<-sapply(1:nrow(cttR),function(g){
          as.matrix(yR[match(as.character(cttR[g,1:2]),rownames(yR)),])->ppTT
          angle.vecs(ppTT[1,],ppTT[2,])
        })
        dtR<-apply(cttR,1,function(w) cop[match(w[1],rownames(cop)),match(w[2],rownames(cop))])
        c(mean(ang.tipR),mean(ang.tipR/dtR))->ang.by.stateR[i,]
      }
      data.frame(state1=stcomb1[1,],state2=stcomb1[2,],ang.state=ang.by.stateR[,1],
                 ang.state.time=ang.by.stateR[,2])
    })


    j=NULL
    ang2stateR=NULL
    if(round((detectCores() * clus), 0)>1)
      eval(parse(text=paste0(cldef,'\nang2stateR<-parLapply(cl=cl,1:(nsim-1),function(j)',
                             core.chunk2,")\n stopCluster(cl)"))) else
                               eval(parse(text=paste0('ang2stateR<-lapply(1:(nsim-1),function(j)',core.chunk2,")")))

    pval<-rlim<-matrix(ncol=2,nrow=nrow(ang2state))
    for(k in 1:nrow(ang2state)){
      apply(rbind(ang2state[k,3:4],do.call(rbind,lapply(ang2stateR,function(x) x[k,3:4]))),2,rank)[1,]/nsim->pval[k,]
      quantile(do.call(rbind,lapply(ang2stateR,function(x) x[k,]))[,3],c(0.05,0.95))->rlim[k,]
    }
    data.frame(ang2state,p.ang.state=pval[,1],p.ang.state.time=pval[,2])->res.tot

    if(length(unique(state.real))>1){
      prnames<-sapply(which(res.tot[,6]<=0.05),function(q)
        paste("Convergent trajectories between",res.tot[q,1],"and",res.tot[q,2]))
      prnames<-paste(prnames[which(!sapply(prnames,is.null))],collapse = "\n")
    }else if(res.tot[,6]<=0.05) paste("Species within group", unique(state.real),"converge")->prnames else NULL->prnames
    if(!is.null(prnames)) cat(prnames,"\n")


    #### Plot preparation ####
    reflim<-rlim[,1]*res.tot[,4]/res.tot[1,4]
    data.frame(ang=res.tot[,3],do.call(rbind,vs))->bbb
    data.frame(bbb,ang.l1=bbb[,1]/2,ang.l2=360-(bbb[,1]/2),CI5=reflim/2,
               CI95=360-reflim/2)->plotData

    if(nrow(res.tot)==1&&res.tot[,1]==res.tot[,2]) res.tot[,2]<-NULL

    list(state.res=res.tot,plotData=plotData)->res.tot

  }

  class(res.tot)<-c("RRphyloList","list")
  attr(res.tot,"hidden")<-"plotData"
  attr(res.tot,"Call")<-funcall

  return(res.tot)
}
