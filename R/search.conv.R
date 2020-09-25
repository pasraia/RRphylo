#' @title Searching for morphological convergence among species and clades
#' @description The function scans a phylogenetic tree looking for morphological
#'   convergence between entire clades or species evolving under specific
#'   states.
#' @usage search.conv(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
#'   min.dim=NULL,max.dim=NULL,min.dist=NULL,PGLSf=FALSE,declust=FALSE,nsim=1000,rsim=1000,
#'    clus=.5,foldername=NULL)
#' @param RR an object produced by \code{\link{RRphylo}}. This is not indicated
#'   if convergence among states is tested.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous. This is not indicated if convergence among clades is
#'   tested.
#' @param y a multivariate phenotype. The object \code{y} should be either a
#'   matrix or dataframe with species names as rownames.
#' @param nodes node pair to be tested. If unspecified, the function
#'   automatically searches for convergence among clades. Notice the node number
#'   must refer to the dichotomic version of the original tree, as produced by
#'   \code{RRphylo}.
#' @param state the named vector of tip states. The function tests for
#'   convergence within a single state or among different states (this latter
#'   case is especially meant to test for iterative evolution as for example the
#'   appearance of repeated morphotypes into different clades). In both cases,
#'   the state for non-focal species (i.e. not belonging to any convergent
#'   group) must be indicated as "nostate".
#' @param aceV phenotypic values at internal nodes. The object \code{aceV}
#'   should be either a matrix or dataframe with nodes (referred to the
#'   dichotomic version of the original tree, as produced by \code{RRphylo}) as
#'   rownames. If \code{aceV} are not indicated, ancestral phenotypes are
#'   estimated via \code{RRphylo}.
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
#'   clades separated by more than 8 nodes has to be tested, the argument
#'   \code{min.dist} should be "node8". If left unspecified, it automatically
#'   searches for convergence between clades separated by a number of nodes
#'   bigger than one tenth of the tree size.
#' @param PGLSf has been deprecated; please see the argument \code{declust}
#'   instead.
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
#' @param foldername the path of the folder where plots are to be found.
#' @export
#' @seealso \href{../doc/search.conv.html}{\code{search.conv} vignette}
#' @importFrom grDevices chull
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
#'   the function returns a dataframe including for each pair of states (or
#'   single state): \itemize{ \item ang.state: the mean theta angle between
#'   species belonging to different states (or within a single state). \item
#'   ang.state.time: the mean of theta angle between species belonging to
#'   different states (or within a single state) divided by time distance. \item
#'   p.ang.state: the p-value computed for ang.state. \item p.ang.state.time:
#'   the p-value computed for ang.state.time. }
#' @details Regardless the case (either 'state' or 'clade'), the function stores
#'   a plot into the folder specified by \code{foldername}. If convergence among
#'   clades is tested, the clade pair plotted corresponds to those clades with
#'   the smallest \code{$average distance from group centroid}. The figure shows
#'   the Euclidean distances computed between the MRCAs of the clades and the
#'   mean Euclidean distance computed between all the tips belonging to the
#'   converging clades, as compared to the distribution of these same figures
#'   across the rest of the tree. Furthermore, the function stores the PC1/PC2
#'   plot obtained by PCA of the species phenotypes. Convergent clades are
#'   indicated by colored convex hulls. Large colored dots represent the mean
#'   phenotypes per clade (i.e. their group centroids). Eventually, a modified
#'   traitgram plot is produced, highlighting the branches of the clades found
#'   to converge. In both PCA and traitgram, asterisks represent the ancestral
#'   phenotypes of the individual clades. If convergence among states is tested,
#'   the function produces a PC plot with colored convex hulls enclosing species
#'   belonging to different states. Furthermore, it generates circular plots of
#'   the mean angle between states (blue lines) and the range of random angles
#'   (gray shaded area). The p-value for the convergence test is printed within
#'   the circular plots.
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
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",
#'             foldername = tempdir(),clus=cc)
#' # by setting min.dist as time distance
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="time38",
#'             foldername = tempdir(),clus=cc)
#'
#' ## Case 2. searching convergence within a single state
#' search.conv(tree=treefel, y=PCscoresfel, state=statefel,declust=TRUE,
#'             foldername = tempdir(),clus=cc)
#'   }

search.conv<-function(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
                      min.dim=NULL,max.dim=NULL,min.dist=NULL,PGLSf=FALSE,
                      declust=FALSE,nsim=1000,rsim=1000,clus=.5,foldername=NULL)
{
  # require(ape)
  # require(geiger)
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

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package \"cluster\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("plotrix", quietly = TRUE)) {
    stop("Package \"plotrix\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!missing(PGLSf)){
    warning("argument PGLSf is deprecated; please use declust instead.",
            call. = FALSE)
    PGLSf->declust
  }

  phylo.run.test<-function(tree,state,st,nsim=100){
    cophenetic.phylo(tree)->cop
    cop[which(state==st),which(state==st)]->subcop
    mean(subcop[upper.tri(subcop)])->mds
    dim(subcop)[1]->sl

    r.mds<-array()
    for(e in 1:nsim){
      sample(tree$tip.label,sl)->test.tip
      cop[test.tip,test.tip]->r.cop
      mean(r.cop[upper.tri(r.cop)])->r.mds[e]
    }
    return(list(p=length(which(r.mds<mds))/nsim))
  }
  declusterize<-function(tree,state,st){
    remT<-c()
    length(which(state==st))->lenst
    while(phylo.run.test(tree=tree,state=state,st=st)$p<0.05){

      names(state)->nam
      as.numeric(as.factor(state))->rt
      #as.numeric(rt)->rt
      #names(rt)<-nam
      data.frame(state,rt)->def
      unique(def[which(def[,1]==st),2])->stN
      data.frame(V=rle(rt)$values,L=rle(rt)$lengths)->VL
      data.frame(VL,pos=cumsum(VL[,2]))->VL
      #subset(VL,V==stN)->vel
      VL[which(VL$V==stN),]->vel
      vel[which.max(vel$L),]->hit
      tree$tip.label[(hit[,3]-hit[,2]+1):hit[,3]]->hitips
      sample(hitips,0.5*length(hitips))
      sample(hitips,ceiling(.5*length(hitips)))->remtips
      drop.tip(tree,remtips)->tree
      state[-which(names(state)%in%remtips)]->state
      if(length(c(remtips,remT))>(lenst-3)) break else c(remtips,remT)->remT
    }
    if(length(remT)==0) remT<-NA
    #return(list(tree=tree,state=state))
    return(remT=remT)
  }


  if(is.null(state)){

    if(is.null(RR)) stop("missing RRphylo object")

    RR$tree->tree1
    RR$aces->RRaces
    if(is.null(aceV)==FALSE){
      if (inherits(aceV,"data.frame")) aceV <- as.matrix(aceV)
      RRaces[match(rownames(aceV),rownames(RRaces)),]<-aceV
    }

    #if (inherits(y,"data.frame"))
    y <- treedata(tree1, y, sort = TRUE)[[2]]

    RR$tip.path->L
    RR$rates->betas
    RR$node.path->L1
    rbind(RRaces,y)->phen

    if(is.null(nodes)){
      if(is.null(min.dim)) min.dim<-Ntip(tree1)*0.1 else min.dim<-min.dim

      if(is.null(min.dist)) {
        min.dist<-Ntip(tree1)*0.1
        dist.type<-"node"
      }else{
        if(any(grep("time",min.dist))){
          min.dist<-as.numeric(gsub("time","",min.dist))
          dist.type<-"time"
        }else{
          min.dist<-as.numeric(gsub("node","",min.dist))
          dist.type<-"node"
        }
      }

      if(is.null(max.dim))  max.dim<-Ntip(tree1)/5 else max.dim<-max.dim

    }else{
      distNodes(tree1,nodes)->matDist
      if(matDist[,1]<Ntip(tree1)*0.15) print("clades are close to each other. Similarity might not represent convergence")
      if(is.null(min.dim)) min.dim<-min(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2])))) else min.dim<-min.dim
      if(is.null(min.dist)) min.dist<-matDist[,1] else min.dist<-min.dist
      if(is.null(max.dim)){
        if (max(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2]))))>min.dim*2) max.dim<-max(c(length(tips(tree1,nodes[1])),length(tips(tree1,nodes[2])))) else max.dim<-min.dim*2
      } else max.dim<-max.dim
    }



    if(is.null(nodes)){

      subtrees(tree1)->subT->subTN
      if(length(sapply(subT,Ntip)>max.dim)>0){
        subT[-c(which(sapply(subT,Ntip)<min.dim),which(sapply(subT,Ntip)>max.dim))]->subT

      }else{
        subT[-which(unlist(lapply(subT,Ntip))<min.dim)]->subT
      }

      sapply(subT,function(x) getMRCA(tree1,x$tip.label))->names(subT)
      as.numeric(names(subT))->nod
      c(nod,Ntip(tree1)+1)->nod
      subT[[length(subT)+1]]<-tree1
      names(subT)[length(subT)]<-Ntip(tree1)+1

      RRaces[match(nod,rownames(RR$aces)),]->aces

      res<-list()
      if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1) else cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      res <- foreach(i = 1:(length(nod)-1),
                     .packages = c("ape", "geiger", "phytools", "doParallel")) %dopar%
        {
          gc()
          nod[i]->sel1
          tips(tree1,sel1)->tt1

          if(length(which(nod%in%getDescendants(tree1,sel1)))>0)
            nod[-which(nod%in%getDescendants(tree1,sel1))]->mean.sel else
              nod->mean.sel

          if(length(which(mean.sel%in%getMommy(tree1,sel1)))>0)
            mean.sel[-which(mean.sel%in%getMommy(tree1,sel1))]->mean.sel else
              mean.sel->mean.sel
          mean.sel[-which(mean.sel==sel1)]->mean.sel

          if(dist.type=="time") {
            distNodes(tree1,sel1)[1:Nnode(tree1),]->distN
            distN[,2]->matDist
            distN[,1]->matN
          }else{
            distNodes(tree1,sel1)[1:Nnode(tree1),1]->matDist
          }

          matDist->matNod
          if (length(which(mean.sel %in% names(matDist[which(matDist<min.dist)])))>0) mean.sel[-which(mean.sel %in% names(matDist[which(matDist<min.dist)]))]->mean.sel else mean.sel->mean.sel


          if(length(mean.sel)==0) {
            c(dir.diff=NULL,diff=NULL,ang=NULL)->diff.p
            nDD<-NULL
            nTT<-NULL
          }else{
            nDD<-array()
            nTT<-array()
            diff.p<-matrix(ncol=6,nrow=length(mean.sel))
            for(k in 1:length(mean.sel)){

              if(dist.type=="time"){
                matN[match(as.numeric(as.character(mean.sel[k])),names(matN))]->nD
                matDist[match(as.numeric(as.character(mean.sel[k])),names(matDist))]->nT
              }else{
                matDist[match(as.numeric(as.character(mean.sel[k])),names(matDist))]->nD
                dist.nodes(tree1)[sel1,as.numeric(mean.sel[k])]->nT
              }

              nD->nDD[k]
              nT->nTT[k]

              tips(tree1,mean.sel[k])->TT
              expand.grid(tt1,TT)->ctt
              aa<-array()
              for(g in 1:dim(ctt)[1]){
                phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),]->ppTT
                as.matrix(ppTT)->ppTT
                aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
              }
              mean(aa)->ang.tip


              aces[which(rownames(aces)%in%c(sel1,mean.sel[k])),]->ac
              rad2deg(acos((ac[1,] %*% ac[2,])/(unitV(ac[1,]) *unitV(ac[2,]))))->ang.ac

              c(dir.diff=ang.tip/nT,diff=(ang.tip+ang.ac)/nT,ang=ang.ac,ang.tip=ang.tip,nD=nD,nT=nT)->diff.p[k,]


            }
            rownames(diff.p)<-mean.sel
            colnames(diff.p)<-c("ang.bydist.tip","ang.conv","ang.ace","ang.tip","nod.dist","time.dist")

          }

          diff.p->mean.diff

          res[[i]]<-list(matNod,mean.diff,mean.sel,nDD,nTT)
        }

      stopCluster(cl)



      lapply(res,"[[",1)->matNod
      lapply(res,"[[",2)->mean.diff
      lapply(res,"[[",3)->mean.sel
      lapply(res,"[[",4)->nD.sel
      lapply(res,"[[",5)->nT.sel


      nod[1:(length(nod)-1)]->names(mean.diff)->names(matNod)->names(nD.sel)->names(nT.sel)

      sapply(mean.diff,is.null)->nulls
      if(length(which(nulls == TRUE))==length(mean.diff)) stop("There are no nodes distant enough to search convergence, consider using a smaller min.dist")

      if(length(which(nulls==TRUE)>0)) mean.diff[-which(nulls==TRUE)]->mean.diff
      if(length(which(nulls==TRUE)>0)) mean.sel[-which(nulls==TRUE)]->mean.sel
      if(length(which(nulls==TRUE)>0)) matNod[-which(nulls==TRUE)]->matNod
      if(length(which(nulls==TRUE)>0)) nD.sel[-which(nulls==TRUE)]->nD.sel
      if(length(which(nulls==TRUE)>0)) nT.sel[-which(nulls==TRUE)]->nT.sel


      mean.aRd<-array()
      AD<-matrix(ncol=5,nrow=rsim)
      suppressWarnings(for(h in 1:rsim){
        tree1$tip.label[sample(seq(1:length(tree1$tip.label)),2)]->tsam

        phen[match(c(as.character(tsam[1]),as.character(tsam[2])),rownames(phen)),]->ppt
        tt <- rad2deg(acos((ppt[1,]%*%ppt[2,])/(unitV(ppt[1,]) *unitV(ppt[2,]))))
        getMommy(tree1,which(tree1$tip.label==as.character(tsam[1])))[1]->n1
        getMommy(tree1,which(tree1$tip.label==as.character(tsam[2])))[1]->n2
        phen[match(c(n1,n2),rownames(phen)),]->pp
        aa <- rad2deg(acos((pp[1,]%*%pp[2,])/(unitV(pp[1,]) *unitV(pp[2,]))))

        dist.nodes(tree1)[n1,n2]->dt


        paste(n1,n2,sep="/")->cb
        c(diff=(aa+tt)/dt,dist=dt,n1=n1,n2=n2,combo=cb)->AD[h,]
        tt/dt->mean.aRd[h]
      })


      as.data.frame(AD)->AD
      as.numeric(as.character(AD[,1]))->AD[,1]
      as.numeric(as.character(AD[,2]))->AD[,2]
      if(length(which(is.na(AD[,1])))>0){
        which(is.na(AD[,1]))->outs
        AD[-outs,]->AD
        mean.aRd[-outs]->mean.aRd
      }

      if(length(which(AD[,1]=="Inf"))>0){
        which(AD[,1]=="Inf")->outs
        AD[-outs,]->AD
        mean.aRd[-outs]->mean.aRd
      }

      colnames(AD)<-c("ang.by.dist","dist","n1","n2","combo")

      res.ran <- list()
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      res.ran <- foreach(k = 1:nsim,
                         .packages = c("ape", "geiger", "phytools", "doParallel")) %dopar%
        {
          gc()

          mean.diffR<-list()
          for(u in 1:length(mean.diff)){
            names(mean.diff)[[u]]->sel1
            mean.sel[[u]]->msel
            nT.sel[[u]]->ntsel
            diff.pR<-list()
            for(j in 1:length(msel)){
              msel[j]->sel2

              as.numeric(names(L1[,which(colnames(L1)==sel1)][which(L1[,which(colnames(L1)==sel1)]!=0)]))[-1]->des
              ldes<-array()
              for(i in 1:length(des)){
                length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
              }
              des[which(ldes<=1)]->des1

              as.numeric(names(L1[,which(colnames(L1)==sel2)][which(L1[,which(colnames(L1)==sel2)]!=0)]))[-1]->des
              ldes<-array()
              for(i in 1:length(des)){
                length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
              }
              des[which(ldes<=1)]->des2

              expand.grid(c(getMommy(tree1,sel1)[1:2],des1,sel1),c(getMommy(tree1,sel2)[1:2],des2,sel2))->ee
              as.numeric(as.character(ee[,2]))->ee[,2]
              as.numeric(as.character(ee[,1]))->ee[,1]
              ee[,c(2,1)]->e2
              colnames(e2)<-colnames(ee)
              rbind(ee,e2)->ex
              apply(ex,1, function(x) paste(x[1],x[2],sep="/"))->exx
              as.data.frame(exx)->exx
              if (length(which(AD$combo%in%exx[,1]))>0) AD[-which(AD$combo%in%exx[,1]),]->ADD else AD->ADD
              sample(seq(1:dim(ADD)[1]),1)->s
              mean.aRd[s]->rdiff
              ADD[s,1]->mdiff
              data.frame(dir.diff=rdiff,diff=mdiff)->diff.pR[[j]]

            }
            do.call(rbind,diff.pR)->mean.diffR[[u]]
            rownames(mean.diffR[[u]])<-names(msel)
          }
          names(mean.diffR)<-names(mean.diff)
          mean.diffR->res.ran[[k]]

        }
      stopCluster(cl)



      diff.rank<-list()
      for(i in 1:length(mean.diff)){
        rnk<-matrix(ncol=2,nrow=dim(mean.diff[[i]])[1])
        for(k in 1:dim(mean.diff[[i]])[1]){
          apply(rbind(mean.diff[[i]][k,c(1,2)],do.call(rbind,lapply(lapply(res.ran,"[[",i),
                                                                    function(x) x[k,c(1,2)]))[1:(nsim-1),]),2, function(x) rank(x)[1]/nsim )->rnk[k,]

        }
        data.frame(mean.diff[[i]],rnk)->diff.rank[[i]]
        colnames(diff.rank[[i]])[c(7,8)]<-c("p.ang.bydist","p.ang.conv")

      }
      names(mean.diff)->names(diff.rank)

      lapply(diff.rank,function(x) abs(x[order(abs(x[,1])),]))->diff.rank


      as.data.frame(do.call(rbind,lapply(diff.rank,function(x) cbind(rownames(x)[1],x[1,]))))->df
      colnames(df)[1]<-"node"
      df[order(df[,8],df[,7]),]->df

      ######################################################################

      df->sc
      sc[which(sc$p.ang.bydist<=0.05 | sc$p.ang.conv<=0.05),]->sc.sel
      if(nrow(sc.sel)==0){
        print("no candidate node pairs selected, no convergence found")
        sc->sc.sel

        i=1
        while(i<nrow(sc.sel)){
          rbind(data.frame(n1=rownames(sc.sel)[i],n2=sc.sel[i,1]),data.frame(n1=sc.sel[,1],n2=rownames(sc.sel)))->dat
          which(duplicated(dat))-1->out
          if(length(out>0))sc.sel[-out,]->sc.sel
          i<-i+1
        }

        ### phenotypic vectors descending from node pairs resulting from search.conv ###
        phen2<-list()
        for(i in 1:nrow(sc.sel)){
          c(rownames(sc.sel)[i],getDescendants(tree1,rownames(sc.sel)[i]),as.numeric(as.character(sc.sel[i,1])),getDescendants(tree1,as.numeric(as.character(sc.sel[i,1]))))->d2
          d2[which(as.numeric(d2)<(Ntip(tree1)+1))]<-tree1$tip.label[as.numeric(d2[which(as.numeric(d2)<(Ntip(tree1)+1))])]
          phen[which(rownames(phen)%in%d2),]->phen2[[i]]
        }

        sapply(phen2,nrow)->len
        gg<-list()
        for(i in 1:length(len)) rep(paste(rownames(sc.sel)[i],sc.sel[i,1],sep="/"),len[i])->gg[[i]]
        unlist(gg)->gr

        ### perform betadisper and TukeyHSD ###
        dist(do.call(rbind,phen2))->dis
        vegan::betadisper(dis,gr)->bd
        data.frame(bd$group,dist=bd$distance)->M
        sc.sel->x
        data.frame(gr=paste(rownames(x),x[,1],sep="/"),ndist=x[,5],dtime=x[,6])->xdist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist

        data.frame(sc.sel,
                   clade.size.n1=sapply(as.numeric(rownames(sc.sel)),function(x) length(tips(tree1,x))),
                   clade.size.n2=sapply(as.numeric(as.character(sc.sel[,1])),function(x) length(tips(tree1,x))))->sc.sel
        list(sc.sel,mean.dist)->res.tot
        names(res.tot)<-c("node pairs","average distance from group centroids")

      }else{

        ### remove duplicated couples
        i=1
        while(i<nrow(sc.sel)){
          rbind(data.frame(n1=rownames(sc.sel)[i],n2=sc.sel[i,1]),data.frame(n1=sc.sel[,1],n2=rownames(sc.sel)))->dat
          which(duplicated(dat))-1->out
          if(length(out>0))sc.sel[-out,]->sc.sel
          i<-i+1
        }

        if(length(which(sc.sel[,9]<=0.05))>0) print("parallel and convergent trajectories") else {
          if(length(which(sc.sel[,4]/sc.sel[,5]>1.1))>0) print("convergent trajectories") else print("parallel and convergent trajectories")
        }


        ### phenotypic vectors descending from node pairs resulting from search.conv ###
        phen2<-list()
        for(i in 1:nrow(sc.sel)){
          c(rownames(sc.sel)[i],getDescendants(tree1,rownames(sc.sel)[i]),as.numeric(as.character(sc.sel[i,1])),getDescendants(tree1,as.numeric(as.character(sc.sel[i,1]))))->d2
          d2[which(as.numeric(d2)<(Ntip(tree1)+1))]<-tree1$tip.label[as.numeric(d2[which(as.numeric(d2)<(Ntip(tree1)+1))])]
          phen[which(rownames(phen)%in%d2),]->phen2[[i]]
        }

        sapply(phen2,nrow)->len
        gg<-list()
        for(i in 1:length(len)) rep(paste(rownames(sc.sel)[i],sc.sel[i,1],sep="/"),len[i])->gg[[i]]
        unlist(gg)->gr

        ### perform betadisper and TukeyHSD ###
        dist(do.call(rbind,phen2))->dis
        vegan::betadisper(dis,gr)->bd
        data.frame(bd$group,dist=bd$distance)->M
        sc.sel->x
        data.frame(gr=paste(rownames(x),x[,1],sep="/"),ndist=x[,5],dtime=x[,6])->xdist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist

        if(dim(x)[1]>1)
        {
          TukeyHSD(aov(M[,3]~M[,1]))[[1]]->BD
        } else {
          paste("no comparison is performed")->BD

        }

        data.frame(sc.sel,
                   clade.size.n1=sapply(as.numeric(rownames(sc.sel)),function(x) length(tips(tree1,x))),
                   clade.size.n2=sapply(as.numeric(as.character(sc.sel[,1])),function(x) length(tips(tree1,x))))->sc.sel


        list(sc.sel,BD,mean.dist)->res.tot
        names(res.tot)<-c("node pairs","node pairs comparison","average distance from group centroids")

      }
      #####################################################################
    }else{

      nodes[1]->sel1
      nodes[2]->sel2

      as.numeric(names(L1[,which(colnames(L1)==sel1)][which(L1[,which(colnames(L1)==sel1)]!=0)]))[-1]->des
      ldes<-array()
      for(i in 1:length(des)){
        length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
      }
      des[which(ldes<=1)]->des1

      as.numeric(names(L1[,which(colnames(L1)==sel2)][which(L1[,which(colnames(L1)==sel2)]!=0)]))[-1]->des
      ldes<-array()
      for(i in 1:length(des)){
        length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
      }
      des[which(ldes<=1)]->des2

      c(getMommy(tree1,sel1)[1:2],des1,sel1)->nod
      if(length(which(is.na(nod)))>0) nod[-which(is.na(nod))]->nod
      if(length(which(sapply(nod,function(x) length(tips(tree1,x)))<(length(tips(tree1,sel1))/2)))>0) nod[-which(sapply(nod,function(x) length(tips(tree1,x)))<(length(tips(tree1,sel1))/2))]->nod
      if(length(which(sapply(nod,function(x) length(tips(tree1,x)))>(length(tips(tree1,sel1))*1.5)))>0) nod[-which(sapply(nod,function(x) length(tips(tree1,x)))>(length(tips(tree1,sel1))*1.5))]->nod

      c(getMommy(tree1,sel2)[1:2],des2,sel2)->mean.sel
      if(length(which(is.na(mean.sel)))>0) mean.sel[-which(is.na(mean.sel))]->mean.sel
      if(length(which(sapply(mean.sel,function(x) length(tips(tree1,x)))<(length(tips(tree1,sel2))/2)))>0) mean.sel[-which(sapply(mean.sel,function(x) length(tips(tree1,x)))<(length(tips(tree1,sel2))/2))]->mean.sel
      if(length(which(sapply(mean.sel,function(x) length(tips(tree1,x)))>(length(tips(tree1,sel2))*1.5)))>0) mean.sel[-which(sapply(mean.sel,function(x) length(tips(tree1,x)))>(length(tips(tree1,sel2))*1.5))]->mean.sel

      RRaces[match(c(nod,mean.sel),rownames(RR$aces)),]->aces

      res<-list()
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      res <- foreach(i = 1:length(nod),
                     .packages = c("ape", "geiger", "phytools", "doParallel")) %dopar%
        {

          gc()
          nod[i]->sel1
          tips(tree1,sel1)->tt1

          distNodes(tree1,sel1)[1:Nnode(tree1),1]->matDist

          matDist->matNod

          if(length(mean.sel)==0) {
            c(dir.diff=NULL,diff=NULL,ang=NULL)->diff.p
            nDD<-NULL
            nTT<-NULL
          }else{
            nDD<-array()
            nTT<-array()
            diff.p<-matrix(ncol=6,nrow=length(mean.sel))
            for(k in 1:length(mean.sel)){
              matDist[match(as.numeric(as.character(mean.sel[k])),names(matDist))]->nD
              dist.nodes(tree1)[sel1,as.numeric(mean.sel[k])]->nT

              nD->nDD[k]
              nT->nTT[k]


              tips(tree1,mean.sel[k])->TT
              expand.grid(tt1,TT)->ctt
              aa<-array()
              for(g in 1:dim(ctt)[1]){
                phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),]->ppTT
                as.matrix(ppTT)->ppTT
                aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
              }
              mean(aa)->ang.tip


              aces[which(rownames(aces)%in%c(sel1,mean.sel[k])),]->ac
              rad2deg(acos((ac[1,] %*% ac[2,])/(unitV(ac[1,]) *unitV(ac[2,]))))->ang.ac

              c(dir.diff=ang.tip/nT,diff=(ang.tip+ang.ac)/nT,ang=ang.ac,ang.tip=ang.tip,nD=nD,nT=nT)->diff.p[k,]


            }
            rownames(diff.p)<-mean.sel
            colnames(diff.p)<-c("ang.bydist.tip","ang.conv","ang.ace","ang.tip","nod.dist","time.dist")

          }

          diff.p->mean.diff

          res[[i]]<-list(matNod,mean.diff,nDD,nTT)
        }

      stopCluster(cl)

      lapply(res,"[[",1)->matNod
      lapply(res,"[[",2)->mean.diff
      lapply(res,"[[",3)->nD.sel
      lapply(res,"[[",4)->nT.sel


      nod[1:length(nod)]->names(mean.diff)->names(matNod)->names(nD.sel)->names(nT.sel)

      sapply(mean.diff,is.null)->nulls

      if(length(which(nulls==TRUE)>0)) mean.diff[-which(nulls==TRUE)]->mean.diff
      if(length(which(nulls==TRUE)>0)) mean.sel[-which(nulls==TRUE)]->mean.sel
      if(length(which(nulls==TRUE)>0)) matNod[-which(nulls==TRUE)]->matNod
      if(length(which(nulls==TRUE)>0)) nD.sel[-which(nulls==TRUE)]->nD.sel
      if(length(which(nulls==TRUE)>0)) nT.sel[-which(nulls==TRUE)]->nT.sel


      mean.aRd<-array()
      AD<-matrix(ncol=5,nrow=rsim)
      suppressWarnings(for(h in 1:rsim){
        tree1$tip.label[sample(seq(1:length(tree1$tip.label)),2)]->tsam

        phen[match(c(as.character(tsam[1]),as.character(tsam[2])),rownames(phen)),]->ppt
        tt <- rad2deg(acos((ppt[1,]%*%ppt[2,])/(unitV(ppt[1,]) *unitV(ppt[2,]))))
        getMommy(tree1,which(tree1$tip.label==as.character(tsam[1])))[1]->n1
        getMommy(tree1,which(tree1$tip.label==as.character(tsam[2])))[1]->n2
        phen[match(c(n1,n2),rownames(phen)),]->pp
        aa <- rad2deg(acos((pp[1,]%*%pp[2,])/(unitV(pp[1,]) *unitV(pp[2,]))))

        dist.nodes(tree1)[n1,n2]->dt


        paste(n1,n2,sep="/")->cb
        c(diff=(aa+tt)/dt,dist=dt,n1=n1,n2=n2,combo=cb)->AD[h,]
        tt/dt->mean.aRd[h]
      })


      as.data.frame(AD)->AD
      as.numeric(as.character(AD[,1]))->AD[,1]
      as.numeric(as.character(AD[,2]))->AD[,2]
      if(length(which(is.na(AD[,1])))>0){
        which(is.na(AD[,1]))->outs
        AD[-outs,]->AD
        mean.aRd[-outs]->mean.aRd
      }

      if(length(which(AD[,1]=="Inf"))>0){
        which(AD[,1]=="Inf")->outs
        AD[-outs,]->AD
        mean.aRd[-outs]->mean.aRd
      }


      colnames(AD)<-c("ang.by.dist","dist","n1","n2","combo")


      res.ran <- list()
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      res.ran <- foreach(k = 1:nsim,
                         .packages = c("ape", "geiger", "phytools", "doParallel")) %dopar%
        {
          gc()

          mean.diffR<-list()
          for(u in 1:length(mean.diff)){
            names(mean.diff)[[u]]->sel1
            mean.sel->msel
            nT.sel[[u]]->ntsel
            diff.pR<-list()
            for(j in 1:length(msel)){
              msel[j]->sel2

              as.numeric(names(L1[,which(colnames(L1)==sel1)][which(L1[,which(colnames(L1)==sel1)]!=0)]))[-1]->des
              ldes<-array()
              for(i in 1:length(des)){
                length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
              }
              des[which(ldes<=1)]->des1

              as.numeric(names(L1[,which(colnames(L1)==sel2)][which(L1[,which(colnames(L1)==sel2)]!=0)]))[-1]->des
              ldes<-array()
              for(i in 1:length(des)){
                length(which(des%in%getMommy(tree1,des[i])))->ldes[i]
              }
              des[which(ldes<=1)]->des2

              expand.grid(c(getMommy(tree1,sel1)[1:2],des1,sel1),c(getMommy(tree1,sel2)[1:2],des2,sel2))->ee
              as.numeric(as.character(ee[,2]))->ee[,2]
              as.numeric(as.character(ee[,1]))->ee[,1]
              ee[,c(2,1)]->e2
              colnames(e2)<-colnames(ee)
              rbind(ee,e2)->ex
              apply(ex,1, function(x) paste(x[1],x[2],sep="/"))->exx
              as.data.frame(exx)->exx
              if (length(which(AD$combo%in%exx[,1]))>0) AD[-which(AD$combo%in%exx[,1]),]->ADD else AD->ADD
              sample(seq(1:dim(ADD)[1]),1)->s
              mean.aRd[s]->rdiff
              ADD[s,1]->mdiff
              data.frame(dir.diff=rdiff,diff=mdiff)->diff.pR[[j]]

            }
            do.call(rbind,diff.pR)->mean.diffR[[u]]
            rownames(mean.diffR[[u]])<-names(msel)
          }
          names(mean.diffR)<-names(mean.diff)
          mean.diffR->res.ran[[k]]

        }
      stopCluster(cl)



      diff.rank<-list()
      for(i in 1:length(mean.diff)){
        rnk<-matrix(ncol=2,nrow=dim(mean.diff[[i]])[1])
        for(k in 1:dim(mean.diff[[i]])[1]){
          apply(rbind(mean.diff[[i]][k,c(1,2)],do.call(rbind,lapply(lapply(res.ran,"[[",i),
                                                                    function(x) x[k,c(1,2)]))[1:(nsim-1),]),2, function(x) rank(x)[1]/nsim )->rnk[k,]

        }
        data.frame(mean.diff[[i]],rnk)->diff.rank[[i]]
        colnames(diff.rank[[i]])[c(7,8)]<-c("p.ang.bydist","p.ang.conv")

      }
      names(mean.diff)->names(diff.rank)

      lapply(diff.rank,function(x) abs(x[order(abs(x[,1])),]))->diff.rank


      as.data.frame(do.call(rbind,lapply(diff.rank,function(x) cbind(rownames(x)[1],x[1,]))))->df
      colnames(df)[1]<-"node"
      df[order(df[,8],df[,7]),]->df


      df->sc
      sc[which(sc$p.ang.bydist<=0.05 | sc$p.ang.conv<=0.05),]->sc.sel

      if(nrow(sc.sel)==0){
        print("no candidate node pairs selected, no convergence found")
        sc->sc.sel

        i=1
        while(i<nrow(sc.sel)){
          rbind(data.frame(n1=rownames(sc.sel)[i],n2=sc.sel[i,1]),data.frame(n1=sc.sel[,1],n2=rownames(sc.sel)))->dat
          which(duplicated(dat))-1->out
          if(length(out>0))sc.sel[-out,]->sc.sel
          i<-i+1
        }

        ### phenotypic vectors descending from node pairs resulting from search.conv ###
        phen2<-list()
        for(i in 1:nrow(sc.sel)){
          c(rownames(sc.sel)[i],getDescendants(tree1,rownames(sc.sel)[i]),as.numeric(as.character(sc.sel[i,1])),getDescendants(tree1,as.numeric(as.character(sc.sel[i,1]))))->d2
          d2[which(as.numeric(d2)<(Ntip(tree1)+1))]<-tree1$tip.label[as.numeric(d2[which(as.numeric(d2)<(Ntip(tree1)+1))])]
          phen[which(rownames(phen)%in%d2),]->phen2[[i]]
        }

        sapply(phen2,nrow)->len
        gg<-list()
        for(i in 1:length(len)) rep(paste(rownames(sc.sel)[i],sc.sel[i,1],sep="/"),len[i])->gg[[i]]
        unlist(gg)->gr

        ### perform betadisper and TukeyHSD ###
        dist(do.call(rbind,phen2))->dis
        vegan::betadisper(dis,gr)->bd
        data.frame(bd$group,dist=bd$distance)->M
        sc.sel->x
        data.frame(gr=paste(rownames(x),x[,1],sep="/"),ndist=x[,5],dtime=x[,6])->xdist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist

        data.frame(sc.sel,
                   clade.size.n1=sapply(as.numeric(rownames(sc.sel)),function(x) length(tips(tree1,x))),
                   clade.size.n2=sapply(as.numeric(as.character(sc.sel[,1])),function(x) length(tips(tree1,x))))->sc.sel
        list(sc.sel,mean.dist)->res.tot
        names(res.tot)<-c("node pairs","average distance from group centroids")

      }else{

        ### remove duplicated couples
        i=1
        while(i<nrow(sc.sel)){
          rbind(data.frame(n1=rownames(sc.sel)[i],n2=sc.sel[i,1]),data.frame(n1=sc.sel[,1],n2=rownames(sc.sel)))->dat
          which(duplicated(dat))-1->out
          if(length(out>0))sc.sel[-out,]->sc.sel
          i<-i+1
        }

        if(length(which(sc.sel[,9]<=0.05))>0) print("parallel and convergent trajectories") else {
          if(length(which(sc.sel[,4]/sc.sel[,5]>1.1))>0) print("convergent trajectories") else print("parallel and convergent trajectories")
        }



        ### phenotypic vectors descending from node pairs resulting from search.conv ###
        phen2<-list()
        for(i in 1:nrow(sc.sel)){
          c(rownames(sc.sel)[i],getDescendants(tree1,rownames(sc.sel)[i]),as.numeric(as.character(sc.sel[i,1])),getDescendants(tree1,as.numeric(as.character(sc.sel[i,1]))))->d2
          d2[which(as.numeric(d2)<(Ntip(tree1)+1))]<-tree1$tip.label[as.numeric(d2[which(as.numeric(d2)<(Ntip(tree1)+1))])]
          phen[which(rownames(phen)%in%d2),]->phen2[[i]]
        }

        sapply(phen2,nrow)->len
        gg<-list()
        for(i in 1:length(len)) rep(paste(rownames(sc.sel)[i],sc.sel[i,1],sep="/"),len[i])->gg[[i]]
        unlist(gg)->gr

        ### perform betadisper and TukeyHSD ###
        dist(do.call(rbind,phen2))->dis
        vegan::betadisper(dis,gr)->bd
        data.frame(bd$group,dist=bd$distance)->M
        sc.sel->x
        data.frame(gr=paste(rownames(x),x[,1],sep="/"),ndist=x[,5],dtime=x[,6])->xdist
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        M[,2]/xdist[match(M[,1],xdist[,1]),3]->M[,3]
        tapply(M[,3],M[,1],mean)->mean.dist

        if(dim(x)[1]>1)
        {
          TukeyHSD(aov(M[,3]~M[,1]))[[1]]->BD
        } else {
          paste("no comparison is performed")->BD

        }

        data.frame(sc.sel,
                   clade.size.n1=sapply(as.numeric(rownames(sc.sel)),function(x) length(tips(tree1,x))),
                   clade.size.n2=sapply(as.numeric(as.character(sc.sel[,1])),function(x) length(tips(tree1,x))))->sc.sel


        list(sc.sel,BD,mean.dist)->res.tot
        names(res.tot)<-c("node pairs","node pairs comparison","average distance from group centroids")

      }

    }


    #### Plot preparation ####
    aceRR <- (L1 %*% betas[1:Nnode(tree1), ])
    y.multi <- (L %*% betas)

    c(aceRR[,1],y.multi[,1])->phen1
    rep("gray",Ntip(tree1)+Nnode(tree1))->colo
    strsplit(names(which.min(res.tot$`average distance from group centroids`)),"/")[[1]]->nn
    as.numeric(nn)->nn
    c(getMommy(tree1,nn[1]),getDescendants(tree1,nn[1]))->path1
    path1[which(path1<(Ntip(tree1)+1))]<-tree1$tip.label[path1[which(path1<(Ntip(tree1)+1))]]

    c(getMommy(tree1,nn[2]),getDescendants(tree1,nn[2]))->path2
    path2[which(path2<(Ntip(tree1)+1))]<-tree1$tip.label[path2[which(path2<(Ntip(tree1)+1))]]


    colo[which(names(phen1)%in%c(nn[1],path1))]<-"#E7298A"
    colo[which(names(phen1)%in%c(nn[2],path2))]<-"#91003F"

    names(colo)<-names(phen1)
    linwd<-rep(2,length(colo))
    linwd[which(colo!="gray")]<-3
    names(linwd)<-names(colo)

    if(ncol(phen)>nrow(phen)) phen[,1:nrow(phen)]->phen
    princomp(phen)->a
    a$scores[order(a$scores[,1]),]->bb
    suppressWarnings(bb[match(c(nn[1],path1[-which(as.numeric(path1)<=nn[1])]),rownames(bb)),]->a1)
    suppressWarnings(bb[match(c(nn[2],path2[-which(as.numeric(path2)<=nn[2])]),rownames(bb)),]->a2)

    dist(RRaces)->dist.ace
    dist(y)->dist.Y
    as.matrix(dist.Y)->dist.Y
    as.matrix(dist.ace)->dist.ace
    mean(dist.Y[match(tips(tree1,nn[1]),rownames(dist.Y)),match(tips(tree1,nn[2]),colnames(dist.Y))])->dist.Ynn
    dist.ace[match(nn[1],colnames(dist.ace)),][which(names(dist.ace[match(nn[1],colnames(dist.ace)),])%in%nn[2])]->dist.nod


    pdf(file = paste(foldername, "convergence plot.pdf",
                     sep = "/"))

    layout(matrix(c(1,3,2,4),ncol=2,nrow=2, byrow = TRUE),widths = c(1,2))

    par(mar = c(3, 1, 2, 1))
    hist(dist.ace[lower.tri(dist.ace)],xlab="",ylab="",axes=FALSE,main="", col=rgb(58/255,137/255,255/255,1))
    title(main="ace distances",line=0.5,cex.main=1.5)
    axis(1,mgp=c(1,0.5,0))
    abline(v=dist.nod,col=rgb(255/255,213/255,0/255,1),lwd=4)


    par(mar = c(3, 1, 2, 1))
    hist(dist.Y[lower.tri(dist.Y)],main="",xlab="",ylab="",axes=FALSE,col=rgb(58/255,137/255,255/255,1))
    title(main="tip distances",line=0.5,cex.main=1.5)
    axis(1,mgp=c(1,0.5,0))
    abline(v=dist.Ynn,col=rgb(255/255,213/255,0/255,1),lwd=3)

    par(mar = c(3, 2.5, 2, 1))
    plot(bb[-match(nn, rownames(bb)),1:2], ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,pch=21,col="black",bg="gray",asp=1)
    Plot_ConvexHull(xcoord = a1[,1], ycoord = a1[,2], lcolor = "#E7298A",lwd=3,lty=2, col.p = rgb(212/255,185/255,218/255,0.5))
    Plot_ConvexHull(xcoord = a2[,1], ycoord = a2[,2], lcolor = "#91003F",lwd=3,lty=2, col.p = rgb(212/255,185/255,218/255,0.5))
    points(cluster::pam(a1, 1)$medoids,xlim=range(bb[,1]),ylim=range(bb[,2]),bg="#E7298A",type = "p",pch=21,col="black", cex=2.5)
    points(cluster::pam(a2, 1)$medoids,xlim=range(bb[,1]),ylim=range(bb[,2]),bg="#91003F",type = "p",pch=21,col="black", cex=2.5)
    points(a1[1,1],a1[1,2],xlim=range(bb[,1]),ylim=range(bb[,2]),col="black",type = "p",pch=8, cex=1.5,lwd=2)
    points(a2[1,1],a2[1,2],xlim=range(bb[,1]),ylim=range(bb[,2]),col="black",type = "p",pch=8, cex=1.5,lwd=2)



    par(mar = c(2, 2.5, 2, 1))
    suppressWarnings(traitgram(y.multi[,1],tree1,col=colo, method = 'ML',show.names = FALSE,lwd=linwd,mgp=c(0,0.5,0)))->traitres
    points(traitres[match(nn,rownames(traitres)),],cex=1.5,col="black",pch=8,lwd=2)

    dev.off()

  }else{

    if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
      data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
      tree$tip.label<-tips(tree,(Ntip(tree)+1))
      data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
      tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
    }

    tree->tree1
    #if (inherits(y,"data.frame"))
    y <- treedata(tree1, y, sort = TRUE)[[2]]

    state[match(rownames(y),names(state))]->state
    if("nostate"%in%state) state[-which(state=="nostate")]->state.real else state->state.real
    if(sort(table(state.real),decreasing = TRUE)[1]>Ntip(tree1)*0.5) warning("one or more states apply to a large portion of the tree, this might be inappropriate for testing convergence")

    if(("nostate"%in%state)&length(unique(state.real))<2){
      y->yOri
      state->stateOri

      unique(state.real)->onestate

      if(isFALSE(declust)&(length(which(state==onestate))>3))
        if(phylo.run.test(tree1,state,onestate)$p<=0.05) print(paste("species within state",onestate,"are phylogenetically clustered; consider running search.conv with declust=TRUE"))

      if(declust&(length(which(state==onestate))>3)){
        declusterize(tree1,state,onestate)->decl
        decl->remtips
        if(!any(is.na(remtips)>0)){
          y[-match(remtips,rownames(y)),]->y
          state[-match(remtips,names(state))]->state
          #drop.tip(tree1,remtips)->tree1
        }
      }

      cophenetic.phylo(tree1)->cop
      mean(apply(y[which(state==onestate),],1,unitV))->vs1->vs2
      t(combn(names(state[which(state==onestate)]),2))->ctt
      aa<-array()
      dt<-array()
      for(g in 1:dim(ctt)[1]){
        y[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(y)),]->ppTT
        as.matrix(ppTT)->ppTT
        aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
        cop[match(as.character(ctt[g,1]),rownames(cop)),match(as.character(ctt[g,2]),rownames(cop))]->dt[g]
      }
      c(mean(aa),mean(aa/dt),vs1,vs2)->ang.by.state


      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      ang.by.stateR<- foreach(j = 1:nsim,.combine = 'rbind') %dopar%
        {
          gc()
          sample(state,length(which(state==onestate)))->sam
          t(combn(names(sam),2))->ctt
          aa<-array()
          dt<-array()
          for(g in 1:dim(ctt)[1]){
            y[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(y)),]->ppTT
            as.matrix(ppTT)->ppTT
            aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
            cop[match(as.character(ctt[g,1]),rownames(cop)),match(as.character(ctt[g,2]),rownames(cop))]->dt[g]
          }
          c(mean(aa),mean(aa/dt))->ang.by.stateR
        }
      stopCluster(cl)


      apply(rbind(ang.by.state[1:2],ang.by.stateR[1:(nsim-1),]),2,rank)[1,]/nsim->pval
      quantile(ang.by.stateR[,1],c(0.05,0.95))->rlim

      c(ang.state=ang.by.state[1],ang.state.time=ang.by.state[2],ang.by.state[c(3,4)],p.ang.state=pval[1],p.ang.state.time=pval[2])->res.tot

      if(res.tot[6]<=0.05) print(paste("species within group", onestate, "converge"))

      #### Plot preparation ####
      yOri->y
      stateOri->state
      c(res.tot[1],rlim=rlim[1],res.tot[3:4])->bbb
      c(bbb,l1=bbb[1]/2,l2=360-(bbb[1]/2),rlim1=bbb[2]/2,
        rlim2=360-bbb[2]/2,p=res.tot[6])->ccc



      pdf(file = paste(foldername, "convergence plot.pdf",
                       sep = "/"))

      mat<-matrix(c(1,2),ncol=1,nrow=2,byrow=TRUE)
      ht<-c(1,1)

      layout(mat,heights = ht)

      c("nostate",names(sort(table(state.real), decreasing = TRUE)))->statetoplot
      if(ncol(y)>nrow(y)) y[,1:nrow(y)]->ypc else y->ypc
      princomp(ypc)->compy
      compy$scores[,1:2]->sco
      #y[,1:2]->sco
      sco[match(names(state),rownames(sco)),]->sco
      RColorBrewer::brewer.pal(3,"Set2")[1]->cols
      par(mar=c(2.5,2.5,1,1))
      plot(sco, ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,xlim=c(range(sco[,1])[1]*1.1,range(sco[,1])[2]),col="white",asp=1)
      for(h in 1:length(statetoplot)){
        if(h==1){
          points(sco[which(state==statetoplot[h]),],pch=21,bg="gray",cex=1.5)
          Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = "#bebebe",lwd=3,lty=2, col.p = paste("#bebebe","4D",sep=""))
        }else{
          points(sco[which(state==statetoplot[h]),],pch=21,bg=cols,cex=1.5)
          Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = cols,lwd=3,lty=2, col.p = paste(cols,"4D",sep=""))
        }
      }
      legend(min(sco[,1])*1.1,max(sco[,2])*1.1, legend = statetoplot,fill = c("#bebebe",cols), bg = rgb(0, 0, 0, 0),
             box.col = rgb(0,0, 0, 0), border = NA, x.intersp = 0.25,y.intersp=0.8)


      lp<-seq(0,340,20)
      lbs<-seq(0,340,20)

      par(cex.axis=.75)
      plotrix::polar.plot(lengths=c(0,mean(unname(as.matrix(ccc)[c(3,4)])),mean(unname(as.matrix(ccc)[c(3,4)]))),
                 polar.pos = c(0,ccc[7],ccc[8]),rp.type="p",line.col=rgb(0,0,0,0.6),
                 poly.col=rgb(127/255,127/255,127/255,0.4),start=90,radial.lim=range(0,max(ccc[c(3,4)])),radial.labels=""
                 ,boxed.radial = FALSE,mar=c(1,1,2,1),label.pos=lp,labels=lbs)
      title(main=onestate,cex.main = 2)
      plotrix::polar.plot(lengths=c(0,ccc[3],ccc[4]),polar.pos = c(0,ccc[5],ccc[6]),
                 line.col="blue",lwd=4,start=90,add=TRUE,radial.labels="",boxed.radial = FALSE)
      text(paste("p-value = ",ccc[9]),x=0,y=-max(ccc[c(3,4)])/2.3,cex=1.5,col="red")

      dev.off()

      res.tot[-c(3,4)]->res.tot
      t(as.data.frame(res.tot))->res.tot
      rownames(res.tot)<-onestate

    }else{

      combn(unique(state.real),2)->stcomb1
      if("nostate"%in%state) combn(unique(state),2)->stcomb else stcomb1->stcomb

      state2decl<-list()
      for(k in 1:ncol(stcomb1)){
        if(any(table(state[which(state%in%stcomb1[,k])])<4)) state2decl[[k]]<-NA else{
          replace(state,which(state%in%stcomb1[,k]),"A")->f
          replace(f,which(!state%in%stcomb1[,k]),"B")->f

          state2decl[[k]]<-f
        }
      }

      if(declust){
        if(ncol(stcomb1)==1&("nostate"%in%state)==FALSE)
          decl<-as.list(rep(NA,ncol(stcomb1))) else {
            decl<-list()
            for(k in 1:ncol(stcomb1)){
              if(any(!is.na(state2decl[[k]]))){
                setTimeLimit(elapsed=60,transient=TRUE)
                try(repeat({
                  declusterize(tree1,state2decl[[k]],"A")->dec
                  if(!any(is.na(dec))){
                    if(all(table(state[-which(names(state)%in%dec)])>3)&
                       length(table(state))==length(table(state[-which(names(state)%in%dec)]))) break
                  }else break
                }),silent=TRUE)->rep.try
                if(class(rep.try)=="try-error") decl[[k]]<-NA else decl[[k]]<-dec
                setTimeLimit(elapsed=Inf)
              }else decl[[k]]<-NA
            }
          }
      }else{
        sapply(state2decl,function(x){
          if(!any(is.na(x))) phylo.run.test(tree1,x,"A")$p->ret else 1->ret
          return(ret)
        })->prt
        if(any(prt<=0.05))
          for(q in which(prt<=0.05) ) print(paste("species in state",
                                                  stcomb1[1,q],"and",
                                                  stcomb1[2,q],
                                                  "are phylogenetically clustered; consider running search.conv with declust=TRUE"))
        decl<-as.list(rep(NA,ncol(stcomb1)))
      }

      y->yOri
      state->stateOri

      cophenetic.phylo(tree1)->cop
      ang.by.state<-matrix(ncol=4,nrow=ncol(stcomb1))
      y.stcomb<-state.stcomb<-list()
      for(i in 1:ncol(stcomb1)){
        if(!any(is.na(decl[[i]]))){
          yOri[-match(decl[[i]],rownames(yOri)),]->y->y.stcomb[[i]]
          stateOri[-match(decl[[i]],names(stateOri))]->state->state.stcomb[[i]]
        }else{
          yOri->y->y.stcomb[[i]]
          stateOri->state->state.stcomb[[i]]
        }

        y[which(state==stcomb1[1,i]),]->tt1
        mean(apply(tt1,1,unitV))->vs1
        y[which(state==stcomb1[2,i]),]->TT
        mean(apply(TT,1,unitV))->vs2
        expand.grid(rownames(tt1),rownames(TT))->ctt
        aa<-array()
        dt<-array()
        for(g in 1:dim(ctt)[1]){
          y[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(y)),]->ppTT
          as.matrix(ppTT)->ppTT
          aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
          cop[match(as.character(ctt[g,1]),rownames(cop)),match(as.character(ctt[g,2]),rownames(cop))]->dt[g]
        }
        c(mean(aa),mean(aa/dt),vs1,vs2)->ang.by.state[i,]
      }
      data.frame(state1=t(stcomb1)[,1],state2=t(stcomb1)[,2],ang.state=ang.by.state[,1],ang.state.time=ang.by.state[,2],size.v1=ang.by.state[,3],size.v2=ang.by.state[,4])->ang2state

      ang2stateR <- list()
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      ang2stateR <- foreach(j = 1:nsim) %dopar%
        {
          gc()
          ang.by.stateR<-matrix(ncol=2,nrow=ncol(stcomb1))
          for(i in 1:ncol(stcomb1)){
            y.stcomb[[i]]->y
            state.stcomb[[i]]->state
            sample(state)->state
            y[which(state==stcomb1[1,i]),]->tt1
            y[which(state==stcomb1[2,i]),]->TT
            expand.grid(rownames(tt1),rownames(TT))->ctt
            aa<-array()
            dt<-array()
            for(g in 1:dim(ctt)[1]){
              y[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(y)),]->ppTT
              as.matrix(ppTT)->ppTT
              aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
              cop[match(as.character(ctt[g,1]),rownames(cop)),match(as.character(ctt[g,2]),rownames(cop))]->dt[g]
            }
            c(mean(aa),mean(aa/dt))->ang.by.stateR[i,]
          }
          data.frame(state1=t(stcomb1)[,1],state2=t(stcomb1)[,2],ang.state=ang.by.stateR[,1],ang.state.time=ang.by.stateR[,2])->ang2stateR[[j]]
        }
      stopCluster(cl)
      pval<-matrix(ncol=2,nrow=nrow(ang2state))
      rlim<-matrix(ncol=2,nrow=nrow(ang2state))
      for(k in 1:nrow(ang2state)){
        apply(rbind(ang2state[k,3:4],do.call(rbind,lapply(ang2stateR,function(x) x[k,3:4]))[1:(nsim-1),]),2,rank)[1,]/nsim->pval[k,]
        quantile(do.call(rbind,lapply(ang2stateR,function(x) x[k,]))[,3],c(0.05,0.95))->rlim[k,]
      }
      data.frame(ang2state,p.ang.state=pval[,1],p.ang.state.time=pval[,2])->res.tot

      for(j in 1:nrow(res.tot)){
        if(res.tot[j,8]<=0.05) print(paste("convergent trajectories between",res.tot[j,1], "and",res.tot[j,2]))
      }

      #### Plot preparation ####
      yOri->y
      stateOri->state
      data.frame(res.tot[,3],rlim=rlim[,1]*res.tot[,4]/res.tot[1,4],res.tot[,5:6])->bbb
      data.frame(bbb,l1=bbb[,1]/2,l2=360-(bbb[,1]/2),rlim1=bbb[,2]/2,
                 rlim2=360-bbb[,2]/2,p=res.tot[,8])->ccc

      pdf(file = paste(foldername, "convergence plot for different states.pdf",
                       sep = "/"))

      if(nrow(res.tot)==1) {
        mat<-matrix(c(1,2),ncol=1,nrow=2,byrow=TRUE)
        ht<-c(1,1)
      }else{
        if(nrow(res.tot)%%2==0) matrix(ncol=nrow(res.tot)/2,nrow=3)->mat else matrix(ncol=(nrow(res.tot)/2+1),nrow=3)->mat
        mat[1,]<-rep(1,ncol(mat))
        if(nrow(res.tot)%%2==0) mat[2:nrow(mat),]<-seq(2,(nrow(res.tot))+1,1) else mat[2:nrow(mat),]<-seq(2,(nrow(res.tot))+2,1)
        mat->mat
        ht<-(c(1.5,1,1))

      }
      layout(mat,heights = ht)

      if("nostate"%in%state) c("nostate",names(sort(table(state.real), decreasing = TRUE)))->statetoplot else names(sort(table(state.real), decreasing = TRUE))->statetoplot
      if(ncol(y)>nrow(y)) y[,1:nrow(y)]->ypc else y->ypc
      princomp(ypc)->compy
      compy$scores[,1:2]->sco
      #y[,1:2]->sco
      sco[match(names(state),rownames(sco)),]->sco
      RColorBrewer::brewer.pal(length(unique(state)),"Set2")->cols
      par(mar=c(2.5,2.5,1,1))
      plot(sco, ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,xlim=c(range(sco[,1])[1]*1.1,range(sco[,1])[2]),col="white",asp=1)
      for(h in 1:length(statetoplot)){
        if("nostate"%in%statetoplot){
          if(h==1){
            points(sco[which(state==statetoplot[h]),],pch=21,bg="gray",cex=1.5)
            Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = "#bebebe",lwd=3,lty=2, col.p = paste("#bebebe","4D",sep=""))
          }else{
            points(sco[which(state==statetoplot[h]),],pch=21,bg=cols[h-1],cex=1.5)
            Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = cols[h-1],lwd=3,lty=2, col.p = paste(cols[h-1],"4D",sep=""))
          }
        }else{
          points(sco[which(state==statetoplot[h]),],pch=21,bg=cols[h],cex=1.5)
          Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = cols[h],lwd=3,lty=2, col.p = paste(cols[h],"4D",sep=""))
        }
      }
      if("nostate"%in%statetoplot){
        legend(min(sco[,1])*1.1,max(sco[,2])*1.1, legend = statetoplot,fill = c("#bebebe",cols), bg = rgb(0, 0, 0, 0),
               box.col = rgb(0,0, 0, 0), border = NA, x.intersp = 0.25,y.intersp=0.8)
      }else{
        legend(min(sco[,1])*1.1,max(sco[,2])*1.1, legend = statetoplot,fill = cols, bg = rgb(0, 0, 0, 0),
               box.col = rgb(0,0, 0, 0), border = NA, x.intersp = 0.25,y.intersp=0.8)
      }

      if(nrow(res.tot)>6) {
        lp<-seq(0,340,40)
        lbs<-seq(0,340,40)
      }else{
        lp<-seq(0,340,20)
        lbs<-seq(0,340,20)
      }
      for(i in 1:nrow(res.tot)){
        par(cex.axis=.75)
        plotrix::polar.plot(lengths=c(0,mean(unname(as.matrix(ccc)[i,c(3,4)])),mean(unname(as.matrix(ccc)[i,c(3,4)]))),
                   polar.pos = c(0,ccc[i,7],ccc[i,8]),rp.type="p",line.col=rgb(0,0,0,0.6),
                   poly.col=rgb(127/255,127/255,127/255,0.4),start=90,radial.lim=range(0,max(ccc[,c(3,4)])),radial.labels=""
                   ,boxed.radial = FALSE,mar=c(1,1,2,1),label.pos=lp,labels=lbs)
        title(main=paste(as.character(res.tot[i,1]),as.character(res.tot[i,2]),sep="-"),cex.main = 2)
        plotrix::polar.plot(lengths=c(0,ccc[i,3],ccc[i,4]),polar.pos = c(0,ccc$l1[i],ccc$l2[i]),
                   line.col="blue",lwd=4,start=90,add=TRUE,radial.labels="",boxed.radial = FALSE)
        text(paste("p-value = ",ccc[i,9]),x=0,y=-max(ccc[i,c(3,4)])/2.3,cex=1.5,col="red")
      }
      dev.off()

      res.tot[,-c(5,6)]->res.tot
    }
  }
  return(res.tot)
}
