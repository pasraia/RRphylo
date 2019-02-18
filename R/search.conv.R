#' @title Searching for morphological convergence among species and clades
#' @description The function scans a phylogenetic tree looking for morphological convergence between entire clades or species evolving under specific states.
#' @usage search.conv(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
#' min.dim=NULL,max.dim=NULL,min.dist=NULL,PGLSf=TRUE,nsim=1000,rsim=1000,
#' clus=.5,foldername=NULL)
#' @param RR an object produced by \code{\link{RRphylo}}. This is not indicated if convergence among states is tested.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous. This is not indicated if convergence among clades is tested.
#' @param y a multivariate phenotype. The object \code{y} should be either a matrix or dataframe with species names as rownames.
#' @param nodes node pair to be tested. If unspecified, the function automatically searches for convergence among clades.
#' @param state the state of the tips. If a single convergent state is indicated species within the state are tested. Otherwise, species belonging to different parts of the tree could be tested for convergence on each other by indicating different states. This latter case is especially meant to test for iterative evolution (i.e. the appearance of repeated morphotypes into different clades). The state for non-focal species (i.e. not belonging to any convergent group) must be indicated as "nostate".
#' @param aceV phenotypic values at internal nodes. The object \code{aceV} should be either a matrix or dataframe with nodes as rownames. If \code{aceV} are not indicated, ancestral phenotypes are estimated via \code{RRphylo}.
#' @param min.dim the minimum size of the clades to be compared. When \code{nodes} is indicated, it indicates the minimum size of the smallest clades in \code{nodes}, otherwise it is set at one tenth of the tree size.
#' @param max.dim the maximum size of the clades to be compared. When \code{nodes} is indicated, it is \code{min.dim}*2 if the largest clade in \code{nodes} is smaller than this value, otherwise it corresponds to the size of the largest clade. Whitout \code{nodes} it is set at one third of the tree size.
#' @param min.dist the minimum distance, in terms of number of nodes, between the clades to be compared. When \code{nodes} is indicated, it indicates the minimum distance between the pair, otherwise it is set at 10 nodes.
#' @param PGLSf if the argument is set to \code{TRUE} (default), the function tests whether states have phylogenetic structure, by running a runs test. If the states are phylogenetically structured, a \code{\link{PGLS_fossil}} is performed using states as the predictor variables. PGLS residuals will be used to test for convergence.
#' @param nsim number of simulations to perform sampling within the theta random distribution. It is set at 1000 by default.
#' @param rsim number of simulations to be performed to produce the random distribution of theta values. It is set at 1000 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @param foldername the path of the folder where plots are to be found.
#' @export
#' @importFrom grDevices chull
#' @importFrom vegan betadisper
#' @importFrom cluster pam
#' @importFrom tseries runs.test
#' @importFrom graphics axis layout lines segments
#' @importFrom stats TukeyHSD aov princomp
#' @importFrom plotrix polar.plot
#' @return If convergence between clades is tested, the function returns a list including:
#' @return \itemize{\item\strong{$node pairs}: a dataframe containing for each pair of nodes:
#' \itemize{\item ang.bydist.tip: the mean theta angle between clades divided by the time distance.
#' \item ang.conv: the mean theta angle between clades plus the angle between aces, divided by the time distance.
#' \item ang.ace: the angle between aces.
#' \item ang.tip: the mean theta angle between clades.
#' \item nod.dist: the distance intervening between clades in terms of number of nodes.
#' \item time.dist: the time distance intervening between the clades.
#' \item p.ang.bydist: the p-value computed for ang.bydist.tip.
#' \item p.ang.conv: the p-value computed for ang.conv.
#' \item clade.size: the size of clades.
#' }
#' \item\strong{$node pairs comparison}: pairwise comparison between significantly convergent pairs (all pairs if no istance of significance was found) performed on the distance from group centroids (the mean phenotype per clade).
#' \item\strong{$average distance from group centroids}: smaller average distances mean less variable phenotypes within the pair.
#' }
#' @return If convergence between (or within a single state) states is tested, the function returns a dataframe including for each pair of states (or single state):
#' \itemize{
#' \item ang.state: the mean theta angle between species belonging to different states (or within a single state).
#' \item ang.state.time: the mean of theta angle between species belonging to different states (or within a single state) divided by time distance.
#' \item p.ang.state: the p-value computed for ang.state.
#' \item p.ang.state.time: the p-value computed for ang.state.time.
#' }
#' @details Regardless the case (either ‘state’ or ‘clade’), the function stores a plot into the folder specified by \code{foldername}. If convergence among clades is tested, the clade pair plotted corresponds to those clades with the smallest \code{$average distance from group centroid}. The figure shows the Euclidean distances computed between the MRCAs of the clades and the mean Euclidean distance computed between all the tips belonging to the converging clades, as compared to the distribution of these same figures across the rest of the tree. Furthermore, the function stores the PC1/PC2 plot obtained by PCA of the species phenotypes. Convergent clades are indicated by colored convex hulls. Large colored dots represent the mean phenotypes per clade (i.e. their group centroids). Eventually, a modified traitgram plot is produced, highlighting the branches of the clades found to converge. In both PCA and traitgram, asterisks represent the ancestral phenotypes of the individual clades. If convergence among states is tested, the function produces a PC plot with colored convex hulls enclosing species belonging to different states. Furthermore, it generates circular plots of the mean angle between states (blue lines) and the range of random angles (gray shaded area). The p- value for the convergence test is printed within the circular plots.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto, Paolo Piras, Davide Tamagnini
#' @examples
#'   data("DataUng")
#'   DataUng$PCscoresung->PCscoresung
#'   DataUng$treeung->treeung
#'   DataUng$stateung->stateung
#'
#'   data("DataFelids")
#'   DataFelids$PCscoresfel->PCscoresfel
#'   DataFelids$treefel->treefel
#'
#'   \donttest{
#'   search.conv(tree=treeung, y=PCscoresung, state=stateung, foldername = tempdir())
#'
#'   RRphylo(treefel,PCscoresfel)->RRfel
#'   search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, foldername = tempdir())
#'   }


search.conv<-function(RR=NULL,tree=NULL,y,nodes=NULL,state=NULL,aceV=NULL,
                      min.dim=NULL,max.dim=NULL,min.dist=NULL,PGLSf=TRUE,
                      nsim=1000,rsim=1000,clus=.5,foldername=NULL)
{
  # require(ape)
  # require(geiger)
  # require(phytools)
  # require(foreach)
  # require(doParallel)
  # require(parallel)
  # require(vegan)
  # require(cluster)
  # require(picante)
  # require(plotrix)
  # require(RColorBrewer)
  # require(tseries)


  traitgram = function(
    x, phy,
    xaxt='s',
    underscore = FALSE,
    show.names = TRUE,
    show.xaxis.values = TRUE,
    method = c('ML','pic'),
    col=NULL,
    lwd=NULL,
    mgp=NULL,...)
  {

    method <- match.arg(method)
    Ntaxa = length(phy$tip.label)
    Ntot = Ntaxa + phy$Nnode
    phy = picante::node.age(phy)
    ages = phy$ages[match(1:Ntot,phy$edge[,2])]
    ages[Ntaxa+1]=0


    if (class(x) %in% c('matrix','array')) {
      xx = as.numeric(x)
      names(xx) = row.names(x)
    } else xx = x

    if (!is.null(names(xx))) {
      umar = 0.1
      if (!all(names(xx) %in% phy$tip.label)) {
        print('trait and phy names do not match')
        return()
      }
      xx = xx[match(phy$tip.label,names(xx))]
    } else umar = 0.1

    lmar = 0.2
    if (xaxt=='s') if (show.xaxis.values) lmar = 1 else lmar = 0.5
    xanc <- ape::ace(xx, phy, method=method)$ace
    xall = c(xx,xanc)

    a0 = ages[phy$edge[,1]]
    a1 = ages[phy$edge[,2]]
    x0 = xall[phy$edge[,1]]
    x1 = xall[phy$edge[,2]]


    if (show.names) {
      maxNameLength = max(nchar(names(xx)))
      ylim = c(min(ages),max(ages)*(1+maxNameLength/50))
      if (!underscore) names(xx) = gsub('_',' ',names(xx))
    } else ylim = range(ages)

    if(is.null(mgp)) magp<-NULL else{
      mgp->magp
    }

    if(is.null(col)) colo<-par("fg") else{
      col->colo
      colo[match(names(x1),names(colo))]->colo
      data.frame(x0,x1,a0,a1,colo)->dato
      dato[order(dato[,5],decreasing=T),]->dato
      dato[,1]->x0
      dato[,2]->x1
      dato[,3]->a0
      dato[,4]->a1
      as.character(dato[,5])->colo
      rownames(dato)->names(x1)->names(x0)->names(a1)->names(a0)->names(colo)
    }

    par(mar = c(3, 2.5, 2, 1))
    plot(range(c(a0,a1)),range(c(x0,x1)),

         type='n',xaxt='n',yaxt='n',
         xlab='',ylab='',bty='n',cex.axis=0.8)
    if (xaxt=='s') if (show.xaxis.values) axis(1,labels=TRUE,mgp=magp) else axis(1,labels=FALSE,mgp=magp)

    if(is.null(lwd)) linwd<-1 else{
      lwd->linwd
      linwd[match(names(x1),names(linwd))]->linwd
    }

    segments(a0,x0,a1,x1,col=colo,lwd=linwd)

    if (show.names) {
      text(max(ages),sort(xx),
           labels = names(xx)[order(xx)],
           adj = -0,
           srt=90,
           cex=.3)
    }

    return(data.frame(a1,x1))
  }


  Plot_ConvexHull<-function(xcoord, ycoord, lcolor,lwd=NULL, lty=NULL,col.p=NULL){
    hpts <- chull(x = xcoord, y = ycoord)
    hpts <- c(hpts, hpts[1])
    lines(xcoord[hpts], ycoord[hpts], col = lcolor,lwd=lwd, lty=lty)
    polygon(xcoord[hpts], ycoord[hpts], col=col.p, border=NA)
  }


  unitV <- function(x) {
    sum(x^2)^0.5
  }
  deg2rad <- function(deg) {
    (deg * pi)/(180)
  }
  rad2deg <- function(rad) {
    (rad * 180)/(pi)
  }


  if(is.null(state)){

    if(is.null(RR)) stop("missing RRphylo object")

    RR$tree->tree1
    RR$aces->RRaces
    if(is.null(aceV)==FALSE){
      if (class(aceV) == "data.frame") aceV <- as.matrix(aceV)
      RRaces[match(rownames(aceV),rownames(RRaces)),]<-aceV
    }
    if (class(y) == "data.frame")
      y <- treedata(tree1, y, sort = TRUE)[[2]]

    RR$tip.path->L
    RR$rates->betas
    RR$node.path->L1
    rbind(RRaces,y)->phen

    if(is.null(nodes)){
      if(is.null(min.dim)) min.dim<-Ntip(tree1)*0.1 else min.dim<-min.dim
      if(is.null(min.dist)) min.dist<-Ntip(tree1)*0.1 else min.dist<-min.dist
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
      unlist(lapply(subT,function(x) getMRCA(tree1,x$tip.label)))->names(subT)
      as.numeric(names(subT))->nodo
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

      res<-list()
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      res <- foreach(i = 1:(length(nod)-1),
                     .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
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

                       distNodes(tree1,sel1)[1:Nnode(tree1),1]->matDist

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
      res.ran <- foreach(i = 1:nsim,
                         .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
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
                           mean.diffR->res.ran[[i]]

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
        betadisper(dis,gr)->bd
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
        betadisper(dis,gr)->bd
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
                     .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
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
      res.ran <- foreach(i = 1:nsim,
                         .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
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
                           mean.diffR->res.ran[[i]]

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
        betadisper(dis,gr)->bd
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
        betadisper(dis,gr)->bd
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
    plot(bb[-match(nn, rownames(bb)),1:2], ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,pch=21,col="black",bg="gray")
    Plot_ConvexHull(xcoord = a1[,1], ycoord = a1[,2], lcolor = "#E7298A",lwd=3,lty=2, col.p = rgb(212/255,185/255,218/255,0.5))
    Plot_ConvexHull(xcoord = a2[,1], ycoord = a2[,2], lcolor = "#91003F",lwd=3,lty=2, col.p = rgb(212/255,185/255,218/255,0.5))
    points(pam(a1, 1)$medoids,xlim=range(bb[,1]),ylim=range(bb[,2]),bg="#E7298A",type = "p",pch=21,col="black", cex=2.5)
    points(pam(a2, 1)$medoids,xlim=range(bb[,1]),ylim=range(bb[,2]),bg="#91003F",type = "p",pch=21,col="black", cex=2.5)
    points(a1[1,1],a1[1,2],xlim=range(bb[,1]),ylim=range(bb[,2]),col="black",type = "p",pch=8, cex=1.5,lwd=2)
    points(a2[1,1],a2[1,2],xlim=range(bb[,1]),ylim=range(bb[,2]),col="black",type = "p",pch=8, cex=1.5,lwd=2)



    par(mar = c(2, 2.5, 2, 1))
    suppressWarnings(traitgram(y.multi[,1],tree1,col=colo, method = 'ML',show.names = FALSE,lwd=linwd,mgp=c(0,0.5,0)))->traitres
    points(traitres[match(nn,rownames(traitres)),],cex=1.5,col="black",pch=8,lwd=2)

    dev.off()

  }else{
    tree->tree1
    if (class(y) == "data.frame")
      y <- treedata(tree1, y, sort = TRUE)[[2]]
    state[match(rownames(y),names(state))]->state
    if(sort(table(state[-which(state=="nostate")]),decreasing = TRUE)[1]>Ntip(tree1)*0.5) warning("one or more states apply to a large portion of the tree, this might be inappropriate for testing convergence")

    if(length(unique(state[-which(state=="nostate")]))<2){
      unique(state[-which(state=="nostate")])->onestate
      if(PGLSf){
        rep("B",length(state))->f
        f[which(state==onestate)]<-"A"
        runs.test(as.factor(f))$p->ranp
      }else{
        ranp<-1
      }


      ape::cophenetic.phylo(tree1)->cop
      if(ranp<=0.05) suppressWarnings(residuals(PGLS_fossil(tree1,state,y))->y)
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

      if(res.tot[4]<=0.05) print(paste("species within group", onestate, "converge"))

      #### Plot preparation ####
      c(res.tot[1],rlim=rlim[1],res.tot[3:4])->bbb
      c(bbb,l1=bbb[1]/2,l2=360-(bbb[1]/2),rlim1=bbb[2]/2,
        rlim2=360-bbb[2]/2,p=res.tot[6])->ccc



      pdf(file = paste(foldername, "convergence plot.pdf",
                       sep = "/"))

      mat<-matrix(c(1,2),ncol=1,nrow=2,byrow=TRUE)
      ht<-c(1,1)

      layout(mat,heights = ht)

      names(sort(table(state),decreasing = TRUE))->statetoplot
      princomp(y)->compy
      compy$scores[,1:2]->sco
      sco[match(names(state),rownames(sco)),]->sco
      brewer.pal(3,"Set2")[1]->cols
      par(mar=c(2.5,2.5,1,1))
      plot(sco, ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,xlim=c(range(sco[,1])[1]*1.1,range(sco[,1])[2]),col="white")
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
      polar.plot(lengths=c(0,mean(unname(as.matrix(ccc)[c(3,4)])),mean(unname(as.matrix(ccc)[c(3,4)]))),
                 polar.pos = c(0,ccc[7],ccc[8]),rp.type="p",line.col=rgb(0,0,0,0.6),
                 poly.col=rgb(127/255,127/255,127/255,0.4),start=90,radial.lim=range(0,max(ccc[c(3,4)])),radial.labels=""
                 ,boxed.radial = F,mar=c(1,1,2,1),label.pos=lp,labels=lbs)
      title(main=onestate,cex.main = 2)
      polar.plot(lengths=c(0,ccc[3],ccc[4]),polar.pos = c(0,ccc[5],ccc[6]),
                 line.col="blue",lwd=4,start=90,add=TRUE,radial.labels="",boxed.radial = F)
      text(paste("p-value = ",ccc[9]),x=0,y=-max(ccc[c(3,4)])/2.3,cex=1.5,col="red")

      dev.off()

      res.tot[-c(3,4)]->res.tot
      t(as.data.frame(res.tot))->res.tot
      rownames(res.tot)<-onestate

    }else{

      combn(unique(state),2)->stcomb
      combn(unique(state)[-which(unique(state)=="nostate")],2)->stcomb1

      if(PGLSf){
        ranp<-array()
        for(k in 1:ncol(stcomb1)){
          f<-array()
          for(i in 1:length(state)) {
            if(state[i]%in%stcomb1[,k]) f[i]<-"A" else f[i]<-"B"
          }
          runs.test(as.factor(f))$p->ranp[k]
        }
      }else{
        ranp<-rep(1,ncol(stcomb1))
      }


      ape::cophenetic.phylo(tree1)->cop
      ang.by.state<-matrix(ncol=4,nrow=ncol(stcomb1))
      for(i in 1:ncol(stcomb1)){
        if(ranp[i]<=0.05) suppressWarnings(residuals(PGLS_fossil(tree1,state,y))->y)
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
        ang.by.stateR<-matrix(ncol=2,nrow=ncol(stcomb))
        for(i in 1:ncol(stcomb)){
          sample(state)->state
          y[which(state==stcomb[1,i]),]->tt1
          y[which(state==stcomb[2,i]),]->TT
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
        data.frame(state1=t(stcomb)[,1],state2=t(stcomb)[,2],ang.state=ang.by.stateR[,1],ang.state.time=ang.by.stateR[,2])->ang2stateR[[j]]
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

      names(sort(table(state),decreasing = TRUE))->statetoplot
      princomp(y)->compy
      compy$scores[,1:2]->sco
      sco[match(names(state),rownames(sco)),]->sco
      brewer.pal(length(unique(state)),"Set2")->cols
      par(mar=c(2.5,2.5,1,1))
      plot(sco, ylab="PC2",xlab="PC1",cex=1.5,mgp=c(1.5,0.5,0),font.lab=2,xlim=c(range(sco[,1])[1]*1.1,range(sco[,1])[2]),col="white")
      for(h in 1:length(statetoplot)){
        if(h==1){
          points(sco[which(state==statetoplot[h]),],pch=21,bg="gray",cex=1.5)
          Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = "#bebebe",lwd=3,lty=2, col.p = paste("#bebebe","4D",sep=""))
        }else{
          points(sco[which(state==statetoplot[h]),],pch=21,bg=cols[h-1],cex=1.5)
          Plot_ConvexHull(xcoord=sco[which(state==statetoplot[h]),1], ycoord=sco[which(state==statetoplot[h]),2], lcolor = cols[h-1],lwd=3,lty=2, col.p = paste(cols[h-1],"4D",sep=""))
        }
      }
      legend(min(sco[,1])*1.1,max(sco[,2])*1.1, legend = statetoplot,fill = c("#bebebe",cols), bg = rgb(0, 0, 0, 0),
             box.col = rgb(0,0, 0, 0), border = NA, x.intersp = 0.25,y.intersp=0.8)

      if(nrow(res.tot)>6) {
        lp<-seq(0,340,40)
        lbs<-seq(0,340,40)
      }else{
        lp<-seq(0,340,20)
        lbs<-seq(0,340,20)
      }
      for(i in 1:nrow(res.tot)){
        par(cex.axis=.75)
        polar.plot(lengths=c(0,mean(unname(as.matrix(ccc)[i,c(3,4)])),mean(unname(as.matrix(ccc)[i,c(3,4)]))),
                   polar.pos = c(0,ccc[i,7],ccc[i,8]),rp.type="p",line.col=rgb(0,0,0,0.6),
                   poly.col=rgb(127/255,127/255,127/255,0.4),start=90,radial.lim=range(0,max(ccc[,c(3,4)])),radial.labels=""
                   ,boxed.radial = F,mar=c(1,1,2,1),label.pos=lp,labels=lbs)
        title(main=paste(as.character(res.tot[i,1]),as.character(res.tot[i,2]),sep="-"),cex.main = 2)
        polar.plot(lengths=c(0,ccc[i,3],ccc[i,4]),polar.pos = c(0,ccc$l1[i],ccc$l2[i]),
                   line.col="blue",lwd=4,start=90,add=TRUE,radial.labels="",boxed.radial = F)
        text(paste("p-value = ",ccc[i,9]),x=0,y=-max(ccc[i,c(3,4)])/2.3,cex=1.5,col="red")
      }
      dev.off()

      res.tot[,-c(5,6)]->res.tot
    }
  }
  return(res.tot)
}


