#'@title Graphical representation of search.conv results
#'@description This function generates customized functions to produce plots out of
#'  \code{\link{search.conv}} results.
#'@usage plotConv(SC,y,variable,RR=NULL,state=NULL,aceV=NULL)
#'@param SC an object produced by \code{\link{search.conv}}.
#'@param y the multivariate phenotype used to perform \code{\link{search.conv}}.
#'@param variable the index of result to plot. If convergence between clades is
#'  inspected, this is the position within the \code{SC$'average distance from
#'  group centroids'} vector of the clade pair to be plotted. In the case of
#'  convergence between states, this is the number of the line of
#'  \code{SC$state.res} where results for the state pair are returned.
#'@param RR the object produced by \code{\link{RRphylo}} used to perform
#'  \code{\link{search.conv}}. This is not indicated if convergence between states is
#'  tested.
#'@param state the named vector of tip states used to perform \code{\link{search.conv}}.
#'  This is not indicated if convergence between clades is tested.
#'@param aceV phenotypic values at internal nodes to be provided if used to
#'  perform \code{\link{search.conv}}.
#'@return If convergence between clades was tested, \code{plotConv} returns a
#'  list of four functions:
#'@return \strong{\code{$plotHistTips}} shows the mean Euclidean distance
#'  computed between phenotypic vectors of all the tips belonging to the
#'  converging clades as compared to the distribution of distances between all
#'  possible pair of tips across the rest of the tree. The usage is:
#'  \code{...$plotHistTips(hist.args=NULL,line.args=NULL)}, where
#'  \code{hist.args} is a list of further arguments passed to the function
#'  \code{hist}, and \code{line.args} is a list of further arguments passed to
#'  the function \code{lines}.
#'@return \strong{\code{$plotHistAces}} shows the Euclidean distance computed
#'  between phenotypic vectors of the MRCAs of the converging clades as compared
#'  to the distribution of distances between all possible pairs of nodes across
#'  the rest of the tree. The usage is identical to \code{$plotHistTips}.
#'@return \strong{\code{$plotPChull}} generates a PC1/PC2 plot obtained by
#'  performing a PCA of the species phenotypes. Convergent clades are indicated
#'  by colored convex hulls. Large dots represent the mean phenotypes per clade
#'  (i.e. their group centroids) and asterisks (customizable) represent the
#'  ancestral phenotypes of the individual clades. The usage is:
#'  \code{...$plotPChull(plot.args=NULL,chull.args=NULL,means.args=NULL,}
#'  \code{ace.args=NULL,legend.args=list()},where \code{plot.args} is a list of
#'  further arguments passed to the function \code{plot}, \code{chull.args} is a
#'  list of further arguments passed to the function \code{polygon},
#'  \code{means.args} and \code{ace.args} are lists of further argument passed
#'  to the function \code{points} to customize the dots representing the
#'  centroids and the ancestral phenotypes respectively, and \code{legend.args}
#'  is a list of additional arguments passed to the function \code{legend} (if
#'  \code{= NULL} the legend is not plotted).
#'@return \strong{\code{$plotTraitgram}} produces a modified traitgram plot (see
#'  package \pkg{picante}) highlighting the branches of the clades found to
#'  converge. The usage is:
#'  \code{plotTraitgram(colTree=NULL,colNodes=NULL,...)}, where \code{colTree}
#'  is the color to represent the traitgram lines not pertaining the converging
#'  clades, \code{colNodes} is the color (or the vector of colors) to represent
#'  the traitgram lines pertaining the converging clades, and \code{...} are
#'  further arguments passed to the function \code{plot} to plot lines.
#'@return If convergence between states was tested, \code{plotConv} returns a
#'  list of two functions:
#'@return \strong{\code{$plotPChull}} generates a PC1/PC2 plot obtained by
#'  performing a PCA of the species phenotypes, with colored convex hulls
#'  enclosing species belonging to different states. The usage is:
#'  \code{...$plotPChull(plot.args=NULL,chull.args=NULL,points.args=NULL,}
#'  \code{legend.args=list()}, where \code{plot.args} is a list of further
#'  arguments passed to the function \code{plot}, \code{chull.args} is a list of
#'  further arguments passed to the function \code{polygon}, \code{points.args}
#'  is a list of further argument passed to the function \code{points}, and
#'  \code{legend.args} is a list of additional arguments passed to the function
#'  \code{legend} (if \code{= NULL} the legend is not plotted).
#'@return \strong{\code{$plotPolar}} produces a polar plot of the mean angle
#'  within/between state/s as compared to the 95% confidence interval of random
#'  angles. The usage is:
#'  \code{...$plotPolar(polar.args=NULL,polygon.args=NULL,line.args=NULL)},
#'  where \code{polar.args} is a list of further arguments passed to the
#'  function \code{\link[plotrix]{polar.plot}} to set the plot basics (i.e.
#'  \code{radial.lim}, \code{start}, and so on), \code{polygon.args} is a list
#'  of further arguments passed to the function \code{polar.plot} under the
#'  condition \code{rp.type="p"} (see \code{plotrix::polar.plot} for details) to
#'  set the angles distribution graphics, and \code{line.args} is a list of
#'  further arguments passed to the function \code{polar.plot} under the
#'  condition \code{rp.type="r"} to set the mean angle graphics.
#'@author Silvia Castiglione
#'@export
#'@seealso \href{../doc/search.conv.html}{\code{search.conv} vignette}
#'@seealso \href{../doc/Plotting-tools.html#plotConv}{\code{plotConv} vignette}
#'@references Castiglione, S., Serio, C., Tamagnini, D., Melchionna, M.,
#'  Mondanaro, A., Di Febbraro, M., Profico, A., Piras, P.,Barattolo, F., &
#'  Raia, P. (2019). A new, fast method to search for morphological convergence
#'  with shape data. \emph{PLoS ONE}, 14, e0226949.
#'  https://doi.org/10.1371/journal.pone.0226949
#' @examples
#' \dontrun{
#' data("DataFelids")
#' DataFelids$PCscoresfel->PCscoresfel
#' DataFelids$treefel->treefel
#' DataFelids$statefel->statefel->state2
#' state2[sample(which(statefel=="nostate"),20)]<-"st2"
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(treefel,PCscoresfel,clus=cc)->RRfel
#'
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",clus=cc)->sc.clade
#' plotConv(sc.clade,PCscoresfel,variable=2,RR=RRfel)->pc.clade
#'
#' pc.clade$plotHistTips(hist.args = list(col="gray80",yaxt="n",cex.axis=0.8,cex.main=1.5),
#'                 line.args = list(lwd=3,lty=4,col="purple"))
#' pc.clade$plotHistAces(hist.args = list(col="gray80",cex.axis=0.8,cex.main=1.5),
#'                 line.args = list(lwd=3,lty=4,col="gold"))
#' pc.clade$plotPChull(chull.args = list(border=c("cyan","magenta"),lty=1),
#'               means.args = list(pch=c(23,22),cex=3,bg=c("cyan2","magenta2")),
#'               ace.args=list(pch=9),legend.args = NULL)
#' pc.clade$plotTraitgram(colTree = "gray70",colNodes = c("cyan","magenta"))
#'
#'
#' search.conv(tree=treefel, y=PCscoresfel, state=statefel,declust=TRUE,clus=cc)->sc.state
#' plotConv(sc.state,PCscoresfel,variable=1,state=statefel)->pc.state
#' pc.state$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
#'               points.args = list(pch=c(23,22),bg="gray"),
#'               legend.args = list(pch=c(23,22),x="top"))
#'
#' pc.state$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
#'              polygon.args = list(line.col="green",poly.col=NA,lwd=2),
#'              line.args = list(line.col="deeppink",lty=2,lwd=3))
#'
#'     }

plotConv<-function(SC,y,variable,RR=NULL,state=NULL,aceV=NULL){

  misspacks<-sapply(c("picante","scales","RColorBrewer","cluster","plotrix"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  variable->i
  if(!is.null(SC$`average distance from group centroids`)){
    if(is.null(RR)) stop("Please provide the RR object used to perform search.conv")

    RR$tree->tree1
    RR$aces->RRaces
    RR$tip.path->L
    RR$rates->betas
    RR$node.path->L1
    if(!is.null(aceV)){
      if (inherits(aceV,"data.frame")) aceV <- as.matrix(aceV)
      RRaces[match(rownames(aceV),rownames(RRaces)),]<-aceV
    }
    rbind(RRaces,y)->phen

    as.numeric(strsplit(names(SC$`average distance from group centroids`)[i],"/")[[1]])->nn
    c(getMommy(tree1,nn[1]),getDescendants(tree1,nn[1]))->path1
    path1[which(path1<(Ntip(tree1)+1))]<-tree1$tip.label[path1[which(path1<(Ntip(tree1)+1))]]

    c(getMommy(tree1,nn[2]),getDescendants(tree1,nn[2]))->path2
    path2[which(path2<(Ntip(tree1)+1))]<-tree1$tip.label[path2[which(path2<(Ntip(tree1)+1))]]
    {# Histograms
      dist(RRaces)->dist.ace
      dist(y)->dist.Y
      as.matrix(dist.Y)->dist.Y
      as.matrix(dist.ace)->dist.ace
      mean(dist.Y[match(tips(tree1,nn[1]),rownames(dist.Y)),match(tips(tree1,nn[2]),colnames(dist.Y))])->dist.Ynn
      dist.ace[match(nn[1],colnames(dist.ace)),][which(names(dist.ace[match(nn[1],colnames(dist.ace)),])%in%nn[2])]->dist.nod

      plotHist<-function(dista,vert){
        if(all(colnames(dista)%in%tree1$tip.label))
          mtitle<-"tip distances" else mtitle<-"ace distances"
          plotHistD<-function(hist.args=NULL,line.args=NULL){

            if(!"xlab"%in%names(hist.args)) hist.args$xlab<-""
            if(!"ylab"%in%names(hist.args)) hist.args$ylab<-""
            if(!"main"%in%names(hist.args)) hist.args$main<-mtitle
            if(!"col"%in%names(hist.args)) hist.args$col<-"#3a89ff"

            if(!"col"%in%names(line.args)) line.args$col<-"#ffd500"
            if(!"lwd"%in%names(line.args)) line.args$lwd<-4

            do.call(hist,c(x=list(dista[lower.tri(dista)]),hist.args))
            do.call(abline,c(v=unname(vert),line.args))
          }
          return(plotHistD)
      }

      plotHistAces<-plotHist(dist.ace,dist.nod)
      plotHistTips<-plotHist(dist.Y,dist.Ynn)
    }

    {# PCplot
      if(ncol(phen)>nrow(phen)) phen[,1:nrow(phen)]->phen
      princomp(phen)->a
      a$scores[order(a$scores[,1]),]->bb
      suppressWarnings(bb[match(c(nn[1],path1[-which(as.numeric(path1)<=nn[1])]),rownames(bb)),]->a1)
      suppressWarnings(bb[match(c(nn[2],path2[-which(as.numeric(path2)<=nn[2])]),rownames(bb)),]->a2)

      plotPChull<-function(plot.args=NULL,chull.args=NULL,means.args=NULL,ace.args=NULL,legend.args=list()){
        if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"PC1"
        if(!"ylab"%in%names(plot.args)) plot.args$ylab<-"PC2"
        if(!"cex"%in%names(plot.args)) plot.args$cex<-1.5
        if(!"main"%in%names(plot.args)) plot.args$font.lab<-2
        if(!"pch"%in%names(plot.args)) plot.args$pch<-21
        if(!"col"%in%names(plot.args)) plot.args$col<-"black"
        if(!"bg"%in%names(plot.args)) plot.args$bg<-"gray"

        if(!"border"%in%names(chull.args)) chull.args$border<-c("#E7298A","#91003F") else if(length(chull.args$border)==1) chull.args$border<-rep(chull.args$border,2)
        if(!"lwd"%in%names(chull.args)) chull.args$lwd<-c(3,3) else if(length(chull.args$lwd)==1) chull.args$lwd<-rep(chull.args$lwd,2)
        if(!"lty"%in%names(chull.args)) chull.args$lty<-c(2,2) else if(length(chull.args$lty)==1) chull.args$lty<-rep(chull.args$lty,2)
        if(!"col"%in%names(chull.args)) chull.args$col<-scales::alpha(chull.args$border,0.5) else if(length(chull.args$col)==1) chull.args$col<-rep(chull.args$col,2)

        if(!"cex"%in%names(means.args)) means.args$cex<-2.5
        if(!"pch"%in%names(means.args)) means.args$pch<-21
        if(any(means.args$pch>20)){
          if(!"col"%in%names(means.args)) means.args$col<-"black" else means.args$col[which(means.args$pch>20)]<-"black"
          if(!"bg"%in%names(means.args)) means.args$bg<-chull.args$border
        }else{
          if(!"col"%in%names(means.args)) means.args$col<-chull.args$border
        }

        if(!"cex"%in%names(ace.args)) ace.args$cex<-1.5
        if(!"pch"%in%names(ace.args)) ace.args$pch<-8
        if(!"lwd"%in%names(ace.args)) ace.args$lwd<-2
        if(!"col"%in%names(ace.args)) ace.args$col<-"black" else
          if(any(ace.args$pch>20)){
            rep(ace.args$col,length(ace.args$pch))->colp
            colp[which(ace.args$pch>20)]<-"black"
            ace.args$col<-colp
          }

        if(!is.null(legend.args)){
          if(!"x"%in%names(legend.args)) legend.args$x<-"topleft"
          if(!"legend"%in%names(legend.args)) legend.args$legend<-nn
          if(!"fill"%in%names(legend.args)&!any(c("lty","lwd")%in%names(legend.args))){
            if(!"pch"%in%names(legend.args)) legend.args$pch<-rep(21,2)
            if(any(legend.args$pch>20)&!"pt.bg"%in%names(legend.args)){
              legend.args$pt.bg<-rep(NA,length(legend.args$pch))
              legend.args$pt.bg[which(legend.args$pch>20)]<-chull.args$border[which(legend.args$pch>20)]
            }
            if(any(legend.args$pch<21)&!"col"%in%names(legend.args)){
              legend.args$col<-rep(NA,length(legend.args$pch))
              legend.args$col[which(legend.args$pch<21)]<-chull.args$border[which(legend.args$pch<21)]
              legend.args$col[which(is.na(legend.args$col))]<-"black"
            }
            if(!"pt.cex"%in%names(legend.args)) legend.args$pt.cex<-1.5
          }
          if(!"bty"%in%names(legend.args)) legend.args$bty<-"n"

        }

        list(c(list(xcoord=a1[,1],ycoord=a1[,2]),lapply(chull.args,"[[",1)),
             c(list(xcoord=a2[,1],ycoord=a2[,2]),lapply(chull.args,"[[",2)))->chull.args
        c(list(x=c(cluster::pam(a1, 1)$medoids[,1],cluster::pam(a2, 1)$medoids[,1]),
               y=c(cluster::pam(a1, 1)$medoids[,2],cluster::pam(a2, 1)$medoids[,2])),means.args)->means.args
        c(list(x=c(a1[1,1],a2[1,1]),y=c(a1[1,2],a2[1,2])),ace.args)->ace.args

        do.call(plot,c(list(x=bb[-match(nn, rownames(bb)),1],
                            y=bb[-match(nn, rownames(bb)),2],asp=1),plot.args))
        lapply(chull.args,function(k) do.call(Plot_ConvexHull,k))
        do.call(points,means.args)
        do.call(points,ace.args)
        if(!is.null(legend.args)) do.call(legend,legend.args)
      }
    }

    {# Traitgram
      aceRR <- (L1 %*% betas[1:Nnode(tree1), ])
      y.multi <- (L %*% betas)
      c(aceRR[,1],y.multi[,1])->phen1

      plotTraitgram<-function(colTree=NULL,colNodes=NULL,...){
        if(is.null(colTree)) colTree<-"gray"
        if(is.null(colNodes)) colNodes<-c("#E7298A","#91003F") else if(length(colNodes)==1) colNodes<-rep(colNodes,2)

        rep(colTree,Ntip(tree1)+Nnode(tree1))->colo
        colo[which(names(phen1)%in%c(nn[1],path1))]<-colNodes[1]
        colo[which(names(phen1)%in%c(nn[2],path2))]<-colNodes[2]
        names(colo)<-names(phen1)

        plot.args<-list(...)
        if(!"lwd"%in%names(plot.args)){
          plot.args$lwd<-rep(2,length(colo))
          plot.args$lwd[which(colo!="gray")]<-3
          names(plot.args$lwd)<-names(colo)
        }

        do.call(traitgram,c(list(x=y.multi[,1],phy=tree1,col=colo,
                                 method="pic",show.names=FALSE),plot.args))->traitres
      }
    }
    return(list(plotHistTips=plotHistTips,plotHistAces=plotHistAces,plotPChull=plotPChull,plotTraitgram=plotTraitgram))
  }else{
    if(is.null(state)) stop("Please provide the vector of states used to perform search.conv")

    unique(state[which(state!="nostate")])->state.real
    if("nostate"%in%state) c("nostate",names(sort(table(state.real), decreasing = TRUE)))->statetoplot else
      names(sort(table(state.real), decreasing = TRUE))->statetoplot
    if(ncol(y)>nrow(y)) y[,1:nrow(y)]->ypc else y->ypc
    princomp(ypc)->compy
    compy$scores[,1:2]->sco
    sco[match(names(state),rownames(sco)),]->sco

    plotPChull_state<-function(plot.args=NULL,chull.args=NULL,points.args=NULL,legend.args=list()){
      if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"PC1"
      if(!"ylab"%in%names(plot.args)) plot.args$ylab<-"PC2"
      if(!"cex"%in%names(plot.args)) plot.args$cex<-1.5
      if(!"main"%in%names(plot.args)) plot.args$font.lab<-2
      # if(!"pch"%in%names(plot.args)) plot.args$pch<-21

      chull.args$xcoord<-lapply(statetoplot,function(h) sco[which(state==h),1])
      chull.args$ycoord<-lapply(statetoplot,function(h) sco[which(state==h),2])
      if(!"border"%in%names(chull.args)){
        if("nostate"%in%statetoplot)
          chull.args$border<-c("#bebebe",suppressWarnings(RColorBrewer::brewer.pal(length(unique(state.real)),"Set2"))) else
            chull.args$border<-suppressWarnings(RColorBrewer::brewer.pal(length(unique(state.real)),"Set2"))
      } else if(length(chull.args$border)<length(statetoplot))
        chull.args$border<-rep(chull.args$border,length.out=length(statetoplot))

      if(!"lwd"%in%names(chull.args)) chull.args$lwd<-rep(3,length(statetoplot)) else if(length(chull.args$lwd)<length(statetoplot)) chull.args$lwd<-rep(chull.args$lwd,length.out=length(statetoplot))
      if(!"lty"%in%names(chull.args)) chull.args$lty<-rep(2,length(statetoplot)) else if(length(chull.args$lty)<length(statetoplot)) chull.args$lty<-rep(chull.args$lty,length.out=length(statetoplot))
      if(!"col"%in%names(chull.args)) chull.args$col<-scales::alpha(chull.args$border,0.5) else if(length(chull.args$col)<length(statetoplot)) chull.args$col<-rep(chull.args$col,length.out=length(statetoplot))

      points.args$x<-chull.args$xcoord
      points.args$y<-chull.args$ycoord
      if(!"cex"%in%names(points.args)) points.args$cex<-rep(1.5,length(statetoplot)) else
        if(length(points.args$cex)<length(statetoplot)) points.args$cex<-rep(points.args$cex,length.out=length(statetoplot))
      if(!"pch"%in%names(points.args)) points.args$pch<-rep(21,length(statetoplot)) else
        if(length(points.args$pch)<length(statetoplot)) points.args$pch<-rep(points.args$pch,length.out=length(statetoplot))
      if(any(points.args$pch>20)){
        if(!"col"%in%names(points.args)) points.args$col<-rep("black",length(statetoplot)) else points.args$col[which(points.args$pch>20)]<-"black"
        if(!"bg"%in%names(points.args)) points.args$bg<-chull.args$border else
          if(length(points.args$bg)<length(statetoplot)) points.args$bg<-rep(points.args$bg,length.out=length(statetoplot))
      }else{
        if("bg"%in%names(points.args)) points.args$bg<-NULL
        if(!"col"%in%names(points.args)) points.args$col<-chull.args$border else
          if(length(points.args$col)<length(statetoplot)) points.args$col<-rep(points.args$col,length.out=length(statetoplot))
      }

      if(!is.null(legend.args)){
        if(!"x"%in%names(legend.args)) legend.args$x<-"topleft"
        if(!"legend"%in%names(legend.args)) legend.args$legend<-statetoplot

        if(!"fill"%in%names(legend.args)&!any(c("lty","lwd")%in%names(legend.args))){
          if(!"pch"%in%names(legend.args)) legend.args$pch<-rep(21,length(statetoplot))
          if(any(legend.args$pch>20)&!"pt.bg"%in%names(legend.args)){
            legend.args$pt.bg<-rep(NA,length(legend.args$pch))
            legend.args$pt.bg[which(legend.args$pch>20)]<-chull.args$border[which(legend.args$pch>20)]
          }
          if(any(legend.args$pch<21)&!"col"%in%names(legend.args)){
            legend.args$col<-rep(NA,length(legend.args$pch))
            legend.args$col[which(legend.args$pch<21)]<-chull.args$border[which(legend.args$pch<21)]
            legend.args$col[which(is.na(legend.args$col))]<-"black"
          }
          if(!"pt.cex"%in%names(legend.args)) legend.args$pt.cex<-1.5
        }
        if(!"bty"%in%names(legend.args)) legend.args$bty<-"n"
      }

      do.call(plot,c(list(x=NA,xlim=c(range(sco[,1])),ylim=range(sco[,2]),asp=1),plot.args))
      lapply(1:length(statetoplot),function(w) do.call(points,lapply(points.args,"[[",w)))
      lapply(1:length(statetoplot),function(w) do.call(Plot_ConvexHull,lapply(chull.args,"[[",w)))
      if(!is.null(legend.args)) do.call(legend,legend.args)
    }

    plotPolar<-function(polar.args=NULL,polygon.args=NULL,line.args=NULL){

      if(!"radial.lim"%in%names(polar.args)) polar.args$radial.lim<-range(0,max(SC$plotData[,2:3]))
      if(!"start"%in%names(polar.args)) polar.args$start<-90
      if(!"radial.labels"%in%names(polar.args)) polar.args$radial.labels<-""
      if(!"boxed.radial"%in%names(polar.args)) polar.args$boxed.radial<-FALSE
      if(!"main"%in%names(polar.args)) {
        if(length(state.real)>1)
          polar.args$main<-paste(as.character(SC$state.res[i,1]),as.character(SC$state.res[i,2]),sep="-") else
            polar.args$main<-rownames(SC$state.res)
      }


      if(!"line.col"%in%names(polygon.args)) polygon.args$line.col<-scales::alpha("black",0.6)
      if(!"poly.col"%in%names(polygon.args)) polygon.args$poly.col<-scales::alpha("gray50",0.6)

      if(!"line.col"%in%names(line.args)) line.args$line.col<-"blue"
      if(!"lwd"%in%names(line.args)) line.args$lwd<-4

      do.call(plotrix::polar.plot,c(list(lengths=c(0,rep(mean(unname(as.matrix(SC$plotData)[i,2:3])),2)),
                                         polar.pos=c(0,SC$plotData[i,6],SC$plotData[i,7]),
                                         rp.type="p"),polar.args,polygon.args))

      do.call(plotrix::polar.plot,c(list(lengths=c(0,SC$plotData[i,2],SC$plotData[i,3]),
                                         polar.pos=c(0,SC$plotData$ang.l1[i],SC$plotData$ang.l2[i]),
                                         add=TRUE),polar.args,line.args))
    }
  }
  return(list(plotPChull=plotPChull_state,plotPolar=plotPolar))
}




