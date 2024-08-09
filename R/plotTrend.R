#'@title Graphical representation of search.trend results
#'@description This function generates customized functions to produce plots of
#'  phenotype versus time and absolute evolutionary rates versus time
#'  regressions for the entire phylogeny and individual clades. Each custom
#'  function takes as first argument the index or name of the variable (as in
#'  the \code{search.trend} output in \code{$trends.data$phenotypeVStime}) to be
#'  plotted.
#'@usage plotTrend(ST)
#'@param ST an object produced by \code{\link{search.trend}}.
#'@return The function returns a list of functions:
#'@return \strong{\code{$plotSTphen}} returns the plot of rescaled phenotype versus age
#'  regression. The 95\% confidence intervals of slopes produced by regressing
#'  phenotypes simulated under the Brownian motion are plotted as a polygon. The
#'  usage is
#'  \code{...$plotSTphen(variable,plot.args=NULL,polygon.args=NULL,line.args=NULL)},
#'   where \code{variable} is the index or name of the variable to be plotted,
#'  \code{plot.args} is a list of further arguments passed to the function
#'  \code{plot}, \code{polygon.args} is a list of further arguments passed to
#'  the function \code{polygon}, and \code{line.args} is a list of further
#'  arguments passed to the function \code{lines}. The functions automatically
#'  plots white points for internal nodes, red points for tips, a blue
#'  regression line, and a gray shaded polygon to represent the 95\% confidence
#'  intervals of Brownian motion slopes.
#'@return \strong{\code{$plotSTrates}} returns the plot of log rescaled rates versus age
#'  regression. The 95\% confidence intervals of slopes produced by regressing
#'  rates simulated under the Brownian motion are plotted as a shaded area. The
#'  arguments are the same as described for \code{$plotSTphen}. In the case of
#'  multivariate \code{y}, the 2-norm vector of multiple rates (see
#'  \code{\link{RRphylo}} for details) can be plotted by setting \code{variable
#'  = "rate"} or \code{variable = } the number of actual variables + 1.
#'@return \strong{\code{$plotSTphenNode}} returns plots of rescaled phenotype versus age
#'  regression for individual clades. The usage is
#'  \code{...$plotSTphenNode(variable,node,plot.args=NULL,lineTree.args=NULL,}
#'  \code{lineNode.args=NULL,node.palette=NULL)}, where \code{variable} is the
#'  same as \code{plotSTphen} and \code{node} is the vector of indices or numbers
#'  (as inputted to \code{search.trend}, to be indicated as character) of nodes
#'  to be plotted. The function allows up to nine nodes at the same time.
#'  \code{plot.args} is a list of further arguments passed to the function
#'  \code{plot}, including \code{pch.node}, a custom argument to set individual
#'  \code{pch} at nodes (its length can be one or as long as the number of
#'  nodes). \code{lineTree.args} is a list of further arguments passed to the
#'  function \code{lines} used to plot the regression line for the entire tree
#'  \code{lineNode.args} is a list of further arguments passed to the function
#'  \code{lines} used to plot the regression line for individual nodes.
#'  \code{node.palette} is the vector of colors specific to nodes points and
#'  lines. Its length can be one or as long as the number of nodes.
#'@return \strong{\code{$plotSTratesNode}} returns plots of absolute rates versus age
#'  regression for individual clades. The arguments are the same as described
#'  for \code{$plotSTphenNode}. In the case of multivariate \code{y}, the 2-norm
#'  vector of multiple rates (see \code{\link{RRphylo}} for details) can be
#'  plotted by setting \code{variable = "rate"} or \code{variable = } the number
#'  of actual variables + 1.
#'@author Silvia Castiglione, Carmela Serio
#'@importFrom graphics points text title polygon pairs plot
#'@export
#'@seealso \href{../doc/search.trend.html}{\code{search.trend} vignette}
#'@seealso \href{../doc/Plotting-tools.html}{\code{plotTrend} vignette}
#'@references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M.,
#'  Profico, A., Girardi, G., & Raia, P. (2019) Simultaneous detection of
#'  macroevolutionary patterns in phenotypic means and rate of change with and
#'  within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14:
#'  e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @examples
#'  \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' cc<- 2/parallel::detectCores()
#'
#' # Extract Pterosaurs tree and data
#' library(ape)
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' RRphylo(tree=treeptero,y=log(massptero),clus=cc)->RRptero
#'
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, node=143, clus=cc,cov=NULL)->ST
#'
#' plotTrend(ST)->plotST
#'
#' plotST$plotSTphen(1) # to plot phenotypic trend through time for entire tree
#' plotST$plotSTphen("log(massptero)",plot.args=list(cex=1.2,col="blue")) # equivalent to above
#'
#' plotST$plotSTrates(1) # to plot rates trend through time for entire tree
#'
#' # to plot phenotypic trend through time for the clade
#' plotST$plotSTphenNode("log(massptero)",plot.args=list(xlab="Age",main="Phenotypic trend"),
#'                     lineTree.args=list(lty=1,col="pink"),lineNode.args=list(lwd=3),
#'                     node.palette=c("chartreuse"))
#'
#' # to plot rates trend through time for the clade
#' plotST$plotSTratesNode("rate")
#'
#'    }

plotTrend<-function(ST){
  ST$trend.data$phenotypeVStime->phen.plot
  ST$trend.data$absrateVStime->absrate.plot
  ST$trend.data$rescaledrateVStime->resrate.plot
  ST$ConfInts$phenotype->CIphen
  ST$ConfInts$rescaled_rate->CIresrate
  ST$phenotypic.regression->p.phen
  if(any(grepl(".multi",colnames(phen.plot)))){
    phen.plot[,-grep(".multi",colnames(phen.plot)),drop=FALSE]->phen.plot
    ncol(phen.plot)-1->ny
    ny+1->nyrate
  } else ncol(phen.plot)-1->ny->nyrate


  apply(phen.plot[,1:ny,drop=FALSE],2,range01)->phen.plot[,1:ny]
  abs(absrate.plot[,1:nyrate,drop=FALSE])->absrate.plot[,1:nyrate]
  phen.plot$age<-max(phen.plot$age)-phen.plot$age

  if(!is.null(absrate.plot$group)){
    phen.plot$group<-absrate.plot$group[match(rownames(phen.plot),rownames(absrate.plot))]
    if("others"%in%unique(absrate.plot$group))
      unique(absrate.plot$group)[-which(unique(absrate.plot$group)=="others")]->groups
    absrate.plot$age<-max(absrate.plot$age)-absrate.plot$age

    car::outlierTest(lm(absrate.plot[,1]~absrate.plot[,2]))->outT
    if(any(outT$p<=0.05))
      max(absrate.plot[-as.numeric(names(outT$p)),1])->maxy else
        max(absrate.plot[,1])->maxy

    plotSTphenNode<-function(variable,node=NULL,plot.args=NULL,lineTree.args=NULL,lineNode.args=NULL,node.palette=NULL){
      if(is.character(variable)){
        if(is.na(match(variable,colnames(phen.plot)))) stop("variable do not match column names")
        match(variable,colnames(phen.plot))->variable
      }else if(variable>ny) stop("variable is out of y bounds")

      if(!is.null(node)){
        if(all(node>length(groups))) as.character(node)->node
        if(is.character(node)){
          if(any(is.na(match(node,gsub("g","",groups))))) stop("node not in tested clades")
          groups[match(node,gsub("g","",groups))]->groups
        }else{
          if(any(node>length(groups))) stop("node are out of tested clades bounds")
          groups[node]->groups
        }
      }

      if(!is.null(node.palette)){
        if(length(node.palette)<length(groups)) rep(node.palette,length.out=length(groups))->node.palette
      }else node.palette<-suppressWarnings(RColorBrewer::brewer.pal(length(groups), "Set2"))

      if("pch.node"%in%names(plot.args)){
        plot.args$pch.node->pch.node
        if(length(pch.node)<length(groups)) rep(pch.node,length.out=length(groups))->pch.node
        plot.args[which(names(plot.args)!="pch.node")]->plot.args
      }else pch.node<-NULL

      plotind<-1:length(groups)
      cbind(rep(1:3,each=2)[-6],rep(1:3,each=2)[-1])->dfmat
      cbind(dfmat,dfmat[,1]*dfmat[,2])->dfmat
      cbind(dfmat,dfmat[,3]-length(groups))->dfmat
      dfmat[which(dfmat[,4]>=0),,drop=FALSE][which.min(dfmat[which(dfmat[,4]>=0),3]),1:2]->nn
      if((nn[1]*nn[2])>length(groups)) c(plotind,rep(0,(nn[1]*nn[2])-length(groups)))->plotind
      matrix(plotind,nrow=nn[1],ncol=nn[2],byrow=TRUE)->mat

      layout(mat)

      variable->i
      lapply(1:length(groups),function(j){
        phen.plot[which(phen.plot$group==groups[j]),c(i,ny+1)]->node.phen
        phen.plot[which(phen.plot$group!=groups[j]),c(i,ny+1)]->tree.phen
        predict(lm(tree.phen[,1]~tree.phen[,2]))->pred.tree
        predict(lm(node.phen[,1]~node.phen[,2]))->pred.node

        if(is.null(plot.args)) plot.args<-list()
        if(is.null(lineTree.args)) lineTree.args<-list()
        if(is.null(lineNode.args)) lineNode.args<-list()


        if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"age"
        if("ylab"%in%names(plot.args)) ylabs<-plot.args$ylab else ylabs<-NULL
        if(is.null(ylabs)) plot.args$ylab<-paste("Clade",gsub("g","",groups[j])) else{
          if(length(ylabs)>1) plot.args$ylab<-ylabs[j]
        }
        if("main"%in%names(plot.args)) maint<-plot.args$main else maint<-NULL
        if(is.null(maint)) plot.args$main<-paste("Phenotypic Trend for Variable",colnames(phen.plot)[i]) else{
          if(length(maint)>1) plot.args$main<-maint[j]
        }
        if(!"xlim"%in%names(plot.args)) plot.args$xlim<-c(max(phen.plot[,ny+1]),min(phen.plot[,ny+1]))
        if(!"ylim"%in%names(plot.args)) plot.args$ylim<-range(phen.plot[,i])
        if(!"col"%in%names(plot.args)) "black"->plot.args$col
        if(!"bg"%in%names(plot.args)) "white"->plot.args$bg
        if(!"pch"%in%names(plot.args)) plot.args$pch<-21
        rep(plot.args$pch,length(phen.plot[, i]))->plot.args$pch
        if(!is.null(pch.node)) plot.args$pch[which(rownames(phen.plot)%in%rownames(node.phen))]<-pch.node[j]

        if(unique(plot.args$pch[which(rownames(phen.plot)%in%rownames(node.phen))])%in%21:25){
          rep(plot.args$bg,length(phen.plot[, i]))->plot.args$bg
          plot.args$bg[which(rownames(phen.plot)%in%rownames(node.phen))]<-node.palette[j]
        }else{
          rep(plot.args$col,length(phen.plot[, i]))->plot.args$col
          plot.args$col[which(rownames(phen.plot)%in%rownames(node.phen))]<-node.palette[j]
        }

        c(x=list(phen.plot[, ny+1]),y=list(phen.plot[, i]),plot.args)->points.args
        c(x=list(NA),y=list(NA),plot.args)->plot.args

        c(x=list(tree.phen$age[order(tree.phen$age)]),
          y=list(pred.tree[order(tree.phen$age)]),type="l",lineTree.args)->lineTree.args
        if(!"col"%in%names(lineTree.args)) lineTree.args$col<-"gray70"
        if(!"lty"%in%names(lineTree.args)) lineTree.args$lty<-3
        if(!"lwd"%in%names(lineTree.args)) lineTree.args$lwd<-3

        c(x=list(node.phen$age[order(node.phen$age)]),
          y=list(pred.node[order(node.phen$age)]),type="l",lineNode.args)->lineNode.args
        if(!"lwd"%in%names(lineNode.args)) lineNode.args$lwd<-4
        lineNode.args$col<-node.palette[j]

        do.call(plot,plot.args)
        do.call(points,points.args)
        do.call(points,lineTree.args)
        do.call(points,lineNode.args)
      })
    }

    plotSTratesNode<-function(variable,node=NULL,plot.args=NULL,lineTree.args=NULL,lineNode.args=NULL,node.palette=NULL){
      if(is.character(variable)){
        if(any(is.na(match(variable,colnames(absrate.plot))))) stop("variable do not match column names")
        match(variable,colnames(absrate.plot))->variable
      }else if(any(variable>nyrate)) stop("variable is out of y bounds")

      if(!is.null(node)){
        if(all(node>length(groups))) as.character(node)->node
        if(is.character(node)){
          if(any(is.na(match(node,gsub("g","",groups))))) stop("node not in tested clades")
          groups[match(node,gsub("g","",groups))]->groups
        }else{
          if(any(node>length(groups))) stop("node are out of tested clades bounds")
          groups[node]->groups
        }
      }

      if(!is.null(node.palette)){
        if(length(node.palette)<length(groups)) rep(node.palette,length.out=length(groups))->node.palette
      }else node.palette<-suppressWarnings(RColorBrewer::brewer.pal(length(groups), "Set2"))

      if("pch.node"%in%names(plot.args)){
        plot.args$pch.node->pch.node
        if(length(pch.node)<length(groups)) rep(pch.node,length.out=length(groups))->pch.node
        plot.args[which(names(plot.args)!="pch.node")]->plot.args
      }else pch.node<-NULL

      plotind<-1:length(groups)
      cbind(rep(1:3,each=2)[-6],rep(1:3,each=2)[-1])->dfmat
      cbind(dfmat,dfmat[,1]*dfmat[,2])->dfmat
      cbind(dfmat,dfmat[,3]-length(groups))->dfmat
      dfmat[which(dfmat[,4]>=0),,drop=FALSE][which.min(dfmat[which(dfmat[,4]>=0),3]),1:2]->nn
      if((nn[1]*nn[2])>length(groups)) c(plotind,rep(0,(nn[1]*nn[2])-length(groups)))->plotind
      matrix(plotind,nrow=nn[1],ncol=nn[2],byrow=TRUE)->mat

      layout(mat)

      variable->i
      lapply(1:length(groups),function(j){
        absrate.plot[which(absrate.plot$group==groups[j]),c(i,nyrate+1)]->node.rates
        absrate.plot[which(absrate.plot$group!=groups[j]),c(i,nyrate+1)]->tree.rates
        predict(lm(tree.rates[,1]~tree.rates[,2]))->pred.tree
        predict(lm(node.rates[,1]~node.rates[,2]))->pred.node


        if(is.null(plot.args)) plot.args<-list()
        if(is.null(lineTree.args)) lineTree.args<-list()
        if(is.null(lineNode.args)) lineNode.args<-list()

        if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"age"
        if("ylab"%in%names(plot.args)) ylabs<-plot.args$ylab else ylabs<-NULL
        if(is.null(ylabs)) plot.args$ylab<-paste("Clade",gsub("g","",groups[j])) else{
          if(length(ylabs)>1) plot.args$ylab<-ylabs[j]
        }
        if("main"%in%names(plot.args)) maint<-plot.args$main else maint<-NULL
        if(is.null(maint)) plot.args$main<-paste("Absolute Rates for Variable",colnames(absrate.plot)[i]) else{
          if(length(maint)>1) plot.args$main<-maint[j]
        }

        if(!"xlim"%in%names(plot.args)) plot.args$xlim<-c(max(absrate.plot[,nyrate+1]),min(absrate.plot[,nyrate+1]))
        if(!"ylim"%in%names(plot.args)) plot.args$ylim<-c(min(absrate.plot[,1]),maxy)
        if(!"col"%in%names(plot.args)) "black"->plot.args$col
        if(!"bg"%in%names(plot.args)) "white"->plot.args$bg
        if(!"pch"%in%names(plot.args)) plot.args$pch<-21
        rep(plot.args$pch,length(absrate.plot[, i]))->plot.args$pch
        if(!is.null(pch.node)) plot.args$pch[which(rownames(absrate.plot)%in%rownames(node.rates))]<-pch.node[j]

        if(unique(plot.args$pch[which(rownames(absrate.plot)%in%rownames(node.rates))])%in%21:25){
          rep(plot.args$bg,length(absrate.plot[, i]))->plot.args$bg
          plot.args$bg[which(rownames(absrate.plot)%in%rownames(node.rates))]<-node.palette[j]
        }else{
          rep(plot.args$col,length(absrate.plot[, i]))->plot.args$col
          plot.args$col[which(rownames(absrate.plot)%in%rownames(node.rates))]<-node.palette[j]
        }

        c(x=list(absrate.plot[, nyrate+1]),y=list(absrate.plot[, i]),plot.args)->points.args
        c(x=list(NA),y=list(NA),plot.args)->plot.args

        c(x=list(tree.rates$age[order(tree.rates$age)]),
          y=list(pred.tree[order(tree.rates$age)]),type="l",lineTree.args)->lineTree.args
        if(!"col"%in%names(lineTree.args)) lineTree.args$col<-"gray70"
        if(!"lty"%in%names(lineTree.args)) lineTree.args$lty<-3
        if(!"lwd"%in%names(lineTree.args)) lineTree.args$lwd<-3

        c(x=list(node.rates$age[order(node.rates$age)]),
          y=list(pred.node[order(node.rates$age)]),type="l",lineNode.args)->lineNode.args
        if(!"lwd"%in%names(lineNode.args)) lineNode.args$lwd<-4
        lineNode.args$col<-node.palette[j]

        do.call(plot,plot.args)
        do.call(points,points.args)
        do.call(points,lineTree.args)
        do.call(points,lineNode.args)
      })
    }


  }else plotSTphenNode<-NULL


  plotSTphen<-function(variable,plot.args=NULL,polygon.args=NULL,line.args=NULL){
    if(is.character(variable)){
      if(any(is.na(match(variable,colnames(phen.plot))))) stop("variable do not match column names")
      match(variable,colnames(phen.plot))->variable
    }else if(any(variable>ny)) stop("variable is out of y bounds")

    variable->i
    lm(phen.plot[, i] ~ I(max(phen.plot$age)-phen.plot$age))->plot.reg
    quantile(CIphen[[i]],c(0.025,0.975))->ci.slope
    sapply(ci.slope,function(k) coef(plot.reg)[1]+k*c(max(phen.plot$age),0))->ci.points


    if(is.null(plot.args)) plot.args<-list()
    if(is.null(polygon.args)) polygon.args<-list()
    if(is.null(line.args)) line.args<-list()

    if(!"main"%in%names(plot.args)) plot.args$main<-"Phenotypic Trend Test"
    if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"age"
    if(!"ylab"%in%names(plot.args)) plot.args$ylab<-paste("Rescaled",colnames(phen.plot)[i])
    if(!"xlim"%in%names(plot.args)) plot.args$xlim<-c(max(phen.plot[,ny+1]),min(phen.plot[,ny+1]))
    if(!"ylim"%in%names(plot.args)) plot.args$ylim<-range(phen.plot[,i])
    if(!"pch"%in%names(plot.args)) plot.args$pch<-21
    if(length(plot.args$pch)<length(phen.plot[, i])) rep(plot.args$pch,length.out=length(phen.plot[, i]))->plot.args$pch
    if(!"col"%in%names(plot.args)) plot.args$col<-"black"
    if(!"bg"%in%names(plot.args)) plot.args$bg<-"white"

    pch.type<-plot.args$pch
    pch.type[which(plot.args$pch>=20)]<-"bg"
    pch.type[which(plot.args$pch<20)]<-"col"

    if(length(plot.args$col)==1&any(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="col")){
      rep(plot.args$col,length(phen.plot[, i]))->plot.args$col
      plot.args$col[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="col")]<-"red"
    }

    if(length(plot.args$bg)==1&any(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="bg")){
      rep(plot.args$bg,length(phen.plot[, i]))->plot.args$bg
      plot.args$bg[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="bg")]<-"red"
    }

    c(x=list(phen.plot[, ny+1]),y=list(phen.plot[, i]),plot.args)->points.args
    c(x=list(NA),y=list(NA),plot.args)->plot.args

    c(x=list(c(c(0,max(phen.plot$age),c(max(phen.plot$age),0)))),
      y=list(c(ci.points[,1], rev(ci.points[, 2]))),polygon.args)->polygon.args
    if(!"col"%in%names(polygon.args)) polygon.args$col<-rgb(0.5, 0.5,0.5, 0.4)
    if(!"border"%in%names(polygon.args)) polygon.args$border<-NA

    c(reg=list(lm(phen.plot[, i] ~ phen.plot$age)),line.args)->line.args
    if(!"col"%in%names(line.args)) line.args$col<-"blue"
    if(!"lwd"%in%names(line.args)) line.args$lwd<-3

    do.call(plot,plot.args)
    do.call(polygon,polygon.args)
    do.call(points,points.args)
    do.call(abline,line.args)
  }

  plotSTrates<-function(variable,plot.args=NULL,polygon.args=NULL,line.args=NULL){
    if(is.character(variable)){
      if(any(is.na(match(variable,colnames(absrate.plot))))) stop("variable do not match column names")
      match(variable,colnames(absrate.plot))->variable
    }else if(any(variable>nyrate)) stop("variable is out of y bounds")

    resrate.plot$age<-max(na.omit(resrate.plot$age))-resrate.plot$age

    variable->i

    lm(resrate.plot[,i] ~ I(max(resrate.plot$age)-resrate.plot$age))->plot.reg
    quantile(CIresrate[[i]],c(0.025,0.975))->ci.slope.scal
    sapply(ci.slope.scal,function(k) coef(plot.reg)[1]+k*c(max(resrate.plot$age),0))->ci.points.scal

    if(is.null(plot.args)) plot.args<-list()
    if(is.null(polygon.args)) polygon.args<-list()
    if(is.null(line.args)) line.args<-list()

    if(!"main"%in%names(plot.args)) plot.args$main<-"Rescaled Rate Trend Test"
    if(!"xlab"%in%names(plot.args)) plot.args$xlab<-"age"
    if(!"ylab"%in%names(plot.args)) plot.args$ylab<-paste("log rescaled rate for variable",colnames(absrate.plot)[i])
    if(!"ylim"%in%names(plot.args)) plot.args$ylim<-range(na.omit(resrate.plot[,i]))
    if(!"xlim"%in%names(plot.args)) plot.args$xlim<-c(max(na.omit(resrate.plot$age)),min(na.omit(resrate.plot$age)))
    if(!"pch"%in%names(plot.args)) plot.args$pch<-21
    if(length(plot.args$pch)<length(phen.plot[, i])) rep(plot.args$pch,length.out=length(phen.plot[, i]))->plot.args$pch
    if(!"col"%in%names(plot.args)) plot.args$col<-"black"
    if(!"bg"%in%names(plot.args)) plot.args$bg<-"white"

    pch.type<-plot.args$pch
    pch.type[which(pch.type>=20)]<-"bg"
    pch.type[which(pch.type<20)]<-"col"

    if(length(plot.args$col)==1&any(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="col")){
      rep(plot.args$col,length(phen.plot[, i]))->plot.args$col
      plot.args$col[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="col")]<-"red"
    }

    if(length(plot.args$bg)==1&any(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="bg")){
      rep(plot.args$bg,length(phen.plot[, i]))->plot.args$bg
      plot.args$bg[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))&pch.type=="bg")]<-"red"
    }

    c(x=list(resrate.plot$age),y=list(resrate.plot[,i]),plot.args)->points.args
    c(x=list(NA),y=list(NA),plot.args)->plot.args

    c(x=list(c(c(0,max(resrate.plot$age),c(max(resrate.plot$age),0)))),
      y=list(c(ci.points.scal[,1], rev(ci.points.scal[, 2]))),polygon.args)->polygon.args
    if(!"col"%in%names(polygon.args)) polygon.args$col<-rgb(0.5, 0.5,0.5, 0.4)
    if(!"border"%in%names(polygon.args)) polygon.args$border<-NA

    c(reg=list(lm(resrate.plot[,i] ~ resrate.plot$age)),line.args)->line.args
    if(!"col"%in%names(line.args)) line.args$col<-"blue"
    if(!"lwd"%in%names(line.args)) line.args$lwd<-3

    do.call(plot,plot.args)
    do.call(polygon,polygon.args)
    do.call(points,points.args)
    do.call(abline,line.args)

  }

  list(plotSTphen=plotSTphen,plotSTrates=plotSTrates)->res
  if(!is.null(plotSTphenNode)) c(res,list(plotSTphenNode=plotSTphenNode,plotSTratesNode=plotSTratesNode))->res
  return(res)
}
