#'@title Plot RRphylo rates at a specified node
#'@description This function generates customized functions to produce
#'  histograms and lollipop charts of the \code{\link{RRphylo}} rates computed
#'  for a given clade as compared to the rates computed for the rest of the
#'  tree.
#'@usage plotRates(RR,node,export.tiff=NULL,filename=NULL)
#'@param RR an object produced by \code{\link{RRphylo}}.
#'@param node the node subtending the clade of interest.
#'@param export.tiff is deprecated.
#'@param filename is deprecated.
#'@export
#'@seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#'@seealso \href{../doc/Plotting-tools.html}{\code{plotRates} vignette}
#'@author Silvia Castiglione, Pasquale Raia
#'@importFrom graphics abline legend hist par
#'@return The function returns a list of functions:
#'@return \strong{$plotHist} returns the histograms of rates (in ln absolute
#'  values) computed for the focal clade against rates computed the rest of the
#'  tree. The usage is:
#'  \code{...$plotHist(hist.args=list(col1,col2),legend.args=list())}, where
#'  \code{legend.args} is a list of additional arguments passed to the function
#'  \code{legend} (if \code{= NULL} the legend is not plotted) and
#'  \code{hist.args} is a list of further arguments passed to the function
#'  \code{plot.histogram}. \code{hist.args} default arguments include histogram
#'  colors for background rates (\code{col1}) and rates of the clade under
#'  inspection (\code{col2}).
#'@return \strong{$plotLollipop} returns the lollipop chart of the rates of
#'  individual branches of the focal clade collated in increasing rate value,
#'  and contrasted to the average rate computed over the rest of the tree
#'  branches (the vertical line). The usage is:
#'  \code{...$plotLollipop(lollipop.args=NULL,line.args=NULL)}, where
#'  \code{lollipop.args} is a list of further arguments passed to the function
#'  \code{\link{lollipoPlot}}  and \code{line.args} is a list of additional
#'  arguments passed to the function \code{line}. This function additionally
#'  returns the vector of rates for the focal clade, collated in increasing
#'  order.
#'@seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#'@seealso \href{../doc/Plotting-tools.html}{\code{plotRates} vignette}
#' @examples
#'\donttest{
#' data("DataApes")
#' DataApes$PCstage->PCstage
#' DataApes$Tstage->Tstage
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=Tstage,y=PCstage,clus=cc)->RR
#'
#' plotRates(RR,node=72)->pR
#' pR$plotHist(hist.args=list(col1="cyan1",col2="blue"),legend.args=list(x="topright"))
#' pR$plotLollipop(lollipop.args=list(col="chartreuse",bg="chartreuse"),
#'                 line.args=list(col="deeppink",lwd=2))
#' }

plotRates<-function(RR,node,export.tiff=NULL,filename=NULL){
  # require(phytools)

  RR$rates->rates
  RR$tree->tree
  getDescendants(tree,node)->des

  c(rates[match(des[des>Ntip(tree)],rownames(rates)),],rates[match(tree$tip.label[des[des<=Ntip(tree)]],rownames(rates)),])->shift.rates
  shift.rates[order(shift.rates,decreasing=TRUE)]->shift.rates

  plotHist<-function(hist.args=list(col1="green",col2="blue"),legend.args=list()){
    hist(log(abs(shift.rates)),plot=FALSE)->H1
    hist(log(abs(rates[-match(names(shift.rates),rownames(rates)),])),plot=FALSE)->H2
    log(abs(rates))[(log(abs(rates))!="-Inf")]->Xa

    if(is.null(hist.args$col1)) hist.args$col1<-"green"
    if(is.null(hist.args$col2)) hist.args$col2<-"blue"
    hist.args$col1->col1
    hist.args$col2->col2
    hist.args[-which(names(hist.args)%in%c("col1","col2"))]->hist.args
    if(!"xlim"%in%names(hist.args)) hist.args$xlim<-c(range(Xa)[1]-diff(range(Xa))*0.1,range(Xa)[2]+diff(range(Xa))*0.1)
    if(!"ylim"%in%names(hist.args)) hist.args$ylim<-c(0,max(H2$count)*1.2)
    if(!"xlab"%in%names(hist.args)) hist.args$xlab<-"log absolute rates"
    if(!"main"%in%names(hist.args)) hist.args$main<-""

    do.call(plot,c(x=list(H2),hist.args,col=col1))
    do.call(plot,c(x=list(H1),hist.args,col=col2,add=TRUE))

    if(!is.null(legend.args)){
      if(!"x"%in%names(legend.args)) legend.args$x<-"topleft"
      if(!"legend"%in%names(legend.args)) legend.args$legend<-c("back rates", "shift node rates")
      if(!"fill"%in%names(legend.args)&!any(c("lty","lwd")%in%names(legend.args))){
        if(!"pch"%in%names(legend.args)) legend.args$pch<-rep(22,2)
        if(any(legend.args$pch>20)&!"pt.bg"%in%names(legend.args)){
          legend.args$pt.bg<-rep(NA,length(legend.args$pch))
          legend.args$pt.bg[which(legend.args$pch>20)]<-c(col1,col2)[which(legend.args$pch>20)]
        }
        if(any(legend.args$pch<21)&!"col"%in%names(legend.args)){
          legend.args$col<-rep(NA,length(legend.args$pch))
          legend.args$col[which(legend.args$pch<21)]<-c(col1,col2)[which(legend.args$pch<21)]
          legend.args$col[which(is.na(legend.args$col))]<-"black"
        }
        if(!"pt.cex"%in%names(legend.args)) legend.args$pt.cex<-1.5
      }

      if(!"bty"%in%names(legend.args)) legend.args$bty<-"n"
      do.call(legend,legend.args)
    }
  }

  plotLollipop<-function(lollipop.args=NULL,line.args=NULL){
    if(is.null(lollipop.args)) lollipop.args<-list()
    if(is.null(line.args)) line.args<-list()

    lollipop.args$values<-shift.rates
    lollipop.args$type<-"h"
    if(!"xlab"%in%names(lollipop.args)) lollipop.args$xlab<-"rates"
    if(!"lwd"%in%names(lollipop.args)) lollipop.args$lwd<-2
    if(!"pch"%in%names(lollipop.args)){
      if("pt.col"%in%names(lollipop.args)&!"bg"%in%names(lollipop.args)) lollipop.args$pch<-1 else lollipop.args$pch<-21
    }
    if(!"pt.col"%in%names(lollipop.args)) lollipop.args$pt.col<-"black"
    if(!"col"%in%names(lollipop.args)) lollipop.args$col<-"blue"
    if(!"bg"%in%names(lollipop.args)) lollipop.args$bg<-"blue"
    if(!"cex"%in%names(lollipop.args)) lollipop.args$cex<-1.5

    line.args$v<-mean(rates)
    if(!"col"%in%names(line.args)) line.args$col<-"green"
    if(!"lwd"%in%names(line.args)) line.args$lwd<-3

    if(identical(par("mar"),c(5.1,4.1,4.1,2.1))){
      pmar<-c(5.1,4.1,4.1,2.1)
      par(mar=c(5.1,15,4.1,2.1))
    }else par("mar")->pmar
    do.call(lollipoPlot,lollipop.args)
    do.call(abline,line.args)

    par(mar=pmar)
    return(shift.rates)
  }

  return(list(plotHist=plotHist,plotLollipop=plotLollipop))

}

