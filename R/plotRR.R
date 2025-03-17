#'@title Plot the RRphylo output onto the phylogenetic tree
#'@description This function generates customized functions to plot the
#'  phylogenetic tree (as returned by \code{\link{RRphylo}}) with branches
#'  colored according to phenotypic values or phenotypic evolutionary rates.
#'@usage plotRR(RR,y,multivariate=NULL)
#'@param RR an object produced by \code{\link{RRphylo}}.
#'@param y the vector/matrix of phenotypic values used to perform
#'  \code{\link{RRphylo}}.
#'@param multivariate if \code{\link{RRphylo}} was performed on multivariate data, this
#'  argument indicates whether individual rates for each variables (\code{=
#'  "multiple.rates"}) or the norm2 vector of multivariate rates (\code{=
#'  "rates"}) should be plotted.
#'@export
#'@seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#'@seealso \href{../doc/Plotting-tools.html#plotRR}{\code{plotRR} vignette}
#'@author Silvia Castiglione, Pasquale Raia
#'@importFrom graphics rect strheight strwidth
#'@importFrom grDevices colorRampPalette
#'@return The function returns a list of functions:
#'@return \strong{$plotRRphen} charts phenotypic values along the tree branches.
#'  Phenotypes at tips are taken as they are from the \code{y} object.
#'  Phenotypic values for internal branches are derived from the \code{RR$aces}
#'  object. The usage is:
#'  \code{...$plotRRphen(variable=NULL,tree.args=NULL,color.pal=NULL,colorbar.args=list())},
#'  where \code{variable} is the index or name of the variable to be plotted in
#'  case of multivariate data, \code{tree.args} is a list of further arguments
#'  passed to the function \code{plot.phylo} plus a logical indicating whether
#'  the tree should be ladderized  before plotting (see examples below),
#'  \code{color.pal} is a function to generate the color palette, and
#'  \code{colorbar.args} is a list of further arguments passed to the function
#'  \code{\link{colorbar}} (if \code{= NULL} the bar is not plotted).
#'@return \strong{$plotRRrates} charts evolutionary rate values along the tree
#'  branches. The usage is identical to \code{$plotRRphen}. In case of
#'  multivariate data and \code{multivariate = "rates"}, the argument
#'  \code{variable} can be left unspecified.
#'@seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#'@seealso \href{../doc/Plotting-tools.html}{\code{plotRR} vignette}
#' @examples
#'\dontrun{
#' data("DataApes")
#' DataApes$PCstage->PCstage
#' DataApes$Tstage->Tstage
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=Tstage,y=PCstage,clus=cc)->RRstage
#'
#' plotRR(RRstage,y=PCstage,multivariate="multiple.rates")->pRR1
#' pRR1$plotRRphen(variable=1,tree.args=list(edge.width=2),color.pal=rainbow,
#'                colorbar.args = list(x="bottomleft",labs.adj=0.7,xpd=TRUE))
#' pRR1$plotRRrates(variable=2,tree.args=list(edge.width=2,direction="leftwards",ladderize=TRUE),
#'                 color.pal=rainbow,colorbar.args = list(x="topright",labs.adj=0.7,xpd=TRUE))
#'
#'
#' plotRR(RRstage,y=PCstage,multivariate="rates")->pRR2
#' pRR2$plotRRrates(tree.args=list(edge.width=2),
#'                 color.pal=hcl.colors,
#'                 colorbar.args = list(x="topleft",labs.adj=0.7,xpd=TRUE,title.pos="bottom"))
#' }



plotRR<-function(RR,y,multivariate=NULL){
  RR$tree->tree
  y <- treedataMatch(tree, y)[[1]]

  rbind(RR$aces,y)->phen
  rownames(phen)[which(rownames(phen)%in%tree$tip.label)]<-
    match(rownames(phen)[which(rownames(phen)%in%tree$tip.label)],tree$tip.label)

  if(!is.null(multivariate)){
    if(is.na(match(multivariate,names(RR)))) stop("The argument multivariate must be one of 'rates' or 'multiple.rates'")
    RR[[match(multivariate,names(RR))]]->rates
  }else RR$rates->rates
  rownames(rates)[which(rownames(rates)%in%tree$tip.label)]<-
    match(rownames(rates)[which(rownames(rates)%in%tree$tip.label)],tree$tip.label)

  plotRRphen<-function(variable=NULL,tree.args=NULL,color.pal=NULL,colorbar.args=list()){
    if(is.null(variable)){
      if(ncol(phen)>1) stop("Please indicate the variable to plot") else variable<-1
    }else{
      if(is.character(variable)){
        if(is.na(match(variable,colnames(phen)))) stop("variable do not match column names")
        match(variable,colnames(phen))->variable
      }else if(variable>ncol(phen)) stop("variable is out of y bounds")

    }

    phen[-1,variable,drop=FALSE]->phenvar

    if(Ntip(tree)>100){
      if(all(!c("show.tip.label","cex")%in%names(tree.args))) tree.args$show.tip.label<-FALSE
    }

    if(!is.null(tree.args$ladderize)&&isTRUE(tree.args$ladderize)){
      ladderize(tree)->tree
      tree.args<-tree.args[which(names(tree.args)!="ladderize")]
    }

    # if(isTRUE(tree.args$no.margin)){
    #     mars <- par("mar")
    #     on.exit(par(mar = mars))
    # }

    if(is.null(color.pal)) colorRampPalette(c("plum1","purple3"))->color.pal
    color.pal(nrow(phenvar))->colpal
    names(colpal)<-rownames(phenvar)[order(phenvar[,1])]
    colpal[match(tree$edge[,2],names(colpal))]->colpal
    tree.args$edge.color<-colpal

    do.call(plot.phylo,c(x=list(tree),tree.args))

    if(!is.null(colorbar.args)){
      colorbar.args$labs<-sapply(range(phenvar),function(x)
        ifelse(x<0.01&x>-0.01,format(x,scientific=TRUE,digits=2),format(x,digits=3)))
      colorbar.args$colors<-colpal[match(rownames(phenvar)[order(phenvar[,1])],names(colpal))]
      if(is.null(colorbar.args$x)) colorbar.args$x<-"topleft"
      if(is.null(colorbar.args$direction)) colorbar.args$direction<-"vertical"
      if(is.null(colorbar.args[["title",exact=TRUE]])) c(colorbar.args,list(title="phenotypes"))->colorbar.args
      do.call(colorbar,colorbar.args)
    }
  }

  plotRRrates<-function(variable=NULL,tree.args=NULL,color.pal=NULL,colorbar.args=list()){
    if(is.null(variable)){
      if(ncol(rates)>1) stop("Please indicate the variable to plot") else variable<-1
    }else{
      if(is.character(variable)){
        if(is.na(match(variable,colnames(rates)))) stop("variable do not match column names")
        match(variable,colnames(rates))->variable
      }else if(variable>ncol(rates)) stop("variable is out of y bounds")
    }

    rates[-1,variable,drop=FALSE]->ratevar

    if(Ntip(tree)>100){
      if(all(!c("show.tip.label","cex")%in%names(tree.args))) tree.args$show.tip.label<-FALSE
    }

    if(!is.null(tree.args$ladderize)&&isTRUE(tree.args$ladderize)){
      ladderize(tree)->tree
      tree.args<-tree.args[which(names(tree.args)!="ladderize")]
    }

    # if(isTRUE(tree.args$no.margin)){
    #   mars <- par("mar")
    #   on.exit(par(mar = mars))
    # }

    if(is.null(color.pal)) colorRampPalette(c("lightblue1","darkblue"))->color.pal
    color.pal(nrow(ratevar))->colpal
    names(colpal)<-rownames(ratevar)[order(ratevar)]
    colpal[match(tree$edge[,2],names(colpal))]->colpal
    tree.args$edge.color<-colpal


    do.call(plot.phylo,c(x=list(tree),tree.args))

    if(!is.null(colorbar.args)){
      colorbar.args$labs<-sapply(range(ratevar),function(x)
        ifelse(x<0.01&x>-0.01,format(x,scientific=TRUE,digits=2),format(x,digits=3)))
      colorbar.args$colors<-colpal[match(rownames(ratevar)[order(ratevar)],names(colpal))]
      if(is.null(colorbar.args$x)) colorbar.args$x<-"topleft"
      if(is.null(colorbar.args$direction)) colorbar.args$direction<-"vertical"
      if(is.null(colorbar.args[["title",exact=TRUE]])) c(colorbar.args,list(title="rates"))->colorbar.args
      do.call(colorbar,colorbar.args)
    }
  }

  return(list(plotRRphen=plotRRphen,plotRRrates=plotRRrates) )
}
