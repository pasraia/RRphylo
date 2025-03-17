#'@title Graphical representation of search.shift results
#'@description \code{plotShift} generates customized functions to produce plots
#'  out of \code{\link{search.shift}} results. \code{addShift} adds circles to
#'  highlight shifting clades onto the currently plotted tree.
#'@aliases plotShift
#'@aliases addShift
#'@usage plotShift(RR,SS,mode,state=NULL)
#'@usage addShift(SS,symbols.args=NULL)
#'@param RR the object produced by \code{\link{RRphylo}} used to perform
#'  \code{\link{search.shift}}.
#'@param SS an object produced by \code{\link{search.shift}}.
#'@param mode if \code{\link{search.shift}} was performed under \code{status.type =
#'  "clade"} by setting the \code{node} argument, \code{mode} is a numeric
#'  indicating whether the output number 1, that is
#'  \code{...$single.clades$singles} or the output number 2
#'  \code{...$single.clades$no.others} should be plotted.
#'@param state if \code{\link{search.shift}} was performed under \code{status.type =
#'  "sparse"}, this is the same \code{state} vector provided to the function.
#'@param symbols.args as described for \code{$plotClades} below.
#'@details Using \code{...$plotClades()} or plotting the tree and applying
#'  \code{addShift()} returns the same plot. The latter function might be useful
#'  in combination with \code{\link{plotRR}} to add the shifts information to
#'  the branch-wise plot of evolutionary rate values.
#'@return The function returns a function to generate a customizable plot of
#'  \code{\link{search.shift}} results.
#'@return If \code{\link{search.shift}} was performed under \code{status.type =
#'  "clade"}, \code{plotShift} returns a \strong{\code{$plotClades}} function
#'  which highlights the shifting clades onto the phylogenetic tree. The usage
#'  is:
#'  \code{...$plotClades(tree.args=NULL,symbols.args=NULL)}, where
#'  \code{tree.args} is a list of further arguments passed to the function
#'  \code{plot.phylo}, and \code{symbols.args} is a list of further arguments
#'  passed to the function \code{symbols} (n.b. the shape of the symbol is not
#'  customizable). The function automatically plots red circles for negative
#'  shifts and blue circles for positive shifts, in each cases with the radium
#'  proportional to the absolute value of rate difference. The user can choose
#'  different color options for positive/negative shifts by setting
#'  \code{symbols.args=list(fg=c(pos="color for positive shift",neg="color for
#'  negative shift"))}, or provide a vector of as many colors as the number of
#'  shifting clades. The same applies to the argument "bg".
#'@return If \code{\link{search.shift}} was performed under \code{status.type =
#'  "sparse"}, \code{plotShift} returns a \strong{\code{$plotStates}} function
#'  which plots the comparison between the real difference and the distributions
#'  of random differences (see
#'  \href{../doc/search.shift.html}{\code{search.shift} vignette#sparse}). The
#'  usage is:
#'  \code{...$plotStates(plot.args=NULL,points.args=NULL,legend.args=list())},
#'  where \code{plot.args} is a list of further arguments passed to the function
#'  \code{plot} (used in the form: \code{plot(y~x)}), \code{points.args} is a
#'  list of further arguments passed to the function \code{points}, and
#'  \code{legend.args} is a list of additional arguments passed to the function
#'  \code{legend} (if \code{= NULL} the legend is not plotted). If as many
#'  colors/pch values as the number of different states are provided in
#'  \code{points.args}, each of them is assigned to each states taken in the
#'  same alphabetical order.
#'@author Silvia Castiglione, Giorgia Girardi
#'@importFrom graphics symbols
#'@importFrom ape plot.phylo
#'@export
#'@seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#'@seealso \href{../doc/Plotting-tools.html#plotShift}{\code{plotShift} vignette}
#'@references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'  Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method
#'  for testing evolutionary rate variation and shifts in phenotypic evolution.
#'  \emph{Methods in Ecology and Evolution}, 9:
#'  974-983.doi:10.1111/2041-210X.12954
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#'
#' search.shift(RR=dinoRates,status.type="clade")->SSauto
#' plotShift(RR=dinoRates,SS=SSauto)->plotSSauto
#' plotSSauto$plotClades()
#'
#' plot(dinoRates$tree)
#' addShift(SS=SSauto)
#'
#' search.shift(RR=dinoRates,status.type="clade",node=c(696,746))->SSnode
#' plotShift(RR=dinoRates,SS=SSnode,mode=2)->plotSSnode
#' plotSSnode$plotClades(tree.args=list(no.margin=TRUE),
#'                   symbols.args=list(fg=NA,bg=c(pos="cyan",neg="magenta")))
#'
#'
#' search.shift(RR=dinoRates,status.type= "sparse",state=statedino)->SSstate
#' plotShift(RR=dinoRates,SS=SSstate,state=statedino)->plotSSstate
#' plotSSstate$plotStates(points.args=list(bg=c("gold","forestgreen","royalblue","white"),
#'                                    col=c("black","black","black","orangered"),
#'                                    pch=c(21,22,24,11)),legend.args=list())
#'     }

plotShift<-function(RR,SS,mode=NULL,state=NULL){

  misspacks<-sapply(c("scales","RColorBrewer"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  RR$tree->tree
  if(is.null(SS$single.clades)&is.null(SS$state.results)){
    stop("No significant result available")
  }
  if(is.null(SS$single.clades$singles)) mode<-NULL
  if(!is.null(SS$single.clades$singles)&is.null(mode))
    stop("Please set the argument `mode` to indicate which output of search.shift you wish to plot")

  if(!is.null(SS$single.clades)){
    if(is.null(mode)) SS$single.clades->single else SS$single.clades[[mode]]->single
    plotClades<-function(tree.args=NULL,symbols.args=NULL){

      if(Ntip(tree)>100){
        if(all(!c("show.tip.label","cex")%in%names(tree.args))) tree.args$show.tip.label<-FALSE
      }

      # if(isTRUE(tree.args$no.margin)){
      #   mars <- par("mar")
      #   on.exit(par(mar = mars))
      # }

      if(any(c("pos","neg")%in%names(symbols.args$fg))|is.null(symbols.args$fg)){
        symbols.args$fg[which(names(symbols.args$fg)=="pos")]->colpos
        symbols.args$fg[which(names(symbols.args$fg)=="neg")]->colneg
        coltot<-NA
        if(!is.null(colneg)) coltot[which(single[,2]<=0.025)]<-colneg else coltot[which(single[,2]<=0.025)]<-"red"
        if(!is.null(colpos)) coltot[which(single[,2]>=0.975)]<-colpos else coltot[which(single[,2]>=0.975)]<-"royalblue"
        symbols.args$fg<-coltot
      }
      if(!"inches"%in%names(symbols.args)) symbols.args$inches<-0.25

      if(any(c("pos","neg")%in%names(symbols.args$bg))|is.null(symbols.args$bg)){
        symbols.args$bg[which(names(symbols.args$bg)=="pos")]->colpos1
        symbols.args$bg[which(names(symbols.args$bg)=="neg")]->colneg1
        coltot1<-NA
        if(!is.null(colneg1)) coltot1[which(single[,2]<=0.025)]<-colneg1 else coltot1[which(single[,2]<=0.025)]<-scales::alpha("red", 0.5)
        if(!is.null(colpos1)) coltot1[which(single[,2]>=0.975)]<-colpos1 else coltot1[which(single[,2]>=0.975)]<-scales::alpha("royalblue", 0.5)
        symbols.args$bg<-coltot1
      }

      if(!"bg"%in%names(symbols.args)) symbols.args$bg<-scales::alpha(symbols.args$fg, 0.5)
      if(any(c("squares","rectangles","stars","thermometers","boxplots")%in%names(symbols.args)))
        warning("The shape of the symbol cannot be modified",.immediate=TRUE)
      symbols.args<-symbols.args[which(!names(symbols.args)%in%c("squares","rectangles","stars","thermometers","boxplots"))]

      do.call(plot.phylo,c(x=list(tree),tree.args))
      xx<-sapply(rownames(single),function(w) get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[22]][as.numeric(w)])
      yy<-sapply(rownames(single),function(w) get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[23]][as.numeric(w)])
      do.call(symbols,c(list(x=xx,y=yy,add=TRUE,circles=abs(single[,1])^0.5),symbols.args))
    }

    return(list(plotClades=plotClades))
  }

  if(!is.null(SS$state.results)){
    if(is.null(state)) stop("Please provide the vector of species states")

    SS$state.results->stateres
    if(any(rownames(stateres)%in%unique(state)))
      stateres[which(rownames(stateres)%in%unique(state)),]->stateres
    stateres[order(rownames(stateres)),]->stateres
    stateres[which(stateres[,2]<=0.025|stateres[,2]>=0.975),]->statesign
    state <- treedataMatch(tree, state)[[1]][,1]
    SS$plotData[,match(rownames(stateres),colnames(SS$plotData)),drop=FALSE]->pldata
    pldata<-do.call(rbind,lapply(1:ncol(pldata),function(j)
      data.frame(v=pldata[,j],state=colnames(pldata)[j])))

    plotStates<-function(plot.args=NULL,points.args=NULL,legend.args=list()){
      # if(is.null(plot.args)) plot.args<-list()
      # if(is.null(points.args)) points.args<-list()

      if(!"xlab"%in%names(plot.args)) plot.args$xlab<-""
      if(!"ylab"%in%names(plot.args)) plot.args$ylab<-"random differences"
      if(!"ylim"%in%names(plot.args)) plot.args$ylim<-range(pldata[,1],stateres[,1])

      if(!"cex"%in%names(points.args)){
        points.args$cex<-rep(2.5,nrow(stateres))
        if((!"pch"%in%names(points.args))&nrow(statesign)>0)
          points.args$cex[which(!rownames(stateres)%in%rownames(statesign))]<-3
      }
      if(!"pch"%in%names(points.args)) {
        points.args$pch<-rep("\U2022",nrow(stateres))
        if(nrow(statesign)>0) points.args$pch[match(rownames(statesign),rownames(stateres))]<-"\U2605"
      }else if(length(points.args$pch)<nrow(stateres)) points.args$pch<-rep(points.args$pch,length.out=nrow(stateres))

      pch.type<-points.args$pch
      if(all(is.numeric(pch.type))){
        pch.type[which(points.args$pch>=20)]<-"bg"
        pch.type[which(points.args$pch<20)]<-"col"
      }
      if(!"col"%in%names(points.args)){
        points.args$col<-suppressWarnings(RColorBrewer::brewer.pal(nrow(stateres), "Set2"))
        if(any(pch.type=="bg")) points.args$col[which(pch.type=="bg")]<-"black"
      } else if(length(points.args$col)<nrow(stateres)) points.args$col<-rep(points.args$col,length.out=nrow(stateres))
      if(!"bg"%in%names(points.args)){
        points.args$bg<-suppressWarnings(RColorBrewer::brewer.pal(nrow(stateres), "Set2"))
        if(any(pch.type=="col")) points.args$bg[which(pch.type=="col")]<-"white"
      }else  if(length(points.args$bg)<nrow(stateres)) points.args$bg<-rep(points.args$bg,length.out=nrow(stateres))

      do.call(plot,c(x=list(as.factor(pldata$state)),y=list(pldata$v),plot.args))
      do.call(points,c(list(x=as.factor(rownames(stateres)),y=stateres[,1]),points.args))

      if(!is.null(legend.args)){
        if(!"x"%in%names(legend.args)) legend.args$x<-"topleft"
        if(!"legend"%in%names(legend.args))
          legend.args$legend<-rownames(stateres)
        if(!"fill"%in%names(legend.args)&!any(c("lty","lwd")%in%names(legend.args))){
          if(!"pch"%in%names(legend.args)){
            legend.args$pch<-points.args$pch
            if(!"pt.bg"%in%names(legend.args)) legend.args$pt.bg<-points.args$bg
            if(!"col"%in%names(legend.args)) legend.args$col<-points.args$col
          }

          if(!"pt.cex"%in%names(legend.args)) legend.args$pt.cex<-1.5
        }

        if(!"bty"%in%names(legend.args)) legend.args$bty<-"n"
        do.call(legend,legend.args)
      }
    }
    return(list(plotStates=plotStates))
  }

}

#'@export
addShift<-function(SS,symbols.args=NULL){
  if(is.null(SS$single.clades)&is.null(SS$state.results)){
    stop("No significant result available")
  }

  SS$single.clades->single
  if(any(c("pos","neg")%in%names(symbols.args$fg))|is.null(symbols.args$fg)){
    symbols.args$fg[which(names(symbols.args$fg)=="pos")]->colpos
    symbols.args$fg[which(names(symbols.args$fg)=="neg")]->colneg
    coltot<-NA
    if(!is.null(colneg)) coltot[which(single[,2]<=0.025)]<-colneg else coltot[which(single[,2]<=0.025)]<-"red"
    if(!is.null(colpos)) coltot[which(single[,2]>=0.975)]<-colpos else coltot[which(single[,2]>=0.975)]<-"royalblue"
    symbols.args$fg<-coltot
  }
  if(!"inches"%in%names(symbols.args)) symbols.args$inches<-0.25

  if(any(c("pos","neg")%in%names(symbols.args$bg))|is.null(symbols.args$bg)){
    symbols.args$bg[which(names(symbols.args$bg)=="pos")]->colpos1
    symbols.args$bg[which(names(symbols.args$bg)=="neg")]->colneg1
    coltot1<-NA
    if(!is.null(colneg1)) coltot1[which(single[,2]<=0.025)]<-colneg1 else coltot1[which(single[,2]<=0.025)]<-scales::alpha("red", 0.5)
    if(!is.null(colpos1)) coltot1[which(single[,2]>=0.975)]<-colpos1 else coltot1[which(single[,2]>=0.975)]<-scales::alpha("royalblue", 0.5)
    symbols.args$bg<-coltot1
  }

  if(!"bg"%in%names(symbols.args)) symbols.args$bg<-scales::alpha(symbols.args$fg, 0.5)
  if(any(c("squares","rectangles","stars","thermometers","boxplots")%in%names(symbols.args)))
    warning("The shape of the symbol cannot be modified",.immediate=TRUE)
  symbols.args<-symbols.args[which(!names(symbols.args)%in%c("squares","rectangles","stars","thermometers","boxplots"))]

  xx<-sapply(rownames(single),function(w) get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[22]][as.numeric(w)])
  yy<-sapply(rownames(single),function(w) get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[23]][as.numeric(w)])
  do.call(symbols,c(list(x=xx,y=yy,add=TRUE,circles=abs(single[,1])^0.5),symbols.args))
}
