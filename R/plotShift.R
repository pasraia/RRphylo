#'@title Graphical representation of search.shift results
#'@description \code{plotShift} generates customized functions to produce plots
#'  out of \code{\link{search.shift}} results.
#'  \code{addShift} adds circles to highlight shifting clades onto the currently
#'  plotted tree.
#'@aliases plotShift
#'@aliases addShift
#'@usage plotShift(RR,SS,state=NULL)
#'@usage addShift(SS,symbols.args=NULL)
#'@param RR the object produced by \code{\link{RRphylo}} used to perform
#'  \code{search.shift}.
#'@param SS an object produced by \code{search.shift}.
#'@param state if \code{search.shift} was performed under \code{status.type =
#'  "sparse"}, this is the same \code{state} vector provided to the function.
#'@param symbols.args as described for \code{$plotClades} below.
#'@details Using \code{...$plotClades()} or plotting the tree and applying
#'  \code{addShift()} returns the same plot. The latter function might be useful
#'  in combination with \code{\link{plotRR}} to add the shifts information to the
#'  branch-wise plot of evolutionary rate values.
#'@return The function returns a function to generate a customizable plot of
#'  \code{search.shift} results.
#'@return If \code{search.shift} was performed under \code{status.type =
#'  "clade"}, \code{plotShift} returns a \strong{\code{$plotClades}} function
#'  which highlights the shifting clades onto the phylogenetic tree. The usage
#'  is: \code{...$plotClades(tree.args=NULL,symbols.args=NULL)}, where
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
#'@return If \code{search.shift} was performed under \code{status.type =
#'  "sparse"}, \code{plotShift} returns a \strong{\code{$plotStates}} function
#'  which plots the states onto the phylogenetic trees. The usage is:
#'  \code{...$plotStates(tree.args=NULL,points.args=NULL,legend.args=list())},
#'  where \code{tree.args} is a list of further arguments passed to the function
#'  \code{plot.phylo}, \code{points.args} is a list of further arguments passed
#'  to the function \code{points}, and \code{legend.args} is a list of
#'  additional arguments passed to the function \code{legend} (if \code{= NULL}
#'  the legend is not plotted). If as many colors/pch values as the number of
#'  different states are provided in \code{points.args}, each of them is
#'  assigned to each states taken in the same order they are returned by
#'  \code{search.shift}.
#'@author Silvia Castiglione
#'@importFrom graphics symbols
#'@importFrom ape plot.phylo
#'@export
#'@seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#'@seealso \href{../doc/Plotting-tools.html}{\code{plotShift} vignette}
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
#' plotShift(RR=dinoRates,SS=SSauto)->plotSS
#' plotSS$plotClades()
#'
#' plot(dinoRates$tree)
#' addShift(SS=SSauto)
#'
#' search.shift(RR=dinoRates,status.type="clade",node=c(696,746))->SSnode
#' plotShift(RR=dinoRates,SS=SSnode)->plotSS
#' plotSS$plotClades(tree.args=list(no.margin=TRUE),
#'                   symbols.args=list(fg=NA,bg=c(pos="cyan",neg="magenta")))
#'
#'
#' search.shift(RR=dinoRates,status.type= "sparse",state=statedino)->SSstate
#' plotShift(RR=dinoRates,SS=SSstate,state=statedino)->plotSS
#' plotSS$plotStates(tree.args=list(no.margin=TRUE),
#'                   points.args=list(bg=c("gold","forestgreen","royalblue","white"),
#'                                    col=c("black","black","black","orangered"),
#'                                    pch=c(21,22,24,11)),legend.args=list())
#'     }

plotShift<-function(RR,SS,state=NULL){
  RR$tree->tree
  if(is.null(SS$single.clades)&is.null(SS$state.results)){
    stop("No significant result available")
  }

  if(!is.null(SS$single.clades)){
    SS$single.clades->single
    plotClades<-function(tree.args=NULL,symbols.args=NULL){

      if(Ntip(tree)>100){
        if(all(!c("show.tip.label","cex")%in%names(tree.args))) tree.args$show.tip.label<-FALSE
      }

      if(isTRUE(tree.args$no.margin)){
        mars <- par("mar")
        on.exit(par(mar = mars))
      }

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
    stateres[which(rownames(stateres)%in%unique(state)),]->stateres
    stateres[which(stateres[,2]<=0.025|stateres[,2]>=0.975),]->statesign
    state <- treedataMatch(tree, state)[[1]][,1]

    plotStates<-function(tree.args=NULL,points.args=NULL,legend.args=list()){
      if(Ntip(tree)>100){
        if(all(!c("show.tip.label","cex")%in%names(tree.args))) tree.args$show.tip.label<-FALSE
        if(!"type"%in%names(tree.args)) tree.args$type<-"fan"
      }

      if(isTRUE(tree.args$no.margin)){
        mars <- par("mar")
        on.exit(par(mar = mars))
      }

      if(!"pch"%in%names(points.args)) {
        points.args$pch<-16
      }else{
        if(length(points.args$pch)==nrow(stateres)){
          pchs<-state
          sapply(1:nrow(stateres),function(x) pchs[which(pchs==rownames(stateres)[x])]<<-points.args$pch[x])
          points.args$pch<-as.numeric(pchs)
        }
      }
      if(length(points.args$pch)<length(state)) rep(points.args$pch,length.out=length(state))->points.args$pch
      pch.type<-points.args$pch
      pch.type[which(points.args$pch>=20)]<-"bg"
      pch.type[which(points.args$pch<20)]<-"col"

      if(!"col"%in%names(points.args)){
        points.args$col<-suppressWarnings(RColorBrewer::brewer.pal(nrow(stateres), "Set2")[as.numeric(as.factor(state))])
      }else{
        if(length(points.args$col)==nrow(stateres)){
          colo<-state
          sapply(1:nrow(stateres),function(x) colo[which(colo==rownames(stateres)[x])]<<-points.args$col[x])
          points.args$col<-colo
        }
      }

      if(!"bg"%in%names(points.args)){
        points.args$bg<-suppressWarnings(RColorBrewer::brewer.pal(nrow(stateres), "Set2")[as.numeric(as.factor(state))])
      }else{
        if(length(points.args$bg)==nrow(stateres)){
          colo<-state
          sapply(1:nrow(stateres),function(x) colo[which(colo==rownames(stateres)[x])]<<-points.args$bg[x])
          points.args$bg<-colo
        }
      }

      if(any(pch.type=="bg")) points.args$col[which(pch.type=="bg")]<-"black"
      if(any(pch.type=="col")) points.args$bg[which(pch.type=="col")]<-"white"

      # leg.col<-sapply(1:length(pch.type),function(j) points.args[[grep(pch.type[j],names(points.args))]][j])

      do.call(plot.phylo,c(x=list(tree),tree.args))
      xx<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[22]][1:Ntip(tree)]
      yy<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)[[23]][1:Ntip(tree)]
      do.call(points,c(list(x=xx,y=yy),points.args))

      if(!is.null(legend.args)){
        if(!"x"%in%names(legend.args)) legend.args$x<-"topleft"
        if(!"legend"%in%names(legend.args))
          legend.args$legend<-sapply(1:nrow(statesign),function(x){
            if(statesign[x,2]<=0.025) paste(rownames(statesign)[x],"- negative shift") else
              paste(rownames(statesign)[x],"- positive shift")
          })
        if(!"fill"%in%names(legend.args)&!any(c("lty","lwd")%in%names(legend.args))){
          if(!"pch"%in%names(legend.args)) legend.args$pch<-points.args$pch[match(rownames(statesign),state)]
          legend.args$pch<-rep(legend.args$pch,length.out=nrow(statesign))
          if(any(legend.args$pch>20)&!"pt.bg"%in%names(legend.args)){
            legend.args$pt.bg<-rep(NA,length(legend.args$pch))
            legend.args$pt.bg[which(legend.args$pch>20)]<-sapply(match(rownames(statesign),state)[which(legend.args$pch>20)],function(k) ifelse(points.args$pch[k]<20,points.args$col[k],points.args$bg[k]))
          }
          if(any(legend.args$pch<21)&!"col"%in%names(legend.args)){
            legend.args$col<-rep(NA,length(legend.args$pch))
            legend.args$col[which(legend.args$pch<21)]<-sapply(match(rownames(statesign),state)[which(legend.args$pch<21)],function(k) ifelse(points.args$pch[k]<20,points.args$col[k],points.args$bg[k]))
            legend.args$col[which(is.na(legend.args$col))]<-"black"
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
