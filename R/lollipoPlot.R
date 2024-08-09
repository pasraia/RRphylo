#' @title Lollipop charts
#' @description The function generates lollipop or dumbbell dots charts.
#' @param values either a vector, matrix, or data.frame of data. If matrix or
#'   data.frame including two columns, a dumbbell dots chart is plotted.
#' @param type plot direction, either vertical ("v", the default) or horizontal
#'   ("h").
#' @param pt.lwd points lwd
#' @param pt.col points color
#' @param ... other arguments passed to the functions \code{plot},
#'   \code{points}, and \code{segments}.
#' @export
#' @details If a dumbbell dots chart is plotted, different parameters (i.e.
#'   col/cex/pch/bg/lwd) for starting and ending points can be supplied. See
#'   example for further details.
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @examples
#' require(emmeans)
#'
#' lollipoPlot(values=feedlot[,4],pt.col="green",pt.lwd=2,lwd=0.8,col="gray20",
#'             ylab="swt",xlab="samples")
#'
#' line.col<-sample(colors()[-1],length(levels(feedlot[,1])))
#' line.col<-rep(line.col,times=table(feedlot[,1]))
#'
#' lollipoPlot(values=feedlot[order(feedlot[,1]),3],ylab="ewt",xlab="samples",
#'             bg=as.numeric(as.factor(feedlot[order(feedlot[,1]),2])),
#'             cex=1.2,pch=21,col=line.col)
#'
#'
#'
#' lollipoPlot(values=feedlot[order(feedlot[,1]),3:4],type="h",ylab="ewt",xlab="samples",
#'             pt.col=c("blue","cyan"),cex=1.2,pch=c(3,4),col=line.col)
#'
#' lollipoPlot(values=feedlot[order(feedlot[,1]),3:4],type="h",ylab="ewt",xlab="samples",
#'             bg=cbind(line.col,line.col),cex=c(1.2,1),pch=c(21,22))


lollipoPlot<-function(values,type="v",pt.lwd=NULL,pt.col=NULL,...){
  if(!type%in%c("v","h")) stop("'type' must be one of 'h' or 'v'")

  deparse1(substitute(values))->namval
  as.matrix(values)->values

  if(ncol(values)>1){
    values[,2]->endp
    values[,1]->startp
    lim<-range(c(startp,endp))
  }else{
    rep(0,nrow(values))->startp
    values[,1]->endp
    lim<-range(endp)
  }

  list(...)->argum
  if(any(names(argum)%in%c("lendp","lwd","lty","col"))){
    argum[which(names(argum)%in%c("lendp","lwd","lty","col"))]->sargs
    argum[-which(names(argum)%in%c("lendp","lwd","lty","col"))]->argum

    sapply(1:length(sargs),function(k){
      if(length(sargs[[k]])<length(endp)) rep(sargs[[k]],length.out=length(endp))->>sargs[[k]]
    })

    lapply(1:length(sargs[[1]]),function(j) lapply(sargs,"[[",j))->seg.args

  }else seg.args<-replicate(length(endp),list())
  argum->argum1->argum2

  if(!is.null(pt.lwd)){
    if(is.null(ncol(pt.lwd))){
      if(ncol(values)>1) matrix(pt.lwd,ncol=2)->pt.lwd else matrix(pt.lwd,ncol=1)->pt.lwd
    }
    if(ncol(pt.lwd)>1) {
      argum1$lwd<-pt.lwd[,1]
      argum2$lwd<-pt.lwd[,2]
    } else argum1$lwd<-argum2$lwd<-pt.lwd[,1]
  }

  if(!is.null(pt.col)){
    if(is.null(ncol(pt.col))){
      if(ncol(values)>1) matrix(pt.col,ncol=2)->pt.col else matrix(pt.col,ncol=1)->pt.col
    }
    if(ncol(pt.col)>1) {
      argum1$col<-pt.col[,1]
      argum2$col<-pt.col[,2]
    } else argum1$col<-argum2$col<-pt.col[,1]
  }

  if(!is.null(argum$cex)){
    if(is.null(ncol(argum$cex))){
      if(ncol(values)>1) matrix(argum$cex,ncol=2)->argum$cex else matrix(argum$cex,ncol=1)->argum$cex
    }
    if(ncol(argum$cex)>1) {
      argum1$cex<-argum$cex[,1]
      argum2$cex<-argum$cex[,2]
    } else argum1$cex<-argum2$cex<-argum$cex[,1]
  }

  if(!is.null(argum$pch)){
    if(is.null(ncol(argum$pch))){
      if(ncol(values)>1) matrix(argum$pch,ncol=2)->argum$pch else matrix(argum$pch,ncol=1)->argum$pch
    }
    if(ncol(argum$pch)>1) {
      argum1$pch<-argum$pch[,1]
      argum2$pch<-argum$pch[,2]
    } else argum1$pch<-argum2$pch<-argum$pch[,1]
  }

  if(!is.null(argum$bg)){
    if(is.null(ncol(argum$bg))){
      if(ncol(values)>1) matrix(argum$bg,ncol=2)->argum$bg else matrix(argum$bg,ncol=1)->argum$bg
    }
    if(ncol(argum$bg)>1){
      argum1$bg<-argum$bg[,1]
      argum2$bg<-argum$bg[,2]
    } else argum1$bg<-argum2$bg<-argum$bg[,1]
  }

  if(type=="v"){
    y<-endp
    y1<-startp
    x1<-x<-1:length(endp)
    xlim<-c(0,length(endp))
    ylim<-lim

    if(!"xaxt"%in%names(argum)||argum$xaxt!="n"){
      argum$xaxt<-"n"
      plot.axis<-TRUE
      ax<-1
    }else{
      if(argum$xaxt=="n") plot.axis<-FALSE
    }

    if(!"xlab"%in%names(argum)) argum$xlab<-""
    if(!"ylab"%in%names(argum)) argum$ylab<-namval
  } else{
    x<-endp
    x1<-startp
    y1<-y<-1:length(endp)
    ylim<-c(0,length(endp))
    xlim<-lim
    if(!"yaxt"%in%names(argum)||argum$yaxt!="n"){
      argum$yaxt<-"n"
      plot.axis<-TRUE
      ax<-2
    }else{
      if(argum$yaxt=="n") plot.axis<-FALSE
    }
    if(!"ylab"%in%names(argum)) argum$ylab<-""
    if(!"xlab"%in%names(argum)) argum$xlab<-namval
  }

  if(!"xlim"%in%names(argum)) argum$xlim<-xlim
  if(!"ylim"%in%names(argum)) argum$ylim<-ylim
  if(!"bty"%in%names(argum)) argum$bty<-"n"

  do.call(plot,c(x=NA,argum))
  if(type=="v")
    sapply(1:length(endp),function(k) do.call(segments,c(x0=k,y0=unname(startp[k]),x1=k,y1=unname(endp[k]),seg.args[[k]]))) else
      sapply(1:length(endp),function(k) do.call(segments,c(x0=unname(startp[k]),y0=k,x1=unname(endp[k]),y1=k,seg.args[[k]])))
  do.call(points,c(x=list(x),y=list(y),argum1))
  if(!all(startp==0)) do.call(points,c(x=list(x1),y=list(y1),argum2))

  if(plot.axis){
    axis(ax,at= 1:length(endp),
         labels = names(endp),las=2,tick=FALSE)
  }
}
