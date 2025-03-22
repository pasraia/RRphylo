#' @title Draw colorbar on a plot
#' @description The function adds a color bar to the current plot.
#' @usage colorbar(colors,x,y=NULL,direction="vertical",
#'   height=1,width=1,border="black",lwd=2,lty=1,
#'   labs=NULL,labs.pos=NULL,title=NULL,title.pos=NULL,
#'   ticks=TRUE,tck.pos=NULL,tck.length=1,xpd=FALSE,...)
#' @param colors vector of colors.
#' @param x,y the x and y coordinates where the bottom left corner of the bar is
#'   positioned. Keywords as in \code{\link[graphics]{legend}} are allowed.
#' @param direction either \code{"vertical"} or \code{"horizontal"}.
#' @param height a number indicating the amount by which the height of the bar
#'   should be scaled relative to the default.
#' @param width a number indicating the amount by which the width of the bar
#'   should be scaled relative to the default.
#' @param border color of the border around the bar. Set \code{NA} to suppress
#'   border drawing.
#' @param lwd border line width.
#' @param lty border line type.
#' @param labs the vector of labels to place next to the bar.
#' @param labs.pos either \code{"left"}/\code{"right"} for
#'   \code{direction="vertical"} or \code{"top"}/\code{"bottom"} for
#'   \code{direction="horizontal"}. Default settings are \code{"right"} and
#'   \code{"bottom"}.
#' @param title the title to be placed next to the bar.
#' @param title.pos either on the \code{"top"} or at the \code{"bottom"} of the
#'   bar. Default setting is \code{"top"}.
#' @param ticks logical indicating whether ticks should be drawn next to each
#'   label.
#' @param tck.pos indicates whether ticks should be plotter \code{"in"}side or
#'   \code{"out"}side the bar border.
#' @param tck.length tick lengths
#' @param xpd a value of the \code{\link[graphics]{par}} \code{xpd}.
#' @param ... further arguments passed to the functions \code{text} (for labels
#'   and title) and \code{segments}. All these arguments must be hooked to the
#'   element they refer to by indicating: \code{labs.}\* for labels,
#'   \code{title.}\* for title, and \code{tck.}\* for ticks. See example for
#'   further details.
#' @export
#' @author Silvia Castiglione
#' @examples
#'
#' rainbow(30)->cols
#' replicate(4,paste(sample(letters,4),collapse=""))->labs
#'
#' plot(rnorm(20),rnorm(20))
#' colorbar(cols,"topleft")
#'
#' plot(rnorm(20),rnorm(20))
#' colorbar(cols,"topright",
#'          height=1.2,width=1.2,lwd=2,
#'          labs=labs,labs.pos="left",labs.cex=1.3,labs.adj=1,
#'          title="Colorbar!",title.cex=1.4,title.font=2,title.adj=c(0,0),
#'          tck.pos="out",tck.lwd=2,xpd=TRUE)

colorbar<-function(colors,x,y=NULL,direction="vertical",
                   height=1,width=1,border="black",lwd=2,lty=1,
                   labs=NULL,labs.pos=NULL,title=NULL,title.pos=NULL,
                   ticks=TRUE,tck.pos=NULL,tck.length=1,xpd=FALSE,...){
  all.args<-list(...)
  labs.args<-all.args[grep("labs.",names(all.args))]
  labs.args<-labs.args[which(names(labs.args)!="labs.pos")]
  names(labs.args)<-gsub("labs.","",names(labs.args))
  title.args<-all.args[grep("title.",names(all.args))]
  title.args<-title.args[which(names(title.args)!="title.pos")]
  names(title.args)<-gsub("title.","",names(title.args))
  tck.args<-all.args[grep("tck.",names(all.args))]
  names(tck.args)<-gsub("tck.","",names(tck.args))
  tck.args<-tck.args[which(names(tck.args)!="tck.pos")]
  all.args<-all.args[-c(grep("title.",names(all.args)),grep("labs.",names(all.args)),grep("tck.",names(all.args)))]

  if(is.null(labs.args$cex)) labs.args$cex<-par("cex")
  if(is.null(title.args$cex)) title.args$cex<-par("cex")

  par("usr")->lims
  if (!is.null(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }

  (lims[4]-lims[3])->yscale
  (lims[2]-lims[1])->xscale

  leny<-ifelse(direction=="vertical",yscale/4,yscale/16)
  lenx<-ifelse(direction=="vertical",xscale/20,xscale/4)

  Htot<-H<-leny*height
  Wtot<-W<-lenx*width

  if(is.character(x)){
    if(!is.null(title)){
      Wtitle<-do.call(strwidth,c(list(s=title),title.args[which(names(title.args)!="adj")]))
      Htitle<-do.call(strheight,c(list(s=title),title.args[which(names(title.args)!="adj")]))
    }else Wtitle<-Htitle<-0

    if(direction=="vertical"){
      if(!is.null(labs)){
        do.call(strheight,c(list(s=labs),labs.args))/2->hlabs
        Htot<-H+sum(hlabs[c(1,length(hlabs))])

        max(do.call(strwidth,c(list(s=labs),labs.args)))/2->maxWlabs
        Wtot<-W+W*3/8+maxWlabs
      }else hlabs<-0

      top.adj<-H+hlabs[1]+Htitle
      bottom.adj<-hlabs[1]
      left.adj<-ifelse(Wtitle>Wtot,Wtitle-Wtot,0)
      right.adj<-Wtot

    }else{
      if(!is.null(labs)){
        do.call(strwidth,c(list(s=labs),labs.args))/2->wlabs
        # Wlabs<-W+sum(wlabs[c(1,length(wlabs))])
        # if(Wlabs>Wtot) Wtot<-Wlabs

        max(do.call(strheight,c(list(s=labs),labs.args))/2)->maxHlabs
        Htot<-H+H*3/8+maxHlabs
      }else wlabs<-0

      top.adj<-H+Htitle
      bottom.adj<-Htot-H
      left.adj<-wlabs[1]
      right.adj<-Wtot

    }

    # if(grepl("bottom",x)) ystart<-par("yaxp")[1]+bottom.adj else if(grepl("top",x)) ystart<-par("yaxp")[2]-top.adj else{
    #   ystart<-lims[3]+yscale/2-Htot/2
    # }
    # if(grepl("left",x)) xstart<-par("xaxp")[1]+left.adj else if(grepl("right",x)) xstart<-par("xaxp")[2]-right.adj else{
    #   xstart<-lims[1]+xscale/2-Wtot/2
    # }

    if(grepl("bottom",x)) ystart<-lims[3]+yscale/80+bottom.adj else
      if(grepl("top",x)) ystart<-lims[4]-yscale/80-top.adj else{
        ystart<-lims[3]+yscale/2-Htot/2
      }
      if(grepl("left",x)) xstart<-lims[1]+xscale/80+left.adj else
        if(grepl("right",x)) xstart<-lims[2]-xscale/80-right.adj else{
          xstart<-lims[1]+xscale/2-Wtot/2
        }

  }else{
    xstart<-x
    if(is.null(y)) ystart<-0 else ystart<-y
    maxWlabs<-maxHlabs<-0
  }

  if(direction=="vertical"){
    yend<-ystart+H
    xend<-xstart+W
    seq(ystart,yend,length.out=length(colors)+1)->yseq
    cbind(yseq[-length(yseq)],yseq[-1])->ymat
    matrix(c(xstart,xend),ncol=2,nrow=nrow(ymat),byrow = TRUE)->xmat

    if(!is.null(labs)){
      if(is.null(labs.pos)||labs.pos%in%c("bottom","top")) labs.pos<-"right"
      seq(ystart,yend,length.out=length(labs))->yleg

      if(labs.pos=="right")
        xend+W*3/8+maxWlabs->xleg else xstart-W*3/8-maxWlabs->xleg
      labs.args<-c(list(x=xleg,y=yleg,labels=labs),labs.args)

      if(ticks){
        tck.len<-W/4*tck.length
        if(is.null(tck.pos)) tck.pos<-"in"
        if(labs.pos=="right")
          ifelse(tck.pos=="out",c(xend,xend+tck.len)->postck,c(xend-tck.len,xend)->postck) else
            ifelse(tck.pos=="out",c(xstart-tck.len,xstart)->postck,c(xstart,xstart+tck.len)->postck)

        matrix(postck,ncol=2,nrow=length(yleg),byrow = TRUE)->xtck
        cbind(yleg,yleg)->ytck
      }
    }
  }else{
    xend<-xstart+lenx*width
    yend<-ystart+leny*height
    seq(xstart,xend,length.out=length(colors)+1)->xseq
    cbind(xseq[-length(xseq)],xseq[-1])->xmat
    matrix(c(ystart,yend),ncol=2,nrow=nrow(xmat),byrow = TRUE)->ymat

    if(!is.null(labs)){
      if(is.null(labs.pos)||labs.pos%in%c("left","right")) labs.pos<-"bottom"
      seq(xstart,xend,length.out=length(labs))->xleg

      if(labs.pos=="bottom") ystart-H*3/8-maxHlabs->yleg else yend+H*3/8+maxHlabs->yleg
      labs.args<-c(list(x=xleg,y=yleg,labels=labs),labs.args)

      if(ticks){
        tck.len<-H/4*tck.length
        if(is.null(tck.pos)) tck.pos<-"in"
        if(labs.pos=="bottom")
          ifelse(tck.pos=="out",c(ystart-tck.len,ystart)->postck,c(ystart,ystart+tck.len)->postck) else
            ifelse(tck.pos=="out",c(yend,yend+tck.len)->postck,c(yend-tck.len,yend)->postck)

        matrix(postck,ncol=2,nrow=length(xleg),byrow = TRUE)->ytck
        cbind(xleg,xleg)->xtck

      }
    }
  }

  if(!is.null(title)){
    if(direction=="vertical"&!is.null(labs.pos)&&labs.pos=="left")
      xtitle<-xstart-Wtot/2 else xtitle<-xstart+Wtot/2
      if(is.null(title.pos)||title.pos=="top")
        ytitle<- ifelse(grepl("top",x),mean(c(lims[4],yend)),yend + yscale/20) else
          ytitle<-ystart-yscale/20
        title.args<-c(list(x=xtitle,y=ytitle,labels=title),title.args)
  }


  if(!is.na(border)){
    all.args<-c(list(xleft=xstart,ybottom=ystart,xright=xend,ytop=yend,border=border,lwd=lwd),all.args)
    do.call(rect,all.args)
  }
  sapply(1:nrow(xmat),function(i) rect(xmat[i,1],ymat[i,1],xmat[i,2],ymat[i,2],
                                       col=colors[i],border=NA))


  if(!is.null(labs)){
    do.call(text,labs.args)
    if(ticks){
      sapply(1:nrow(xtck),function(i)
        do.call(segments,c(list(x0=xtck[i,1],y0=ytck[i,1],x1=xtck[i,2],y1=ytck[i,2]),tck.args)))
    }
  }
  if(!is.null(title)) do.call(text,title.args)

}
