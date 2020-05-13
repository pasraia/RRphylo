range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unitV <- function(x) sum(x^2)^0.5

deg2rad <- function(deg) (deg * pi)/(180)

rad2deg <- function(rad)  (rad * 180)/(pi)

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
    dato[order(dato[,5],decreasing=TRUE),]->dato
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
