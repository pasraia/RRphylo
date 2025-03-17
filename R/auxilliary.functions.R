ftipTree<-function(tr,N,tar.tips,ftiplen){
  i = 1
  while (i <= length(N)) {
    nn <- getMRCA(tr, tar.tips[[i]])
    tr$edge.length[which(tr$edge[,2]==nn)]->edlen
    # if(edlen<=0.001) edlen/10->pp else 0.001->pp
    pp<-1e-8
    tr <- bind.tip(tr, tip.label = paste0("fnd",N[i]),
                   edge.length = ftiplen,
                   where = nn,position = pp)
    i = i + 1
  }
  return(tr)
}

rootVfun<-function(tree,toriginal,yoriginal){
  if (!is.binary(tree)) u <- data.frame(yoriginal, (1/diag(vcv(tree))^2)) else
    u <- data.frame(yoriginal,(1/diag(vcv(toriginal))^2))
  u <- u[order(u[, ncol(u)], decreasing = TRUE),]
  u1 <- u[1:(nrow(u) * 0.1), ,drop=FALSE]
  rv <- apply(u1[, 1:(ncol(u1)-1),drop=FALSE],2,function(x)
    weighted.mean(x,u1[, dim(u1)[2]]))
  return(rv)
}

RRcore<-function(lambda,y,rootV,L,L1,Lprod,tr,y1=NULL){
  if(!is.null(y1)){
    LX<-L
    LX1<-L1
    betas <- (solve(Lprod + lambda * diag(ncol(LX))) %*%
                t(LX)) %*% (as.matrix(y)-rootV)

    aceRR <- ((LX1 %*% betas[c(1:Nnode(tr),(length(betas)+1-ncol(y1)):length(betas)), ]))+rootV
    y.hat <- (LX %*% betas)+rootV
    betas[(length(betas)+1-ncol(y1)):length(betas),,drop=FALSE]->x1.rate
    betas[1:(length(betas)-ncol(y1)),,drop=FALSE]->betas
    colnames(betas)<-NULL
    list(aceRR,betas,y.hat,x1.rate)
  }else{
    betas <- (solve(Lprod + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(y) - rootV)
    aceRR <- (L1 %*% betas[1:Nnode(tr), ]) + rootV
    y.hat <- (L %*% betas) + rootV
    list(aceRR,betas,y.hat)
  }
}

covRates<-function(cov,betas){
  cov[match(rownames(betas),names(cov))]->cov
  Y <- abs(cov)
  R <- log(abs(betas))
  #### Covariate multi ####
  if (length(which(apply(betas, 1, sum) == 0)) > 0) {
    zeroes <- which(apply(betas, 1, sum) == 0)
    R <- R[-zeroes, ]
    Y <- Y[-zeroes]
    res <- residuals(lm(R ~ Y))
    factOut <- which(apply(betas, 1, sum) != 0)
    betas[factOut, ] <- res
    betas[zeroes, ] <- 0
  } else { #### Covariate uni ####
    res <- residuals(lm(R ~ Y))
    betas <- as.matrix(res)
  }
  return(betas)
}

phylo.run.test<-function(tree,state,st,nsim=100){
  cophenetic.phylo(tree)->cop
  cop[which(state==st),which(state==st)]->subcop
  mean(subcop[upper.tri(subcop)])->mds
  nrow(subcop)->sl
  r.mds<-replicate(nsim,{
    sample(tree$tip.label,sl)->test.tip
    cop[test.tip,test.tip]->r.cop
    mean(r.cop[upper.tri(r.cop)])
  })
  return(list(p=length(which(r.mds<mds))/nsim))
}

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

unitV <- function(x,na.rm=FALSE) sum(x^2,na.rm=na.rm)^0.5

deg2rad <- function(deg) (deg * pi)/(180)

rad2deg <- function(rad)  (rad * 180)/(pi)

angle.vecs<-function(vec1,vec2){
  ((vec1%*%vec2)/(unitV(vec1) *unitV(vec2)))->ppA
  if((ppA-1)>0) ppA<-1 else if((1+ppA)<0) ppA<-(-1)
  rad2deg(acos(ppA))
}

traitgram <- function(x, phy,
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

  xanc <- ape::ace(xx, phy, method=method)$ace
  xall = c(xx,xanc)

  a0 = max(ages)-ages[phy$edge[,1]]
  a1 = max(ages)-ages[phy$edge[,2]]
  x0 = xall[phy$edge[,1]]
  x1 = xall[phy$edge[,2]]


  if (show.names) {
    maxNameLength = max(nchar(names(xx)))
    ylim = c(min(ages),max(ages)*(1+maxNameLength/50))
    if (!underscore) names(xx) = gsub('_',' ',names(xx))
  } else ylim = range(ages)

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

  plot.args<-list(...)

  if(!"xlab"%in%names(plot.args)) plot.args$xlab<-""
  if(!"ylab"%in%names(plot.args)) plot.args$ylab<-""
  if(!"bty"%in%names(plot.args)) plot.args$bty<-"n"
  if(!"cex.axis"%in%names(plot.args)) plot.args$cex.axis<-0.8
  if(!"yaxt"%in%names(plot.args)) plot.args$yaxt<-"n"
  if(!"lwd"%in%names(plot.args)) plot.args$lwd<-2

  do.call(plot,c(list(x=NA,xlim=rev(range(c(a0,a1))),ylim=range(c(x0,x1))),plot.args))
  do.call(segments,c(list(x0=a0,y0=x0,x1=a1,y1=x1,col=colo),
                     plot.args[which(names(plot.args)%in%c("lwd","lty","lend","ljoin","lmitre"))]))


  if (show.names) {
    text(max(ages),sort(xx),
         labels = names(xx)[order(xx)],
         adj = -0,
         srt=90,
         cex=0.3)
  }

  return(data.frame(a1,x1))
}

Plot_ConvexHull<-function(xcoord, ycoord,...){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])

  poly.args<-list(...)
  do.call(polygon,c(list(x=xcoord[hpts],y=ycoord[hpts]),poly.args))
}

dosur <- function(scores,pcs,sel=NULL,mshape,radius=0){
  if(is.null(sel)==TRUE) {
    temp<-Morpho::restoreShapes(scores,pcs,mshape)
  } else {
    temp<-Morpho::restoreShapes(scores[sel],pcs[,sel],mshape)
  }
  mshape<-Rvcg::vcgBallPivoting(mshape, radius = radius)
  sur<-mshape
  sur$vb[1:3,]<-t(temp)
  return(sur)

}

areadiff<-function(mesh1,mesh2,out.rem=FALSE,scale01=TRUE,fact=1.5){

  area_shape1<-Rvcg::vcgArea(mesh1,perface=T)$pertriangle
  area_shape2<-Rvcg::vcgArea(mesh2,perface=T)$pertriangle
  diff_areas<-(area_shape1-area_shape2)/area_shape1
  sel<-which(is.na(diff_areas))

  if(length(sel)>0){
    mesh1$it<-mesh1$it[,-sel]
    mesh2$it<-mesh2$it[,-sel]
    mesh1<-Morpho::rmUnrefVertex(mesh1)
    mesh2<-Morpho::rmUnrefVertex(mesh2)
    area_shape1<-Rvcg::vcgArea(mesh1,perface=T)$pertriangle
    area_shape2<-Rvcg::vcgArea(mesh2,perface=T)$pertriangle
    diff_areas<-(area_shape1-area_shape2)/area_shape1
  }

  if(out.rem==TRUE){
    x=diff_areas
    qq <- quantile(x, c(1,3)/4, names=FALSE)
    r <- diff(qq) * fact
    tst <- x < qq[1] - r | x > qq[2] + r
    tstp<-qq[2] + r
    tstn<-qq[1] - r
    diff_areas[x>tstp]<-tstp
    diff_areas[x<tstn]<-tstn
  }else diff_areas=diff_areas

  if(scale01==TRUE) diff_areas<-range01(diff_areas)

  return(list("ash1"=area_shape1,"ash2"=area_shape2,"dareas"=diff_areas))
}

localmeshdiff <- function(mesh1, mesh2, ploton, paltot = rainbow(200),
                          from = NULL, to = NULL, n.int = 200, out.rem = TRUE, fact = 1.5,
                          visual = 1, scale01 = TRUE, vec = NULL, plot = FALSE,
                          densityplot = TRUE) {


  if(is.null(vec)) {
    areadiff(mesh1,mesh2,out.rem=out.rem,scale01=scale01,fact=fact)->adiff
    area_shape1<-adiff$ash1
    area_shape2<-adiff$ash2
    diff_areas<-adiff$dareas
  } else diff_areas<-vec

  cat("the range of diff_areas is ", range(diff_areas),
      sep = "\n")
  if (is.null(to) == TRUE) {
    to <- max(diff_areas) * 1.01
  }
  if (is.null(from) == TRUE) {
    from <- min(diff_areas) * 1.01
  }
  selfromto <- which(diff_areas < to & diff_areas >= from)

  colmap_tot <- colorRampPalette(paltot)
  breaks_tot <- cut(c(from, diff_areas, to), n.int)
  cols_tot <- colmap_tot(n.int)[breaks_tot]
  cols_tot <- cols_tot[-c(1, length(cols_tot))]
  selfromto <- which(diff_areas < to & diff_areas >= from)
  cols_tot[-selfromto] <- "#FFFFFF"

  if (isTRUE(densityplot)) {
    plot(density(c(from, diff_areas, to)), main = "",
         xlab = "", ylab = "")
    abline(v = seq(from, to, length.out = n.int), col = colmap_tot(n.int),
           lwd = 5)
    points(density(diff_areas), type = "l", lwd = 2)
  }
  if (ploton == 1) {
    meshtobeplotted <- mesh1
  }
  if (ploton == 2) {
    meshtobeplotted <- mesh2
  }
  if (plot == TRUE) {
    if (visual == 1) {
      rgl::triangles3d(t(meshtobeplotted$vb[, meshtobeplotted$it]),
                       col = rep(cols_tot, each = 3), alpha = 1, lit = TRUE,
                       specular = "black")
    }
    if (visual == 2) {
      rgl::triangles3d(t(meshtobeplotted$vb[, meshtobeplotted$it]),
                       col = rep(cols_tot, each = 3), alpha = 1, lit = TRUE,
                       specular = "black")
      rgl::wire3d(meshtobeplotted, col = "blue",
                  lit = FALSE)
    }
  }
  diffvert <- NULL
  track_v <- Rvcg::vcgVFadj(meshtobeplotted)
  for (i in 1:length(track_v)) {
    diffvert[i] <- mean(diff_areas[track_v[[i]]])
  }
  selfromto <- which(diffvert < to & diffvert >= from)
  colmap_tot <- colorRampPalette(paltot)
  breaks_tot <- cut(c(from, diffvert, to), n.int)
  cols_tot <- colmap_tot(n.int)[breaks_tot]
  cols_tot <- cols_tot[-c(1, length(cols_tot))]
  selfromto <- which(diffvert < to & diffvert >= from)
  cols_tot[-selfromto] <- "#FFFFFF"
  mesh <- meshtobeplotted
  mesh$material$color <- cols_tot
  if (!is.null(vec)) list(mesh = mesh) else
    return(list(ash1 = area_shape1, ash2 = area_shape2,
                dareas = diff_areas, mesh = mesh))
}

