#' @title Mapping morphological convergence on 3D surfaces
#' @description Given vectors of RW (or PC) scores, the function selects the
#'   RW(PC) axes which best account for convergence and maps convergent areas on
#'   the corresponding 3D surfaces.
#' @usage conv.map(dataset,pcs,mshape,conv=NULL, exclude=NULL,out.rem=TRUE,
#'   show.consensus=FALSE, plot=TRUE,col="blue",names = TRUE)
#' @param dataset data frame (or matrix) with the RW (or PC) scores of the group
#'   or species to be compared.
#' @param pcs RW (or PC) vectors (eigenvectors of the covariance matrix) of all
#'   the samples.
#' @param mshape the Consensus configuration.
#' @param conv a named character vector indicating convergent species as
#'   (indicated as "conv" in \code{dataset}) and not convergent species
#'   (indicated as "noconv").
#' @param exclude integer: the index number of the RW (or PC) to be excluded
#'   from the comparison.
#' @param out.rem logical: if \code{TRUE} triangles with outlying area
#'   difference are removed.
#' @param show.consensus logical: if \code{TRUE}, the Consensus configuration is
#'   included in the comparison.
#' @param plot logical: if \code{TRUE}, the pairwise comparisons are be plotted.
#'   For more than 5 pairwise comparisons, the plot is not shown.
#' @param col character: the colour for the plot.
#' @param names logical: if \code{TRUE}, the names of the groups or species are
#'   displayed in the 3d plot.
#' @details \code{conv.map} automatically builds a 3D mesh on the mean shape
#'   calculated from the Relative Warp Analysis (RWA) or Principal Component
#'   Analysis (PCA) (\cite{Schlager 2017}) by applying the function
#'   \code{\link[Rvcg]{vcgBallPivoting}} (\pkg{Rvcg}). \code{conv.map} further gives
#'   the opportunity to exclude some RW (or PC) axes from the analysis because,
#'   for example, in most cases the first axes are mainly related to high-order
#'   morphological differences driven by phylogeny and size variations.
#'   \code{conv.map} finds and plots the strength of convergence on 3D surfaces.
#'   An output of \code{conv.map} (if the dataset contains a number equal or
#'   lower then 5 items) is an interactive plot mapping the convergence on the
#'   3D models. In the upper triangle of the 3D multiple layouts the rows
#'   representing the reference models and the columns the target models. On the
#'   contrary, on the lower triangle the rows correspond to the target models
#'   and the columns the reference models. In the calculation of the differences
#'   of areas we supply the possibility to find and remove outliers from the
#'   vectors of areas calculated on the reference and target surfaces. We
#'   suggest considering this possibility if the mesh may contain degenerate
#'   facets.
#' @export
#' @seealso \href{../doc/search.conv.html}{\code{search.conv} vignette} ; \code{\link[Morpho]{relWarps}} ;
#'   \code{\link[Morpho]{procSym}}
#' @importFrom grDevices colorRampPalette rainbow
#' @return The function returns a list including:
#'   \itemize{\item\strong{$angle.compare} data frame including the real angles
#'   between the given shape vectors, the angles conv computed between vectors
#'   of the selected RWs (or PCs), the angles between vectors of the
#'   non-selected RWs (or PCs), the difference conv, and its p values.
#'   \item\strong{$selected.pcs} RWs (or PCs) axes selected for convergence.
#'   \item\strong{$average.dist} symmetric matrix of pairwise distances between
#'   3D surfaces. \item\strong{$suface1} list of coloured surfaces, if two
#'   meshes are given, it represents convergence between mesh A and B charted on
#'   mesh A. \item\strong{$suface2} list of coloured surfaces, if two meshes are
#'   given, it represents convergence between mesh A and B charted on mesh B.
#'   \item \strong{$scale} the value used to set the colour gradient, computed
#'   as the maximum of all differences between each surface and the mean shape.}
#' @author Marina Melchionna, Antonio Profico, Silvia Castiglione, Carmela
#'   Serio, Gabriele Sansalone, Pasquale Raia
#' @references Schlager, S. (2017). \emph{Morpho and Rvcg–Shape Analysis in R:
#'   R-Packages for geometric morphometrics, shape analysis and surface
#'   manipulations.} In: Statistical shape and deformation analysis. Academic
#'   Press.
#' @examples
#'   \dontrun{
#'   data(DataSimians)
#'   DataSimians$pca->pca
#'
#'   ## Case 1. Convergent species only
#'      dato<-pca$PCscores[c(1,4),]
#'
#'      CM<-conv.map(dataset = dato,
#'                  pcs = pca$PCs,
#'                  mshape = pca$mshape,
#'                  show.consensus = TRUE)
#'
#'   ## Case 2. Convergent and non-convergent species
#'      dato<-pca$PCscores[c(1,4,7),]
#'      conv<-c("conv","conv","noconv")
#'      names(conv)<-rownames(dato)
#'
#'      CM<-conv.map(dataset = dato,
#'                   pcs = pca$PCs,
#'                   mshape = pca$mshape,
#'                   conv = conv,
#'                   show.consensus = TRUE,
#'                   col = "orange")
#'   }


conv.map<-function(dataset,pcs,mshape,
                   conv=NULL,exclude=NULL,
                   out.rem=TRUE,show.consensus=FALSE,
                   plot=TRUE,col="blue",names = TRUE){
  # require(inflection)
  # require(ddpcr)
  # require(rgl)
  # require(Rvcg)
  # require(Morpho)


  if (!requireNamespace("inflection", quietly = TRUE)) {
    stop("Package \"inflection\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("ddpcr", quietly = TRUE)) {
    stop("Package \"ddpcr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package \"rgl\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("Rvcg", quietly = TRUE)) {
    stop("Package \"Rvcg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("Morpho", quietly = TRUE)) {
    stop("Package \"Morpho\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(is.null(conv)) {
    conv<-rep("conv",nrow(dataset))
    names(conv)<-rownames(dataset)

  } else  conv<-conv[match(rownames(dataset),names(conv))]


  if(isTRUE(show.consensus)) {
    dataset<-rbind(dataset,rep(10^-6,ncol(dataset)))
    rownames(dataset)[nrow(dataset)]<-"consensus"
    conv<-c(conv,"noconv")
    names(conv)[length(conv)]<-"consensus"
  }

  df<-dataset
  colnames(df)<-sapply(1:ncol(df),function(x) paste("S",x,sep=""))
  comb.df<-combn(rownames(df),2)
  comb.c.df<-combn(names(which(conv=="conv")),2)

  selected<-list()

  for(k in 1:ncol(comb.c.df)){

    vec1<-df[match(comb.c.df[1,k],rownames(df)),]
    vec2<-df[match(comb.c.df[2,k],rownames(df)),]

    # if(is.null(names(vec1)))  names(vec1)<-paste("S",1:length(vec1),sep="");  names(vec2)<-paste("S",1:length(vec2),sep="")
    # if(is.null(names(vec2)))  names(vec1)<-paste("S",1:length(vec1),sep="");  names(vec2)<-paste("S",1:length(vec2),sep="")

    if(is.null(exclude)==FALSE) {
      v.vec1<-vec1[-exclude]
      v.vec2<-vec2[-exclude]
    } else {
      v.vec1<-vec1
      v.vec2<-vec2}
    rad2deg<-function(rad) (rad * 180)/(pi)
    unitV<-function (x)  sum(x^2)^0.5

    j.ang<-array()
    for(i in 1:length(v.vec1)){
      v.vec1[-i]->v1
      v.vec2[-i]->v2
      rad2deg(acos((v1 %*% v2)/(unitV(v1) * unitV(v2))))->j.ang[i]
    }

    names(j.ang)<-names(v.vec1)
    inflection::ede(seq(1:length(j.ang)),j.ang[order(j.ang,decreasing=TRUE)],0)[1]->cutter
    if(cutter==1){
      inflection::ede(seq(1:length(j.ang[-1])),j.ang[order(j.ang,decreasing=TRUE)][-1],0)[1]->cutter2 ### CHANGED IN V.2.4.10
      cutter2+1->cutter
    }
    if(cutter>.5*length(j.ang)) cutter=round(0.5*(length(j.ang)),0)

    j.ang[order(j.ang,decreasing=TRUE)][seq(1:cutter)]->main
    names(main)->sele

    selected[[k]] <- as.numeric(unlist(lapply(strsplit(sele,"S"),"[[",2)))

  }

  sele0<-table(unlist(selected))
  if(length(which(sele0>1))>1) selected<-as.numeric(names(sele0[which(sele0>1)])) else selected<-as.numeric(names(sele0))

  surfs<-list()
  surfs1<-list()
  surfsd<-list()
  surfsd1<-list()
  res<-list()
  res2<-list()

  selected.S<-sapply(selected,function(x) paste("S",x,sep = ""))

  for(k in 1:ncol(comb.df)){

    vec1<-df[match(comb.df[1,k],rownames(df)),]
    vec2<-df[match(comb.df[2,k],rownames(df)),]

    # if(is.null(names(vec1)))  names(vec1)<-paste("S",1:length(vec1),sep="");  names(vec2)<-paste("S",1:length(vec2),sep="")
    # if(is.null(names(vec2)))  names(vec1)<-paste("S",1:length(vec1),sep="");  names(vec2)<-paste("S",1:length(vec2),sep="")

    if(is.null(exclude)==FALSE) {
      v.vec1<-vec1[-exclude]
      v.vec2<-vec2[-exclude]
    } else {
      v.vec1<-vec1
      v.vec2<-vec2}

    ang<- rad2deg(acos((v.vec1 %*% v.vec2)/(unitV(v.vec1) * unitV(v.vec2))))
    a.sel<-rad2deg(acos((v.vec1[match(selected.S,names(v.vec1))] %*% v.vec2[match(selected.S,names(v.vec2))])/(unitV(v.vec1[match(selected.S,names(v.vec1))]) * unitV(v.vec2[match(selected.S,names(v.vec2))]))))
    a.others<-rad2deg(acos((v.vec1[-match(selected.S,names(v.vec1))] %*% v.vec2[-match(selected.S,names(v.vec2))])/(unitV(v.vec1[-match(selected.S,names(v.vec1))]) * unitV(v.vec2[-match(selected.S,names(v.vec2))]))))
    (a.sel-a.others)[1]->ang.diff

    length(selected.S)->nn
    inn<-array()
    sell<-list()
    for(i in 1:10000){
      v.vec1[sample(seq(1:length(v.vec1)),nn)]->x1
      v.vec2[which(names(v.vec2)%in%names(x1))]->x2
      v.vec1[-which(names(v.vec1)%in%names(x1))]->xn1
      v.vec2[-which(names(v.vec2)%in%names(x1))]->xn2
      rad2deg(acos((x1 %*% x2)/(unitV(x1) * unitV(x2))))->a1
      rad2deg(acos((xn1 %*% xn2)/(unitV(xn1) * unitV(xn2))))->a2
      list(selected=names(x1),angs=c(a1,a2),d=a1-a2)->sell[[i]]
      if(length(which(sell[[i]]$selected%in%selected.S))==nn) inn[i]<-0 else inn[i]<-1
    }

    if(length(which(inn==0))!=0) sell[-which(inn==0)]->sell
    1-length(which(unlist(lapply(sell, "[[", 3))>ang.diff))/length(sell)->p.sell

    sur.mshape<-Rvcg::vcgBallPivoting(mshape, radius = 0)
    sur <- dosur(scores = vec1, pcs = pcs, sel = selected, mshape = mshape, radius = 0)
    sur1 <- dosur(scores = vec2, pcs = pcs, sel = selected, mshape = mshape, radius = 0)

    mapD.sur<-areadiff(sur,sur.mshape,out.rem = out.rem,fact = 1.5)
    mapD.sur1<-areadiff(sur1,sur.mshape,out.rem = out.rem,fact = 1.5)

    msurD<-mapD.sur[[3]]-mapD.sur1[[3]]
    msurD1<-mapD.sur1[[3]]-mapD.sur[[3]]
    thr<-max(abs(c(msurD,msurD1)))

    res[[k]]<-list(angle.compare=c(real.angle=ang[1],selected=a.sel,others=a.others,ang.diff=ang.diff,p.value=p.sell))
    res2[[k]]<-list(average.dist=mean(abs(c(msurD,msurD1))),thr=thr)

    surfs[[k]]<-sur
    surfs1[[k]]<-sur1
    surfsd[[k]]<-msurD
    surfsd1[[k]]<-msurD1

  }

  thrs<-sapply(res2,"[[",2)


  meshes1<-list()
  meshes2<-list()
  for(k in 1:ncol(comb.df)){


    pale=colorRampPalette(c(col,"white","white","white","white"))
    v.pale=pale(1000)
    v.pale=v.pale[seq(1,1000,by=100)-1]
    v.pale=c(rev(v.pale[-1]),pale(1),v.pale[-1])


    ddpcr::quiet(meshes1[[k]]<-localmeshdiff(surfs[[k]],surfs1[[k]],ploton = 1,n.int = 100,vec=surfsd[[k]],
                                             paltot=v.pale,densityplot = FALSE,
                                             from = -max(thrs),to=max(thrs),out.rem = out.rem,fact = 1.5,visual = 1,scale01 = FALSE, plot = FALSE))
    ddpcr::quiet(meshes2[[k]]<-localmeshdiff(surfs1[[k]],surfs[[k]],ploton = 1,n.int = 100,vec=surfsd1[[k]],
                                             paltot=v.pale,densityplot = FALSE,
                                             from = -max(thrs),to=max(thrs),out.rem = out.rem,fact = 1.5,visual = 1,scale01 = FALSE, plot = FALSE))

  }

  if(isTRUE(plot) & nrow(df)<=5){

    matplot<-matrix(0,nrow(df),nrow(df))
    matplot[lower.tri(matplot)]<-seq(1:ncol(comb.df))
    matplot<-t(matplot)
    matplot[lower.tri(matplot)]<-seq((ncol(comb.df)+1),ncol(comb.df)*2)
    diagmat<-diag(matplot)
    diagmat[1]<-1
    matplot<-matplot+1
    diag(matplot)<-diagmat


    rgl::open3d()
    rgl::layout3d(matplot,sharedMouse = TRUE)
    rgl::spheres3d(mshape,radius = 0.000000001)
    if(names == TRUE) rgl::title3d(paste(rownames(df)[1]),font=2)
    if(names == TRUE) rgl::mtext3d(text=paste(rownames(df)[1]),font=2,edge="y",line=3)
    rgl::next3d()
    for (i in 1:length(meshes1)){
      rgl::shade3d(meshes1[[i]]$mesh,specular="black")
      if (i<ncol(matplot) & names == TRUE) rgl::title3d(paste(rownames(df)[i+1]),font=2)
      rgl::next3d()
    }
    for (i in 1:length(meshes2)){
      rgl::shade3d(meshes2[[i]]$mesh,specular="black")
      if (i<ncol(matplot) & names == TRUE) rgl::mtext3d(text=paste(rownames(df)[i+1]),font=2,edge="y",line=3)
      if(i<length(meshes2)) rgl::next3d()
    }

    colmap_tot<-colorRampPalette(v.pale)
    plot(seq(-max(thrs),max(thrs),length.out = 200),rep(0,200),main="area differences legend",xlab="",ylab="",col="white")
    abline(v=seq(-max(thrs),max(thrs),length.out = 200),col=colmap_tot(200),lwd=5)

  }

  angle.compare<-do.call(rbind,lapply(res,"[[",1))
  rownames(angle.compare)<-apply(comb.df,2,function(x) paste(x[1],x[2],sep = "-"))
  as.data.frame(angle.compare[order(angle.compare[,2]),])->angle.compare


  selected.pcs<-selected

  average.dist<-sapply(res2,"[[",1)
  mat<-matrix(0,nrow(df),nrow(df))
  mat[lower.tri(mat)]<-average.dist
  mat<-t(mat)
  mat[lower.tri(mat)]<-average.dist
  colnames(mat)<-rownames(mat)<-rownames(df)

  meshes1<-lapply(meshes1,"[[",1)
  meshes2<-lapply(meshes2,"[[",1)
  names(meshes1)<-apply(comb.df,2,function(x) paste(x[1],x[2],sep = "-"))
  names(meshes2)<-apply(comb.df,2,function(x) paste(x[2],x[1],sep = "-"))

  return(list(angle.compare=angle.compare,
              selected.pcs=selected.pcs,
              average.dist=mat,
              surfaces1=meshes1,
              surfaces2=meshes2,
              scale= max(thrs)))

}

