## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
require(rgl)

options(knitr.kable.NA = '',mc.cores=2,rgl.useNULL=TRUE,
        rmarkdown.html_vignette.check_title = FALSE)


## ----include=FALSE------------------------------------------------------------
require(RRphylo)
require(RColorBrewer)
require(rgl)
require(ape)
require(phytools)
require(geiger)
require(mvMORPH)

theta<-function(a,b){
  unitV <- function(x) {
    sum(x^2)^0.5
  }
  deg2rad <- function(deg) {
    (deg * pi)/(180)
  }
  rad2deg <- function(rad) {
    (rad * 180)/(pi)
  }
  rad2deg(acos((a%*%b)/(unitV(a) *unitV(b))))->theta
  return(theta)
}

veclen <- function(v) sqrt(sum(v^2))
xprod <- function(v, w) c( v[2]*w[3] - v[3]*w[2],
                           v[3]*w[1] - v[1]*w[3],
                           v[1]*w[2] - v[2]*w[1] )

arc3d <- function(from, to, center, radius, n, circle = 50, base = 0, plot = TRUE, ...) {
  fixarg <- function(arg) {
    if (is.matrix(arg))
      arg[, 1:3, drop = FALSE]
    else
      matrix(arg, 1, 3)
  }
  normalize <- function(v)
    v / veclen(v)
  getrow <- function(arg, i) {
    arg[1 + (i - 1) %% nrow(arg),]
  }
  from <- fixarg(from)
  to <- fixarg(to)
  center <- fixarg(center)

  m <- max(nrow(from), nrow(to), nrow(center), length(base))
  base <- rep_len(base, m)

  result <- matrix(NA_real_, nrow = 1, ncol = 3)

  for (j in seq_len(m)) {
    from1 <- getrow(from, j)
    to1 <- getrow(to, j)
    center1 <- getrow(center, j)
    base1 <- base[j]
    logr1 <- log(veclen(from1 - center1))
    logr2 <- log(veclen(to1 - center1))
    A <- normalize(from1 - center1)
    B <- normalize(to1 - center1)
    steps <- if (base1 <= 0) 4*abs(base1) + 1 else 4*base1 - 1
    for (k in seq_len(steps)) {
      if (k %% 2) {
        A1 <- A * (-1)^(k %/% 2)
        B1 <- B * (-1)^(k %/% 2 + (base1 > 0))
      } else {
        A1 <- B * (-1)^(k %/% 2 + (base1 <= 0))
        B1 <- A * (-1)^(k %/% 2)
      }
      theta <- acos(sum(A1*B1))
      if (isTRUE(all.equal(theta, pi)))
        warning("Arc ", j, " points are opposite each other!  Arc is not well defined.")
      if (missing(n))
        n1 <- ceiling(circle*theta/(2*pi))
      else
        n1 <- n

      if (missing(radius)) {
        pretheta <- (k %/% 2)*pi - (k %% 2 == 0)*theta
        if (k == 1)
          totaltheta <- (steps %/% 2)*pi - (steps %% 2 == 0)*theta + theta
        p1 <- pretheta/totaltheta
        p2 <- (pretheta + theta)/totaltheta
        radius1 <- exp(seq(from = (1 - p1)*logr1 + p1*logr2,
                           to   = (1 - p2)*logr1 + p2*logr2,
                           length.out = n1 + 1))
      } else
        radius1 <- rep_len(radius, n1)
      arc <- matrix(NA_real_, nrow = n1 + 1, ncol = 3)
      p <- seq(0, 1, length.out = n1 + 1)
      arc[1,] <- center1 + radius1[1]*A1
      arc[n1 + 1,] <- center1 + radius1[n1 + 1]*B1
      AB <- veclen(A1 - B1)
      for (i in seq_len(n1)[-1]) {
        ptheta <- p[i]*theta
        phi <- pi/2 + (0.5 - p[i])*theta
        q <- (sin(ptheta) / sin(phi))/AB
        D <- (1-q)*A1 + q*B1
        arc[i,] <- center1 + radius1[i] * normalize(D)
      }
      if (k == 1)
        result <- rbind(result, arc)
      else
        result <- rbind(result[-nrow(result), ,drop = FALSE], arc)
    }
    result <- rbind(result, result[1,])
  }
  if (plot)
    lines3d(result[c(-1, -nrow(result)), , drop = FALSE], ...)
  else
    result[c(-1, -nrow(result)), , drop = FALSE]
}

### creating a phylogenetic tree with 100 species
set.seed(14)
rtree(50)->tree

### select, extract and then modify the clade to be duplicated
dist.nodes(tree)[(Ntip(tree)+1),]->cd
cd[cd<0.8*max(cd)]->cd
cd[which(as.numeric(names(cd))>Ntip(tree))]->cd
cdt<-array()
for(e in 1:length(cd)) length(tips(tree,names(cd)[e]))->cdt[e]
names(cd[-which(cdt<Ntip(tree)*.1|cdt>Ntip(tree)/4)])->cd
cdd<-array()
for(i in 1:length(cd)){
  if((length(tips(tree,cd[i]))-2)>(length(tips(tree,cd[i]))+Ntip(tree)-2)*0.1) cdd[i]<-"ok" else cdd[i]<-"no"
}
if(length(which(cdd=="ok"))>0) cd[which(cdd=="ok")]->cd else cd->cd
as.numeric(sample(cd,1))->n

extract.clade(tree,n)->t1
max(nodeHeights(t1))->nH1
suppressWarnings(swapONE(t1)[[1]]->t1)
drop.tip(t1,t1$tip.label[c(1,length(t1$tip.label))])->t1

if(max(nodeHeights(t1))!=nH1) geiger:::rescale(t1,nH1,model = "depth")->t1

t1$root.edge<-data.frame(tree$edge,tree$edge.length)[which(
  data.frame(tree$edge,tree$edge.length)[,2]==n),3]


### selecting the node where the new clade is to be binded
distNodes(tree,n)[1:Nnode(tree),]->dfN
dfN[-which(rownames(dfN)%in%c(n,getDescendants(tree,n))),]->dfN
distNodes(tree,(Ntip(tree)+1))[1:Nnode(tree),][-1,]->dR
dR[match(rownames(dfN),rownames(dR)),]->dR1
dfN[match(rownames(dR1)[which(dR1[,2]<max(dR[,2])*0.8)],rownames(dfN)),]->dfN
rownames(dfN)[which(dfN[,1]<=Ntip(tree)/10)]->bar
dfN[-which(dfN[,1]<=Ntip(tree)/10),]->dfn2


if(class(dfn2)!="matrix"){
  t(as.matrix(dfn2))->dfn2
  rownames(dfN)[-which(dfN[,1]<=Ntip(tree)/10)]->rownames(dfn2)
}

if(dim(dfn2)[1]==0){
  rownames(dfN)[which.max(dfN[,1])]->nodN
  dfN[which.max(dfN[,1]),1]->minD
} else{
  dR[match(rownames(dfn2),rownames(dR)),]->dfn3
  if(class(dfn3)!="matrix"){
    t(as.matrix(dfn3))->dfn3
    rownames(dR)[match(rownames(dfn2),rownames(dR))]->rownames(dfn3)
  }
  rownames(dfn3)[which.min(dfn3[,2])]->nodN
  dfn2[rownames(dfn2)==nodN,1]->minD
}
nodN->tar

at<-tar
data.frame(tree$edge,tree$edge.length)[which(data.frame(tree$edge,
                                                        tree$edge.length)[,2]==at),3]->pos
diff(dist.nodes(tree)[(Ntip(tree)+1),c(n,tar)])->dH
if(dH>0) {
  geiger:::rescale(t1,(abs(diff(c(max(nodeHeights(t1)),dH)))),model = "depth")->t1
  t1$root.edge+pos/2->t1$root.edge
}else{
  geiger:::rescale(t1,(max(nodeHeights(t1))+abs(dH)),model = "depth")->t1
  t1$root.edge+pos/2->t1$root.edge
}

tree->treeO
t1->t1O


## ----noconv, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(93)
tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
bind.tree(tree,t1,where=at,position=pos/2)->tree1
ntraits<-3
thetaX=c(0,0,0)
matrix(c(1,0.5,0.5,0.5,1,0.5,.5,.5,1),3,3)->sigma
mvSIM(tree1,param=list(ntraits=ntraits,theta=thetaX,sigma=sigma))->y
c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec
match(rownames(y.vec),tree1$tip.label)->aa

apply(y[tips(tree1,nod.par[1]),],2,mean)->y.vec[1,]
apply(y[tips(tree1,nod.par[2]),],2,mean)->y.vec[2,]
matrix(c(0,0,0),1,3)->rV

rep("gray80",length(y))->coltips
coltips[match(tips(tree1,nod.par[1]),rownames(y))]<-"firebrick1"
coltips[match(tips(tree1,nod.par[2]),rownames(y))]<-"deepskyblue1"
plot3d(y,col=coltips,size=1,bbox=F,axes=F,box=F,type="s")
box3d()
title3d(xlab="y[,1]",ylab="y[,2]",zlab="y[,3]")
plot3d(matrix(rV,1,3),col="green",add=T,size=2,bbox=F,axes=F,box=F,type="s")
plot3d(ant,col="gray16",size=1,bbox=F,axes=T,box=F,add=T,type="s")
text3d(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),texts="mrca1",add=T,adj=1,cex=1,col="red4")
text3d(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),texts="mrca2",add=T,adj=1,cex=1,col="blue3")


segments3d(rbind(rV,matrix(ant[match(nod.par[1],rownames(ant)),],1,3)),col="red4",lwd=4)
segments3d(rbind(rV,matrix(ant[match(nod.par[2],rownames(ant)),],1,3)),col="blue3",lwd=4)
segments3d(rbind(rV,y.vec[2,]),col="deepskyblue1",lwd=4)
segments3d(rbind(rV,y.vec[1,]),col="firebrick1",lwd=4)

matrix(ant[match(nod.par[1],rownames(ant)),],1,3)->a
matrix(ant[match(nod.par[2],rownames(ant)),],1,3)->b
theta(c(a),c(b))->theta.ace
as.numeric(round(theta.ace,1))->theta.ace
veclen(a)->la
veclen(b)->lb
if(la>lb) (la/lb)*b->b else (lb/la)*a->a
if(la<.5*lb) segments3d(rbind(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),(lb/la)*matrix(ant[match(nod.par[1],rownames(ant)),],1,3)))
if(la>.5*lb) segments3d(rbind(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),(la/lb)*matrix(ant[match(nod.par[2],rownames(ant)),],1,3)))


veclen(a)->la
arc3d(a,b,matrix(rV,1,3),add=T,col="red",radius=.5*la)
arc3d(b,a,matrix(rV,1,3),add=T,col="red",radius=.5*la)

text3d(rV,texts=paste("theta.ace",theta.ace),adj=c(-.5,.5),col="red")

matrix(y.vec[1,],1,3)->at1
matrix(y.vec[2,],1,3)->at2
theta(c(at1),c(at2))->theta.real
as.numeric(round(theta.real,1))->theta.real
veclen(at1)->lat1
veclen(at2)->lat2
if(lat1>lat2) (lat1/lat2)*at2->at2 else (lat2/lat1)*at1->at1
if(lat1<.5*lat2) segments3d(rbind(matrix(y.vec[1,],1,3),(lat2/lat1)*matrix(y.vec[1,],1,3))) 
if(lat2<.5*lat1) segments3d(rbind(matrix(y.vec[2,],1,3),(lat1/lat2)*matrix(y.vec[2,],1,3)))
veclen(at1)->lat1
arc3d(at1,at2,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
arc3d(at2,at1,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
text3d(rV,texts=paste("theta.tips",theta.real),adj=c(-.5,lat1),col="blue")

range(y[,3])->ran
if(ran[1]<0) ran[1]*0.9->ran[1] else ran[1]*1.1->ran[1]
if(ran[2]<0) ran[2]*1.1->ran[2] else ran[2]*0.9->ran[2]


points3d(matrix(c(rep(range(y[,1])[2]*1.5,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
         col=c("gray80","gray16","green","firebrick1","deepskyblue1","blue3","red4"),size=8)
text3d(matrix(c(rep(range(y[,1])[2]*1.6,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
       texts=c("phenotypes at tips","phenotypes at nodes","origin","clade 2","clade 1"),adj=c(0,0.5),cex=1.2)

rglwidget(elementId = "plot3drgl")

## ----simpleconv, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(14)
treeO->tree
t1O->t1
ntraits<-3
thetaX=c(0,0,0)
matrix(c(1,0.5,0.5,0.5,1,0.5,.5,.5,1),3,3)->sigma
mvSIM(tree,param=list(ntraits=ntraits,theta=thetaX,sigma=sigma))->y
y[match(tree$tip.label,rownames(y)),]->y
y[match(tips(tree,n),rownames(y)),]->a


apply(y,2,range)->m.a
m.a[2,]*1.2->m.a

a1<-matrix(ncol=dim(y)[2],nrow=dim(a)[1])
for(m in 1:dim(a)[1])
{
  v<-array()
  for(i in 1:length(m.a)) jitter(m.a[i],amount=(sd(a[,i])*1))->v[i]
  v->a1[m,]
}

rownames(a1)<-rownames(a)
y[match(tips(tree,n),rownames(y)),]<-a1

apply(y[match(t1$tip.label,rownames(y)),],2,jitter)->y.t1

tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
rownames(y)<-tree$tip.label
rownames(y.t1)<-t1$tip.label
bind.tree(tree,t1,where=at,position=pos/2)->tree1
rbind(y,y.t1)->y

### nod.par is the node pair set to converge
c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec
match(rownames(y.vec),tree1$tip.label)->aa


apply(y[tips(tree1,nod.par[1]),],2,mean)->y.vec[1,]
apply(y[tips(tree1,nod.par[2]),],2,mean)->y.vec[2,]
matrix(c(0,0,0),1,3)->rV

rep("gray80",length(y))->coltips
coltips[match(tips(tree1,nod.par[1]),rownames(y))]<-"firebrick1"
coltips[match(tips(tree1,nod.par[2]),rownames(y))]<-"deepskyblue1"
plot3d(y,col=coltips,size=1,bbox=F,axes=F,box=F,type="s")
box3d()
title3d(xlab="y[,1]",ylab="y[,2]",zlab="y[,3]")
plot3d(matrix(rV,1,3),col="green",add=T,size=2,bbox=F,axes=F,box=F,type="s")
plot3d(ant,col="gray16",size=1,bbox=F,axes=T,box=F,add=T,type="s")
text3d(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),texts="mrca1",add=T,adj=1,cex=1,col="red4")
text3d(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),texts="mrca2",add=T,adj=1,cex=1,col="blue3")


segments3d(rbind(rV,matrix(ant[match(nod.par[1],rownames(ant)),],1,3)),col="red4",lwd=4)
segments3d(rbind(rV,matrix(ant[match(nod.par[2],rownames(ant)),],1,3)),col="blue3",lwd=4)
segments3d(rbind(rV,y.vec[2,]),col="deepskyblue1",lwd=4)
segments3d(rbind(rV,y.vec[1,]),col="firebrick1",lwd=4)

matrix(ant[match(nod.par[1],rownames(ant)),],1,3)->a
matrix(ant[match(nod.par[2],rownames(ant)),],1,3)->b
theta(c(a),c(b))->theta.ace
as.numeric(round(theta.ace,1))->theta.ace
veclen(a)->la
veclen(b)->lb
if(la>lb) (la/lb)*b->b else (lb/la)*a->a
if(la<.5*lb) segments3d(rbind(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),(lb/la)*matrix(ant[match(nod.par[1],rownames(ant)),],1,3)))
if(la>.5*lb) segments3d(rbind(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),(la/lb)*matrix(ant[match(nod.par[2],rownames(ant)),],1,3)))


veclen(a)->la
arc3d(a,b,matrix(rV,1,3),add=T,col="red",radius=.5*la)
arc3d(b,a,matrix(rV,1,3),add=T,col="red",radius=.5*la)

text3d(rV,texts=paste("theta.ace",theta.ace),adj=c(-.5,.5),col="red")

matrix(y.vec[1,],1,3)->at1
matrix(y.vec[2,],1,3)->at2
theta(c(at1),c(at2))->theta.real
as.numeric(round(theta.real,1))->theta.real
veclen(at1)->lat1
veclen(at2)->lat2
if(lat1>lat2) (lat1/lat2)*at2->at2 else (lat2/lat1)*at1->at1
if(lat1<.5*lat2) segments3d(rbind(matrix(y.vec[1,],1,3),(lat2/lat1)*matrix(y.vec[1,],1,3))) 
if(lat2<.5*lat1) segments3d(rbind(matrix(y.vec[2,],1,3),(lat1/lat2)*matrix(y.vec[2,],1,3)))
veclen(at1)->lat1
arc3d(at1,at2,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
arc3d(at2,at1,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
text3d(rV,texts=paste("theta.tips",theta.real),adj=c(-.5,lat1),col="blue")

range(y[,3])->ran
if(ran[1]<0) ran[1]*0.9->ran[1] else ran[1]*1.1->ran[1]
if(ran[2]<0) ran[2]*1.1->ran[2] else ran[2]*0.9->ran[2]


points3d(matrix(c(rep(range(y[,1])[2]*1.5,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
         col=c("gray80","gray16","green","firebrick1","deepskyblue1","blue3","red4"),size=8)
text3d(matrix(c(rep(range(y[,1])[2]*1.6,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
       texts=c("phenotypes at tips","phenotypes at nodes","origin","clade 2","clade 1"),adj=c(0,0.5),cex=1.2)
rglwidget(elementId = "plot3drgl1")

## ----conv&par, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(14)
treeO->tree
t1O->t1

ntraits<-3
thetaX=c(0,0,0)
matrix(c(1,0.5,0.5,0.5,1,0.5,.5,.5,1),3,3)->sigma
mvSIM(tree,param=list(ntraits=ntraits,theta=thetaX,sigma=sigma))->y
y[match(tree$tip.label,rownames(y)),]->y
y[match(tips(tree,n),rownames(y)),]->a


apply(y,2,range)->m.a
sample(seq(.8,1.2,.1),1)->ff
m.a[2,]*ff->m.a

a1<-matrix(ncol=dim(y)[2],nrow=dim(a)[1])
for(m in 1:dim(a)[1])
{
  v<-array()
  for(i in 1:length(m.a)) jitter(m.a[i],amount=(sd(a[,i])*1))->v[i]
  v->a1[m,]
}

rownames(a1)<-rownames(a)
y[match(tips(tree,n),rownames(y)),]<-a1

apply(y[match(t1$tip.label,rownames(y)),],2,jitter)->y.t1

tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
rownames(y)<-tree$tip.label
rownames(y.t1)<-t1$tip.label
bind.tree(tree,t1,where=at,position=pos/2)->tree1
rbind(y,y.t1)->y

c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant

### sample one tip descending from each node set to converge
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec
match(rownames(y.vec),tree1$tip.label)->aa

### modify the phenotypes at the mrcas (nod.par) as the median of the distribution of phenotypes of descending tips   
apply(y[tips(tree1,nod.par[1]),],2,median)->mrca.1
ant[match(nod.par[1],rownames(ant)),]<-mrca.1
apply(y[tips(tree1,nod.par[2]),],2,median)->mrca.2
ant[match(nod.par[2],rownames(ant)),]<-mrca.2


apply(y[tips(tree1,nod.par[1]),],2,mean)->y.vec[1,]
apply(y[tips(tree1,nod.par[2]),],2,mean)->y.vec[2,]
matrix(c(0,0,0),1,3)->rV

rep("gray80",length(y))->coltips
coltips[match(tips(tree1,nod.par[1]),rownames(y))]<-"firebrick1"
coltips[match(tips(tree1,nod.par[2]),rownames(y))]<-"deepskyblue1"
plot3d(y,col=coltips,size=1
       ,bbox=F,axes=F,box=F,type="s")
box3d()
title3d(xlab="y[,1]",ylab="y[,2]",zlab="y[,3]")
plot3d(matrix(rV,1,3),col="green",add=T,size=2,bbox=F,axes=F,box=F,type="s")
plot3d(ant,col="gray16",size=1,bbox=F,axes=T,box=F,add=T,type="s")
text3d(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),texts="mrca1",add=T,adj=1,cex=1,col="red4")
text3d(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),texts="mrca2",add=T,adj=1,cex=1,col="blue3")


segments3d(rbind(rV,matrix(ant[match(nod.par[1],rownames(ant)),],1,3)),col="red4",lwd=4)
segments3d(rbind(rV,matrix(ant[match(nod.par[2],rownames(ant)),],1,3)),col="blue3",lwd=4)
segments3d(rbind(rV,y.vec[2,]),col="deepskyblue1",lwd=4)
segments3d(rbind(rV,y.vec[1,]),col="firebrick1",lwd=4)

matrix(ant[match(nod.par[1],rownames(ant)),],1,3)->a
matrix(ant[match(nod.par[2],rownames(ant)),],1,3)->b
theta(c(a),c(b))->theta.ace
as.numeric(round(theta.ace,1))->theta.ace
veclen(a)->la
veclen(b)->lb
if(la>lb) (la/lb)*b->b else (lb/la)*a->a
if(la<.5*lb) segments3d(rbind(matrix(ant[match(nod.par[1],rownames(ant)),],1,3),(lb/la)*matrix(ant[match(nod.par[1],rownames(ant)),],1,3)))
if(la>.5*lb) segments3d(rbind(matrix(ant[match(nod.par[2],rownames(ant)),],1,3),(la/lb)*matrix(ant[match(nod.par[2],rownames(ant)),],1,3)))


veclen(a)->la
arc3d(a,b,matrix(rV,1,3),add=T,col="red",radius=.5*la)
arc3d(b,a,matrix(rV,1,3),add=T,col="red",radius=.5*la)

text3d(rV,texts=paste("theta.ace",theta.ace),adj=c(-.5,.5),col="red")

matrix(y.vec[1,],1,3)->at1
matrix(y.vec[2,],1,3)->at2
theta(c(at1),c(at2))->theta.real
as.numeric(round(theta.real,1))->theta.real
veclen(at1)->lat1
veclen(at2)->lat2
if(lat1>lat2) (lat1/lat2)*at2->at2 else (lat2/lat1)*at1->at1
if(lat1<.5*lat2) segments3d(rbind(matrix(y.vec[1,],1,3),(lat2/lat1)*matrix(y.vec[1,],1,3))) 
if(lat2<.5*lat1) segments3d(rbind(matrix(y.vec[2,],1,3),(lat1/lat2)*matrix(y.vec[2,],1,3)))
veclen(at1)->lat1
arc3d(at1,at2,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
arc3d(at2,at1,matrix(rV,1,3),add=T,col="blue",radius=.5*lat1)
text3d(rV,texts=paste("theta.tips",theta.real),adj=c(-.5,lat1),col="blue")

range(y[,3])->ran
if(ran[1]<0) ran[1]*0.9->ran[1] else ran[1]*1.1->ran[1]
if(ran[2]<0) ran[2]*1.1->ran[2] else ran[2]*0.9->ran[2]


points3d(matrix(c(rep(range(y[,1])[2]*1.5,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
         col=c("gray80","gray16","green","firebrick1","deepskyblue1","blue3","red4"),size=8)
text3d(matrix(c(rep(range(y[,1])[2]*1.6,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=F),
       texts=c("phenotypes at tips","phenotypes at nodes","origin","clade 2","clade 1"),adj=c(0,0.5),cex=1.2)
rglwidget(elementId = "plot3drgl2")

## ----echo=FALSE,fig.width=5,fig.height=6,fig.align="center",message=FALSE,warning=FALSE,out.width='60%',dpi=200----
rad2deg <- function(rad) (rad * 180)/(pi)
unitV <- function(x) sum(x^2)^0.5

set.seed(22)
rtree(14)->tree
getDescendants(tree,18)->des1
getDescendants(tree,24)->des2
c(18,getMommy(tree,18))->mom1
c(24,getMommy(tree,24))->mom2
mom1[-length(mom1)]->mom1
mom2[-length(mom2)]->mom2
colo<-rep("gray50",length(tree$edge.length))
colo[which(tree$edge[,2]%in%des1)]<-"deepskyblue1"
colo[which(tree$edge[,2]%in%des2)]<-"firebrick1"
colo[which(tree$edge[,2]%in%c(mom1,mom2))]<-"black"

wid<-rep(3,length(tree$edge.length))
wid[which(tree$edge[,2]%in%mom1)]<-4
wid[which(tree$edge[,2]%in%mom2)]<-4

set.seed(14)
fastBM(tree,nsim=3)->y

apply(y[tips(tree,18)[1:3],],2,function(x) jitter(x,2))->y[tips(tree,24),]

RRphylo(tree,y,clus=2/parallel::detectCores())->RR
rbind(RR$aces,y)->phen
tree->tree1

c(18,24)->nod
dist.nodes(tree1)[nod[1],nod[2]]->nT
tips(tree1,nod[1])->tt1
tips(tree1,nod[2])->TT
expand.grid(tt1,TT)->ctt
aa<-array()
for(g in 1:dim(ctt)[1]){
  phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),]->ppTT
  as.matrix(ppTT)->ppTT
  aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
}
mean(aa)->ang.tip

RR$aces[which(rownames(RR$aces)%in%nod),]->ac
rad2deg(acos((ac[1,] %*% ac[2,])/(unitV(ac[1,]) *unitV(ac[2,]))))->ang.ac

apply(ctt,2,as.character)->ctt
rbind(c("$\\theta_{real}$", "=",round(ang.tip,3)),c(paste("$\\theta_{real}$","+","$\\theta_{ace}$",sep=""),"=",round(ang.tip+ang.ac,3)),c("$distance_{mrcas}$","=",round(nT,3)))->abc
as.data.frame(rbind(cbind(ctt,round(aa,3)),c("mrcaC-18","mrca-24",round(ang.ac,3)),abc))->angs
colnames(angs)<-c("clade18","clade24","angle")

col2hex <- function(col, alpha) rgb(t(col2rgb(col)), alpha=alpha, maxColorValue=255)
require(kableExtra)
knitr::kable(angs,digits=3,align="c") %>%
  kable_styling(full_width = F, position = "float_right")  %>%
  column_spec(1, color = col2hex("deepskyblue1"),bold=TRUE) %>%
  column_spec(2, color = col2hex("firebrick1"),bold=TRUE) %>%
  pack_rows("", 14, 16) %>%
  row_spec(14:16, color = "black",bold=TRUE)

plot(tree,no.margin = TRUE,direction="downward",srt=90,adj=0.5,edge.color = colo,edge.width = wid)
nodelabels(node=c(18,24),text=c("mrcaC1","mrcaC2"),bg="w",frame="n",adj=c(-0.1,-0.2),font=2)
nodelabels(bg="w",frame="n",adj=c(0.5,1.2),font=2)
legend("bottomright",legend=c("Clade 18","Clade 24","mrca-18 to mrca-24 distance"),
       lty=1,col=c("deepskyblue1","firebrick1","black"),lwd=c(3,3,4))


## ----message=FALSE, warning=FALSE,echo=2--------------------------------------
pdf.options(useDingbats=FALSE)
search.conv(RR=RR,y=y,min.dim=3,max.dim=4,nsim=100,rsim=100,
            clus=2/parallel::detectCores(),foldername=getwd())->SC
paste(rownames(SC$`node pairs`),SC$`node pairs`[,1],sep="-")->SC$`node pairs`[,1]
colnames(SC$`node pairs`)[c(1,6:11)]<-c("node.pair","node","time","ang.bydist","ang.conv","n1","n2")
rownames(SC$`node pairs`)<-NULL

knitr::kable(SC$`node pairs`,digits=3,align="c") %>%
  column_spec(1, bold=TRUE) %>%
  add_header_above( c(" "=5,"distance"=2,"p-value"=2,"Clade size" =2))


## ----message=FALSE, warning=FALSE,echo=FALSE----------------------------------
knitr::kable(as.data.frame(cbind(SC$`node pairs comparison`,
                                 t(as.matrix(SC$`average distance from group centroids`)))),digits=3,align="c")%>%
  kable_styling(full_width = F, position = "center")%>%
  add_header_above( c("node pairs comparison"=5,"average distance from group centroids"=2))


## ----echo=FALSE,fig.align="center",message=FALSE,warning=FALSE, out.width='98%'----
require(pdftools)
ddpcr::quiet(pdf_convert(paste(getwd(),"/convergence plot.pdf",sep=""), format = "png",dpi=300))
ddpcr::quiet(file.remove(paste(getwd(),"/convergence plot.pdf",sep="")))
knitr::include_graphics("convergence plot_1.png")

## ----echo=FALSE,fig.width=5,fig.height=6,fig.align="center",message=FALSE,warning=FALSE,out.width='60%',dpi=200----

set.seed(22)
fastBM(tree,nsim=3)->y
apply(y,2,range)[2,]->m.a
sample(seq(1:Ntip(tree)),3)->ff

p1<-matrix(ncol=length(m.a),nrow=length(ff))
for(e in 1:length(ff)){
  p0<-array()
  for(i in 1:length(m.a)) jitter(m.a[i],amount=sd(y[,i]))->p0[i]
  p0->p1[e,]
}

### simulating the phenotypic state
y[ff,]<-p1
rep("nostate",Ntip(tree))->state
names(state)<-rownames(y)
state[ff]<-"a"

sample(seq(1:Ntip(tree))[-ff],4)->f2
p2<-matrix(ncol=length(m.a),nrow=length(f2))
for(e in 1:length(f2)){
  p0<-array()
  for(i in 1:length(m.a)) jitter(m.a[i],amount=2*sd(y[,i]))->p0[i]
  p0->p2[e,]
}
y[f2,]<-p2
state[f2]<-"b"

tree->tree1
state[match(rownames(y),names(state))]->state
if("nostate"%in%state) state[-which(state=="nostate")]->state.real else state->state.real

combn(unique(state.real),2)->stcomb1
if("nostate"%in%state) combn(unique(state),2)->stcomb else stcomb1->stcomb

cophenetic.phylo(tree1)->cop
i=1
y[which(state==stcomb1[1,i]),]->tt1
mean(apply(tt1,1,unitV))->vs1
y[which(state==stcomb1[2,i]),]->TT
mean(apply(TT,1,unitV))->vs2
expand.grid(rownames(tt1),rownames(TT))->ctt
aa<-array()
dt<-array()
for(g in 1:dim(ctt)[1]){
  y[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(y)),]->ppTT
  as.matrix(ppTT)->ppTT
  aa[g] <- rad2deg(acos((ppTT[1,]%*%ppTT[2,])/(unitV(ppTT[1,]) *unitV(ppTT[2,]))))
  cop[match(as.character(ctt[g,1]),rownames(cop)),match(as.character(ctt[g,2]),rownames(cop))]->dt[g]
}


rbind(c(paste("mean","$\\theta_{real}$",sep=" "), "=",round(mean(aa),3)),
      c(paste("mean", "$\\frac{\\theta_{real}}{distance}$",sep=" "),"=",round(mean(aa/dt),3)))->abc
as.data.frame(rbind(cbind(as.matrix(ctt)[,2:1],round(aa,3),round(dt,3)),c("","","",""),c("","","",""),c("","","","")))->angs
colnames(angs)<-c("state a","state b","angle","distance")

col2hex <- function(col, alpha) rgb(t(col2rgb(col)), alpha=alpha, maxColorValue=255)
require(kableExtra)
knitr::kable(angs,digits=3,align="c") %>%
  kable_styling(full_width = F, position = "float_right")  %>%
  column_spec(1, color = col2hex("deepskyblue1"),bold=TRUE) %>%
  column_spec(2, color = col2hex("firebrick1"),bold=TRUE)%>%  
  pack_rows(paste(paste("mean","$\\theta_{real}$",sep=" "), round(mean(aa),3),sep=" = "), 14,14,
            label_row_css = "",latex_gap_space = "0em") %>% 
  pack_rows(paste(paste("mean", "$\\frac{\\theta_{real}}{distance}$",sep=" "),round(mean(aa/dt),3),sep=" = "), 15,15,
            label_row_css = "",latex_gap_space = "0em") 

colo<-rep("gray50",length(tree$tip.label))
colo[which(tree$tip.label%in%names(which(state=="a")))]<-"deepskyblue1"
colo[which(tree$tip.label%in%names(which(state=="b")))]<-"firebrick1"

plot(tree,no.margin = TRUE,direction="downward",srt=90,adj=0.5,label.offset = 0.08)
tiplabels(bg=colo,text=rep("",length(tree$tip.label)),frame="circle",adj=0.5,offset=0.05)
legend("bottomright",legend=c("State a","State b","nostate"),pch=21,pt.cex=2,
       pt.bg=c("deepskyblue1","firebrick1","gray50"))


## ----message=FALSE, warning=FALSE,echo=2--------------------------------------
pdf.options(useDingbats=FALSE)
search.conv(tree=tree,y=y,state=state,nsim=100,
            clus=2/parallel::detectCores(),foldername=getwd())->SC

knitr::kable(SC,digits=3,align="c") %>%
  column_spec(1:2, bold=TRUE)


## ----echo=FALSE,fig.align="center",message=FALSE,warning=FALSE, out.width='98%'----
require(pdftools)
ddpcr::quiet(pdf_convert(paste(getwd(),"/convergence plot for different states.pdf",sep=""), format = "png",dpi=300))
ddpcr::quiet(file.remove(paste(getwd(),"/convergence plot for different states.pdf",sep="")))
knitr::include_graphics("convergence plot for different states_1.png")

## ----fig.width=5,fig.height=5,fig.align="center",message=FALSE,warning=FALSE,out.width='90%',dpi=200----
# load the RRphylo example dataset including Felids tree and data
data("DataFelids")
DataFelids$PCscoresfel->PCscoresfel # mandible shape data
DataFelids$treefel->treefel # phylogenetic tree
DataFelids$statefel->statefel # conical-toothed ("nostate") or saber-toothed condition

library(ape)
plot(ladderize(treefel),show.tip.label = F,no.margin = T)
colo<-rep("gray50",length(treefel$tip.label))
colo[match(names(which(statefel=="saber")),treefel$tip.label)]<-"firebrick1"
tiplabels(text=rep("",Ntip(treefel)),bg=colo,frame="circle",cex=.4)
legend("bottomleft",legend=c("Sabertooths","nostate"),pch=21,pt.cex=1.5,
       pt.bg=c("firebrick1","gray50"))

## ----eval=FALSE---------------------------------------------------------------
#  # perform RRphylo on Felids tree and data
#  RRphylo(tree=treefel,y=PCscoresfel)->RRfel
#  
#  ## Example 1: search for morphological convergence between clades (automatic mode)
#  ## by setting 9 nodes as minimum distance between the clades to be tested
#  search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",
#              foldername = tempdir())->SC.clade
#  
#  ## Example 2: search for morphological convergence within sabertoothed species
#  search.conv(tree=treefel, y=PCscoresfel, state=statefel,
#              foldername = tempdir())->SC.state
#  

