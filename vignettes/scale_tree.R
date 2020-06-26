## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----echo=c(14:15,32:33), fig.dim=c(6,6), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
library(ape)
library(phytools)
library(geiger)

set.seed(14)
DataFelids$treefel->tree
tree$tip.label<-paste("t",seq(1:Ntip(tree)),sep="")

max(nodeHeights(tree))->H
sample(H-diag(vcv(tree)),8)->sp.ages
sp.ages+tree$edge.length[match(match(names(sp.ages),tree$tip.label),tree$edge[,2])]/2->sp.ages
sp.ages[c(3,7)]<-0

sp.ages
scaleTree(tree,tip.ages=sp.ages)->treeS1

edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%match(names(sp.ages),tree$tip.label))]<-"blue"

par(mfcol=c(2,2),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("species ages rescaled",cex.main=1.2)

set.seed(0)
sample(seq((Ntip(tree)+2),(Nnode(tree)+Ntip(tree))),8)->nods
H-dist.nodes(tree)[(Ntip(tree)+1),nods]->nod.ages

sapply(1:length(nods),function(x) {
H-dist.nodes(tree)[(Ntip(tree)+1),c(getMommy(tree,nods[x])[1],getDescendants(tree,nods[x])[1:2])]->par.ages
nod.ages[x]+((min(abs(nod.ages[x]-par.ages))-0.2)*sample(c(-1,1),1))
})->nod.ages

nod.ages
scaleTree(tree,node.ages=nod.ages)->treeS2
treeS2->treeS1
unlist(lapply(1:length(nods), function(x) c(nods[x],getDescendants(tree,nods[x])[1:2])))->brcol
edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%brcol)]<-"red"

#par(mfrow=c(2,1),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("node ages rescaled",cex.main=1.2)

## ----echo=7:8, fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
H-dist.nodes(tree)[(Nnode(tree)+1),91]->sp.ages
names(sp.ages)<-tree$tip.label[1]

H-dist.nodes(tree)[(Nnode(tree)+1),164]->nod.ages
names(nod.ages)<-96

c(sp.ages,nod.ages)
scaleTree(tree,tip.ages = sp.ages,node.ages = nod.ages,min.branch = 1)->treeS

par(mfrow=c(1,2))
plot(tree,edge.color = "gray40",show.tip.label=F,no.margin = TRUE,edge.width=1.5)
plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(plotinfo$xx[1],plotinfo$yy[1],pch=16,col="blue",cex=1.2)
points(plotinfo$xx[91],plotinfo$yy[1],pch=4,col="blue",cex=1.5,lwd=2)
text("target age",x= plotinfo$xx[91]-1,y=plotinfo$yy[1],adj=1,col="blue",cex=0.8)
points(plotinfo$xx[96],plotinfo$yy[96],pch=16,col="red",cex=1.2)
points(plotinfo$xx[164],plotinfo$yy[96],pch=4,col="red",cex=1.5,lwd=2)
text("target age",x= plotinfo$xx[164]-1,y=plotinfo$yy[96],adj=1,col="red",cex=0.8)

edge.col<-rep("gray40",nrow(tree$edge))
edge.col[which(treeS$edge[,2]%in%c(1,getMommy(tree,1)))]<-"blue"
edge.col[which(treeS$edge[,2]%in%c(94,95,96))]<-"red"

plot(treeS,edge.color = edge.col,show.tip.label=F,no.margin = TRUE,edge.width=1.5)
plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(plotinfo$xx[1],plotinfo$yy[1],pch=16,col="blue",cex=1.2)
points(plotinfo$xx[96],plotinfo$yy[96],pch=16,col="red",cex=1.2)


## ----echo=c(1:16,26:37,50:56), fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
# load the RRphylo example dataset including Felids tree
data("DataFelids")
DataFelids$treefel->tree

# get species and nodes ages 
# (meant as distance from the youngest species, that is the Recent in this case)
max(nodeHeights(tree))->H
H-dist.nodes(tree)[(Ntip(tree)+1),(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))]->age.nodes
H-diag(vcv(tree))->age.tips

# apply Pagel's lambda transformation to change node ages only 
rescale(tree,"lambda",0.8)->tree1

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at nodes as node.ages
scaleTree(tree1,node.ages=age.nodes)->treeS1

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree1,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("lambda rescaled",cex.main=1.2)
plot(treeS1,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

# change leaf length of 10 sampled species
tree->tree2
set.seed(14)
sample(tree2$tip.label,10)->sam.sp
age.tips[sam.sp]->age.sam
age.sam[which(age.sam>0.1)]<-age.sam[which(age.sam>0.1)]-1.5
age.sam[which(age.sam<0.1)]<-age.sam[which(age.sam<0.1)]+0.2
tree2$edge.length[match(match(sam.sp,tree$tip.label),tree$edge[,2])]<-age.sam

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at sampled tips as tip.ages
scaleTree(tree2,tip.ages=age.tips[sam.sp])->treeS2

edge.col<-rep("black",nrow(tree$edge))
edge.col[which(treeS2$edge[,2]%in%match(sam.sp,tree$tip.label))]<-"red"

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree2,edge.color = edge.col,show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("leaves cut",cex.main=1.2)
plot(treeS2,edge.color = edge.col,show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

# apply Pagel's kappa transformation to change both species and node ages, 
# including the age at the tree root
rescale(tree,"kappa",0.5)->tree3

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at nodes as node.ages
scaleTree(tree1,tip.ages = age.tips,node.ages=age.nodes)->treeS3

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree3,edge.color = "black",show.tip.label=F)
axisPhylo(tck=-0.02,cex.axis=0.8)
title("kappa rescaled",cex.main=1.2)
plot(treeS3,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

