## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}


misspacks<-sapply(c("manipulate","scales","kableExtra"),requireNamespace,quietly=TRUE)
if(any(!misspacks)){
  warning(call. = FALSE,paste(names(misspacks)[which(!misspacks)],collapse=", "), "not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)
source("functions4vignettes.R")

## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,6),out.width="98%",dpi=220----
require(ape)
require(phytools)

set.seed(14)
rtree(10,tip.label=c("genusB_1","genusD_1","genusB_2","genusA_1","genusC_1",
                     "genusE_1","genusC_2","genusD_2","genusB_3","genusC_3"))->tree.back
tree.back$node.label[6]<-"DC"

par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
plot(tree.back,edge.width=1.2)
nodelabels(text=tree.back$node.label[6],node=Ntip(tree.back)+6,col="forestgreen",frame="n",adj=0)
title("backbone tree")
axisPhylo()

data.frame(bind=c("genusE_1a","genusC_3a","genusB_1a","genusC_4","genusF_1",
                  "genusG_1","genusH_1","genusC_5","genusC_6","genusB_3","genusE_1b","genusI_1"),
           reference=c("genusE_1","genusC_3","genusB_1","Genus genusC","Clade DC",
                       "genusF_1-genusB_2","genusE_1a-genusC_3","Genus genusC",
                       "genusC_5-genusC_1",
                       "genusB_2-genusB_1","genusE_1","Genus genusC"),
           poly=c(FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))->dato
par(mar=c(2,0,1,0.5))
plotTab(tab=dato,
        text.cex=0.9,text.box=c(box="n"),
        text.grid=data.frame(rown=1:nrow(dato),coln=rep(0,nrow(dato)),col="gray80"),
        colN.highlight=c(col="gray97"),colN.box=list(box="c",col=c("black","gray97","gray80"),lwd=1),
        colN.grid = "n",colN.font=2,colN.cex=0.9,
        main="dato",main.font=2,main.cex=0.9,main.highlight =c(col="gray97"),main.box = list(box="7",col=c("black","gray97")))

require(kableExtra)

## ----echo=1,results = "hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato,plot=FALSE)

suppressWarnings(tm(backbone=tree.back,data=dato,title="merged tree"))

## ----echo=7:8,message=FALSE---------------------------------------------------

c(1,2,1.7,1.5,0.8,1.5,0.3,1.2,0.2)->ages.tip
c("genusH_1","genusE_1a","genusE_1","genusE_1b","genusF_1","genusC_5","genusC_3a","genusG_1","genusB_1a")->names(ages.tip)
c(2.2,2.9,3.5)->ages.node
c("genusB_1-genusF_1","genusE_1a-genusE_1b","genusH_1-genusB_1")->names(ages.node)

ages.tip
ages.node

## ----echo=1,results="hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato,tip.ages=ages.tip,node.ages = ages.node,plot=FALSE)

suppressWarnings(tm(backbone=tree.back,data=dato,tip.ages=ages.tip,node.ages = ages.node,title="merged and calibrated tree"))

## ----echo=FALSE,warnings=FALSE,message=FALSE,fig.dim=c(8,4),out.width="98%",dpi=220----
 set.seed(1)
 rtree(13,tip.label=c("genusK_1","genusA_2","genusA_4","genusM_1","genusH_1","genusF_2","genusF_1",
 "genusJ_1","genusI_1","genusB_5","genusA_3","genusN_1","genusG_1") )->tree.source
 tree.source$node.labels[7]<-"HI"

par(mfrow=c(1,2),mar=c(2,0,1,0))
plot(tree.back,edge.width=1.8)
nodelabels(text=tree.back$node.label[6],node=Ntip(tree.back)+6,col="forestgreen",frame="n",adj=0)
title("backbone tree")
axisPhylo()
plot(tree.source,edge.width=1.8)
nodelabels(text=tree.source$node.label[7],node=Ntip(tree.source)+7,col="forestgreen",frame="n",adj=0)
title("source tree")
axisPhylo()

data.frame(bind=c("Genus genusA","genusG_1-genusF_2","Clade HI","genusL_1","genusB_4","genusM_1-genusN_1"),
           reference=c("genusA_1","Clade DC","Genus genusB","Genus genusA","genusB_3","Genus genusB"),
           poly=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))->dato.clade

knitr::kable(dato.clade,align="c") %>%
  kable_styling(full_width = TRUE, position = "center") %>%
  add_header_above(c("dato.clade" = 3)) 


## ----echo=1,results = "hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato.clade,source.tree=tree.source,plot=FALSE)

tm(backbone=tree.back,data=dato.clade,source.tree=tree.source,title="merged tree")


## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
dato.new<-data.frame(bind=c("sp_1","sp_3","sp_4","sp_5","sp_6","sp_7","sp_8"),
           reference=c("sp_2","sp_1-sp_2","sp_3","sp_3-sp_4","sp_1-sp_4","sp_6","sp_1-sp_6"),
           poly=rep(FALSE,7))

## ----echo=1,results = "hide",message=FALSE------------------------------------
tree.merger(data=dato.new,plot=FALSE)

## ----echo=FALSE,results = "hide",message=FALSE,fig.dim=c(6,4),out.width="99%",dpi=220,fig.align='center'----
layout(matrix(c(1,0,2,2),ncol=2,nrow=2),width=c(3,7),height=c(6,4))
par(mar=c(0,1,0,0))
plotTab(tab=dato.new,
        text.cex=0.9,text.box=c(box="n"),
        text.grid=data.frame(rown=1:nrow(dato),coln=rep(0,nrow(dato)),col="gray80"),
        colN.highlight=c(col="gray97"),colN.box=list(box="c",col=c("black","gray97","gray80"),lwd=1),
        colN.grid = "n",colN.font=2,colN.cex=0.9,
        main="dato.new",main.font=2,main.cex=0.9,main.highlight =c(col="gray97"),main.box = list(box="7",col=c("black","gray97")))
suppressWarnings(tm(data=dato.new,title="new tree"))

## ----echo=c(1:19),message=FALSE,fig.dim=c(6,6),out.width="98%",dpi=220,fig.align='center'----
#### Merging phylogenetic information 
### load the RRphylo example dataset including Cetaceans tree 
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree
treecet$node.label[(131-Ntip(treecet))]<-"Crown Mysticeti" # assigning node labels

### Select two clades and some species to be removed
tips(treecet,131)->crown.Mysticetes
tips(treecet,193)->Delphininae
c("Aetiocetus_weltoni","Saghacetus_osiris",
  "Zygorhiza_kochii","Ambulocetus_natans",
  "Kentriodon_pernix","Kentriodon_schneideri","Kentriodon_obscurus",
  "Eurhinodelphis_cristatus","Eurhinodelphis_bossi")->extinct

plot(treecet,show.tip.label = FALSE,no.margin=TRUE)
nodelabels(frame="n",col="blue",font=2,node=c(131,193),text=c("crown\nMysticetes","Delphininae"))
tiplabels(frame="circle",bg="red",cex=.3,text=rep("",length(c(crown.Mysticetes,Delphininae,extinct))),
          tip=which(treecet$tip.label%in%c(crown.Mysticetes,Delphininae,extinct)))

### Create the backbone and source trees
drop.tip(treecet,c(crown.Mysticetes[-which(tips(treecet,131)%in%
                                             c("Caperea_marginata","Eubalaena_australis"))],
                   Delphininae[-which(tips(treecet,193)=="Tursiops_aduncus")],extinct))->backtree
keep.tip(treecet,c(crown.Mysticetes,Delphininae,extinct))->sourcetree


### Create the data object
data.frame(bind=c("Clade Crown Mysticeti",
                  "Aetiocetus_weltoni",
                  "Saghacetus_osiris",
                  "Zygorhiza_kochii",
                  "Ambulocetus_natans",
                  "Genus Kentriodon",
                  "Sousa_chinensis-Delphinus_delphis",
                  "Kogia_sima",
                  "Eurhinodelphis_cristatus",
                  "Grampus_griseus",
                  "Eurhinodelphis_bossi"),
           reference=c("Fucaia_buelli-Aetiocetus_weltoni",
                       "Aetiocetus_cotylalveus",
                       "Fucaia_buelli-Tursiops_truncatus",
                       "Saghacetus_osiris-Fucaia_buelli",
                       "Dalanistes_ahmedi-Fucaia_buelli",
                       "Phocoena_phocoena-Delphinus_delphis",
                       "Sotalia_fluviatilis",
                       "Kogia_breviceps",
                       "Eurhinodelphis_longirostris",
                       "Globicephala_melas-Pseudorca_crassidens",
                       "Eurhinodelphis_longirostris"),
           poly=c(FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE))->dato

knitr::kable(dato,align="c") %>%
  kable_styling(full_width = TRUE, position = "center") %>%
  add_header_above(c("dato" = 3)) 

## ----results="hide",message=FALSE---------------------------------------------
### Merge the backbone and the source trees according to dato without calibrating tip and node ages
tree.merger(backbone = backtree,data=dato,source.tree = sourcetree,plot=FALSE)

### Set tips and nodes calibration ages
c(Aetiocetus_weltoni=28.0,
  Saghacetus_osiris=33.9,
  Zygorhiza_kochii=34.0,
  Ambulocetus_natans=40.4,
  Kentriodon_pernix=15.9,
  Kentriodon_schneideri=11.61,
  Kentriodon_obscurus=13.65,
  Eurhinodelphis_bossi=13.65,
  Eurhinodelphis_cristatus=5.33)->tipages
c("Ambulocetus_natans-Fucaia_buelli"=52.6,
  "Balaena_mysticetus-Caperea_marginata"=21.5)->nodeages

### Merge the backbone and the source trees and calibrate tips and nodes ages
tree.merger(backbone = backtree,data=dato,source.tree = sourcetree,
            tip.ages=tipages,node.ages=nodeages,plot=FALSE)


## ----message=FALSE,results="hide"---------------------------------------------
#### Building a new phylogenetic tree: build the phylogenetic tree shown in 
#### Pandolfi et al. 2020 - Figure 2 (see reference)

### Create the data object
data.frame(bind=c("Hippopotamus_lemerlei",
                  "Hippopotamus_pentlandi",
                  "Hippopotamus_amphibius",
                  "Hippopotamus_antiquus",
                  "Hippopotamus_gorgops",
                  "Hippopotamus_afarensis",
                  "Hexaprotodon_sivalensis",
                  "Hexaprotodon_palaeindicus",
                  "Archaeopotamus_harvardi",
                  "Saotherium_mingoz",
                  "Choeropsis_liberiensis"),
           reference=c("Hippopotamus_madagascariensis",
                       "Hippopotamus_madagascariensis-Hippopotamus_lemerlei",
                       "Hippopotamus_pentlandi-Hippopotamus_madagascariensis",
                       "Hippopotamus_amphibius-Hippopotamus_madagascariensis",
                       "Hippopotamus_antiquus-Hippopotamus_madagascariensis",
                       "Hippopotamus_gorgops-Hippopotamus_madagascariensis",
                       "Genus Hippopotamus",
                       "Hexaprotodon_sivalensis",
                       "Hexaprotodon_sivalensis-Hippopotamus_madagascariensis",
                       "Archaeopotamus_harvardi-Hippopotamus_madagascariensis",
                       "Saotherium_mingoz-Hippopotamus_madagascariensis"),
           poly=c(FALSE,
                  TRUE,
                  FALSE,
                  FALSE,
                  TRUE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE,
                  FALSE))->dato

## ----echo=FALSE,message=FALSE-------------------------------------------------
knitr::kable(dato,align="c") %>%
  kable_styling(full_width = TRUE, position = "center")

## ----message=FALSE,results="hide"---------------------------------------------
### Build an uncalibrated version of the tree
tree.merger(data=dato,plot=FALSE)->tree.uncal

### Set tips and nodes calibration ages
## Please note: the following ages are only used to show how to use the function
## they are not assumed to be correct.
c("Hippopotamus_lemerlei"=0.001,
  "Hippopotamus_pentlandi"=0.45,
  "Hippopotamus_amphibius"=0,
  "Hippopotamus_antiquus"=0.5,
  "Hippopotamus_gorgops"=0.4,
  "Hippopotamus_afarensis"=0.75,
  "Hexaprotodon_sivalensis"=1,
  "Hexaprotodon_palaeindicus"=0.4,
  "Archaeopotamus_harvardi"=5.2,
  "Saotherium_mingoz"=4,
  "Choeropsis_liberiensis"=0)->tip.ages
c("Choeropsis_liberiensis-Hippopotamus_amphibius"=13,
  "Archaeopotamus_harvardi-Hippopotamus_amphibius"=8.5,
  "Hexaprotodon_sivalensis-Hexaprotodon_palaeindicus"=6)->node.ages


### Build a calibrated version of the tree
tree.merger(data=dato,tip.ages=tip.ages,node.ages=node.ages,plot=FALSE)->tree.cal

## ----echo=13, message=FALSE, warning=FALSE------------------------------------
library(ape)
library(phytools)

set.seed(14)
DataFelids$treefel->tree
tree$tip.label<-paste("t",seq(1:Ntip(tree)),sep="")

max(nodeHeights(tree))->H
sample(H-diag(vcv(tree)),8)->sp.ages
sp.ages+tree$edge.length[match(match(names(sp.ages),tree$tip.label),tree$edge[,2])]/2->sp.ages
sp.ages[c(3,7)]<-0

sp.ages

## ----echo=c(1),results="hide"-------------------------------------------------
scaleTree(tree,tip.ages=sp.ages)->treeS1

## ----echo=c(19), fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%match(names(sp.ages),tree$tip.label))]<-"blue"

par(mfcol=c(1,2),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("species ages rescaled",cex.main=1.2)

## ----echo=10------------------------------------------------------------------
set.seed(0)
sample(seq((Ntip(tree)+2),(Nnode(tree)+Ntip(tree))),8)->nods
H-dist.nodes(tree)[(Ntip(tree)+1),nods]->nod.ages

sapply(1:length(nods),function(x) {
H-dist.nodes(tree)[(Ntip(tree)+1),c(getMommy(tree,nods[x])[1],getDescendants(tree,nods[x])[1:2])]->par.ages
nod.ages[x]+((min(abs(nod.ages[x]-par.ages))-0.2)*sample(c(-1,1),1))
})->nod.ages

nod.ages

## ----echo=c(1),results="hide", fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
scaleTree(tree,node.ages=nod.ages)->treeS2
treeS2->treeS1
unlist(lapply(1:length(nods), function(x) c(nods[x],getDescendants(tree,nods[x])[1:2])))->brcol
edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%brcol)]<-"red"

#par(mfrow=c(2,1),mar=c(0.1,0.1,1,0.1))
par(mfcol=c(1,2),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("node ages rescaled",cex.main=1.2)

## ----echo=7, fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
H-dist.nodes(tree)[(Nnode(tree)+1),91]->sp.ages
names(sp.ages)<-tree$tip.label[1]

H-dist.nodes(tree)[(Nnode(tree)+1),164]->nod.ages
names(nod.ages)<-96

c(sp.ages,nod.ages)

## ----echo=1,results="hide", fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
scaleTree(tree,tip.ages = sp.ages,node.ages = nod.ages)->treeS

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


## ----echo=c(1:16,26:37,50:56), fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center",results="hide"----
# load the RRphylo example dataset including Felids tree
data("DataFelids")
DataFelids$treefel->tree

# get species and nodes ages 
# (meant as distance from the youngest species, that is the Recent in this case)
max(nodeHeights(tree))->H
H-dist.nodes(tree)[(Ntip(tree)+1),(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))]->age.nodes
H-diag(vcv(tree))->age.tips

# apply Pagel's lambda transformation to change node ages only
rescaleRR(tree,lambda=0.8)->tree1

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
rescaleRR(tree,kappa=0.5)->tree3

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

## ----results="hide"-----------------------------------------------------------
require(phytools)
DataCetaceans$tree->treecet
###  Moving a single tip
## sister to a tip
move.lineage(treecet,focal="Orcinus_orca",sister="Balaenoptera_musculus")->mol1
## sister to a clade
move.lineage(treecet,focal="Orcinus_orca",sister=131)->mol2

### Moving a clade
## sister to a tip
move.lineage(treecet,focal="Genus Mesoplodon",sister="Balaenoptera_musculus")->mol7
## sister to a clade
move.lineage(treecet,focal="Clade Delphinida",sister=131)->mol11
## sister to a clade by using treecet$node.label
move.lineage(treecet,focal="Clade Delphinida",sister="Clade Plicogulae")->mol14
## sister to the tree root with and without rootage
move.lineage(treecet,focal="Genus Mesoplodon",sister=117)->mol19
move.lineage(treecet,focal="Clade Delphinida",
            sister=117,rootage=max(diag(vcv(treecet))))->mol23


## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
DataFelids$treefel->tree
max(nodeHeights(tree))->H

par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
plot(tree,show.tip.label=F)
axis(side=1,at=seq(2,32,5),labels=rev(c(0,5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
abline(v=H-5,col="blue")
title("original",cex.main=1.2)

plot(tree,show.tip.label=F)
axis(side=1,at=seq(2,32,5),labels=rev(c(0,5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
title("original",cex.main=1.2)
points(plotinfo$xx[129],plotinfo$yy[129],pch=16,col="red",cex=1.2)

## ----eval=FALSE---------------------------------------------------------------
# cutPhylo(tree,age=5,keep.lineage = TRUE)
# cutPhylo(tree,age=5,keep.lineage = FALSE)

## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
age<-5
{
  distNodes(tree,(Ntip(tree)+1))->dN
  max(nodeHeights(tree))-age->cutT
  dN[,2]->dd
  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter
  
  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
    #}
  }else{
    ### Tips and nodes ###
    as.numeric(cutter[which(suppressWarnings(as.numeric(cutter))>Ntip(tree))])->cutn
    i=1
    while(i<=length(cutn)){
      if(any(cutn%in%getDescendants(tree,cutn[i]))) cutn[-which(cutn%in%getDescendants(tree,cutn[i]))]->cutn
      i=i+1
    }
    
    tree->tt
    i=1
    while(i<=length(cutn)){
      getMRCA(tt,tips(tree,cutn[i]))->nn
      drop.clade(tt,tips(tt,nn))->tt
      tt$tip.label[which(tt$tip.label=="NA")]<-paste("l",i,sep="")
      i=i+1
    }
    
    diag(vcv(tt))->times
    if(any(times>cutT)){
      times[which(times>cutT)]->times
      tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]<-
        tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]-(times-cutT)
    }
    
  }
  
  tt->tt.age
  plot(tt,show.tip.label=FALSE)
  tiplabels(frame="n",col="blue",cex=.6,offset=0.5,
            tip=which(!tt$tip.label%in%tree$tip.label),
            text=tt$tip.label[which(!tt$tip.label%in%tree$tip.label)])
  
  axis(side=1,at=seq(2,27,5),labels=rev(c(5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
  title("cut at 5: keeping lineages",cex.main=1.2)
  
   drop.tip(tt.age,which(!tt.age$tip.label%in%tree$tip.label))->tt.age
   plot(tt.age,show.tip.label=FALSE)
   axis(side=1,at=seq(2,27,5),labels=rev(c(5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
   title("cut at 5: removing lineages",cex.main=1.2)
  
}

## ----eval=FALSE---------------------------------------------------------------
# cutPhylo(tree,node=129,keep.lineage = TRUE)
# cutPhylo(tree,node=129,keep.lineage = FALSE)

## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
node<-129
{
  distNodes(tree,(Ntip(tree)+1))->dN
  dN[match(node,rownames(dN)),2]->cutT
  dN[,2]->dd
  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter
  
  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
    #}
  }else{
    ### Tips and nodes ###
    as.numeric(cutter[which(suppressWarnings(as.numeric(cutter))>Ntip(tree))])->cutn
    i=1
    while(i<=length(cutn)){
      if(any(cutn%in%getDescendants(tree,cutn[i]))) cutn[-which(cutn%in%getDescendants(tree,cutn[i]))]->cutn
      i=i+1
    }
    
    tree->tt
    i=1
    while(i<=length(cutn)){
      getMRCA(tt,tips(tree,cutn[i]))->nn
      drop.clade(tt,tips(tt,nn))->tt
      tt$tip.label[which(tt$tip.label=="NA")]<-paste("l",i,sep="")
      i=i+1
    }
    
    diag(vcv(tt))->times
    if(any(times>cutT)){
      times[which(times>cutT)]->times
      tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]<-
        tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]-(times-cutT)
    }
    
  }
  
  tt->tt.node
  par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
  plot(tt,show.tip.label=FALSE)
  tiplabels(frame="n",col="red",cex=.6,offset=0.5,
            tip=which(!tt$tip.label%in%tree$tip.label),
            text=tt$tip.label[which(!tt$tip.label%in%tree$tip.label)])
  
  axis(side=1,at=seq(2,23.5,5),labels=rev(c(10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
  title("cut at node: keeping lineages",cex.main=1.2)
  
   drop.tip(tt.node,which(!tt.node$tip.label%in%tree$tip.label))->tt.node
   plot(tt.node,show.tip.label=FALSE)
   axis(side=1,at=seq(2,23.5,5),labels=rev(c(10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
   title("cut at node: removing lineages",cex.main=1.2)
  
}


## ----warnings=FALSE,message=FALSE,fig.dim=c(6,3),out.width="98%",dpi=220,results="hide"----
 ### load the RRphylo example dataset including Cetaceans tree 
 data("DataCetaceans")
 DataCetaceans$treecet->treecet

 ### Resolve all the polytomies within Cetaceans phylogeny
 fix.poly(treecet,type="resolve")->treecet.fixed
 
 ## Set branch colors
 unlist(sapply(names(which(table(treecet$edge[,1])>2)),function(x) 
   c(x,getDescendants(treecet,as.numeric(x)))))->tocolo
 unlist(sapply(names(which(table(treecet$edge[,1])>2)),function(x) 
   c(getMRCA(treecet.fixed,tips(treecet,x)),
     getDescendants(treecet.fixed,as.numeric(getMRCA(treecet.fixed,tips(treecet,x)))))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.fixed$edge))
 names(colo2)<-treecet.fixed$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.fixed,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)

 ### Resolve the polytomies pertaining the genus Kentriodon
 fix.poly(treecet,type="resolve",node=221)->treecet.fixed2
 
 ## Set branch colors
 c(221,getDescendants(treecet,as.numeric(221)))->tocolo
 c(getMRCA(treecet.fixed2,tips(treecet,221)),
   getDescendants(treecet.fixed2,as.numeric(getMRCA(treecet.fixed2,tips(treecet,221)))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.fixed2$edge))
 names(colo2)<-treecet.fixed2$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.fixed2,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)

 ### Collapse Delphinidae into a polytomous clade
 fix.poly(treecet,type="collapse",node=179)->treecet.collapsed
 
 # Set branch colors
 c(179,getDescendants(treecet,as.numeric(179)))->tocolo
 c(getMRCA(treecet.collapsed,tips(treecet,179)),
   getDescendants(treecet.collapsed,as.numeric(getMRCA(treecet.collapsed,tips(treecet,179)))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.collapsed$edge))
 names(colo2)<-treecet.collapsed$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.collapsed,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)


