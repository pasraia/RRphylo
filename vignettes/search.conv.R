## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !rmarkdown::pandoc_available()) {
  warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
  knitr::knit_exit()
}

misspacks<-sapply(c("rgl","mvMORPH","RColorBrewer","kableExtra","phangorn"),requireNamespace,quietly=TRUE)
if(any(!misspacks)){
  warning(call. = FALSE,paste(names(misspacks)[which(!misspacks)],collapse=", "), "not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
require(rgl)
source("functions4vignettes.R")
options(knitr.kable.NA = '',mc.cores=2,rgl.useNULL=TRUE,
        rmarkdown.html_vignette.check_title = FALSE)

## ----include=FALSE------------------------------------------------------------
require(RRphylo)
require(RColorBrewer)
require(rgl)
require(ape)
require(phytools)
require(mvMORPH)

### creating a phylogenetic tree with 100 species
set.seed(14)
rtree(50)->tree

### select, extract and then modify the clade to be duplicated
dist.nodes(tree)[(Ntip(tree)+1),]->cd
cd[cd<0.8*max(cd)]->cd
cd[which(as.numeric(names(cd))>Ntip(tree))]->cd
cdt<-sapply(names(cd),function(e) length(tips(tree,e)))
names(cd)[which(cdt>=Ntip(tree)*0.1&cdt<=Ntip(tree)/4)]->cd
cdd<-sapply(cd,function(i) ifelse((length(tips(tree,i))-2)>(length(tips(tree,i))+Ntip(tree)-2)*0.1,"ok","no"))
if(sum(cdd=="ok")>0) cd[which(cdd=="ok")]->cd
as.numeric(sample(cd,1))->n

extract.clade(tree,n)->t1
max(nodeHeights(t1))->nH1
suppressWarnings(swapONE(t1)[[1]]->t1)
drop.tip(t1,t1$tip.label[c(1,length(t1$tip.label))])->t1

if(max(nodeHeights(t1))!=nH1) rescaleRR(t1,height=nH1)->t1
t1$root.edge<-data.frame(tree$edge,tree$edge.length)[which(tree$edge[,2]==n),3]


### selecting the node where the new clade is to be binded
distNodes(tree,n)[1:Nnode(tree),]->dfN
dfN[which(!rownames(dfN)%in%c(n,getDescendants(tree,n))),]->dfN
distNodes(tree,(Ntip(tree)+1))[2:Nnode(tree),]->dR
dR[match(rownames(dfN),rownames(dR)),]->dR1
dfN[match(rownames(dR1)[which(dR1[,2]<max(dR[,2])*0.8)],rownames(dfN)),]->dfN
rownames(dfN)[which(dfN[,1]<=Ntip(tree)/10)]->bar
dfN[which(dfN[,1]>Ntip(tree)/10),,drop=FALSE]->dfn2

if(nrow(dfn2)==0){
  rownames(dfN)[which.max(dfN[,1])]->nodN
  dfN[which.max(dfN[,1]),1]->minD
} else{
  dR[match(rownames(dfn2),rownames(dR)),,drop=FALSE]->dfn3
  rownames(dfn3)[which.min(dfn3[,2])]->nodN
  dfn2[rownames(dfn2)==nodN,1]->minD
}
nodN->tar

tree$edge.length[which(tree$edge[,2]==tar)]->pos
diff(dist.nodes(tree)[(Ntip(tree)+1),c(n,tar)])->dH
if(dH>0)
  rescaleRR(t1,height=(abs(diff(c(max(nodeHeights(t1)),dH)))))->t1 else 
    rescaleRR(t1,height=(max(nodeHeights(t1))+abs(dH)))->t1
t1$root.edge+pos/2->t1$root.edge

tree->treeO
t1->t1O


## ----noconv, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(93)
tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
bind.tree(tree,t1,where=tar,position=pos/2)->tree1
matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3)->sigma
mvSIM(tree1,param=list(ntraits=3,theta=c(0,0,0),sigma=sigma))->y
c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec

conv3dplot()

rglwidget(elementId = "plot3drgl")

## ----simpleconv, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(14)
treeO->tree
t1O->t1
matrix(c(1,0.5,0.5,0.5,1,0.5,.5,.5,1),3,3)->sigma
mvSIM(tree,param=list(ntraits=3,theta=c(0,0,0),sigma=sigma))->y
y[match(tree$tip.label,rownames(y)),]->y
y[match(tips(tree,n),rownames(y)),]->a

apply(y,2,range)->m.a
m.a[2,]*1.2->m.a

t(sapply(1:nrow(a),function(e) sapply(1:length(m.a),function(i) jitter(m.a[i],amount=(sd(a[,i])*1)))))->a1
rownames(a1)<-rownames(a)
y[match(tips(tree,n),rownames(y)),]<-a1

apply(y[match(t1$tip.label,rownames(y)),],2,jitter)->y.t1

tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
rownames(y)<-tree$tip.label
rownames(y.t1)<-t1$tip.label
bind.tree(tree,t1,where=tar,position=pos/2)->tree1
rbind(y,y.t1)->y

### nod.par is the node pair set to converge
c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec

conv3dplot()
rglwidget(elementId = "plot3drgl1")

## ----conv&par, webgl=TRUE,echo=FALSE,fig.width=6,fig.height=6,fig.align="center",message=FALSE,warning=FALSE----
set.seed(14)
treeO->tree
t1O->t1

matrix(c(1,0.5,0.5,0.5,1,0.5,.5,.5,1),3,3)->sigma
mvSIM(tree,param=list(ntraits=3,theta=c(0,0,0),sigma=sigma))->y
y[match(tree$tip.label,rownames(y)),]->y
y[match(tips(tree,n),rownames(y)),]->a


apply(y,2,range)->m.a
m.a[2,]*sample(seq(0.8,1.2,0.1),1)->m.a

t(sapply(1:nrow(a),function(e) sapply(1:length(m.a),function(i) jitter(m.a[i],amount=sd(a[,i])))))->a1
rownames(a1)<-rownames(a)
y[match(tips(tree,n),rownames(y)),]<-a1

apply(y[match(t1$tip.label,rownames(y)),],2,jitter)->y.t1

tree$tip.label<-paste("a",seq(1,Ntip(tree)),sep="")
t1$tip.label<-paste("a",seq((Ntip(tree)+1),(Ntip(tree)+Ntip(t1))),sep="")
rownames(y)<-tree$tip.label
rownames(y.t1)<-t1$tip.label
bind.tree(tree,t1,where=tar,position=pos/2)->tree1
rbind(y,y.t1)->y

c(getMRCA(tree1,tree$tip.label[which(tree$tip.label%in%tips(tree,n))]),
  getMRCA(tree1,t1$tip.label))->nod.par
apply(y,2,function(x) fastAnc(tree1,x))->ant

### sample one tip descending from each node set to converge
y[c(sample(tips(tree1,nod.par[1]),1),sample(tips(tree1,nod.par[2]),1)),]->y.vec

### modify the phenotypes at the mrcas (nod.par) as the median of the distribution of phenotypes of descending tips   
invisible(sapply(nod.par,function(v) ant[match(v,rownames(ant)),]<<-apply(y[tips(tree1,v),],2,median)))

conv3dplot()
rglwidget(elementId = "plot3drgl2")

## ----echo=FALSE,fig.width=5,fig.height=6,fig.align="center",message=FALSE,warning=FALSE,out.width='60%',dpi=200----
load("sc-data.Rda")

set.seed(22)
# rtree(14)->tree
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
wid[which(tree$edge[,2]%in%c(mom1,mom2))]<-4

set.seed(14)
rbind(RR$aces,y)->phen
tree->tree1

c(18,24)->nod
dist.nodes(tree1)[nod[1],nod[2]]->nT
tips(tree1,nod[1])->tt1
tips(tree1,nod[2])->TT
expand.grid(tt1,TT)->ctt

aa<-sapply(1:nrow(ctt),function(g){
  as.matrix(phen[match(c(as.character(ctt[g,1]),as.character(ctt[g,2])),rownames(phen)),])->ppTT
  angle.vecs(ppTT[1,],ppTT[2,])
})
mean(aa)->ang.tip

RR$aces[which(rownames(RR$aces)%in%nod),]->ac
angle.vecs(ac[1,],ac[2,])->ang.ac

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


## ----message=FALSE, warning=FALSE,eval=FALSE----------------------------------
#  search.conv(RR=RR,y=y,min.dim=3,max.dim=4,nsim=100,rsim=100,clus=2/parallel::detectCores())->SC

## ----message=FALSE, warning=FALSE,echo=FALSE----------------------------------
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


## ----echo=FALSE,fig.width=5,fig.height=6,fig.align="center",message=FALSE,warning=FALSE,out.width='60%',dpi=200----
tree->tree1
state[match(rownames(y2),names(state))]->state
if("nostate"%in%state) state[which(state!="nostate")]->state.real else state->state.real

combn(unique(state.real),2)->stcomb1
if("nostate"%in%state) combn(unique(state),2)->stcomb else stcomb1->stcomb

cophenetic.phylo(tree)->cop
tt1<-lapply(stcomb1[,1], function(w) y[which(state==w),])
expand.grid(rownames(tt1[[1]]),rownames(tt1[[2]]), stringsAsFactors = FALSE)->ctt
vs<-sapply(tt1,function(w) mean(apply(w,1,unitV)))

ang.tip<-sapply(1:nrow(ctt),function(g){
  as.matrix(y[match(as.character(ctt[g,1:2]),rownames(y)),])->ppTT
  angle.vecs(ppTT[1,],ppTT[2,])
})
dt<-apply(ctt,1,function(w) cop[match(w[1],rownames(cop)),match(w[2],rownames(cop))])

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


## ----message=FALSE, warning=FALSE,eval=FALSE----------------------------------
#  search.conv(tree=tree,y=y,state=state,nsim=100,clus=2/parallel::detectCores())->SC

## ----message=FALSE, warning=FALSE,echo=FALSE----------------------------------
knitr::kable(SCstate$state.res,digits=3,align="c") %>%
  column_spec(1:2, bold=TRUE)


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
#  search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9")->SC.clade
#  
#  ## Example 2: search for morphological convergence within sabertoothed species
#  search.conv(tree=treefel, y=PCscoresfel, state=statefel)->SC.state
#  

