## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(knitr.kable.NA = '',rmarkdown.html_vignette.check_title = FALSE)

## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
require(ape)
require(geiger)
require(phytools)
require(RRphylo)
require(ddpcr)

RRphylo:::range01->range01

sim.bdtree(b = .5, d = 0.2,seed=14)->tree
rt<-0
s2m<-1
es=2

set.seed(14)
fastBM(tree,a=rt,sig2=s2m)->y

### TREND ###
yt <- (diag(vcv(tree))^es)/(diag(vcv(tree))) * y
data.frame(tree$edge[,2],nodeHeights(tree)[,2])->hei
tree$tip.label[hei[which(hei[,1]<(Ntip(tree)+1)),1]]->hei[which(hei[,1]<(Ntip(tree)+1)),1]
hei[,2]->rootD
c(0.0001,rootD)->rootD
names(rootD)<-c((Ntip(tree)+1),hei[,1])

makeL(tree)->L
makeL1(tree)->L1
RRphylo(tree,yt)->rr

quiet(search.trend(rr,yt,clus=2/parallel::detectCores(),foldername=tempdir())->st.rates)

c(rr$aces,yt)->phen

rr$rates->rts
abs(rts)->rts
as.matrix(as.data.frame(rts[match(names(rootD),rownames(rts)),]))->rts
range01(rts)->rts
rootD->rootC
tend<-2

rr$ace[1]->a
rootV<-a
rr$lambda->lambda
BTS<-list()
for(i in 1:100){
  rootC->rootB
  fastBM(tree,sig2=1,a=a,bounds=c(min(yt),max(yt)))->yb
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
              t(L)) %*% (as.matrix(yb)-rootV)
  bts <- abs(betas)
  range01(bts)->BTS[[i]]
  as.matrix(as.data.frame(bts[match(names(rootB),rownames(bts)),]))->bts
  bts
}

sapply(BTS,function(x) x[match(names(rootD),rownames(x)),])->fu

apply(fu,1,function(x) quantile(x,.975))->ful
apply(fu,1,function(x) quantile(x,.025))->mins
predict(loess(ful~rootD))->fk
predict(loess(mins~rootD))->fk.min
names(fk)<-names(ful)
names(fk.min)<-names(ful)

##abline(lm(rts~rootD))
data.frame(rootD,fk,fk.min)->dfk
dfk[order(dfk[,1]),]->dfk

data.frame(rts=rts,col=rep("red",nrow(rts)))->rts
as.character(rts[,2])->rts[,2]
rts[which(rownames(rts)%in%tree$tip.label),2]<-"green"

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(NA,ylim=range(rts[,1]),xlim=range(rootD),
     mgp=c(1.8,0.5,0),
     xlab="distance from the root",ylab="rescaled rates",main="Simulation for trend\nin absolute rates")
polygon(c(dfk$rootD,rev(dfk$rootD)),c(dfk$fk,rev(dfk$fk.min)),col = rgb(0.5, 0.5, 0.5,
                                                                        0.3), border = NA)
points(rootD,rts[,1],cex=1.5,bg=rts$col,
       col="black",pch=21)
abline(lm(rts[,1]~rootD),lwd=2,col="#ff00ff")

legend("topleft",legend=c("nodes","tips","brownian range"),fill=c("red","green","#dadad9"),bty="n")

#### DRIFT
ds=1
yt1<-y+diag(vcv(tree))*ds

range01(yt1)->yd


data.frame(tree$edge[,2],nodeHeights(tree)[,2])->hei
tree$tip.label[hei[which(hei[,1]<(Ntip(tree)+1)),1]]->hei[which(hei[,1]<(Ntip(tree)+1)),1]
hei[,2]->rootD
c(0.0001,rootD)->rootD
names(rootD)<-c((Ntip(tree)+1),hei[,1])

makeL(tree)->L
makeL1(tree)->L1
RRphylo(tree,yd)->rr
quiet(search.trend(rr,yd,clus=2/parallel::detectCores(),foldername=tempdir())->st.phen)

c(rr$aces,yd)->phen
names(phen)<-c(rownames(rr$aces),names(yd))
rootD[match(names(phen),names(rootD))]->rootD

rootD->rootP

phen[[1]]->a->rootV
lm(phen~rootD)->regr

#ratematrix(tree,yt1)[1]->S2[m] 
phenD<-list()
yl<-list()
for(i in 1:100)
{
  fastBM(tree,sig2=s2m,a=rootV,bounds=c(min(yd),max(yd)))->yc
  yc->yl[[i]]
  range01(yc)->yc
  lambda <- rr$lambda
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
              t(L)) %*% (as.matrix(yc) - rootV)
  aceRR <- (L1 %*% betas[1:Nnode(tree), ]) + rootV
  c(aceRR,yc)->phenC
  names(phenC)<-c(rownames(aceRR),names(yc))
  phenC->phenD[[i]]
  #coef(summary(lm(phenC~rootP)))[2]->slopec[i]
  
}

sapply(phenD,function(x) x[match(names(rootD),names(x))])->fu

apply(fu,1,function(x) quantile(x,.975))->ful
apply(fu,1,function(x) quantile(x,.025))->mins
predict(loess(ful~rootD))->fk
predict(loess(mins~rootD))->fk.min
names(fk)<-names(ful)
names(fk.min)<-names(ful)

data.frame(rootD,fk,fk.min)->dfk
dfk[order(dfk[,1]),]->dfk

data.frame(phen=phen,col=rep("red",length(phen)))->phen
as.character(phen[,2])->phen[,2]
phen[which(rownames(phen)%in%tree$tip.label),2]<-"green"

plot(NA,ylim=range(phen[,1]),xlim=range(rootD),
     mgp=c(1.8,0.5,0),
     xlab="distance from the root",ylab="rescaled phenotypes",main="Simulation for trend\nin mean phenotypes")
polygon(c(dfk$rootD,rev(dfk$rootD)),c(dfk$fk,rev(dfk$fk.min)),col = rgb(0.5, 0.5, 0.5,
                                                                        0.3), border = NA)
points(rootD,phen[,1],cex=1.5,bg=phen$col,
       col="black",pch=21)
abline(lm(phen[,1]~rootD),lwd=2,col="#ff00ff")



as.data.frame(rbind(c(st.rates[[4]],NA),c(st.phen[[3]][1:3],NA,st.phen[[3]][4])))->res

colnames(res)[4:5]<-c("spread","dev")
rownames(res)<-c("absolute rate regression","phenotypic regression")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y,foldername=getwd())->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
require(kableExtra)

knitr::kable(res,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center")  


## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
require(plotrix)

sim.bdtree(b = .5, d = 0.2,seed=14)->T

set.seed(93)
sizedsubtree(T,20)->n
repeat({
  c(getMommy(T,n),getDescendants(T,n))->outs
  sizedsubtree(T,20)->n2
  if(n2%in%outs==FALSE & getSis(T,n2,printZoom = FALSE)!=n & n!=n2) break
})

max(nodeHeights(T))/nodeHeights(T)[,2][which(T$edge[,2]==n)]->dN1
max(nodeHeights(T))/nodeHeights(T)[,2][which(T$edge[,2]==n2)]->dN2
ds1<-.5
ds2<- -.5
es1=.01
es2=1.9
runif(1,-10,10)->rt
S2<-sample(seq(0.1,10,0.01),1)
sample(c(es1,es2),1)->esx1
sample(c(es1,es2),1)->esx2

fastBM(T,a=rt,sig2=S2)->yB
yB->y

### Phenotypes ###
y[match(tips(T,n),names(y))]->yDR
diag(vcv(T)[names(yDR),names(yDR)])->timeDR
coef(lm(yDR~timeDR))[2]->tend
if(tend*ds1>0) ds1->ds1 else ds1*(-1)->ds1
ds1*dN1->dsx1
yD <- (timeDR*dsx1)+ yDR
y[match(tips(T,n),names(y))]<-yD

y[match(tips(T,n2),names(y))]->yDR
diag(vcv(T)[names(yDR),names(yDR)])->timeDR
coef(lm(yDR~timeDR))[2]->tend
if(tend*ds2>0) ds2->ds2 else ds2*(-1)->ds2
ds2*dN2->dsx2
yD <- (timeDR*dsx2)+ yDR
y[match(tips(T,n2),names(y))]<-yD

RRphylo(T,y)->RR
quiet(search.trend(RR,y,node=c(n,n2),clus=2/parallel::detectCores(),foldername=tempdir())->STphen)
c(RR$aces[,1],y)->phen

c(n,getDescendants(T,n))->desn
T$tip.label[desn[which(desn<=Ntip(T))]]->desn[which(desn<=Ntip(T))]

c(n2,getDescendants(T,n2))->desn2
T$tip.label[desn2[which(desn2<=Ntip(T))]]->desn2[which(desn2<=Ntip(T))]

dist.nodes(T)[Ntip(T)+1,]->ages
T$tip.label[as.numeric(names(ages)[which(as.numeric(names(ages))<=Ntip(T))])]->names(ages)[which(as.numeric(names(ages))<=Ntip(T))]
ages[match(names(phen),names(ages))]->ages


### Rates ###
es1=.01
es2=1.9
sample(c(es1,es2),1)->esx1
sample(c(es1,es2),1)->esx2

yB->y
y[match(tips(T,n),names(y))]->yT
diag(vcv(T)[names(yT),names(yT)])->timeT
yT <- ((timeT^esx1)/timeT)*yT

y[match(tips(T,n2),names(y))]->yT2
diag(vcv(T)[names(yT2),names(yT2)])->timeT2
yT2 <- ((timeT2^esx2)/timeT2)*yT2
y[match(tips(T,n),names(y))]<-yT
y[match(tips(T,n2),names(y))]<-yT2

RRphylo(T,y)->RR
quiet(search.trend(RR,y,node=c(n,n2),clus=2/parallel::detectCores(),foldername=tempdir())->STrates)
abs(RR$rates[,1])->rats


par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],rats[-which(names(rats)%in%c(desn,desn2))],
     ylim=range(rats),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="absolute rates",main="Simulation for trend\nin absolute rates")
abline(lm(rats~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],rats[which(names(rats)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="slateblue2",cex=1.5)
ablineclip(lm(rats[which(names(rats)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(rats[which(names(rats)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="slateblue2",lwd=3)
points(max(ages)-ages[which(names(ages)%in%desn2)],rats[which(names(rats)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="#ae9eff",cex=1.5)
ablineclip(lm(rats[which(names(rats)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(rats[which(names(rats)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="#ae9eff",lwd=3)
legend("topleft",legend=c("node 275","node 198","entire tree"),fill=c("slateblue2","#ae9eff","gray20"),bty="n",x.intersp = .5)

plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],phen[-which(names(phen)%in%c(desn,desn2))],
     ylim=range(phen),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="phenotype",main="Simulation for trend\nin mean phenotypes")
abline(lm(phen~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],phen[which(names(phen)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine3",cex=1.5)
ablineclip(lm(phen[which(names(phen)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(phen[which(names(phen)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="aquamarine3",lwd=3)

points(max(ages)-ages[which(names(ages)%in%desn2)],phen[which(names(phen)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine",cex=1.5)
ablineclip(lm(phen[which(names(phen)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(phen[which(names(phen)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="aquamarine",lwd=3)
legend("topleft",legend=c("node 275","node 198","entire tree"),fill=c("aquamarine3","aquamarine","gray20"),bty="n",x.intersp = .5)

cbind(data.frame(node=names(STrates[[6]]),do.call(rbind,STrates[[6]])),
      data.frame(node=names(STphen[[5]]),do.call(rbind,STphen[[5]])))->res

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y,node=c(275,198),foldername=getwd())->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
rownames(res)<-NULL
knitr::kable(res,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
column_spec(5, border_right = TRUE) %>%
add_header_above(c("Trend in absolute rates" = 5, "Trend in phenotypic means" = 5))


## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
knitr::kable(STrates[[7]][[2]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in absolute rates" = 6))

knitr::kable(STphen[[7]][[1]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in phenotypic means" = 6))


## ----out.width='100%',fig.dim=c(8,7),message=FALSE,dpi=220--------------------
# load the RRphylo example dataset including Cetaceans tree and data
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree
DataCetaceans$masscet->masscet # logged body mass data
DataCetaceans$brainmasscet->brainmasscet # logged brain mass data
DataCetaceans$aceMyst->aceMyst # known phenotypic value for the most recent common ancestor of Mysticeti

require(ggtree)
ggtree(ape::ladderize(treecet))+
  geom_cladelabel(node=128,label="Mysticeti", align=TRUE, angle=270, hjust='center',
                  fontsize=6, offset.text=.6, barsize=1.5,color="red")+
  geom_cladelabel(node=142,label="Odontoceti", align=TRUE, angle=270, hjust='center',
                  fontsize=6, offset.text=.6, barsize=1.5,color="blue")->p
  p+geom_point2(aes(subset=(node==128)), size=3, color='red')+
  geom_text(aes(x=p$data[which(p$data$node==128),]$x,
                y=p$data[which(p$data$node==128),]$y,
                label="Mystacodon"),hjust=-0.08)

## ----eval=FALSE,message=FALSE,dpi=220-----------------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  masscet[match(treecet$tip.label,names(masscet))]->masscet
#  
#  # Set the body mass of Mysticetes ancestor (Mystacodon selenensis)
#  # as known value at node and perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treecet,y=masscet,aces=aceMyst)->RR
#  
#  # Perform search.trend on the RR object and (log) body mass by indicating Mysticeti as focal clade
#  search.trend(RR=RR,y=masscet,node=as.numeric(names(aceMyst)),foldername=tempdir())
#  

