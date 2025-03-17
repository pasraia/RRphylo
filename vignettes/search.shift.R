## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

misspacks<-sapply(c("plotrix","scales","RColorBrewer","kableExtra"),requireNamespace,quietly=TRUE)
if(any(!misspacks)){
  warning(call. = FALSE,paste(names(misspacks)[which(!misspacks)],collapse=", "), "not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
require(kableExtra)
options(rmarkdown.html_vignette.check_title = FALSE)
options(knitr.kable.NA = '')
source("functions4vignettes.R")

## ----echo=13,warning=FALSE,message=FALSE,fig.dim=c(6,4),out.width="99%",dpi=220----
require(ape)
require(phytools)
require(scales)
cc<- 2/parallel::detectCores()

set.seed(22)
rtree(100)->Tree
fastBM(Tree)->y

y[tips(Tree,161)]*5->y[tips(Tree,161)]

RRphylo(Tree,y,clus = cc)->RR
search.shift(RR,status.type = "clade")->SSauto

layout(matrix(c(1,2),ncol=2),width=c(0.64,0.36))
plot(Tree,no.margin = TRUE, show.tip.label = FALSE)
nodelabels(bg="white",frame="none",node=as.numeric(rownames(SSauto$all.clades)),col="red",cex=0.8)

fontm<-cbind(rep(2,nrow(SSauto$all.clades)),c(rep(1,nrow(SSauto$all.clades))),c(rep(1,nrow(SSauto$all.clades))))
plotTab(tab=cbind(rownames(SSauto$all.clades),round(SSauto$all.clades,3)),
        text.cex=0.7,text.font=fontm,text.box=c(box="n"),
        text.grid=data.frame(rown=1:(nrow(SSauto$all.clades)-1),coln=rep(0,nrow(SSauto$all.clades)-1),col="gray80"),
        colN.highlight=c(col="gray97"),colN.box=list(box="o",col=c("black","gray97"),lwd=1),colN.grid = "n",colN.font=2,colN.cex=0.7,
        main="all clades differences",main.font=1,main.cex=0.7)


## ----echo=FALSE,warning=FALSE,message=FALSE, fig.dim=c(4,4),out.width="70%",dpi=220----
SSauto$single.clades->tab2
knitr::kable(tab2,digits=3,align="c", caption="single clades differences") %>%
  kable_styling(full_width = FALSE, position = "center")  %>%
  column_spec(1, bold = TRUE)

## ----echo=1:2,warning=FALSE,message=FALSE,out.width="100%"--------------------
search.shift(RR,status.type = "clade",node=c(162,134,179))->SSnode

SSnode$single.clades$singles->tab3
data.frame(rownames(tab3),tab3)->tab3
SSnode$single.clades$no.others->tab4
data.frame(rownames(tab4),tab4)->tab4
rownames(tab3)<-rownames(tab4)<-NULL
SSnode$all.clades.together->tab5
data.frame("all.clades",tab5)->tab5
cbind(tab5,rep(NA,3),tab3,rep(NA,3),tab4)->tabtot
tabtot[2:3,1:3]<-NA
# cbind(rbind(tab5,rep(NA,ncol(tab5)),rep(NA,ncol(tab5))),
#       rep(NA,3),tab3,rep(NA,3),tab4)->tabtot
colnames(tabtot)<-rep(c("","rate.difference","p.value",""),length.out=11)

knitr::kable(tabtot,digits=3,align="c") %>%
  kable_styling(full_width = FALSE, position = "center")  %>%
  column_spec(c(1,5,9), bold = TRUE) %>%
  column_spec(c(3,7),border_right = "1px solid lightgray") %>%
  add_header_above( c("\\$all.clades.together"=3," "=1,"\\$single.clades\\$singles"=3," "=1,"\\$single.clades\\$no.others"=3))



## ----echo=13:14,warning=FALSE,message=FALSE-----------------------------------
set.seed(14)
rep("a",100)->categ
categ[sample(1:100,30)]<-"b"
names(categ)<-Tree$tip.label
y[which(categ=="b")]*5->y[which(categ=="b")]
categ->two_categ

categ[sample(1:100,20)]<-"c"
y[which(categ=="c")]*1.5->y[which(categ=="c")]
categ->three_categ
RRphylo(Tree,y,clus=cc)->RR

search.shift(RR,status.type = "sparse",state=two_categ)->SSstate2
search.shift(RR,status.type = "sparse",state=three_categ)->SSstate3

cbind(rownames(SSstate2$state.results),SSstate2$state.results)->tab2
cbind(rownames(SSstate3$state.results),SSstate3$state.results)->tab3
cbind(tab2,rep(NA,6),tab3)->tabtot
tabtot[2:6,1:3]<-NA
colnames(tabtot)<-rep(c("","rate.difference","p.value",""),length.out=7)

knitr::kable(tabtot,digits=3,align="c",) %>%
  kable_styling(full_width = FALSE, position = "center")  %>%
  column_spec(c(1,5), bold = TRUE) %>%
  column_spec(c(3,7),border_right = "1px solid lightgray") %>%
  add_header_above(c("SSstate2\\$state.results"=3," "=1,"SSstate3\\$state.results"=3))


## ----out.width='99%',fig.dim=c(7,6),message=FALSE,dpi=220,warning=FALSE,echo=1:5----
# load the RRphylo example dataset including Ornithodirans tree and data
DataOrnithodirans$treedino->treedino # phylogenetic tree
DataOrnithodirans$massdino->massdino # body mass data
DataOrnithodirans$statedino->statedino # locomotory type data
log(massdino)->lmass


require(plotrix)
require(RColorBrewer)
statedino->colo
as.character(colo)->colo
names(colo)<-names(statedino)
colo[which(colo=="B")]<-brewer.pal(n = 6, name = "Set1")[6]
colo[which(colo=="BQ")]<-brewer.pal(n = 6, name = "Set1")[5]
colo[which(colo=="F")]<-brewer.pal(n = 6, name = "Set1")[3]
colo[which(colo=="Q")]<-brewer.pal(n = 6, name = "Set1")[1]
ladderize(treedino,FALSE)->treed
#par(mar=c(9,4,1,3))
plot(treed,type="fan",show.tip.label = FALSE,edge.color = "gray40",no.margin = TRUE)
tiplabels(text=rep("",419),bg=colo,frame="circle",cex=.4)
plotinfo<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)
sqrt(plotinfo$xx[416]^2+plotinfo$yy[416]^2)->rad

c(422,516,583,746,689)->clads
c("Ornithischia","Sauropodomorpha","Theropoda","Pterosauria","Avialae")->nams
#c(2,4,3,3,1)->poss
list(c(0.5,0),c(0,1),c(0.5,-0.5),c(0.3,0),c(1,0))->poss
for(k in 1:length(clads)){
  clads[k]->aa
  plotinfo$yy[which(treed$tip.label%in%tips(treed,aa)[c(1,length(tips(treed,aa)))])]->yc
  plotinfo$xx[which(treed$tip.label%in%tips(treed,aa)[c(1,length(tips(treed,aa)))])]->xc
  sqrt(xc[1]^2+yc[1]^2)->r1
  if(yc[1]>=0) acos(xc[1]/r1)->ang1 else 2*pi-acos(xc[1]/r1)->ang1
  sqrt(xc[2]^2+yc[2]^2)->r2
  if(yc[2]>=0) acos(xc[2]/r2)->ang2 else 2*pi-acos(xc[2]/r2)->ang2
  
  if(k<5){
    draw.arc(x=0,y=0,radius=rad+10,angle1 =ang1+0.05,angle2 =ang2,lwd=3,col=brewer.pal(n = 6, name = "Set1")[2]) 
    c(plotinfo$xx[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])],
      plotinfo$yy[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])])->cc
    xcc=sqrt((cc[1]^2*(180+25)^2)/(cc[1]^2+cc[2]^2))
    if(cc[1]>=0) text(x=xcc,y=cc[2]/cc[1]*xcc,labels=nams[k],font=2,adj=poss[[k]]) else text(x=(-xcc),y=cc[2]/cc[1]*-xcc,labels=nams[k],font=2,adj=poss[[k]])
  }else{
    draw.arc(x=0,y=0,radius=rad+20,angle1 =ang1+0.05,angle2 =ang2,lwd=3,col=brewer.pal(n = 6, name = "Set1")[4])
    c(plotinfo$xx[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])],
      plotinfo$yy[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])])->cc
    xcc=sqrt((cc[1]^2*(180+35)^2)/(cc[1]^2+cc[2]^2))
    if(cc[1]>=0) text(x=xcc,y=cc[2]/cc[1]*xcc,labels=nams[k],font=2,adj=poss[[k]]) else text(x=(-xcc),y=cc[2]/cc[1]*-xcc,labels=nams[k],font=2,adj=poss[[k]])
  }
}

legend("topleft",legend=c("B","BQ","F","Q"),pt.bg =unique(colo)[c(1,3,4,2)],pch=21,col="black",pt.cex=1.5,bty="n")


## ----message=FALSE,eval=FALSE-------------------------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  lmass[match(treedino$tip.label,names(lmass))]->lmass
#  statedino[match(treedino$tip.label,names(statedino))]->statedino
#  
#  # perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treedino,y=lmass)->RRdinomass
#  
#  # search for clades showing significant shifts in mass specific evolutionary rates
#  # (i.e. using the log body mass itself as a covariate)
#  search.shift(RRdinomass, status.type= "clade",cov=lmass)->SSauto
#  
#  # search for shifts in mass specific evolutionary rates pertaining different locomotory types.
#  search.shift(RRdinomass, status.type= "sparse", state=statedino,cov=lmass)->SSstate

