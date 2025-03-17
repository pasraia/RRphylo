## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

misspacks<-sapply(c("plotrix","kableExtra"),requireNamespace,quietly=TRUE)
if(any(!misspacks)){
  warning(call. = FALSE,paste(names(misspacks)[which(!misspacks)],collapse=", "), "not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(knitr.kable.NA = '',rmarkdown.html_vignette.check_title = FALSE)
load("st-data.Rda")

## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
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
plot(NA,ylim=range(phen[,1]),xlim=range(rootP),
     mgp=c(1.8,0.5,0),
     xlab="distance from the root",ylab="rescaled phenotypes",main="Simulation for trend\nin mean phenotypes")
polygon(c(dfk2$rootP,rev(dfk2$rootP)),c(dfk2$fk,rev(dfk2$fk.min)),col = rgb(0.5, 0.5, 0.5,
                                                                        0.3), border = NA)
points(rootP,phen[,1],cex=1.5,bg=phen$col,
       col="black",pch=21)
abline(lm(phen[,1]~rootP),lwd=2,col="#ff00ff")

# rbind(data.frame(st.rates[[3]],dev=NA),data.frame(st.phen[[2]][,1:3,drop=FALSE],spread=NA,st.phen[[2]][,4,drop=FALSE]))->res
# 
# colnames(res)[4:5]<-c("spread","dev")
# rownames(res)<-c("rescaled absolute rate regression","phenotypic regression")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y)->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
require(kableExtra)

knitr::kable(res,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center")  


## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
require(plotrix)

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],rates[-which(names(rates)%in%c(desn,desn2))],
     ylim=range(rates),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="absolute rates",main="Simulation for trend\nin absolute rates")
abline(lm(rates~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],rates[which(names(rates)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="slateblue2",cex=1.5)
ablineclip(lm(rates[which(names(rates)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(rates[which(names(rates)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="slateblue2",lwd=3)
points(max(ages)-ages[which(names(ages)%in%desn2)],rates[which(names(rates)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="#ae9eff",cex=1.5)
ablineclip(lm(rates[which(names(rates)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(rates[which(names(rates)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="#ae9eff",lwd=3)
legend("topleft",legend=c("node 225","node 319","entire tree"),fill=c("slateblue2","#ae9eff","gray20"),bty="n",x.intersp = .5)

plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],phen2[-which(names(phen2)%in%c(desn,desn2))],
     ylim=range(phen2),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="phenotype",main="Simulation for trend\nin mean phenotypes")
abline(lm(phen2~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],phen2[which(names(phen2)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine3",cex=1.5)
ablineclip(lm(phen2[which(names(phen2)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(phen2[which(names(phen2)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="aquamarine3",lwd=3)

points(max(ages)-ages[which(names(ages)%in%desn2)],phen2[which(names(phen2)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine",cex=1.5)
ablineclip(lm(phen2[which(names(phen2)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(phen2[which(names(phen2)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="aquamarine",lwd=3)
legend("topleft",legend=c("node 225","node 319","entire tree"),fill=c("aquamarine3","aquamarine","gray20"),bty="n",x.intersp = .5)

cbind(data.frame(node=names(STrates[[5]]),do.call(rbind,STrates[[5]])),
      data.frame(node=names(STphen[[4]]),do.call(rbind,STphen[[4]])))->res2

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y,node=c(225,319))->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
rownames(res2)<-NULL
knitr::kable(res2,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
column_spec(6, border_right = TRUE) %>%
add_header_above(c("Trend in absolute rates" = 6, "Trend in phenotypic means" = 5))


## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
knitr::kable(STrates[[6]][[2]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in absolute rates" = 7))

knitr::kable(STphen[[6]][[1]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in phenotypic means" = 8))


## ----out.width='100%',fig.dim=c(8,7),message=FALSE,dpi=220--------------------
require(ape)
# load the RRphylo example dataset including Cetaceans tree and data
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree
DataCetaceans$masscet->masscet # logged body mass data
DataCetaceans$brainmasscet->brainmasscet # logged brain mass data
DataCetaceans$aceMyst->aceMyst # known phenotypic value for the most recent 
                               # common ancestor of Mysticeti


par(mar=c(0,0,0,1))
plot(ladderize(treecet,FALSE),show.tip.label = FALSE,edge.color = "gray40")
plotinfo<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)
nodelabels(text="",node=128,frame="circle",bg="red",cex=0.5)
nodelabels(text="Mystacodon",node=128,frame="n",bg="w",cex=1,adj=c(-0.1,0.5),font=2)
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,128))])->yran128
rect(plotinfo$x.lim[2]+0.4,yran128[1],plotinfo$x.lim[2]+0.7,yran128[2],col="red",border="red")
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,142))])->yran142
rect(plotinfo$x.lim[2]+0.4,yran142[1],plotinfo$x.lim[2]+0.7,yran142[2],col="blue",border="blue")
mtext(c("Mysticeti","Odontoceti"), side = 4,line=-0.5,at=c(sum(yran128)/2,sum(yran142)/2),
      cex=1.5,adj=0.5,col=c("red","blue"))

## ----eval=FALSE,message=FALSE,dpi=220-----------------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  masscet[match(treecet$tip.label,names(masscet))]->masscet
#  
#  ## Example 1: search.trend by setting values at internal nodes
#  # Set the body mass of Mysticetes ancestor (Mystacodon selenensis)
#  # as known value at node and perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treecet,y=masscet,aces=aceMyst)->RR
#  
#  # Perform search.trend on the RR object and (log) body mass by indicating Mysticeti as focal clade
#  search.trend(RR=RR,y=masscet,node=as.numeric(names(aceMyst)))
#  
#  
#  ## Example 2: search.trend on multiple regression version of RRphylo
#  # cross-reference the phylogenetic tree and body and brain mass data. Remove from both the tree and
#  # vector of body sizes the species whose brain size is missing
#  drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#  masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#  
#  # check the order of your data: best if
#  # data vectors (i.e. masscet and brainmasscet) are sorted
#  # in the same order of the species on the phylogeny
#  masscet.multi[match(treecet.multi$tip.label,names(masscet.multi))]->masscet.multi
#  brainmasscet[match(treecet.multi$tip.label,names(brainmasscet))]->brainmasscet
#  
#  # perform RRphylo on tree and body mass data
#  RRphylo(tree=treecet.multi,y=masscet.multi)->RRmass.multi
#  
#  # create the predictor vector: retrieve the ancestral character estimates
#  # of body size at internal nodes from the RR object ($aces) and collate them
#  # to the vector of species' body sizes to create
#  c(RRmass.multi$aces[,1],masscet.multi)->x1.mass
#  
#  # perform the multiple regression version of RRphylo by using
#  # the brain size as variable and the body size as predictor
#  RRphylo(treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#  
#  # Perform search.trend on the multiple RR object to inspect the effect of body
#  # size on absolute rates temporal trend only
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc)
#  
#  # Perform search.trend on the multiple RR object to inspect the effect of body
#  # size on trends in both absolute evolutionary rates and phenotypic values
#  # (by using brain versus body mass regression residuals)
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,
#               clus=cc)
#  

