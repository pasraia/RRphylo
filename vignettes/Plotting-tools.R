## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

if (!requireNamespace("scales", quietly = TRUE)) {
   warning(call. = FALSE, "scales not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)
# load("plotting-data.Rda")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ### load the RRphylo example dataset including Apes tree and shape data (multivariate)
#  data("DataApes")
#  DataApes$PCstage->PCstage
#  DataApes$Tstage->Tstage
#  
#  # Perform RRphylo
#  RRphylo(tree=Tstage,y=PCstage)->RR
#  
#  # Visualize phenotypic values and "multiple" evolutionary rates onto the
#  # phylogenetic tree.
#  plotRR(RR,y=PCstage,multivariate="multiple.rates")->pRRmulti
#  
#  # Visualize phenotypic values and the norm2 vector of multiple rates onto the
#  # phylogenetic tree.
#  plotRR(RR,y=PCstage,multivariate="rates")->pRR
#  
#  # Using default options
#  pRRmulti$plotRRphen(variable=1,colorbar.args = list(title.pos="bottom"))
#  pRRmulti$plotRRrates(variable=2,colorbar.args = list(title.pos="bottom"))

## ----plotRR1,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotRR1.png")

## ----eval=FALSE ,message=FALSE,warning=FALSE----------------------------------
#  # Using default options
#  pRR$plotRRphen(variable=1,colorbar.args = list(title.pos="bottom")) # this is the same as before of course
#  pRR$plotRRrates(colorbar.args = list(title.pos="bottom"))

## ----plotRR2,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotRR2.png")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  # Customizing parameters
#  pRRmulti$plotRRphen(variable=3,tree.args=list(edge.width=2,cex=0.7),color.pal=heat.colors,
#                 colorbar.args = list(x="topleft",direction="h",labs.adj=0.7,xpd=TRUE,title.pos="bottom"))
#  pRRmulti$plotRRrates(variable=4,tree.args=list(edge.width=2,direction="leftwards",cex=0.7),
#                  color.pal=terrain.colors,colorbar.args = list(x="topright",labs.adj=0.7,xpd=TRUE,tck.pos="out"))
#  

## ----plotRR3,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotRR3.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ### load the RRphylo example dataset including Cetaceans tree and ln body masses
#  data("DataCetaceans")
#  DataCetaceans$treecet->treecet
#  DataCetaceans$masscet->masscet
#  
#  # Perform RRphylo
#  RRphylo(treecet,masscet)->RRcet
#  
#  # Visualize evolutionary rates computed for the clade including extant Mysticetes
#  plotRates(RRcet,node=131)->pRates
#  
#  # Using default options
#  par(mar=c(3,3,1,1),mgp=c(1.7,0.5,0))
#  pRates$plotHist()
#  par(mar=c(3,12,1,1),mgp=c(1.7,0.5,0))
#  pRates$plotLollipop()->pLol

## ----plotRates1, echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotRates1.png")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  # Customizing parameters
#  par(mar=c(3,3,1,1),mgp=c(1.7,0.5,0))
#  pRates$plotHist(hist.args=list(yaxt="n",ylab=NULL,col1="blue",col2="cyan"),
#                  legend.args = list(pch=21,pt.cex=2))
#  
#  # Setting different pch and color for nodes and species
#  colo<-rep("gray70",length(pLol))
#  colo[names(pLol)%in%treecet$tip.label]<-"cyan"
#  pchc<-rep(22,length(pLol))
#  pchc[names(pLol)%in%treecet$tip.label]<-21
#  
#  par(mar=c(3,12,1,1),mgp=c(1.7,0.5,0))
#  pRates$plotLollipop(lollipop.args = list(pch=pchc,bg=colo,col=colo,cex=2,lwd=2),
#                      line.args = list(col="deeppink2",lwd=4,lty=2))->pLol2
#  
#  

## ----plotRates2, echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotRates2.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ### load the RRphylo example dataset including Felids tree and 4 mandible shape variables (PCs)
#  data("DataFelids")
#  DataFelids$treefel->treefel
#  DataFelids$PCscoresfel->PCscoresfel
#  
#  # Perform RRphylo and search.trend
#  RRphylo(treefel,PCscoresfel)->RRfel
#  search.trend(RRfel,PCscoresfel)->ST
#  
#  # Visualize search.trend results
#  plotTrend(ST)->pTrend
#  
#  # Using default options
#  par(mfrow=c(1,2),mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
#  pTrend$plotSTphen("PC1")
#  pTrend$plotSTrates(1) # This is for PC1 as well
#  

## ----plotTrend1,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotTrend1.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  ## Customizing parameters
#  # Extracting nodes and tips descending from the most recent common ancestor of Felinae
#  library(phytools)
#  Felinae<-getDescendants(treefel,94)
#  Felinae[which(Felinae<=Ntip(treefel))]<-treefel$tip.label[Felinae[which(Felinae<=Ntip(treefel))]]
#  
#  # Setting different pch for Felinae only
#  pchc<-rep(21,nrow(ST$trend.data$phenotypeVStime))
#  pchc[rownames(ST$trend.data$phenotypeVStime)%in%Felinae]<-6
#  
#  
#  par(mfrow=c(1,2),mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
#  pTrend$plotSTphen("PC2",plot.args = list(pch=pchc,cex=1.3),
#                  polygon.args = list(col="white",border="black",lwd=2,lty=3),
#                  line.args = list(lty=4,col="purple3"))
#  
#  # When multivariate data are used, a "rate" vector is calculated as the 2-norm
#  # vector of rates computed for each individual variable. The "total rate" can be plotted:
#  pTrend$plotSTrates("rate") # Equivalent to setting variable = 5

## ----plotTrend2,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotTrend2.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  # Perform search.trend setting Smilodontini and Pantherini as individual clades
#  cc<-2/parallel::detectCores()
#  search.trend(RRfel,PCscoresfel,node=c(129,154),clus=cc)->STclades
#  
#  # Visualize search.trend results
#  plotTrend(STclades)->pTrend2
#  # Using default options
#  par(mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
#  pTrend2$plotSTphenNode("PC1",node=1:2)
#  pTrend2$plotSTratesNode("PC1",node=c(154,129)) # This is the same as indicating node= 2:1

## ----plotTrend3,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotTrend3.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  ## Customizing parameters
#  par(mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
#  pTrend2$plotSTphenNode("PC2",node=1:2,
#                      plot.args = list(pch.node=c(23,24),pch=1,col="gray70",cex=1.2),
#                      lineTree.args = list(col="black",lwd=3),lineNode.args = list(lwd=5),
#                      node.palette = c("orangered","chartreuse"))
#  pTrend2$plotSTratesNode("rate",node=c(154,129),
#                      plot.args = list(pch.node=c(5,6),pch=16,col="gray70",cex=1.2,lwd=2),
#                      lineTree.args = list(col="gold",lwd=3,lty=4),
#                      lineNode.args = list(lwd=5),
#                      node.palette = c("deeppink","cyan2"))
#  

## ----plotTrend4,echo=FALSE,message=FALSE,warning=FALSE,out.width="98%",dpi=300,fig.show="hold"----
knitr::include_graphics("figs/plotTrend4.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ### load the RRphylo example dataset including Ornithodirans tree, body mass and locomotory type
#  data("DataOrnithodirans")
#  DataOrnithodirans$treedino->treedino
#  DataOrnithodirans$massdino->massdino
#  DataOrnithodirans$statedino->statedino
#  
#  # Perform RRphylo and search.shift under status.type="clade"
#  RRphylo(tree=treedino,y=massdino)->dinoRates
#  search.shift(RR=dinoRates,status.type="clade")->SSauto
#  
#  # Visualize search.shift results
#  plotShift(RR=dinoRates,SS=SSauto)->plotSS
#  # Using default options
#  plotSS$plotClades()

## ----plotShift1,echo=FALSE,message=FALSE,warning=FALSE,out.width="60%",fig.align='center',dpi=300----
knitr::include_graphics("figs/plotShift1.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  ## Customizing parameters
#  # Setting different colors for positive and negative shift
#  plotSS$plotClades(tree.args=list(no.margin=TRUE,type="fan"),
#                    symbols.args=list(lwd=2,fg=c(pos="gold",neg="green"),
#                                      bg=scales::alpha("grey70",0.3)))
#  
#  # Setting different colors for each shifting clade
#  plotSS$plotClades(tree.args=list(no.margin=TRUE),
#                    symbols.args=list(lwd=2,fg=NA,bg=scales::alpha(c("red","green","blue"),0.3)))

## ----plotShift2,echo=FALSE,message=FALSE,warning=FALSE,out.width="60%",fig.align='center',dpi=300----
knitr::include_graphics("figs/plotShift2.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  # Perform RRphylo and search.shift under status.type="sparse"
#  search.shift(RR=dinoRates,status.type= "sparse",state=statedino)->SSstate
#  
#  # Visualize search.shift results
#  plotShift(RR=dinoRates,SS=SSstate,state=statedino)->plotSS2
#  
#  # Using default options
#  plotSS2$plotStates()

## ----plotShift3,echo=FALSE,message=FALSE,warning=FALSE,out.width="60%",fig.align='center',dpi=600----
knitr::include_graphics("figs/plotShift3.png")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  ## Customizing parameters
#  # Setting customized colors and pch for different states and suppressing the legend
#  plotSS2$plotStates(plot.args=list(pch=2,main="Dinosaur stances"),
#                    points.args=list(bg=c("gold","forestgreen","royalblue","white"),
#                                     col=c("black","black","black","orangered"),
#                                     pch=c(21,22,24,13)),legend.args=NULL)

## ----plotShift4,echo=FALSE,message=FALSE,warning=FALSE,out.width="60%",fig.align='center',dpi=300----
knitr::include_graphics("figs/plotShift4.png")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  # Setting customized colors and pch for different states and not suppressing the legend
#  plotSS2$plotStates(plot.args=list(pch=2,main="Dinosaur stances"),
#                    points.args=list(bg=c("gold","forestgreen","royalblue","white"),
#                                     col=c("black","black","black","orangered"),
#                                     pch=c(21,22,24,13)),legend.args=list())
#  
#  

## ----plotShift5,echo=FALSE,message=FALSE,warning=FALSE,out.width="60%",fig.align='center',dpi=300----
knitr::include_graphics("figs/plotShift5.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ### load the RRphylo example dataset including Felids tree, 4 mandible shape variables (PCs),
#  ### and a category indicating whether or not the species shows sabertooth morphology
#  data("DataFelids")
#  DataFelids$PCscoresfel->PCscoresfel
#  DataFelids$treefel->treefel
#  DataFelids$statefel->statefel
#  
#  
#  # Perform RRphylo and search.conv between a pair of clades
#  cc<-2/parallel::detectCores()
#  RRphylo(treefel,PCscoresfel,clus=cc)->RRfel
#  search.conv(RR=RRfel, y=PCscoresfel, nodes=c(85,155),clus=cc)->sc_clade
#  
#  ## Visualize search.conv results
#  # Set variable = 1 to see results for the pair "85/155" (the first element in
#  # sc_clade$`average distance from group centroids`)
#  plotConv(SC=sc_clade, y=PCscoresfel, variable=1, RR = RRfel)->plotSC
#  
#  # Using default options
#  par(mfrow=c(1,2))
#  plotSC$plotHistTips()
#  plotSC$plotHistAces()

## ----plotConv1,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv1.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  plotSC$plotPChull()
#  plotSC$plotTraitgram()

## ----plotConv2,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv2.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  ## Customizing parameters
#  # Set variable = 2 to see results for the pair "86/156" (the second element in
#  # sc_clade$`average distance from group centroids`)
#  plotConv(SC=sc_clade, y=PCscoresfel, variable=2, RR = RRfel)->plotSC
#  
#  par(mfrow=c(1,2))
#  plotSC$plotHistTips(hist.args = list(col="gray80",yaxt="n",cex.axis=0.8,cex.main=1.5),
#                  line.args = list(lwd=3,lty=4,col="purple"))
#  plotSC$plotHistAces(hist.args = list(col="gray80",cex.axis=0.8,cex.main=1.5),
#                  line.args = list(lwd=3,lty=4,col="gold"))

## ----plotConv3,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv3.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  plotSC$plotPChull(chull.args = list(border=c("cyan","magenta"),lty=1),
#                means.args = list(pch=c(21,22),cex=3,bg=c("cyan2","magenta2")),
#                ace.args=list(pch=c(7,10)),legend.args = list(pch=c(24,11),bty="o",x="top"))
#  plotSC$plotTraitgram(colTree = "gray70",colNodes = c("cyan","magenta"),yaxt="s")
#  

## ----plotConv4,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv4.png")

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  # Perform search.conv within "saber" category
#  search.conv(tree=treefel, y=PCscoresfel, state=statefel,declust=TRUE,clus=cc)->sc_state
#  
#  ## Visualize search.conv results
#  # variable = 1 must be indicated also when a single row is output
#  plotConv(SC=sc_state, y=PCscoresfel, variable=1, state=statefel)->plotSC_state
#  
#  # Using default options
#  plotSC_state$plotPChull()
#  plotSC_state$plotPolar()

## ----plotConv5,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv5.png")

## ----eval=FALSE,message=FALSE,warning=FALSE,results=FALSE---------------------
#  ## Customizing parameters
#  plotSC_state$plotPChull(chull.args=list(border=c("gold2","blue"),lty=3),
#               points.args=list(pch=c(23,21),bg="gray"),
#               legend.args=list(pch=c(23,21),x="top"))
#  
#  plotSC_state$plotPolar(polar.args=list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
#              polygon.args=list(line.col="green",poly.col=NA,lwd=2),
#              line.args=list(line.col="deeppink",lty=2,lwd=3))
#  

## ----plotConv6,message=FALSE,warning=FALSE,echo=FALSE,out.width="98%",dpi=300----
knitr::include_graphics("figs/plotConv6.png")

