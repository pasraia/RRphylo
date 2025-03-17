## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

if (!requireNamespace("kableExtra", quietly = TRUE)) {
   warning(call. = FALSE, "kableExtra not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% t(L)) %*% as.matrix(y)

## ----eval = FALSE-------------------------------------------------------------
#  y.pred <- L %*% betas
#  ace <- L1 %*% betas[1:Nnode(tree),]

## ----message = FALSE,fig.show='asis',fig.align="center",fig.dim=c(4,4),echo=c(1:8),out.width='60%',dpi=220----
require(ape)
set.seed(76)
rtree(5)->tree
makeL(tree)->L
makeL1(tree)->L1

plot(tree,no.margin=TRUE,edge.width=1.7)
nodelabels(bg="w",frame="n",col="red",adj=-0.01)

require(kableExtra)

knitr::kable(L1,digits=3,align="c", caption="L1 matrix") %>%
  kable_styling(full_width = F, position = "float_left")  %>%
  column_spec(1, bold = T)

knitr::kable(L,digits=3,align="c", caption="L matrix") %>%
  kable_styling(full_width = F, position = "right")  %>%
  column_spec(1, bold = T)


## ----fig.width=4, fig.height=4,fig.align="center",message = FALSE,out.width='60%',dpi=220----
require(phytools)

rtree(100)->tree
makeL(tree)->L

# produce a phenotypic vector for both tips and nodes
fastBM(tree,internal=T)->phen
phen[1:100]->y
phen[101:199]->acey

# set λ close to 0
lambda <-1e-10

# compute evolutionary rates and estimate phenotypic values at tips as within RRphylo
betas<-(solve(t(L) %*% L +  lambda *diag(ncol(L))) %*% t(L)) %*% as.matrix(y)
y.pred <- L %*% betas

plot(y,y.pred,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="simulated",ylab="estimated",main="phenotype at tips")


## ----fig.show="hold",fig.width=8, fig.height=4,fig.align='center',out.width='100%',dpi=220----
par(mfrow=c(1,2),mar=c(4,3,3,1))

plot(c(acey,y),betas,bg=c(rep("red",99),rep("green",100)),pch=21,cex=1.5,mgp=c(1.8,0.5,0),
     xlab="phenotypes at nodes and tips",ylab="rates",main="rates vs simulated phenotypes")
legend("topright",legend=c("nodes","tips"),fill=c("red","green"),bty="n")

# compute ancestral character estimates at nodes as within RRphylo
makeL1(tree)->L1
ace <- L1 %*% betas[1:Nnode(tree),]
plot(acey,ace,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="simulated",ylab="estimated",main="phenotype at nodes (aces)")


## ----fig.show="hold",fig.width=8, fig.height=4,fig.align='center',out.width='100%',dpi=220----
cc<-2/parallel::detectCores()
RRphylo(tree,y,clus=cc)->RR
RR$rates[,1]->betas
RR$predicted.phenotype[,1]->y.pred
RR$aces[,1]->ace

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(y,y.pred,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="simulated",ylab="estimated by RRphylo",main="phenotype at tips")
plot(acey,ace,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="simulated",ylab="estimated by RRphylo",main="phenotype at nodes (aces)")

## ----fig.width=8, fig.height=4,fig.align='center',out.width='100%',dpi=220----
set.seed(76)
rtree(100)->tree
fastBM(tree,sig2=2)->y

# perform RRphylo on the covariate to retrieve ancestral character values 
# and collate them to tip values to create the covariate vector
RRphylo(tree,y,clus=cc)->RR
RR$rates[,1]->betas
c(RR$aces[,1],y)->cov

# within RRphylo the logged absolute rates are regressed against 
# the absolute covariate values, and regression residuals are used 
# as rate values
R <- log(abs(betas))
Y <- abs(cov)
res <- residuals(lm(R ~ Y))
betas <- as.matrix(res)

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(Y,R,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="absolute covariate values",ylab="log absolute rates")
abline(lm(R~Y),lwd=2,col="red")
plot(Y,res,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="absolute covariate values",ylab="regression residuals")
abline(lm(res~Y),lwd=2,col="red")


## ----message = FALSE,fig.show='asis',fig.align="center",fig.dim=c(4,4),out.width='60%',echo=FALSE,dpi=220----
require(ape)
set.seed(76)
rtree(5)->tree
fastBM(tree)->pred
RRphylo(tree,pred,clus=cc)->RRpred

makeL(tree)->L
makeL1(tree)->L1
cbind(L,pred=pred)->LL
cbind(L1,pred=RRpred$aces[,1])->LL1

plot(tree,no.margin=TRUE,edge.width=1.7)
nodelabels(bg="w",frame="n",col="red",adj=-0.01)

require(kableExtra)

knitr::kable(LL1,digits=3,align="c", caption="L1' matrix") %>%
  kable_styling(full_width = F, position = "float_left")  %>%
  column_spec(1, bold = T) %>%
  column_spec(6, background = "#FE7979") 

knitr::kable(LL,digits=3,align="c", caption="L' matrix") %>%
  kable_styling(full_width = F, position = "right")  %>%
  column_spec(1, bold = T) %>%
  column_spec(11, background = "#FE7979")


## ----message = FALSE,echo=FALSE,out.width='47%',dpi=220-----------------------
require(ape)
require(phytools)
set.seed(76)
rtree(8)->tree
tree$node.label<-9:15

N <- c(12,14)
tar.tips <- lapply(N, function(x) tips(tree, x))
names(tar.tips) <- N
treeN <- tree

i = 1
while (i <= length(N)) {
  nn <- getMRCA(treeN, tar.tips[[i]])
  treeN <- bind.tip(treeN,
                    tip.label = paste("nod",N[i], sep = ""),
                    edge.length = 0, where = nn,position = 0.001)
  i = i + 1
}


plot(tree,no.margin=TRUE,edge.width=1.7)
nodelabels(bg="w",frame="n",col="red",adj=-0.01)

edge.col<-rep("black",nrow(treeN$edge))
edge.col[which(treeN$edge[,2]%in%grep("nod",treeN$tip.label))]<-"green"
gsub("NA","",treeN$node.label)->treeN$node.label

plot(treeN,no.margin=TRUE,edge.color = edge.col,edge.width=1.7)
nodelabels(text=na.omit(treeN$node.label),bg="w",frame="n",col="red",adj=-0.01)


## ----out.width='100%',fig.dim=c(8,7),message=FALSE,dpi=220,warning=FALSE------
# load the RRphylo example dataset including Cetaceans tree and data
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree
DataCetaceans$masscet->masscet # logged body mass data
DataCetaceans$brainmasscet->brainmasscet # logged brain mass data
DataCetaceans$aceMyst->aceMyst # known phenotypic value for the most recent 
                               # common ancestor of Mysticeti

par(mar=c(0,0,0,1))
plot(ladderize(treecet,FALSE),show.tip.label = FALSE,edge.color = "gray40",edge.width = 1.5)
plotinfo<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)
nodelabels(text="",node=128,frame="circle",bg="red",cex=0.5)
nodelabels(text="Mystacodon",node=128,frame="n",bg="w",cex=1,adj=c(-0.1,0.5),font=2)
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,128))])->yran128
rect(plotinfo$x.lim[2]+0.4,yran128[1],plotinfo$x.lim[2]+0.7,yran128[2],col="red",border="red")
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,142))])->yran142
rect(plotinfo$x.lim[2]+0.4,yran142[1],plotinfo$x.lim[2]+0.7,yran142[2],col="blue",border="blue")
mtext(c("Mysticeti","Odontoceti"), side = 4,line=-0.5,at=c(sum(yran128)/2,sum(yran142)/2),
      cex=1.5,adj=0.5,col=c("red","blue"))

## ----message=FALSE,warning=FALSE,eval=FALSE-----------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  masscet[match(treecet$tip.label,names(masscet))]->masscet
#  
#  
#  ## Example 1: RRphylo by accounting for the effect of a coviariate
#  # perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treecet,y=masscet)->RRmasscet
#  
#  # create the covariate vector: extract phenotypic character (i.e. log body mass)
#  # estimates at nodes from the RR object ($aces) and collate them
#  # to the vector of (log) body mass
#  c(RRmasscet$aces[,1],masscet)->masscov
#  
#  # perform RRphylo on the vector of (log) body mass by including
#  # the body mass itslef as covariate
#  RRphylo(tree=treecet,y=masscet,cov=masscov)->RRcov
#  
#  
#  ## Example 2: RRphylo by setting values at internal nodes
#  # Set the body mass of Mysticetes ancestor (Mystacodon selenensis)
#  # as known value at node
#  RRphylo(tree=treecet,y=masscet,aces=aceMyst)->RRace
#  
#  
#  ## Example 3: multiple regression version of RRphylo
#  # cross-reference the phylogenetic tree and body and brain mass data.
#  # Remove from both the tree and vector of body sizes the species
#  # whose brain size is missing
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
#  
#  ## Example 4: categorical and multiple regression version of RRphylo with
#  ## 2 additional predictors performed by setting values at internal nodes
#  require(phytools)
#  set.seed(1458)
#  
#  # generate a random tree and a BM phenotypic vector on it
#  rtree(50)->tree
#  fastBM(tree)->y
#  
#  # produce two variables to be used as additional predictors into the multiple
#  # regression version of  RRphylo. One variable is continuous, the other is discrete.
#  jitter(y)*10->y1
#  rep(1,length(y))->y2
#  y2[sample(1:50,20)]<-2
#  names(y2)<-names(y)
#  
#  # perform RRphylo on y1 and y2 to retrieve ancestral state estimates at nodes
#  # and create the x1 matrix
#  c(RRphylo(tree,y1)$aces[,1],y1)->x1
#  RRphylo(tree,y2)->RRcat # this is categorical RRphylo
#  c(RRcat$aces[,1],y2)->x2
#  cbind(x1,x2)->x1mat
#  
#  # create the phenotypes for y, y1, and y2 to be set as known values at internal nodes
#  cbind(c(jitter(mean(y1[tips(tree,83)])),1),
#       c(jitter(mean(y1[tips(tree,53)])),2))->acex
#  c(jitter(mean(y[tips(tree,83)])),jitter(mean(y[tips(tree,53)])))->acesy
#  names(acesy)<-rownames(acex)<-c(83,53)
#  
#  # perform RRphylo by specifying aces, x1, and aces.x1 arguments
#  RRphylo(tree,y,aces=acesy,x1=x1mat,aces.x1 = acex)->RRcat.multi
#  

