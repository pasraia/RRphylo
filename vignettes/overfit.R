## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=FALSE
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  ## overfitRR routine
#  # load the RRphylo example dataset including Ornithodirans tree and data
#  library(ape)
#  data("DataOrnithodirans")
#  DataOrnithodirans$treedino->treedino
#  log(DataOrnithodirans$massdino)->massdino
#  DataOrnithodirans$statedino->statedino
#  
#  # extract Pterosaurs tree and data
#  extract.clade(treedino,746)->treeptero
#  massdino[match(treeptero$tip.label,names(massdino))]->massptero
#  massptero[match(treeptero$tip.label,names(massptero))]->massptero
#  
#  # peform RRphylo on body mass
#  RRphylo(tree=treeptero,y=massptero,clus=cc)->RRptero
#  
#  # generate a list of subsampled and swapped phylogenies to test
#  tree.list<-resampleTree(RRptero$tree,s = 0.25,swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#  
#  # test the robustness of RRphylo
#  ofRRptero<-overfitRR(RR = RRptero,y=massptero,phylo.list=tree.list,clus=cc)
#  
#  
#  ## overfitRR routine on multiple RRphylo
#  # load the RRphylo example dataset including Cetaceans tree and data
#  data("DataCetaceans")
#  DataCetaceans$treecet->treecet
#  DataCetaceans$masscet->masscet
#  DataCetaceans$brainmasscet->brainmasscet
#  DataCetaceans$aceMyst->aceMyst
#  
#  # cross-reference the phylogenetic tree and body and brain mass data. Remove from
#  # both the tree and vector of body sizes the species whose brain size is missing
#  drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#  masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#  
#  # peform RRphylo on the variable (body mass) to be used as additional predictor
#  RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#  RRmass.multi$aces[,1]->acemass.multi
#  
#  # create the predictor vector: retrieve the ancestral character estimates
#  # of body size at internal nodes from the RR object ($aces) and collate them
#  # to the vector of species' body sizes to create
#  c(acemass.multi,masscet.multi)->x1.mass
#  
#  # peform RRphylo on brain mass by using body mass as additional predictor
#  RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#  
#  # generate a list of subsampled and swapped phylogenies to test
#  treecet.list<-resampleTree(RRmulti$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#  
#  # test the robustness of multiple RRphylo
#  ofRRcet<-overfitRR(RR = RRmulti,y=brainmasscet,phylo.list=treecet.list,clus=cc,x1 =x1.mass)
#  

## -----------------------------------------------------------------------------
#  
#  # load the RRphylo example dataset including Ornithodirans tree and data
#  data("DataOrnithodirans")
#  DataOrnithodirans$treedino->treedino
#  log(DataOrnithodirans$massdino)->massdino
#  DataOrnithodirans$statedino->statedino
#  
#  # peform RRphylo on Ornithodirans tree and data
#  RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#  
#  # perform search.shift under both "clade" and "sparse" condition
#  search.shift(RR=dinoRates, status.type= "clade")->SSauto
#  search.shift(RR=dinoRates, status.type= "sparse", state=statedino)->SSstate
#  
#  ## overfitSS routine
#  # generate a list of subsampled and swapped phylogenies, setting as categories/node
#  # the state/node under testing
#  tree.list<-resampleTree(dinoRates$tree,s = 0.25,categories=statedino,
#  node=rownames(SSauto$single.clades),swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#  
#  # test the robustness of search.shift
#  ofRRdino<-overfitRR(RR = dinoRates,y=massdino,phylo.list=tree.list,clus=cc)
#  ofSS<-overfitSS(RR = dinoRates,oveRR = ofRRdino,state=statedino,
#                  node=rownames(SSauto$single.clades))
#  

## -----------------------------------------------------------------------------
#  library(ape)
#  
#  ## Case 1
#  # load the RRphylo example dataset including Ornithodirans tree and data
#  data("DataOrnithodirans")
#  DataOrnithodirans$treedino->treedino
#  log(DataOrnithodirans$massdino)->massdino
#  DataOrnithodirans$statedino->statedino
#  
#  # extract Pterosaurs tree and data
#  extract.clade(treedino,746)->treeptero
#  massdino[match(treeptero$tip.label,names(massdino))]->massptero
#  massptero[match(treeptero$tip.label,names(massptero))]->massptero
#  
#  # perform RRphylo and search.trend on body mass data
#  RRphylo(tree=treeptero,y=massptero,clus=cc)->RRptero
#  search.trend(RR=RRptero, y=massptero,node=143,clus=cc)->ST
#  
#  ## overfitST routine
#  # generate a list of subsampled and swapped phylogenies setting as node
#  # the clade under testing
#  treeptero.list<-resampleTree(RRptero$tree,s = 0.25,node=143,
#                               swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#  
#  # test the robustness of search.trend
#  ofRRptero<-overfitRR(RR = RRptero,y=massptero,phylo.list=treeptero.list,clus=cc)
#  ofSTptero<-overfitST(RR=RRptero,y=massptero,oveRR=ofRRptero,node=143,clus=cc)
#  
#  
#  ## Case 2
#  # load the RRphylo example dataset including Cetaceans tree and data
#  data("DataCetaceans")
#  DataCetaceans$treecet->treecet
#  DataCetaceans$masscet->masscet
#  DataCetaceans$brainmasscet->brainmasscet
#  
#  # cross-reference the phylogenetic tree and body and brain mass data. Remove from
#  # both the tree and vector of body sizes the species whose brain size is missing
#  drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),
#                                                 treecet$tip.label)])->treecet.multi
#  masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#  
#  # peform RRphylo on the variable (body mass) to be used as additional predictor
#  RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#  RRmass.multi$aces[,1]->acemass.multi
#  
#  # create the predictor vector: retrieve the ancestral character estimates
#  # of body size at internal nodes from the RR object ($aces) and collate them
#  # to the vector of species' body sizes to create
#  c(acemass.multi,masscet.multi)->x1.mass
#  
#  # peform RRphylo and search.trend on brain mass by using body mass as additional predictor
#  RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc)->STcet
#  
#  ## overfitST routine
#  # generate a list of subsampled and swapped phylogenies to test
#  treecet.list<-resampleTree(RRmulti$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#  
#  # test the robustness of search.trend with and without x1.residuals
#  ofRRcet<-overfitRR(RR = RRmulti,y=brainmasscet,phylo.list=treecet.list,clus=cc,x1 =x1.mass)
#  ofSTcet1<-overfitST(RR=RRmulti,y=brainmasscet,oveRR=ofRRcet,x1 =x1.mass,clus=cc)
#  ofSTcet2<-overfitST(RR=RRmulti,y=brainmasscet,oveRR=ofRRcet,x1 =x1.mass,x1.residuals = TRUE,clus=cc)
#  

## -----------------------------------------------------------------------------
#  library(Morpho)
#  
#  ## Testing search.conv
#  # load the RRphylo example dataset including Felids tree and data
#  data("DataFelids")
#  DataFelids$treefel->treefel
#  DataFelids$statefel->statefel
#  DataFelids$landfel->feldata
#  
#  # perform data reduction via Procrustes superimposition (in this case)
#  procSym(feldata)->pcafel
#  pcafel$PCscores->PCscoresfel
#  
#  # perform RRphylo on Felids tree and data
#  RRphylo(treefel,PCscoresfel)->RRfel
#  
#  # search for morphologicl convergence between clades (automatic mode)
#  # and within the category
#  search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="time38")->sc1
#  search.conv(tree=treefel, y=PCscoresfel, state=statefel, declust=TRUE)->sc2
#  
#  # select converging clades returned in sc1
#  felnods<-c(85,155)
#  
#  ## overfitSC routine
#  
#  # generate a list of subsampled and swapped phylogenies to test for search.conv
#  # robustness. Use as reference tree the phylogeny returned by RRphylo.
#  # Set the nodes and the categories under testing as arguments of
#  # resampleTree so that it maintains no less than 5 species in each clade/state.
#  tree.list<-resampleTree(RRfel$tree,s=0.15,nodes=felnods,categories=statefel,
#                          nsim=15,swap.si=0.1,swap.si2=0.1)
#  
#  # match the original data with each subsampled-swapped phylogeny in tree.list
#  # and repeat data reduction
#  y.list<-lapply(tree.list,function(k){
#    treedataMatch(k,feldata)[[1]]->ynew
#    procSym(ynew)$PCscores
#  })
#  
#  # test for robustness of search.conv results by overfitSC
#  oSC<-overfitSC(RR=RRfel,phylo.list=tree.list,y.list=y.list,
#                 nodes = felnods,state=statefel)
#  

## -----------------------------------------------------------------------------
#  library(phytools)
#  library(ape)
#  
#  # generate fictional data to test the function
#  rtree(100)->tree
#  fastBM(tree)->resp
#  fastBM(tree,nsim=3)->resp.multi
#  fastBM(tree)->pred1
#  fastBM(tree)->pred2
#  data.frame(y1=resp,x2=pred1,x1=pred2)->dato
#  
#  # perform RRphylo and PGLS_fossil with univariate/multivariate phenotypic data
#  PGLS_fossil(modform=y1~x1+x2,data=dato,tree=tree)->pgls_noRR
#  RRphylo(tree,resp,clus=cc)->RR
#  PGLS_fossil(modform=resp~pred1+pred2,RR=RR)->pgls_RR
#  
#  PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree)->pgls2_noRR
#  RRphylo(tree,resp.multi,clus=cc)->RR2
#  PGLS_fossil(modform=resp.multi~pred1+pred2,tree=tree,RR=RR)->pgls2_RR
#  
#  ## overfitPGLS routine
#  # generate a list of subsampled and swapped phylogenies to test
#  tree.list<-resampleTree(RR$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#  
#  # test the robustnes of PGLS_fossil with univariate/multivariate phenotypic data
#  ofRR<-overfitRR(RR = RR,y=resp,phylo.list=tree.list,clus=cc)
#  ofPGLS<-overfitPGLS(oveRR = ofRR,phylo.list=tree.list,modform = y1~x1+x2,data=dato)
#  
#  ofRR2<-overfitRR(RR = RR2,y=resp.multi,phylo.list=tree.list,clus=cc)
#  ofPGLS2<-overfitPGLS(oveRR = ofRR2,phylo.list=tree.list,modform = y1~x1+x2,
#                       data=list(y1=resp.multi,x2=pred1,x1=pred2))

