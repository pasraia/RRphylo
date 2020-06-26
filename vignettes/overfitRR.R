## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,echo=FALSE,  message = FALSE---------------------------------------
require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  library(ape)
#  
#  # load the RRphylo example dataset including Ornithodirans tree and data
#  DataOrnithodirans$treedino->treedino # phylogenetic tree
#  DataOrnithodirans$massdino->massdino # body mass data
#  DataOrnithodirans$statedino->statedino # locomotory type data
#  
#  ### Testing search.shift
#  # perform RRphylo Ornithodirans tree and data
#  RRphylo(tree=treedino,y=massdino)->dinoRates
#  
#  # perform search.shift under both "clade" and "sparse" condition
#  search.shift(RR=dinoRates, status.type= "clade",foldername=tempdir())->SSnode
#  search.shift(RR=dinoRates, status.type= "sparse", state=statedino,
#               foldername=tempdir())->SSstate
#  
#  # test the robustness of search.shift results
#  overfitRR(RR=dinoRates,y=massdino,swap.args =list(si=0.2,si2=0.2),
#            shift.args = list(node=rownames(SSnode$single.clades),state=statedino),
#            nsim=10)
#  
#  
#  ### Testing search.trend
#  # Extract Pterosaurs tree and data
#  extract.clade(treedino,748)->treeptero # phylogenetic tree
#  massdino[match(treeptero$tip.label,names(massdino))]->massptero # body mass data
#  massptero[match(treeptero$tip.label,names(massptero))]->massptero
#  
#  # perform RRphylo and search.trend on Pterosaurs tree and data
#  # by specifying a clade to be tested
#  RRphylo(tree=treeptero,y=log(massptero))->RRptero
#  
#  search.trend(RR=RRptero, y=log(massptero),node=143,foldername=tempdir(),
#               cov=NULL,ConfInt=FALSE)->STnode
#  
#  # test the robustness of search.trend results
#  overfitRR(RR=RRptero,y=log(massptero),trend.args = list(node=143),nsim=10)
#  
#  ### Applying overfitRR on multiple RRphylo
#  # load the RRphylo example dataset including Cetaceans tree and data
#  data("DataCetaceans")
#  DataCetaceans$treecet->treecet # phylogenetic tree
#  DataCetaceans$masscet->masscet # logged body mass data
#  DataCetaceans$brainmasscet->brainmasscet # logged brain mass data
#  DataCetaceans$aceMyst->aceMyst # known phenotypic value for the most recent common ancestor of Mysticeti
#  
#  # cross-reference the phylogenetic tree and body and brain mass data. Remove from both the tree and
#  # vector of body sizes the species whose brain size is missing
#  drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),
#                                                 treecet$tip.label)])->treecet.multi
#  masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#  
#  # peform RRphylo on the variable (body mass) to be used as additional predictor
#  RRphylo(tree=treecet.multi,y=masscet.multi)->RRmass.multi
#  RRmass.multi$aces[,1]->acemass.multi
#  
#  # create the predictor vector: retrieve the ancestral character estimates
#  # of body size at internal nodes from the RR object ($aces) and collate them
#  # to the vector of species' body sizes to create
#  c(acemass.multi,masscet.multi)->x1.mass
#  
#  # peform RRphylo and search.trend on the brain mass
#  # by using the body mass as additional predictor
#  RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#  
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,foldername=tempdir())->STcet
#  
#  # test the robustness of search.trend results
#  overfitRR(RR=RRmulti,y=brainmasscet,trend.args = list(),x1=x1.mass,nsim=10)
#  
#  
#  ### Testing search.conv
#  # load the RRphylo example dataset including Felids tree and data
#  data("DataFelids")
#  DataFelids$PCscoresfel->PCscoresfel # mandible shape data
#  DataFelids$treefel->treefel # phylogenetic tree
#  DataFelids$statefel->statefel # conical-toothed or saber-toothed condition
#  
#  # perform RRphylo on Felids tree and data
#  RRphylo(tree=treefel,y=PCscoresfel)->RRfel
#  
#  # search for morphologicl convergence between clades (automatic mode) and within the category
#  search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",
#              foldername = tempdir())->SC.clade
#  as.numeric(c(rownames(SC.clade[[1]])[1],as.numeric(as.character(SC.clade[[1]][1,1]))))->conv.nodes
#  
#  search.conv(tree=treefel, y=PCscoresfel, state=statefel,
#              foldername = tempdir())->SC.state
#  
#  # test the robustness of seach.conv results
#  overfitRR(RR=RRfel, y=PCscoresfel,conv.args=
#              list(node=conv.nodes,state=statefel,declust=TRUE),nsim=10)

