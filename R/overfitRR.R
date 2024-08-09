#' @title Testing RRphylo methods overfit
#' @description Testing the robustness of \code{\link{search.trend}}
#'   (\cite{Castiglione et al. 2019a}), \code{\link{search.shift}}
#'   (\cite{Castiglione et al. 2018}),  \code{\link{search.conv}}
#'   (\cite{Castiglione et al. 2019b}), and \code{\link{PGLS_fossil}} results to
#'   sampling effects and phylogenetic uncertainty.
#' @usage
#' overfitRR(RR,y,phylo.list=NULL,s=0.25,swap.args=NULL,trend.args=NULL,shift.args=NULL,
#' conv.args=NULL, pgls.args=NULL,aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,
#' rootV=NULL,nsim=100,clus=0.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a named vector of phenotypes.
#' @param phylo.list a list (or multiPhylo) of alternative phylogenies to be
#'   tested.
#' @param s the percentage of tips to be cut off. It is set at 25\% by default.
#'   If \code{phylo.list} is provided, this argument is ignored.
#' @param swap.args a list of arguments to be passed to the function
#'   \code{\link{swapONE}}, including \code{list(si=NULL,si2=NULL,}
#'   \code{node=NULL)}. If \code{swap.arg} is unspecified, the function
#'   automatically sets both \code{si} and \code{si2} to 0.1. If
#'   \code{phylo.list} is provided, swapping is not performed.
#' @param trend.args a list of arguments specific to the function
#'   \code{search.trend}, including \code{list(node=NULL,x1.residuals=FALSE)}.
#'   If a trend for the whole tree is to be tested, type \code{trend.args =
#'   list()}. No trend is tested if left unspecified.
#' @param shift.args a list of arguments specific to the function
#'   \code{search.shift}, including \code{list(node=NULL,} \code{state=NULL)}.
#'   Arguments \code{node} and \code{state} can be specified at the same time.
#' @param conv.args a list of arguments specific to the function
#'   \code{search.conv}, including \code{list(node=NULL,} \code{state=NULL,
#'   declust=FALSE)}. Arguments \code{node} and \code{state} can be specified at
#'   the same time.
#' @param pgls.args a list of arguments specific to the function
#'   \code{PGLS_fossil}, including \code{list(modform,} \code{data,
#'   tree=FALSE,RR=TRUE,...)}. If \code{tree=TRUE}, \code{PGLS_fossil} is
#'   performed by using the RRphylo output tree as \code{tree} argument. If
#'   \code{RR=TRUE}, \code{PGLS_fossil} is performed by using the RRphylo output
#'   as \code{RR} argument. Arguments \code{tree} and \code{RR} can be
#'   \code{TRUE} at the same time. \code{...} are further argument passed to
#'   \code{PGLS_fossil}.
#' @param aces if used to produce the \code{RR} object, the vector of those
#'   ancestral character values at nodes known in advance must be specified.
#'   Names correspond to the nodes in the tree.
#' @param x1 the additional predictor to be specified if the RR object has been
#'   created using an additional predictor (i.e. multiple version of
#'   \code{RRphylo}). \code{'x1'} vector must be as long as the number of nodes
#'   plus the number of tips of the tree, which can be obtained by running
#'   \code{RRphylo} on the predictor as well, and taking the vector of ancestral
#'   states and tip values to form the \code{x1}.
#' @param aces.x1 a named vector of ancestral character values at nodes for
#'   \code{x1}. It must be indicated if the RR object has been created using
#'   both \code{aces} and \code{x1}. Names correspond to the nodes in the tree.
#' @param cov if used to produce the \code{RR} object, the covariate must be
#'   specified. As in \code{RRphylo}, the covariate vector must be as long as
#'   the number of nodes plus the number of tips of the tree, which can be
#'   obtained by running \code{RRphylo} on the covariate as well, and taking the
#'   vector of ancestral states and tip values to form the covariate.
#' @param rootV if used to produce the \code{RR} object, the phenotypic value at
#'   the tree root must be specified.
#' @param nsim number of simulations to be performed. It is set at 100 by
#'   default.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{overfitRR} set \code{clus} = 0.
#' @return The function returns a 'RRphyloList' object containing:
#' @return \strong{$mean.sampling} the mean proportion of species actually
#'   removed from the tree over the iterations.
#' @return \strong{$tree.list} a 'multiPhylo' list including the trees generated
#'   within \code{overfitRR}
#' @return \strong{$RR.list} a 'RRphyloList' including the results of each
#'   \code{RRphylo} performed within \code{overfitRR}
#' @return \strong{$rootCI} the 95\% confidence interval around the root value.
#' @return \strong{$ace.regressions} a 'RRphyloList' including the results of
#'   linear regression between ancestral state estimates before and after the
#'   subsampling.
#' @return \strong{$conv.results} a list including results for
#'   \code{search.conv} performed under \code{clade} and \code{state}
#'   conditions. If a node pair is specified within \code{conv.args}, the
#'   \code{$clade} object contains the percentage of simulations producing
#'   significant p-values for convergence between the clades, and the proportion
#'   of tested trees (i.e. where the clades identity was preserved; always 1 if
#'   no \code{phylo.list} is supplied). If a state vector is supplied within
#'   \code{conv.args}, the object \code{$state} contains the percentage of
#'   simulations producing significant p-values for convergence within (single
#'   state) or between states (multiple states).
#' @return \strong{$shift.results} a list including results for
#'   \code{search.shift} performed under \code{clade} and \code{sparse}
#'   conditions. If one or more nodes are specified within \code{shift.args},
#'   the \code{$clade} object contains for each node the percentage of
#'   simulations producing significant p-value separated by shift sign, and the
#'   same figures by considering all the specified nodes as evolving under a
#'   single rate (all.clades). For each node the proportion of tested trees
#'   (i.e. where the clade identity was preserved; always 1 if no
#'   \code{phylo.list} is supplied) is also indicated. If a state vector is
#'   supplied within \code{shift.args}, the object \code{$sparse} contains the
#'   percentage of simulations producing significant p-value separated by shift
#'   sign ($p.states).
#' @return \strong{$trend.results} a list including the percentage of
#'   simulations showing significant p-values for phenotypes versus age and
#'   absolute rates versus age regressions for the entire tree separated by
#'   slope sign ($tree). If one or more nodes are specified within
#'   \code{trend.args}, the list also includes the same results at nodes ($node)
#'   and the results for comparison between nodes ($comparison). For each node the proportion
#'   of tested trees (i.e. where the clade identity was preserved; always 1 if
#'   no \code{phylo.list} is supplied) is also indicated.
#' @return \strong{$pgls.results} two 'RRphyloList' objects including results of
#'   \code{PGLS_fossil} performed by using the phylogeny as it is (\code{$tree})
#'   or rescaled according to the \code{RRphylo} rates (\code{$RR}).
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @details Methods using a large number of parameters risk being overfit. This
#'   usually translates in poor fitting with data and trees other than the those
#'   originally used. With \code{RRphylo} methods this risk is usually very low.
#'   However, the user can assess how robust the results got by applying
#'   \code{search.shift}, \code{search.trend}, \code{search.conv} or
#'   \code{PGLS_fossil} are by running \code{overfitRR}. With the latter, the
#'   original tree and data are subsampled by specifying a \code{s} parameter,
#'   that is the proportion of tips to be removed from the tree. In some cases,
#'   though, removing as many tips as imposed by \code{s} would delete too many
#'   tips right in clades and/or states under testing. In these cases, the
#'   function maintains no less than 5 species at least in each clade/state
#'   under testing (or all species if there is less), reducing the sampling
#'   parameter \code{s} if necessary. Internally, \code{overfitRR} further
#'   shuffles the tree by using the function \code{\link{swapONE}}. Thereby,
#'   both the potential for overfit and phylogenetic uncertainty are accounted
#'   for straight away.
#'
#'   Otherwise, a list of alternative phylogenies can be supplied to
#'   \code{overfitRR}. In this case subsampling and swapping arguments are
#'   ignored, and robustness testing is performed on the alternative topologies
#'   as they are. If a clade has to be tested either in \code{search.shift},
#'   \code{search.trend}, or \code{search.conv}, the function scans each
#'   alternative topology searching for the corresponding clade. If the species
#'   within such clade on the alternative topology differ more than 10% from the
#'   species within the clade in the original tree, the identity of the clade is
#'   considered disrupted and the test is not performed.
#' @export
#' @seealso \href{../doc/overfitRR.html}{\code{overfitRR} vignette} ;
#'   \href{../doc/search.trend.html}{\code{search.trend} vignette} ;
#'   \href{../doc/search.shift.html}{\code{search.shift} vignette} ;
#'   \href{../doc/search.conv.html}{\code{search.conv} vignette} ;
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P. (2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M.,
#'   Profico, A., Girardi, G., & Raia, P. (2019a) Simultaneous detection of
#'   macroevolutionary patterns in phenotypic means and rate of change with and
#'   within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14:
#'   e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @references Castiglione, S., Serio, C., Tamagnini, D., Melchionna, M.,
#'   Mondanaro, A., Di Febbraro, M., Profico, A., Piras, P.,Barattolo, F., &
#'   Raia, P. (2019b). A new, fast method to search for morphological
#'   convergence with shape data. \emph{PLoS ONE}, 14, e0226949.
#'   https://doi.org/10.1371/journal.pone.0226949
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' cc<- 2/parallel::detectCores()
#'
#' # Extract Pterosaurs tree and data
#' library(ape)
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#'
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#' RRphylo(tree=treeptero,y=log(massptero),clus=cc)->RRptero
#'
#' # Case 1 search.shift under both "clade" and "sparse" condition
#' search.shift(RR=dinoRates, status.type= "clade")->SSnode
#' search.shift(RR=dinoRates, status.type= "sparse", state=statedino)->SSstate
#'
#' overfitRR(RR=dinoRates,y=massdino,swap.args =list(si=0.2,si2=0.2),
#'           shift.args = list(node=rownames(SSnode$single.clades),state=statedino),
#'           nsim=10,clus=cc)->orr.ss
#'
#' # Case 2 search.trend on the entire tree
#' search.trend(RR=RRptero, y=log(massptero),nsim=100,clus=cc,cov=NULL,node=NULL)->STtree
#'
#' overfitRR(RR=RRptero,y=log(massptero),swap.args =list(si=0.2,si2=0.2),
#'           trend.args = list(),nsim=10,clus=cc)->orr.st1
#'
#' # Case 3 search.trend at specified nodescov=NULL,
#' search.trend(RR=RRptero, y=log(massptero),node=143,clus=cc)->STnode
#'
#' overfitRR(RR=RRptero,y=log(massptero),
#'           trend.args = list(node=143),nsim=10,clus=cc)->orr.st2
#'
#' # Case 4 overfitRR on multiple RRphylo
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#' DataCetaceans$aceMyst->aceMyst
#'
#' ape::drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),
#'                                                treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc)->STcet
#' overfitRR(RR=RRmulti,y=brainmasscet,trend.args = list(),
#'           x1=x1.mass,nsim=10,clus=cc)->orr.st3
#'
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,
#'              clus=cc)->STcet.resi
#' overfitRR(RR=RRmulti,y=brainmasscet,trend.args = list(x1.residuals=TRUE),
#'           x1=x1.mass,nsim=10,clus=cc)->orr.st4
#'
#' # Case 5 searching convergence between clades and within a single state
#' data("DataFelids")
#' DataFelids$PCscoresfel->PCscoresfel
#' DataFelids$treefel->treefel
#' DataFelids$statefel->statefel
#'
#' RRphylo(tree=treefel,y=PCscoresfel,clus=cc)->RRfel
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",clus=cc)->SC.clade
#' as.numeric(c(rownames(SC.clade[[1]])[1],as.numeric(as.character(SC.clade[[1]][1,1]))))->conv.nodes
#'
#' overfitRR(RR=RRfel, y=PCscoresfel,conv.args =
#' list(node=conv.nodes,state=statefel,declust=TRUE),nsim=10,clus=cc)->orr.sc
#'
#' # Case 6 overfitRR on PGLS_fossil
#' library(phytools)
#' rtree(100)->tree
#' fastBM(tree)->resp
#' fastBM(tree,nsim=3)->resp.multi
#' fastBM(tree)->pred1
#' fastBM(tree)->pred2
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree)->pgls_noRR
#'
#' RRphylo(tree,resp,clus=cc)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree,RR=RR)->pgls_RR
#'
#' overfitRR(RR=RR,y=resp,
#'           pgls.args=list(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),
#'                          tree=TRUE,RR=TRUE),nsim=10,clus=cc)->orr.pgls1
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree)->pgls2_noRR
#'
#' RRphylo(tree,resp.multi,clus=cc)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree,RR=RR)->pgls2_RR
#'
#' overfitRR(RR=RR,y=resp.multi,
#'           pgls.args=list(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),
#'                          tree=TRUE,RR=TRUE),nsim=10,clus=cc)->orr.pgls2
#'
#'
#' }
overfitRR<-function(RR,y,
                    phylo.list=NULL,
                    s=0.25,
                    swap.args=NULL,
                    trend.args=NULL,
                    shift.args=NULL,
                    conv.args=NULL,
                    pgls.args=NULL,
                    aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,rootV=NULL,nsim=100,
                    clus=0.5)
{
  # require(phytools)
  # require(ddpcr)
  # require(rlist)

  if (!requireNamespace("ddpcr", quietly = TRUE)) {
    stop("Package \"ddpcr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  '%ni%' <- Negate('%in%')

  RR$tree->tree
  y <- treedataMatch(tree, y)[[1]]
  RR$aces->y.ace
  tree$node.label<-rownames(y.ace)

  if(!is.null(phylo.list)){
    if(!is.null(swap.args)) warning("Swapping is not performed if a list of alternative phylogenies is provided",immediate.=TRUE)
    s<-0
    si<-0
    si2<-0
    swap.node<-NULL
    nsim<-length(phylo.list)
  }else{
    if(!is.null(swap.args)){
      if(any(is.null(names(swap.args)))) stop("All swap.args must be named")
      if(is.null(swap.args$si)) si<-0.1 else si<-swap.args$si
      if(is.null(swap.args$si2)) si2<-0.1 else si2<-swap.args$si2
      if(is.null(swap.args$node)) swap.node<-NULL else swap.node<-swap.args$node
    }else{
      si<-0.1
      si2<-0.1
      swap.node<-NULL
    }
  }

  if(!is.null(trend.args)){
    if(length(trend.args)>0&any(is.null(names(trend.args)))) stop("All trend.args must be named")
    trend<-TRUE
    if(!is.null(trend.args$node)) trend.node<-trend.args$node else trend.node<-NULL
    if(!is.null(trend.args$x1.residuals)) trend.x1.residuals<-trend.args$x1.residuals else trend.x1.residuals<-FALSE
  } else {
    trend<-FALSE
    trend.node<-NULL
    trend.x1.residuals<-FALSE
  }

  if(!is.null(shift.args)){
    if(any(is.null(names(shift.args)))) stop("All shift.args must be named")
    if(!is.null(shift.args$node)) shift.node<-shift.args$node else shift.node<-NULL
    if(!is.null(shift.args$state)) {
      shift.state<-shift.args$state
      shift.state<-treedataMatch(tree,shift.state)[[1]][,1]
    }else shift.state<-NULL
  }else{
    shift.node<-NULL
    shift.state<-NULL
  }

  if(!is.null(conv.args)){
    if(any(is.null(names(conv.args)))) stop("All conv.args must be named")
    if(!is.null(conv.args$node)) conv.node<-conv.args$node else conv.node<-NULL
    if(!is.null(conv.args$state)){
      conv.state<-conv.args$state
      conv.state<-treedataMatch(tree,conv.state)[[1]][,1]
    }else conv.state<-NULL
    if(!is.null(conv.args$declust)) conv.declust<-conv.args$declust else conv.declust<-FALSE
  }else{
    conv.node<-NULL
    conv.state<-NULL
    conv.declust<-NULL
  }

  if(!is.null(pgls.args)){
    if(any(is.null(names(pgls.args)))) stop("All pgls.args must be named")
    modform<-pgls.args$modform
    pgls.data<-pgls.args$data
    if(pgls.args$tree) pgls.tree<-pgls.args$tree else pgls.tree<-NULL
    if(pgls.args$RR) pgls.RR<-pgls.args$RR else pgls.RR<-NULL
    pgls.args<-pgls.args[which(!names(pgls.args)%in%c("modform","data","tree","RR"))]
  }else{
    modform<-NULL
    pgls.data<-NULL
    pgls.tree<-NULL
    pgls.RR<-NULL
  }

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  rootlist<-list()
  RR.list<-tree.list<-list()
  acefit<-STcut<-SScut<-SScutS<-SCcut<-SCcutS<-PGLScut<-PGLScutRR<-list()
  trend.node.match<-shift.node.match<-conv.node.match<-list()
  real.s<-array()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)
    if(s>0){
      unlist(lapply(trend.node,function(x) {
        length(tips(tree,x))->lenx
        if(lenx<=5) tips(tree,x) else sample(tips(tree,x),5)
      }))->out.st
      unlist(lapply(shift.node,function(x) {
        length(tips(tree,x))->lenx
        if(lenx<=5) tips(tree,x) else sample(tips(tree,x),5)
      }))->out.ss
      unlist(lapply(conv.node,function(x) {
        length(tips(tree,x))->lenx
        if(lenx<=5) tips(tree,x) else sample(tips(tree,x),5)
      }))->out.sc

      if(!is.null(shift.state)){
        table(shift.state)->tab.ss
        unlist(lapply(1:length(tab.ss),function(x) {
          if(tab.ss[x]<=5) names(which(shift.state==names(tab.ss)[x])) else
            sample(names(which(shift.state==names(tab.ss)[x])),5)
        }))->out.st.ss
      } else out.st.ss<-NULL

      if(!is.null(conv.state)){
        table(conv.state)->tab.cs
        unlist(lapply(1:length(tab.cs),function(x) {
          if(tab.cs[x]<=5) names(which(conv.state==names(tab.cs)[x])) else
            sample(names(which(conv.state==names(tab.cs)[x])),5)
        }))->out.st.sc
      }else out.st.sc<-NULL

      unique(c(out.st,out.ss,out.sc,out.st.ss,out.st.sc))->outs
      if(length(outs>0)) tree$tip.label[-match(outs,tree$tip.label)]->samtips else tree$tip.label->samtips

      sx<-s
      repeat({
        if(length(samtips)>Ntip(tree)*sx) break else s*.9->sx
      })
      sample(samtips,round(Ntip(tree)*sx,0))->offs
    }

    if(!is.null(phylo.list)) phylo.list[[k]]->tree.swap else
      suppressWarnings(swapONE(tree,si=si,si2=si2,node=swap.node,plot.swap=FALSE)[[1]])->tree.swap

    y[match(tree.swap$tip.label,rownames(y)),,drop=FALSE]->y

    if(s>0){
      tree.swap$edge[tree.swap$edge[,1]==(Ntip(tree.swap)+1),2]->rootdesc
      if(length(which(rootdesc<(Ntip(tree.swap)+1)))>0) tree.swap$tip.label[rootdesc[which(rootdesc<Ntip(tree.swap)+1)]]->saver else saver="xx"
      if(saver%in%offs) offs[-which(offs==saver)]->offs

      y[-which(rownames(y)%in%offs),,drop=FALSE]->ycut
      drop.tip(tree.swap,which(rownames(y)%ni%rownames(ycut)))->treecut
      y.ace[which(rownames(y.ace)%in%treecut$node.label),,drop=FALSE]->y.acecut
    }else{
      y->ycut
      tree.swap->treecut
      y.ace->y.acecut
    }

    treecut->tree.list[[k]]

    1-(Ntip(treecut)/Ntip(tree))->real.s[k]

    if(!is.null(cov)){
      treedataMatch(treecut,cov)$y->covcut
      c(RRphylo(treecut,covcut,clus=clus)$aces[,1],covcut)->covcut
      # cov[match(c(rownames(y.acecut),rownames(ycut)),names(cov))]->covcut
      # names(covcut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else covcut<-NULL

    if(!is.null(x1)) {
      as.matrix(x1)->x1
      treedataMatch(treecut,x1)$y->x1cut
      rbind(RRphylo(treecut,x1cut,clus=clus)$aces,x1cut)->x1cut
      # x1[match(c(rownames(y.acecut),rownames(ycut)),rownames(x1)),,drop=FALSE]->x1cut
      # rownames(x1cut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else x1cut<-NULL

    if(!is.null(aces)){
      if(is.vector(aces)) as.matrix(aces)->aces
      aces->acescut

      drop<-c()
      for(i in 1:nrow(aces)) {
        if(length(which(tips(tree,rownames(aces)[i])%in%treecut$tip.label))>1){
          getMRCA(treecut,tips(tree,rownames(aces)[i])[which(tips(tree,rownames(aces)[i])%in%treecut$tip.label)])->newN
          if(!is.null(phylo.list)){
            length(tips(treecut,newN))/length(tips(tree,rownames(aces)[i]))->sh.tips
            if(sh.tips<=1.1&sh.tips>=0.9) newN->rownames(acescut)[i] else c(drop,i)->drop
          }else newN->rownames(acescut)[i]
          # getMRCA(treecut,tips(tree,rownames(aces)[i])[which(tips(tree,rownames(aces)[i])%in%treecut$tip.label)])->rownames(acescut)[i]
        }else c(drop,i)->drop
      }
      if(length(drop>0)) acescut[-drop,]->acescut
      if(is.null(nrow(acescut))) acescut<-NULL
    }else acescut<-NULL

    if(!is.null(aces.x1)){
      if(is.vector(aces.x1)) as.matrix(aces.x1)->aces.x1
      aces.x1->aces.x1cut
      drop<-c()
      for(i in 1:nrow(aces.x1)) {
        if(length(which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label))>1){
          getMRCA(treecut,tips(tree,rownames(aces.x1)[i])[which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label)])->newN1
          if(!is.null(phylo.list)){
            length(tips(treecut,newN1))/length(tips(tree,rownames(aces.x1)[i]))->sh.tips
            if(sh.tips<=1.1&sh.tips>=0.9) newN1->rownames(aces.x1cut)[i] else c(drop,i)->drop
          }else newN1->rownames(aces.x1cut)[i]
          # getMRCA(treecut,tips(tree,rownames(aces.x1)[i])[which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label)])->rownames(aces.x1cut)[i]
        }else c(drop,i)->drop
      }
      if(length(drop>0)) aces.x1cut[-drop,,drop=FALSE]->aces.x1cut
      if(is.null(nrow(aces.x1cut))) aces.x1cut<-NULL
    }else aces.x1cut<-NULL

    if(!is.null(trend.node)){
      trend.node.cut<-array()
      for(i in 1:length(trend.node)) {
        getMRCA(treecut,tips(tree,trend.node[i])[which(tips(tree,trend.node[i])%in%treecut$tip.label)])->trN
        if(!is.null(phylo.list)){
          length(tips(treecut,trN))/length(tips(tree,trend.node[i]))->sh.tips
          if(sh.tips<=1.1&sh.tips>=0.9) trN->trend.node.cut[i] else NA->trend.node.cut[i]
        } else trN->trend.node.cut[i]
        # getMRCA(treecut,tips(tree,trend.node[i])[which(tips(tree,trend.node[i])%in%treecut$tip.label)])->trend.node.cut[i]
      }
      data.frame(trend.node,trend.node.cut)->trend.node.match[[k]]

      trend.node.cut[which(!is.na(trend.node.cut))]->trend.node.cut
      if(length(trend.node.cut)==0) trend.node.cut<-NULL
    }else trend.node.cut<-NULL

    if(!is.null(shift.node)){
      shift.node.cut<-array()
      for(i in 1:length(shift.node)){
        getMRCA(treecut,tips(tree,shift.node[i])[which(tips(tree,shift.node[i])%in%treecut$tip.label)])->shN
        if(!is.null(phylo.list)){
          length(tips(treecut,shN))/length(tips(tree,shift.node[i]))->sh.tips
          if(sh.tips<=1.1&sh.tips>=0.9) shN->shift.node.cut[i] else NA->shift.node.cut[i]
        } else shN->shift.node.cut[i]
      }
      # getMRCA(treecut,tips(tree,shift.node[i])[which(tips(tree,shift.node[i])%in%treecut$tip.label)])->shift.node.cut[i]
      data.frame(shift.node,shift.node.cut)->shift.node.match[[k]]
      shift.node.cut[which(!is.na(shift.node.cut))]->shift.node.cut
      if(length(shift.node.cut)==0) shift.node.cut<-NULL

    }

    if(!is.null(shift.state)) {
      shift.state[match(c(tree.swap$node.label,tree.swap$tip.label), names(shift.state))]->shift.state
      shift.state[match(rownames(ycut),names(shift.state))]->shift.state.cut
    }

    if(!is.null(conv.node)){
      conv.node.cut<-array()
      for(i in 1:length(conv.node)){
        getMRCA(treecut,tips(tree,conv.node[i])[which(tips(tree,conv.node[i])%in%treecut$tip.label)])->scN
        if(!is.null(phylo.list)){
          length(tips(treecut,scN))/length(tips(tree,conv.node[i]))->sh.tips
          if(sh.tips<=1.1&sh.tips>=0.9) scN->conv.node.cut[i] else NA->conv.node.cut[i]
        } else scN->conv.node.cut[i]
        if(any(is.na(conv.node.cut))) conv.node.cut<-NULL
      }
      # getMRCA(treecut,tips(tree,conv.node[i])[which(tips(tree,conv.node[i])%in%treecut$tip.label)])->conv.node.cut[i]
      if(!is.null(conv.node.cut)) data.frame(conv.node,conv.node.cut)->conv.node.match[[k]] else NULL->conv.node.match[[k]]
    }

    if(!is.null(conv.state)) {
      conv.state[match(c(tree.swap$node.label,tree.swap$tip.label), names(conv.state))]->conv.state
      conv.state[match(rownames(ycut),names(conv.state))]->conv.state.cut
    }

    if(!is.null(pgls.tree)|!is.null(pgls.RR)) {
      ddpcr::quiet(lapply(pgls.data,function(x){
        if(is.null(nrow(x))) treedataMatch(treecut, x)[[1]][,1] else treedataMatch(treecut, x)[[1]]
      })->pgls.datacut)
    }

    if(!is.null(rootV)) rootV->rootVcut else rootVcut<-NULL

    RRphylo(treecut,ycut,aces=acescut,x1=x1cut,aces.x1=aces.x1cut,cov=covcut,rootV = rootVcut,clus=clus)->RRcut->RR.list[[k]]
    if(trend|!is.null(trend.node)) ddpcr::quiet(search.trend(RRcut,ycut,x1=x1cut,x1.residuals = trend.x1.residuals,node=trend.node.cut,cov=covcut,clus=clus)->stcut->STcut[[k]],all=TRUE)
    if(!is.null(shift.node)&&!is.null(shift.node.cut)) ddpcr::quiet(search.shift(RRcut,status.type="clade",node=shift.node.cut)->sscut->SScut[[k]],all=TRUE)
    if(!is.null(shift.state)) ddpcr::quiet(search.shift(RRcut,status.type="sparse",state=shift.state.cut)->sscut->SScutS[[k]],all=TRUE)
    if(!is.null(conv.node)&&any(!is.na(conv.node.cut))) ddpcr::quiet(search.conv(RR=RRcut,y=ycut,nodes=na.omit(conv.node.cut),aceV=acescut,clus=clus)->sccut->SCcut[[k]],all=TRUE)
    if(!is.null(conv.state)) ddpcr::quiet(search.conv(tree=treecut,y=ycut,state=conv.state.cut,aceV=acescut,declust=conv.declust,clus=clus)->sccut->SCcutS[[k]],all=TRUE)
    # if(!is.null(pgls.tree)) ddpcr::quiet(PGLS_fossil(modform,data=pgls.datacut,tree=treecut)->PGLScut[[k]],all=TRUE)
    # if(!is.null(pgls.RR)) ddpcr::quiet(PGLS_fossil(modform,data=pgls.datacut,tree=RRcut$tree,RR=RRcut)->PGLScutRR[[k]],all=TRUE)
    if(!is.null(pgls.tree)) ddpcr::quiet(do.call(PGLS_fossil,c(list(modform=modform,data=pgls.datacut,tree=treecut),pgls.args))->PGLScut[[k]],all=TRUE)
    if(!is.null(pgls.RR)) ddpcr::quiet(do.call(PGLS_fossil,c(list(modform=modform,data=pgls.datacut,RR=RRcut),pgls.args))->PGLScutRR[[k]],all=TRUE)

    RRcut$aces[1,]->rootlist[[k]]
    summary(lm(y.acecut~RRcut$aces))->acefit[[k]]
    do.call(rbind,lapply(seq(1:ncol(y.acecut)),function(x) summary(lm(y.acecut[,x]~RRcut$aces[,x]))$coef[c(1,2,7,8)]))->acefit[[k]]
    if(!is.null(colnames(y))) rownames(acefit[[k]])<-colnames(y) else rownames(acefit[[k]])<-sapply(1:ncol(y),function(x) paste("y",x,sep=""))

    colnames(acefit[[k]])<-c("intercept","slope","p.intercept","p.slope")
  }

  if(length(unlist(rootlist))>length(rootlist)){
    do.call(rbind,rootlist)->rootlist
    apply(rootlist,2,function(x) quantile(x,c(0.025,0.975)))->CIroot
    data.frame(root=t(y.ace)[,1],"CI 2.5"=t(CIroot)[,1],"CI 97.5"=t(CIroot)[,2])->root.conf.int
    if(!is.null(colnames(y))) rownames(root.conf.int)<-colnames(y) else rownames(root.conf.int)<-sapply(1:ncol(y),function(x) paste("y",x,sep=""))
  }else{
    unlist(rootlist)->rootlist
    quantile(rootlist,c(0.025,0.975))->CIroot
    data.frame(root=y.ace[1,,drop=FALSE],"CI 2.5"=CIroot[1],"CI 97.5"=CIroot[2])->root.conf.int
  }

  if(!is.null(shift.node)){
    mapply(a=shift.node.match,b=SScut,function(a,b){
      a[match(rownames(b$single.clade),a[,2]),1]->rownames(b$single.clade)
      b$single.clade
    },SIMPLIFY = FALSE)->singles

    t(sapply(shift.node,function(k){
      t(sapply(singles,function(j) {
        if(any(as.numeric(rownames(j))==k)) j[which(as.numeric(rownames(j))==k),] else c(NA,NA)
      }))->pran
      pran[which(!is.na(pran[,1])),]->pran
      cbind(length(which(pran[,2]>=0.975))/nrow(pran),length(which(pran[,2]<=0.025))/nrow(pran),nrow(pran)/nsim)
    }))->shift.res.clade
    rownames(shift.res.clade)<-shift.node
    colnames(shift.res.clade)<-c("p.shift+","p.shift-","tested.trees")

    lapply(SScut,function(j) j$all.clades)->allcla
    if(!all(is.null(allcla))){
      do.call(rbind, allcla)->allcla
      rbind(cbind(length(which(allcla$p.value>=0.975))/nrow(allcla),
                  length(which(allcla$p.value<=0.025))/nrow(allcla),nrow(allcla)/nsim),
            shift.res.clade)->shift.res.clade
      rownames(shift.res.clade)[1]<-"all.clades"
    }
  }else shift.res.clade<-NULL

  if(!is.null(shift.state)){
    p.shift<-matrix(ncol=2,nrow=nrow(SScutS[[1]][[1]]))
    for(i in 1:nrow(SScutS[[1]][[1]])){
      unlist(lapply(lapply(SScutS,"[[",1),function(x) x[i,2]))->pr
      c(length(which(pr>=0.975))/nsim,length(which(pr<=0.025))/nsim)->p.shift[i,]
    }
    rownames(p.shift)<-rownames(SScutS[[1]][[1]])
    colnames(p.shift)<-c("p.shift+","p.shift-")
    p.shift->shift.res.state

  }else shift.res.state<-NULL

  list(shift.res.clade,shift.res.state)->shift.res
  names(shift.res)<-c("clade","sparse")

  if(trend|!is.null(trend.node)){
    #### Whole tree ####
    if(ncol(y)==1) iter<-1 else iter<-ncol(y)+1
    phen.trend<-rate.trend<-list()
    for(j in 1:iter){
      as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",2),function(x) x[j,]))[,c(1,3)])->pr#->phen.ran[[j]]
      as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",3),function(x) x[j,]))[,c(1,3)])->rr#->rat.ran[[j]]

      c(sum(pr$slope>0&pr$p.random>=0.975)/nsim,
        sum(pr$slope>0&pr$p.random<=0.025)/nsim,
        sum(pr$slope<0&pr$p.random>=0.975)/nsim,
        sum(pr$slope<0&pr$p.random<=0.025)/nsim)->phen.trend[[j]]

      c(sum(rr$slope>0&rr$p.random>=0.975)/nsim,
        sum(rr$slope>0&rr$p.random<=0.025)/nsim,
        sum(rr$slope<0&rr$p.random>=0.975)/nsim,
        sum(rr$slope<0&rr$p.random<=0.025)/nsim)->rate.trend[[j]]

      names(phen.trend[[j]])<-names(rate.trend[[j]])<-c("slope+p.up","slope+p.down","slope-p.up","slope-p.down")
    }
    do.call(rbind,phen.trend)->phen.trend
    do.call(rbind,rate.trend)->rate.trend
    if(!is.null(colnames(y))){
      if(ncol(y)==1) colnam<-colnames(y) else colnam<-c(colnames(y),"multiple")
    }else{
      if(ncol(y)==1) colnam<-"y" else
        colnam<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
    }
    rownames(phen.trend)<-rownames(rate.trend)<-colnam
    list(phen.trend,rate.trend)->p.trend
    names(p.trend)<-c("phenotype","rates")
    p.trend->whole.tree.res

    if(!is.null(trend.node)){
      mapply(a=trend.node.match,b=lapply(STcut,"[[",4),function(a,b){
        a[match(names(b),a[,2]),1]->names(b)
        b
      },SIMPLIFY = FALSE)->phen.node

      mapply(a=trend.node.match,b=lapply(STcut,"[[",5),function(a,b){
        a[match(names(b),a[,2]),1]->names(b)
        b
      },SIMPLIFY = FALSE)->rat.node

      p.phen.node<-list()
      p.rate.node<-list()
      for(k in 1:length(trend.node)){
        lapply(phen.node,function(j) {
          if(any(as.numeric(names(j))==trend.node[k])) j[[which(as.numeric(names(j))==trend.node[k])]] else NA
        })->pran
        pran[which(!sapply(pran,function(w) all(is.na(w))))]->phen.pran

        lapply(rat.node,function(j) {
          if(any(as.numeric(names(j))==trend.node[k])) j[[which(as.numeric(names(j))==trend.node[k])]] else NA
        })->pran
        pran[which(!sapply(pran,function(w) all(is.na(w))))]->rat.pran

        p.phen.node.y<-matrix(ncol=7,nrow=iter)
        p.rate.node.y<-matrix(ncol=5,nrow=iter)
        for(w in 1:iter){
          as.data.frame(do.call(rbind,lapply(phen.pran,function(x) x[w,])))->pnod
          as.data.frame(do.call(rbind,lapply(rat.pran,function(x) x[w,])))->rnod

          c(sum(pnod$slope>0&pnod$p.slope>=0.975)/nrow(pnod),
            sum(pnod$slope>0&pnod$p.slope<=0.025)/nrow(pnod),
            sum(pnod$slope<0&pnod$p.slope>=0.975)/nrow(pnod),
            sum(pnod$slope<0&pnod$p.slope<=0.025)/nrow(pnod),
            sum(pnod$emm.difference>0&pnod$p.emm<=0.05)/nrow(pnod),
            sum(pnod$emm.difference<0&pnod$p.emm<=0.05)/nrow(pnod),
            nrow(pnod)/nsim)->p.phen.node.y[w,]

          c(sum(rnod$emm.difference>0&rnod$p.emm<=0.05)/nrow(rnod),
            sum(rnod$emm.difference<0&rnod$p.emm<=0.05)/nrow(rnod),
            sum((rnod$slope.node-rnod$slope.others)>0&rnod$p.slope<=0.05)/nrow(rnod),
            sum((rnod$slope.node-rnod$slope.others)<0&rnod$p.slope<=0.05)/nrow(rnod),
            nrow(rnod)/nsim)->p.rate.node.y[w,]
        }
        colnames(p.phen.node.y)<-c("slope+p.up","slope+p.down","slope-p.up","slope-p.down","p.emm+","p.emm-","tested.trees")
        colnames(p.rate.node.y)<-c("p.emm+","p.emm-","p.slope+","p.slope-","tested.trees")
        if(!is.null(colnames(y))){
          if(ncol(y)==1) colnam<-colnames(y) else colnam<-c(colnames(y),"multiple")
        }else{
          if(ncol(y)==1) colnam<-"y" else
            colnam<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
        }
        rownames(p.phen.node.y)<-rownames(p.rate.node.y)<-colnam

        p.phen.node.y->p.phen.node[[k]]
        p.rate.node.y->p.rate.node[[k]]
      }
      names(p.phen.node)<-names(p.rate.node)<-trend.node
      list(p.phen.node,p.rate.node)->p.trend.node
      names(p.trend.node)<-c("phenotype","rates")
      node.res<-p.trend.node

      if(length(trend.node)>1){ #### Node comparison ####

        apply(combn(trend.node,2),2,function(j) c(paste(j[1],j[2],sep="-"),paste(j[2],j[1],sep="-")))->tn.pairs

        lapply(STcut,function(j) j$group.comparison)->comptot

        mapply(x=lapply(comptot,"[[",1)[!sapply(comptot,is.null)],
               xx=trend.node.match[!sapply(comptot,is.null)],function(x,xx){
                 if(ncol(y)>1){
                   t(apply(x[[1]],1,function(fx) xx[match(gsub("g","",fx[1:2]),xx[,2]),1]))->x[[1]][,1:2]
                   lapply(2:length(x), function(xw) x[[xw]][,1:2]<<-x[[1]][,1:2])
                   apply(x[[1]][,1:2],1,function(fx) paste(fx,collapse="-"))->realn
                   unlist(apply(tn.pairs,2,function(jj) which(realn%in%jj)))->roword
                   lapply(x,function(fx) fx[roword,])->x
                   realn[roword]->realn
                   apply(tn.pairs,2,function(jj) sum(match(realn,jj,nomatch = 0)))->revcols
                   revcols[which(revcols>0)]->revcols

                   lapply(x,function(kk){
                     if(any(revcols==2)){
                       data.frame(kk[which(revcols==2),c(2,1,4,3)],1-kk[which(revcols==2),5],
                                  -1*kk[which(revcols==2),6],kk[which(revcols==2),7])->kk[which(revcols==2),]
                     }
                     data.frame(kk,pair=apply(kk[,1:2],1,function(jk) paste(jk,collapse = "-")))
                   })
                 }else{
                   t(apply(x,1,function(fx) xx[match(gsub("g","",fx[1:2]),xx[,2]),1]))->x[,1:2]
                   apply(x[,1:2],1,function(fx) paste(fx,collapse="-"))->realn
                   unlist(apply(tn.pairs,2,function(jj) which(realn%in%jj)))->roword
                   x[roword,]->x
                   realn[roword]->realn
                   apply(tn.pairs,2,function(jj) sum(match(realn,jj,nomatch = 0)))->revcols

                   revcols[which(revcols>0)]->revcols
                   if(any(revcols==2)){
                     data.frame(x[which(revcols==2),c(2,1)],-1*x[which(revcols==2),3],
                                x[which(revcols==2),c(4,6,5,7)])->x[which(revcols==2),]
                   }
                   data.frame(x,pair=apply(x[,1:2],1,function(jk) paste(jk,collapse = "-")))
                 }
               },SIMPLIFY = FALSE)->pcomptot

        mapply(x=lapply(comptot,"[[",2)[!sapply(comptot,is.null)],
               xx=trend.node.match[!sapply(comptot,is.null)],function(x,xx){
                 if(ncol(y)>1){
                   t(apply(x[[1]],1,function(fx) xx[match(gsub("g","",fx[1:2]),xx[,2]),1]))->x[[1]][,1:2]
                   lapply(2:length(x), function(xw) x[[xw]][,1:2]<<-x[[1]][,1:2])
                   apply(x[[1]][,1:2],1,function(fx) paste(fx,collapse="-"))->realn
                   unlist(apply(tn.pairs,2,function(jj) which(realn%in%jj)))->roword
                   lapply(x,function(fx) fx[roword,])->x
                   realn[roword]->realn
                   apply(tn.pairs,2,function(jj) sum(match(realn,jj,nomatch = 0)))->revcols
                   revcols[which(revcols>0)]->revcols

                   lapply(x,function(kk){
                     if(any(revcols==2)){
                       data.frame(kk[which(revcols==2),c(2,1)],-1*kk[which(revcols==2),3],kk[which(revcols==2),c(4,6,5,7)])->kk[which(revcols==2),]
                     }
                     data.frame(kk,pair=apply(kk[,1:2],1,function(jk) paste(jk,collapse = "-")))
                   })
                 }else{
                   t(apply(x,1,function(fx) xx[match(gsub("g","",fx[1:2]),xx[,2]),1]))->x[,1:2]
                   apply(x[,1:2],1,function(fx) paste(fx,collapse="-"))->realn
                   unlist(apply(tn.pairs,2,function(jj) which(realn%in%jj)))->roword
                   x[roword,]->x
                   realn[roword]->realn
                   apply(tn.pairs,2,function(jj) sum(match(realn,jj,nomatch = 0)))->revcols
                   revcols[which(revcols>0)]->revcols

                   if(any(revcols==2)){
                     data.frame(x[which(revcols==2),c(2,1)],-1*x[which(revcols==2),3],
                                x[which(revcols==2),c(4,6,5,7)])->x[which(revcols==2),]
                   }
                   data.frame(x,pair=apply(x[,1:2],1,function(jk) paste(jk,collapse = "-")))
                 }
               },SIMPLIFY = FALSE)->rcomptot


        comp.phen.y<-comp.rat.y<-list()
        for(w in 1:iter){
          nod.nam<-list()
          p.comp.phen<-p.comp.rat<-matrix(ncol=5,nrow=ncol(tn.pairs))
          for(k in 1:ncol(tn.pairs)){
            if(ncol(y)>1)
              do.call(rbind,lapply(lapply(pcomptot,"[[",w),function(x) x[which(x$pair==tn.pairs[1,k]),]))->pcomp else
                do.call(rbind,lapply(pcomptot,function(x) x[which(x$pair==tn.pairs[1,k]),]))->pcomp

            if(w==1) pcomp[nrow(pcomp),1:2]->nod.nam[[k]]
            as.data.frame(pcomp[,3:7])->pcomp#->phen.comp[[k]]
            if(ncol(y)>1)
              do.call(rbind,lapply(lapply(rcomptot,"[[",w),function(x) x[which(x$pair==tn.pairs[1,k]),]))[,3:7,drop=FALSE]->rcomp else
                do.call(rbind,lapply(rcomptot,function(x) x[which(x$pair==tn.pairs[1,k]),]))[,3:7,drop=FALSE]->rcomp

            c(sum((pcomp$slope.group_1-pcomp$slope.group_2)>0&pcomp$p.slope>=0.95)/nrow(pcomp),
              sum((pcomp$slope.group_1-pcomp$slope.group_2)<0&pcomp$p.slope<=0.05)/nrow(pcomp),
              sum(pcomp$emm.difference>0&pcomp$p.emm<=0.05)/nrow(pcomp),
              sum(pcomp$emm.difference<0&pcomp$p.emm<=0.05)/nrow(pcomp),
              nrow(pcomp)/nsim)->p.comp.phen[k,]

            c(sum(rcomp$emm.difference>0&rcomp$p.emm<=0.05)/nrow(rcomp),
              sum(rcomp$emm.difference<0&rcomp$p.emm<=0.05)/nrow(rcomp),
              sum((rcomp$slope.group_1-rcomp$slope.group_2)>0&rcomp$p.slope<=0.05)/nrow(rcomp),
              sum((rcomp$slope.group_1-rcomp$slope.group_2)<0&rcomp$p.slope<=0.05)/nrow(rcomp),
              nrow(rcomp)/nsim)->p.comp.rat[k,]

          }
          colnames(p.comp.phen)<-c("p.slope+","p.slope-","p.emm+","p.emm-","tested.trees")
          colnames(p.comp.rat)<-c("p.emm+","p.emm-","p.slope+","p.slope-","tested.trees")
          if(w==1) do.call(rbind, nod.nam)->nam.pair

          rownames(p.comp.phen)<-rownames(p.comp.rat)<-apply(nam.pair,1, function(x) paste(x[1], x[2], sep="-"))

          p.comp.phen->comp.phen.y[[w]]
          p.comp.rat->comp.rat.y[[w]]

        }

        p.comp.phenN<-p.comp.ratN<-list()
        for(q in 1:ncol(tn.pairs)){
          do.call(rbind,lapply(comp.phen.y,function(x) x[q,]))->p.comp.phenN[[q]]
          do.call(rbind,lapply(comp.rat.y,function(x) x[q,]))->p.comp.ratN[[q]]

          if(!is.null(colnames(y))){
            if(ncol(y)==1) colnam<-colnames(y) else colnam<-c(colnames(y),"multiple")
            rownames(p.comp.phenN[[q]])<-rownames(p.comp.ratN[[q]])<-colnam
          }else{
            if(ncol(y)==1) colnam<-"y" else
              colnam<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
          }
          rownames(p.comp.phenN[[q]])<-rownames(p.comp.ratN[[q]])<-colnam
        }
        names(p.comp.phenN)<-names(p.comp.ratN)<-rownames(comp.phen.y[[1]])
        list(p.comp.phenN,p.comp.ratN)->p.comp
        names(p.comp)<-c("phenotype","rates")
      }else{
        p.comp<-NULL
      }

      if(length(trend.node)>1) node.res<-list(node=node.res,comparison=p.comp) else node.res<-list(node=node.res)
      trend.res<-do.call(c,list(tree=list(whole.tree.res),node.res))
    }else trend.res<-whole.tree.res

  }else trend.res<-NULL

  if(!is.null(conv.node)){
    lapply(SCcut,"[[",1)->scres
    scres[which(!sapply(SCcut,is.null))]->scres
    matrix(c(length(which(lapply(scres,function(x) x[1,8])<=0.05))/length(scres),
             length(which(lapply(scres,function(x) x[1,9])<=0.05))/length(scres),
             length(scres)/nsim),ncol=3)->p.convC
    colnames(p.convC)<-c("p.ang.bydist","p.ang.conv","tested.trees")
    rownames(p.convC)<-paste(conv.node,collapse="-")
  }else p.convC<-NULL

  if(!is.null(conv.state)){
    lapply(SCcutS,"[[",1)->scresS
    p.convS<-matrix(ncol=2,nrow=nrow(scresS[[1]]))
    if("nostate"%in%conv.state&length(unique(conv.state)[-which(unique(conv.state)=="nostate")])){
      c(length(which(sapply(scresS,function(x) x[1,3])<=0.05))/nsim,
        length(which(sapply(scresS,function(x) x[1,4])<=0.05))/nsim)->p.convS[1,]
      rownames(p.convS)<-rownames(scresS[[1]])
      colnames(p.convS)<-colnames(scresS[[1]])[3:4]
    }else{
      for(i in 1:nrow(scresS[[1]])){
        c(length(which(sapply(scresS,function(x) x[i,5])<=0.05))/nsim,
          length(which(sapply(scresS,function(x) x[i,6])<=0.05))/nsim)->p.convS[i,]
      }
      rownames(p.convS)<-apply(scresS[[1]][,1:2],1,function(x) paste(x[1],x[2],sep="-"))
      colnames(p.convS)<-colnames(scresS[[1]])[5:6]
    }

  }else p.convS<-NULL

  list(p.convC,p.convS)->conv.res
  names(conv.res)<-c("clade","state")

  if(is.null(pgls.tree)) PGLScut<-NULL else class(PGLScut)<-"RRphyloList"
  if(is.null(pgls.RR)) PGLScutRR<-NULL else class(PGLScutRR)<-"RRphyloList"

  list(PGLScut,PGLScutRR)->pgls.res
  names(pgls.res)<-c("tree","RR")

  # res<-list(mean(real.s),root.conf.int,acefit,conv.res,shift.res,trend.res,pgls.res)
  # names(res)<-c("mean.sampling","rootCI","ace.regressions","conv.results","shift.results","trend.results","pgls.results")

  class(RR.list)<-"RRphyloList"
  class(acefit)<-"RRphyloList"
  class(tree.list)<-"multiPhylo"

  res<-structure(list(mean.sampling = mean(real.s),
                      tree.list=tree.list,
                      RR.list=RR.list,
                      rootCI=root.conf.int,
                      ace.regressions=acefit,
                      conv.results=conv.res,
                      shift.results=shift.res,
                      trend.results=trend.res,
                      pgls.results=pgls.res),
                 class = "RRphyloList")
  res
}


#' @export
print.RRphyloList<-function(x,...){
  if("mean.sampling"%in%attributes(x)[[1]]) cat(paste(length(x[[2]]),"overfitRR simulations",sep=" ")) else
    if("lambda"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"RRphylo outputs",sep=" ")) else
      if(class(x[[1]])[1]%in%c("gls","procD.lm")) cat("List of",paste(length(x),"PGLS_fossil outputs",sep=" ")) else
        cat("List of",paste(length(x),"outputs",sep=" "))

}

