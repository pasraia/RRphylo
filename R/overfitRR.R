#' @title Testing RRphylo methods overfit
#' @description Testing the robustness of \code{\link{search.trend}}
#'   (\cite{Castiglione et al. 2019a}), \code{\link{search.shift}}
#'   (\cite{Castiglione et al. 2018}), and  \code{\link{search.conv}}
#'   (\cite{Castiglione et al. 2019b}) results to sampling effects and
#'   phylogenetic uncertainty.
#' @usage
#' overfitRR(RR,y,s=0.25,swap.args=NULL,trend.args=NULL,shift.args=NULL,conv.args=NULL,
#' aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,rootV=NULL,nsim=100,clus=.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a named vector of phenotypes.
#' @param s the percentage of tips to be cut off. It is set at 25\% by default.
#' @param swap.args a list of arguments to be passed to the function
#'   \code{\link{swapONE}}, including \code{list(si=NULL,si2=NULL,}
#'   \code{node=NULL)}. If \code{swap.arg} is unspecified, the function
#'   automatically sets both \code{si} and \code{si2} to 0.1.
#' @param trend.args a list of arguments specific to the function
#'   \code{search.trend}, including \code{list(node=NULL,x1.residuals=FALSE)}. If a trend for the
#'   whole tree is to be tested, type \code{trend.args = list()}. No trend is
#'   tested if left unspecified.
#' @param shift.args a list of arguments specific to the function
#'   \code{search.shift}, including \code{list(node=NULL,} \code{state=NULL)}.
#'   Arguments \code{node} and \code{state} can be specified at the same time.
#' @param conv.args a list of arguments specific to the function
#'   \code{search.conv}, including \code{list(node=NULL,} \code{state=NULL,
#'   declust=FALSE)}. Arguments \code{node} and \code{state} can be specified at
#'   the same time.
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
#' @return The function returns a 'list' containing:
#' @return \strong{$mean.sampling} the mean proportion of species actually
#'   removed from the tree over the iterations.
#' @return \strong{$rootCI} the 95\% confidence interval around the root value.
#' @return \strong{$ace.regressions} the results of linear regression between
#'   ancestral state estimates before and after the subsampling.
#' @return \strong{$conv.results} a list including results for
#'   \code{search.conv} performed under \code{clade} and \code{state}
#'   conditions. If a node pair is specified within \code{conv.args}, the
#'   \code{$clade} object contains the percentage of simulations producing
#'   significant p-values for convergence between the clades. If a state vector
#'   is supplied within \code{conv.args}, the object \code{$state} contains the
#'   percentage of simulations producing significant p-values for convergence
#'   within (single state) or between states (multiple states).
#' @return \strong{$shift.results} a list including results for
#'   \code{search.shift} performed under \code{clade} and \code{sparse}
#'   conditions. If one or more nodes are specified within \code{shift.args},
#'   the \code{$clade} object contains for each node the percentage of
#'   simulations producing significant p-value separated by shift sign, and the
#'   same figures by considering all the specified nodes as evolving under a
#'   single rate (all.clades).If a state vector is supplied within
#'   \code{shift.args}, the object \code{$sparse} contains the percentage of
#'   simulations producing significant p-value separated by shift sign
#'   ($p.states).
#' @return \strong{$trend.results} a list including the percentage of
#'   simulations showing significant p-values for phenotypes versus age and
#'   absolute rates versus age regressions for the entire tree separated by
#'   slope sign ($tree). If one or more nodes are specified within
#'   \code{trend.args}, the list also includes the same results at nodes ($node)
#'   and the results for comparison between nodes ($comparison).
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @details Methods using a large number of parameters risk being overfit. This
#'   usually translates in poor fitting with data and trees other than the those
#'   originally used. With \code{RRphylo} methods this risk is usually very low.
#'   However, the user can assess how robust the results got by applying
#'   \code{search.shift}, \code{search.trend} or \code{search.conv} are by
#'   running \code{overfitRR}. With the latter, the original tree and data are
#'   subsampled by specifying a \code{s} parameter, that is the proportion of
#'   tips to be removed from the tree. In some cases, though, removing as many
#'   tips as imposed by \code{s} would delete too many tips right in clades
#'   and/or states under testing. In these cases, the function maintains no less
#'   than 5 species at least in each clade/state under testing (or all species
#'   if there is less), reducing the sampling parameter \code{s} if necessary.
#'   Internally, \code{overfitRR} further shuffles the tree by using the
#'   function \code{\link{swapONE}}. Thereby, both the potential for overfit and
#'   phylogenetic uncertainty are accounted for straight away.
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
#' RRphylo(tree=treedino,y=massdino)->dinoRates
#' RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'
#' # Case 1 search.shift under both "clade" and "sparse" condition
#' search.shift(RR=dinoRates, status.type= "clade",foldername=tempdir())->SSnode
#' search.shift(RR=dinoRates, status.type= "sparse", state=statedino,
#'              foldername=tempdir())->SSstate
#'
#' overfitRR(RR=dinoRates,y=massdino,swap.args =list(si=0.2,si2=0.2),
#'           shift.args = list(node=rownames(SSnode$single.clades),state=statedino),
#'           nsim=10,clus=cc)
#'
#' # Case 2 search.trend on the entire tree
#' search.trend(RR=RRptero, y=log(massptero),nsim=100,clus=cc,
#'              foldername=tempdir(),cov=NULL,ConfInt=FALSE,node=NULL)->STtree
#'
#' overfitRR(RR=RRptero,y=log(massptero),swap.args =list(si=0.2,si2=0.2),
#'           trend.args = list(),nsim=10,clus=cc)
#'
#' # Case 3 search.trend at specified nodes
#' search.trend(RR=RRptero, y=log(massptero),node=143,clus=cc,foldername=tempdir(),
#'              cov=NULL,ConfInt=FALSE)->STnode
#'
#' overfitRR(RR=RRptero,y=log(massptero),trend.args = list(node=143),nsim=10,clus=cc)
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
#' RRphylo(tree=treecet.multi,y=masscet.multi)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc,foldername=tempdir())->STcet
#' overfitRR(RR=RRmulti,y=brainmasscet,trend.args = list(),x1=x1.mass,nsim=10,clus=cc)
#'
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,
#'              clus=cc,foldername=tempdir())->STcet.resi
#' overfitRR(RR=RRmulti,y=brainmasscet,trend.args = list(x1.residuals=TRUE),
#'           x1=x1.mass,nsim=10,clus=cc)
#'
#' # Case 5 searching convergence between clades and within a single state
#' data("DataFelids")
#' DataFelids$PCscoresfel->PCscoresfel
#' DataFelids$treefel->treefel
#' DataFelids$statefel->statefel
#'
#' RRphylo(tree=treefel,y=PCscoresfel,clus=cc)->RRfel
#' search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9",
#'             foldername = tempdir(),clus=cc)->SC.clade
#' as.numeric(c(rownames(SC.clade[[1]])[1],as.numeric(as.character(SC.clade[[1]][1,1]))))->conv.nodes
#'
#' overfitRR(RR=RRfel, y=PCscoresfel,conv.args =
#' list(node=conv.nodes,state=statefel,declust=TRUE),nsim=10,clus=cc)
#'
#'
#' }
overfitRR<-function(RR,y,
                    s=0.25,
                    swap.args=NULL,
                    trend.args=NULL,
                    shift.args=NULL,
                    conv.args=NULL,
                    aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,rootV=NULL,nsim=100,
                    clus=.5)
{
  # require(phytools)
  # require(geiger)
  # require(ddpcr)
  # require(rlist)

  if (!requireNamespace("ddpcr", quietly = TRUE)) {
    stop("Package \"ddpcr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  RR$tree->tree

  # if (inherits(y,"data.frame"))
  #   y <- treedata(tree, y, sort = TRUE)[[2]]
  if(is.null(nrow(y))) y <- treedata(tree, y, sort = TRUE)[[2]][,1] else y <- treedata(tree, y, sort = TRUE)[[2]]

  if (length(y) > Ntip(tree)) {
    # if((Ntip(tree)-round(Ntip(tree)*s))<ncol(y))
    # stop(paste("After the sampling there are more variables than observations.
    #            Please consider running your preliminary analyses with",
    #            (Ntip(tree)-round(Ntip(tree)*s)),"variables.",sep=" "))
    RR$aces->y.ace
    tree$node.label<-rownames(y.ace)
    y[match(tree$tip.label,rownames(y)),]->y
  }else{
    RR$aces[,1]->y.ace
    tree$node.label<-names(y.ace)
    y[match(tree$tip.label,names(y))]->y
  }

  if(is.null(swap.args)==FALSE){
    if(is.null(swap.args$si)) si<-0.1 else si<-swap.args$si
    if(is.null(swap.args$si2)) si2<-0.1 else si2<-swap.args$si2
    if(is.null(swap.args$node)) swap.node<-NULL else swap.node<-swap.args$node
  }else{
    si<-0.1
    si2<-0.1
    swap.node<-NULL
  }

  if(is.null(trend.args)==FALSE){
    trend<-TRUE
    if(!is.null(trend.args$node)) trend.node<-trend.args$node else trend.node<-NULL
    if(!is.null(trend.args$x1.residuals)) trend.x1.residuals<-trend.args$x1.residuals else trend.x1.residuals<-FALSE
  } else {
    trend<-FALSE
    trend.node<-NULL
    trend.x1.residuals<-FALSE
  }

  if(is.null(shift.args)==FALSE){
    if(is.null(shift.args$node)==FALSE) shift.node<-shift.args$node else shift.node<-NULL
    if(is.null(shift.args$state)==FALSE) {
      shift.state<-shift.args$state
      shift.state<-treedata(tree,shift.state, sort = TRUE)[[2]][,1]
      }else shift.state<-NULL
  }else{
    shift.node<-NULL
    shift.state<-NULL
  }

  if(is.null(conv.args)==FALSE){
    if(is.null(conv.args$node)==FALSE) conv.node<-conv.args$node else conv.node<-NULL
    if(is.null(conv.args$state)==FALSE){
      conv.state<-conv.args$state
      conv.state<-treedata(tree,conv.state, sort = TRUE)[[2]][,1]
    }else conv.state<-NULL
    if(is.null(conv.args$declust)==FALSE) conv.declust<-conv.args$declust else conv.declust<-FALSE
  }else{
    conv.node<-NULL
    conv.state<-NULL
    conv.declust<-NULL
  }


  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  rootlist<-list()
  acefit<-STcut<-SScut<-SScutS<-SCcut<-SCcutS<-list()
  real.s<-array()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)
    '%ni%' <- Negate('%in%')
    # repeat({
    #   #suppressWarnings(swapONE(tree,0.1,0.1)[[1]])->tree.swap
    #   suppressWarnings(swapONE(tree,si=si,si2=si2,node=swap.node,plot.swap=FALSE)[[1]])->tree.swap
    #   #sample(y,round((Ntip(tree)*s),0))->offs
    #   #sapply(trend.node,function(x) length(tips(tree,x))-length(which(names(offs)%in%tips(tree,x))))->lt
    #   #sapply(shift.node,function(x) length(tips(tree,x))-length(which(names(offs)%in%tips(tree,x))))->lss
    #   sample(tree$tip.label,round((Ntip(tree)*s),0))->offs
    #   sapply(trend.node,function(x) length(tips(tree,x))-length(which(offs%in%tips(tree,x))))->lt
    #   sapply(shift.node,function(x) length(tips(tree,x))-length(which(offs%in%tips(tree,x))))->lss
    #   sapply(conv.node,function(x) length(tips(tree,x))-length(which(offs%in%tips(tree,x))))->lsc
    #   if(length(which(lt<5))==0&length(which(lss<5))==0&length(which(lsc<5))==0) break
    # })

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

    suppressWarnings(swapONE(tree,si=si,si2=si2,node=swap.node,plot.swap=FALSE)[[1]])->tree.swap

    if(length(y)>Ntip(tree)) y[match(tree.swap$tip.label,rownames(y)),]->y else y[match(tree.swap$tip.label,names(y))]->y

    if(s>0){
      tree.swap$edge[tree.swap$edge[,1]==(Ntip(tree.swap)+1),2]->rootdesc
      if(length(which(rootdesc<(Ntip(tree.swap)+1)))>0) tree.swap$tip.label[rootdesc[which(rootdesc<Ntip(tree.swap)+1)]]->saver else saver="xx"
      if(saver%in%offs) offs[-which(offs==saver)]->offs

      if(length(y)>Ntip(tree)) {
        y[-which(rownames(y)%in%offs),]->ycut
        drop.tip(tree.swap,which(rownames(y)%ni%rownames(ycut)))->treecut
        y.ace[which(rownames(y.ace)%in%treecut$node.label),]->y.acecut
      }else{
        y[-which(names(y)%in%offs)]->ycut
        drop.tip(tree.swap,which(names(y)%ni%names(ycut)))->treecut
        y.ace[which(names(y.ace)%in%treecut$node.label)]->y.acecut
      }

    }else{
      y->ycut
      tree.swap->treecut
      y.ace->y.acecut
    }

    1-(Ntip(treecut)/Ntip(tree))->real.s[k]

    if(is.null(cov)==FALSE) {
      cov[match(c(tree.swap$node.label,tree.swap$tip.label), names(cov))]->cov
      if(length(y)>Ntip(tree))
        cov[match(c(rownames(y.acecut),rownames(ycut)),names(cov))]->covcut else
          cov[match(c(names(y.acecut),names(ycut)),names(cov))]->covcut
      names(covcut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else covcut<-NULL

    if(is.null(x1)==FALSE) {
      x1[match(c(tree.swap$node.label,tree.swap$tip.label), names(x1))]->x1
      if(length(y)>Ntip(tree))
        x1[match(c(rownames(y.acecut),rownames(ycut)),names(x1))]->x1cut else
          x1[match(c(names(y.acecut),names(ycut)),names(x1))]->x1cut
      names(x1cut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else x1cut<-NULL

    if(is.null(aces)==FALSE){
      aces->acescut

      if(length(y)>Ntip(tree)){
        drop<-c()
        for(i in 1:nrow(aces)) {
          if(length(which(tips(tree,rownames(aces)[i])%in%treecut$tip.label))>1)
            getMRCA(treecut,tips(tree,rownames(aces)[i])[which(tips(tree,rownames(aces)[i])%in%treecut$tip.label)])->rownames(acescut)[i] else
              c(drop,i)->drop
        }
        if(length(drop>0)) acescut[-drop,]->acescut
      }else{
        drop<-c()
        for(i in 1:length(aces)) {
          if(length(which(tips(tree,names(aces[i]))%in%treecut$tip.label))>1)
            getMRCA(treecut,tips(tree,names(aces[i]))[which(tips(tree,names(aces[i]))%in%treecut$tip.label)])->names(acescut)[i] else
              c(drop,i)->drop
        }
        if(length(drop>0)) acescut[-drop]->acescut
      }

    }else acescut<-NULL

    if(is.null(aces.x1)==FALSE){
      aces.x1->aces.x1cut
      drop<-c()
      for(i in 1:length(aces.x1)) {
        if(length(which(tips(tree,names(aces.x1[i]))%in%treecut$tip.label))>1)
          getMRCA(treecut,tips(tree,names(aces.x1[i]))[which(tips(tree,names(aces.x1[i]))%in%treecut$tip.label)])->names(aces.x1cut)[i] else
            c(drop,i)->drop
      }
      if(length(drop>0)) aces.x1cut[-drop]->aces.x1cut

    }else aces.x1cut<-NULL

    if(is.null(trend.node)==FALSE){
      trend.node.cut<-array()
      for(i in 1:length(trend.node)) getMRCA(treecut,tips(tree,trend.node[i])[which(tips(tree,trend.node[i])%in%treecut$tip.label)])->trend.node.cut[i]
    }else trend.node.cut<-NULL

    if(is.null(shift.node)==FALSE){
      shift.node.cut<-array()
      for(i in 1:length(shift.node)) getMRCA(treecut,tips(tree,shift.node[i])[which(tips(tree,shift.node[i])%in%treecut$tip.label)])->shift.node.cut[i]
    }

    if(is.null(shift.state)==FALSE) {
      shift.state[match(c(tree.swap$node.label,tree.swap$tip.label), names(shift.state))]->shift.state
      if(length(y)>Ntip(tree))
        shift.state[match(rownames(ycut),names(shift.state))]->shift.state.cut else
          shift.state[match(names(ycut),names(shift.state))]->shift.state.cut
    }

    if(is.null(conv.node)==FALSE){
      conv.node.cut<-array()
      for(i in 1:length(conv.node)) getMRCA(treecut,tips(tree,conv.node[i])[which(tips(tree,conv.node[i])%in%treecut$tip.label)])->conv.node.cut[i]
    }

    if(is.null(conv.state)==FALSE) {
      conv.state[match(c(tree.swap$node.label,tree.swap$tip.label), names(conv.state))]->conv.state
      if(length(y)>Ntip(tree))
        conv.state[match(rownames(ycut),names(conv.state))]->conv.state.cut else
          conv.state[match(names(ycut),names(conv.state))]->conv.state.cut
    }

    if(is.null(rootV)==FALSE) rootV->rootVcut else rootVcut<-NULL

    RRphylo(treecut,ycut,aces=acescut,x1=x1cut,aces.x1=aces.x1cut,cov=covcut,rootV = rootVcut,clus=clus)->RRcut
    if(trend|is.null(trend.node)==FALSE) ddpcr::quiet(search.trend(RRcut,ycut,x1=x1cut,x1.residuals = trend.x1.residuals,node=trend.node.cut,foldername=tempdir(),cov=covcut,clus=clus)->stcut->STcut[[k]],all=TRUE)
    if(is.null(shift.node)==FALSE) ddpcr::quiet(search.shift(RRcut,status.type="clade",node=shift.node.cut,foldername=tempdir())->sscut->SScut[[k]],all=TRUE)
    if(is.null(shift.state)==FALSE) ddpcr::quiet(search.shift(RRcut,status.type="sparse",state=shift.state.cut,foldername=tempdir())->sscut->SScutS[[k]],all=TRUE)
    if(is.null(conv.node)==FALSE) ddpcr::quiet(search.conv(RR=RRcut,y=ycut,nodes=conv.node.cut,aceV=acescut,clus=clus,foldername=tempdir())->sccut->SCcut[[k]],all=TRUE)
    if(is.null(conv.state)==FALSE) ddpcr::quiet(search.conv(tree=treecut,y=ycut,state=conv.state.cut,aceV=acescut,declust=conv.declust,clus=clus,foldername=tempdir())->sccut->SCcutS[[k]],all=TRUE)

    RRcut$aces[1,]->rootlist[[k]]
    summary(lm(y.acecut~RRcut$aces))->acefit[[k]]
    if(length(y)>Ntip(tree)){
      do.call(rbind,lapply(seq(1:ncol(y.acecut)),function(x) summary(lm(y.acecut[,x]~RRcut$aces[,x]))$coef[c(1,2,7,8)]))->acefit[[k]]
      if(!is.null(colnames(y))) rownames(acefit[[k]])<-colnames(y) else rownames(acefit[[k]])<-sapply(1:ncol(y),function(x) paste("y",x,sep=""))
    }else matrix(summary(lm(y.acecut~RRcut$aces))$coef[c(1,2,7,8)],nrow=1)->acefit[[k]]
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
    data.frame(root=y.ace[1],"CI 2.5"=CIroot[1],"CI 97.5"=CIroot[2])->root.conf.int
  }

  if(is.null(shift.node)==FALSE){
    if(length(SScut[[1]])>=2){
      unlist(lapply(lapply(SScut,"[[",1),function(x) x[,2]))-> p.ran.whole
      c(length(which(p.ran.whole>=0.975))/nsim,length(which(p.ran.whole<=0.025))/nsim)->p.shift.whole
      do.call(rbind,lapply(lapply(SScut,"[[",2),function(x) x[,2]))-> p.ran
      cbind(apply(p.ran,2,function(x) length(which(x>=0.975)))/nsim,
            apply(p.ran,2,function(x) length(which(x<=0.025)))/nsim)->p.shift
      rbind(p.shift.whole,p.shift)->shift.res.clade
      rownames(shift.res.clade)<-c("all.clades",shift.node)
      colnames(shift.res.clade)<-c("p.shift+","p.shift-")

    }else{
      sapply(lapply(SScut,"[[",1),"[[",2)-> p.ran
      matrix(c(length(which(p.ran>=0.975))/nsim,length(which(p.ran<=0.025))/nsim),ncol=2)->shift.res.clade
      rownames(shift.res.clade)<-shift.node
      colnames(shift.res.clade)<-c("p.shift+","p.shift-")
    }
  }else shift.res.clade<-NULL

  if(is.null(shift.state)==FALSE){
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

  if(trend|is.null(trend.node)==FALSE){
    #### Whole tree ####
    if(length(y)>Ntip(tree)){ #### Multivariate ####
      phen.trend<-rate.trend<-list()
      for(j in 1:(ncol(y)+1)){
        as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",3),function(x) x[j,]))[,c(1,3)])->pr#->phen.ran[[j]]
        as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",4),function(x) x[j,]))[,c(1,3)])->rr#->rat.ran[[j]]

        c(length(which(pr$slope>0&pr$p.random<=0.05))/nsim,
          length(which(pr$slope<0&pr$p.random<=0.05))/nsim)->phen.trend[[j]]

        c(length(which(rr$slope>0&rr$p.random<=0.05))/nsim,
          length(which(rr$slope<0&rr$p.random<=0.05))/nsim)->rate.trend[[j]]

        names(phen.trend[[j]])<-names(rate.trend[[j]])<-c("p.slope+","p.slope-")
      }
      do.call(rbind,phen.trend)->phen.trend
      do.call(rbind,rate.trend)->rate.trend
      if(!is.null(colnames(y))){
        rownames(phen.trend)<-c(colnames(y),"multiple")
        rownames(rate.trend)<-c(colnames(y),"multiple")
      }else{
        rownames(phen.trend)<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
        rownames(rate.trend)<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
      }
      list(phen.trend,rate.trend)->p.trend
      names(p.trend)<-c("phenotype","rates")
    }else{ #### Univariate ####
      as.data.frame(do.call(rbind,lapply(STcut,"[[",3))[,c(1,3)])->phen.ran
      as.data.frame(do.call(rbind,lapply(STcut,"[[",4))[,c(1,3)])->rat.ran

      rbind(c(length(which(phen.ran$slope>0&phen.ran$p.random<=0.05))/nsim,
              length(which(phen.ran$slope<0&phen.ran$p.random<=0.05))/nsim),
            c(length(which(rat.ran$slope>0&rat.ran$p.random<=0.05))/nsim,
              length(which(rat.ran$slope<0&rat.ran$p.random<=0.05))/nsim))->p.trend

      colnames(p.trend)<-c("p.slope+","p.slope-")
      rownames(p.trend)<-c("phenotype","rates")
    }

    p.trend->whole.tree.res

    if(is.null(trend.node)==FALSE){
      lapply(STcut,"[[",5)->phen.node
      lapply(STcut,"[[",6)->rat.node

      if(length(STcut[[1]])==7){ #### Node comparison ####
        if(length(trend.node)>2) { #### More than 2 nodes ####
          if(length(y)>Ntip(tree)){ #### Multivariate ####
            comp.phen.y<-comp.rat.y<-nod.nam<-list()
            for(w in 1:(ncol(y)+1)){
              nod.nam<-list()
              p.comp.phen<-p.comp.rat<-matrix(ncol=4,nrow=nrow(STcut[[1]][[7]][[1]][[1]]))
              for(k in 1:nrow(STcut[[1]][[7]][[1]][[1]])){
                do.call(rbind,lapply(lapply(lapply(lapply(STcut,"[[",7),"[[",1),"[[",w),function(x) x[k,]))->pcomp
                if(w==1) pcomp[nsim,1:2]->nod.nam[[k]]
                as.data.frame(pcomp[,3:6])->pcomp#->phen.comp[[k]]
                as.data.frame(do.call(rbind,lapply(lapply(lapply(lapply(STcut,"[[",7),"[[",2),"[[",w),function(x) x[k,]))[,3:6])->rcomp

                c(length(which(pcomp$slope.difference>0&pcomp$p.slope<=0.05))/nsim,
                  length(which(pcomp$slope.difference<0&pcomp$p.slope<=0.05))/nsim,
                  length(which(pcomp$emm.difference>0&pcomp$p.emm<=0.05))/nsim,
                  length(which(pcomp$emm.difference<0&pcomp$p.emm<=0.05))/nsim)->p.comp.phen[k,]

                c(length(which(rcomp$emm.difference>0&rcomp$p.emm<=0.05))/nsim,
                  length(which(rcomp$emm.difference<0&rcomp$p.emm<=0.05))/nsim,
                  length(which(rcomp$slope.difference>0&rcomp$p.slope<=0.05))/nsim,
                  length(which(rcomp$slope.difference<0&rcomp$p.slope<=0.05))/nsim)->p.comp.rat[k,]

              }
              colnames(p.comp.phen)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
              colnames(p.comp.rat)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
              if(w==1){
                do.call(rbind, nod.nam)->nod.nam
                t(apply(nod.nam, 1, function(x) unlist(lapply(strsplit(x, "g"), "[[", 2))))->nam.pair
              }
              rownames(p.comp.phen)<-rownames(p.comp.rat)<-apply(nam.pair,1, function(x) paste(x[1], x[2], sep="-"))

              p.comp.phen->comp.phen.y[[w]]
              p.comp.rat->comp.rat.y[[w]]

            }
            p.comp.phenN<-p.comp.ratN<-list()
            for(q in 1:length(trend.node)) {
              do.call(rbind,lapply(comp.phen.y,function(x) x[q,]))->p.comp.phenN[[q]]
              do.call(rbind,lapply(comp.rat.y,function(x) x[q,]))->p.comp.ratN[[q]]

              if(!is.null(colnames(y))){
                rownames(p.comp.phenN[[q]])<-c(colnames(y),"multiple")
                rownames(p.comp.ratN[[q]])<-c(colnames(y),"multiple")
              }else{
                rownames(p.comp.phenN[[q]])<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
                rownames(p.comp.ratN[[q]])<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
              }

            }

            names(p.comp.phenN)<-names(p.comp.ratN)<-rownames(comp.phen.y[[1]])
            list(p.comp.phenN,p.comp.ratN)->p.comp
            names(p.comp)<-c("phenotype","rates")
          }else{ #### Univariate ####
            nod.nam<-list()
            p.comp.phen<-p.comp.rat<-matrix(ncol=4,nrow=nrow(STcut[[1]][[7]][[1]]))
            for(k in 1:nrow(STcut[[1]][[7]][[1]])){
              do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",1),function(w) w[k,]))->pcomp
              pcomp[nsim,1:2]->nod.nam[[k]]
              as.data.frame(pcomp[,3:6])->pcomp
              as.data.frame(do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",2),function(w) w[k,]))[,3:6])->rcomp

              c(length(which(pcomp$slope.difference>0&pcomp$p.slope<=0.05))/nsim,
                length(which(pcomp$slope.difference<0&pcomp$p.slope<=0.05))/nsim,
                length(which(pcomp$emm.difference>0&pcomp$p.emm<=0.05))/nsim,
                length(which(pcomp$emm.difference<0&pcomp$p.emm<=0.05))/nsim)->p.comp.phen[k,]

              c(length(which(rcomp$emm.difference>0&rcomp$p.emm<=0.05))/nsim,
                length(which(rcomp$emm.difference<0&rcomp$p.emm<=0.05))/nsim,
                length(which(rcomp$slope.difference>0&rcomp$p.slope<=0.05))/nsim,
                length(which(rcomp$slope.difference<0&rcomp$p.slope<=0.05))/nsim)->p.comp.rat[k,]

            }
            colnames(p.comp.phen)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
            colnames(p.comp.rat)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
            do.call(rbind, nod.nam)->nod.nam
            t(apply(nod.nam, 1, function(x) unlist(lapply(strsplit(x, "g"), "[[", 2))))->nam.pair
            t(apply(nam.pair,1,function(x) trend.node[match(x,trend.node.cut)]))->nam.pair
            rownames(p.comp.phen)<-rownames(p.comp.rat)<-apply(nam.pair,1, function(x) paste(x[1], x[2], sep="-"))

            list(p.comp.phen,p.comp.rat)->p.comp
            names(p.comp)<-c("phenotype","rates")
          }
        }else{ #### 2 nodes ####
          if(length(y)>Ntip(tree)){ #### Multivariate ####
            p.comp.phen<-p.comp.rat<-matrix(ncol=4,nrow=ncol(y)+1)
            for(w in 1:(ncol(y)+1)){
              as.data.frame(do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",1),"[[",w)))[,3:6]->phen.comp
              as.data.frame(do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",2),"[[",w)))[,3:6]->rat.comp

              c(length(which(phen.comp$slope.difference>0&phen.comp$p.slope<=0.05))/nsim,
                length(which(phen.comp$slope.difference<0&phen.comp$p.slope<=0.05))/nsim,
                length(which(phen.comp$emm.difference>0&phen.comp$p.emm<=0.05))/nsim,
                length(which(phen.comp$emm.difference<0&phen.comp$p.emm<=0.05))/nsim)->p.comp.phen[w,]

              c(length(which(rat.comp$emm.difference>0&rat.comp$p.emm<=0.05))/nsim,
                length(which(rat.comp$emm.difference<0&rat.comp$p.emm<=0.05))/nsim,
                length(which(rat.comp$slope.difference>0&rat.comp$p.slope<=0.05))/nsim,
                length(which(rat.comp$slope.difference<0&rat.comp$p.slope<=0.05))/nsim)->p.comp.rat[w,]
            }
            colnames(p.comp.phen)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
            colnames(p.comp.rat)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
            rownames(p.comp.phen)<-names(STcut[[1]][[7]][[1]])
            rownames(p.comp.rat)<-names(STcut[[1]][[7]][[2]])
          }else{ #### Univariate ####
            as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",7),"[[",1))[,3:6])->phen.comp
            as.data.frame(do.call(rbind,lapply(lapply(STcut,"[[",7),"[[",2))[,3:6])->rat.comp

            c(length(which(phen.comp$slope.difference>0&phen.comp$p.slope<=0.05))/nsim,
              length(which(phen.comp$slope.difference<0&phen.comp$p.slope<=0.05))/nsim,
              length(which(phen.comp$emm.difference>0&phen.comp$p.emm<=0.05))/nsim,
              length(which(phen.comp$emm.difference<0&phen.comp$p.emm<=0.05))/nsim)->p.comp.phen

            c(length(which(rat.comp$emm.difference>0&rat.comp$p.emm<=0.05))/nsim,
              length(which(rat.comp$emm.difference<0&rat.comp$p.emm<=0.05))/nsim,
              length(which(rat.comp$slope.difference>0&rat.comp$p.slope<=0.05))/nsim,
              length(which(rat.comp$slope.difference<0&rat.comp$p.slope<=0.05))/nsim)->p.comp.rat

            names(p.comp.phen)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
            names(p.comp.rat)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
          }
          list(p.comp.phen,p.comp.rat)->p.comp
          names(p.comp)<-c("phenotype","rates")

        }
      }else{
        p.comp<-NULL
      }

      if(length(y)>Ntip(tree)){
        p.phen.node<-list()
        p.rate.node<-list()
      }else{
        p.phen.node<-matrix(ncol=4,nrow=length(trend.node))
        p.rate.node<-matrix(ncol=4,nrow=length(trend.node))
      }

      for(i in 1:length(trend.node)){ ### Results nodes ####
        if(length(y)>Ntip(tree)){#### Multivariate ####
          p.phen.node.y<-matrix(ncol=4,nrow=(ncol(y)+1))
          p.rate.node.y<-matrix(ncol=4,nrow=(ncol(y)+1))
          for(j in 1:(ncol(y)+1)){
            as.data.frame(do.call(rbind,lapply(lapply(phen.node,"[[",i),function(x) x[j,])))->pnod
            as.data.frame(do.call(rbind,lapply(lapply(rat.node,"[[",i),function(x) x[j,])))->rnod

            c(length(which(pnod$slope>0&pnod$p.slope<=0.05))/nsim,
              length(which(pnod$slope<0&pnod$p.slope<=0.05))/nsim,
              length(which(pnod$emm.difference>0&pnod$p.emm<=0.05))/nsim,
              length(which(pnod$emm.difference<0&pnod$p.emm<=0.05))/nsim)->p.phen.node.y[j,]

            c(length(which(rnod$emm.difference>0&rnod$p.emm<=0.05))/nsim,
              length(which(rnod$emm.difference<0&rnod$p.emm<=0.05))/nsim,
              length(which(rnod$slope.difference>0&rnod$p.slope<=0.05))/nsim,
              length(which(rnod$slope.difference<0&rnod$p.slope<=0.05))/nsim)->p.rate.node.y[j,]
          }

          colnames(p.phen.node.y)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
          colnames(p.rate.node.y)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
          if(!is.null(colnames(y))) rownames(p.phen.node.y)<-rownames(p.rate.node.y)<-c(colnames(y),"multiple") else
            rownames(p.phen.node.y)<-rownames(p.rate.node.y)<-c(sapply(1:ncol(y),function(x) paste("y",x,sep="")),"multiple")
          p.phen.node.y->p.phen.node[[i]]
          p.rate.node.y->p.rate.node[[i]]
        }else{#### Univariate ####
          as.data.frame(do.call(rbind,lapply(phen.node,"[[",i)))->pnod
          rownames(pnod)<-seq(1:nsim)
          as.data.frame(do.call(rbind,lapply(rat.node,"[[",i)))->rnod

          c(length(which(pnod$slope>0&pnod$p.slope<=0.05))/nsim,
            length(which(pnod$slope<0&pnod$p.slope<=0.05))/nsim,
            length(which(pnod$emm.difference>0&pnod$p.emm<=0.05))/nsim,
            length(which(pnod$emm.difference<0&pnod$p.emm<=0.05))/nsim)->p.phen.node[i,]

          c(length(which(rnod$emm.difference>0&rnod$p.emm<=0.05))/nsim,
            length(which(rnod$emm.difference<0&rnod$p.emm<=0.05))/nsim,
            length(which(rnod$slope.difference>0&rnod$p.slope<=0.05))/nsim,
            length(which(rnod$slope.difference<0&rnod$p.slope<=0.05))/nsim)->p.rate.node[i,]

        }
      }

      if(length(y)>Ntip(tree)){
        names(p.phen.node)<-names(p.rate.node)<-trend.node
      }else{
        colnames(p.phen.node)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
        colnames(p.rate.node)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
        rownames(p.phen.node)<-rownames(p.rate.node)<-trend.node
      }
      list(p.phen.node,p.rate.node)->p.trend.node
      names(p.trend.node)<-c("phenotype","rates")
      node.res<-p.trend.node
      if(length(trend.node)>1) node.res<-list(node=node.res,comparison=p.comp) else node.res<-list(node=node.res)
      trend.res<-do.call(c,list(tree=list(whole.tree.res),node.res))


    }else trend.res<-whole.tree.res

  }else trend.res<-NULL

  if(is.null(conv.node)==FALSE){
    matrix(c(length(which(lapply(lapply(SCcut,"[[",1),function(x) x[1,8])<=0.05))/nsim,
             length(which(lapply(lapply(SCcut,"[[",1),function(x) x[1,9])<=0.05))/nsim),ncol=2)->p.convC
    colnames(p.convC)<-colnames(SCcut[[1]][[1]])[8:9]
    rownames(p.convC)<-paste(conv.node,collapse="-")
  }else p.convC<-NULL

  if(is.null(conv.state)==FALSE){
    p.convS<-matrix(ncol=2,nrow=nrow(SCcutS[[1]]))
    if("nostate"%in%conv.state&length(unique(conv.state)[-which(unique(conv.state)=="nostate")])){
      c(length(which(sapply(SCcutS,function(x) x[1,3])<=0.05))/nsim,
        length(which(sapply(SCcutS,function(x) x[1,4])<=0.05))/nsim)->p.convS[1,]
      rownames(p.convS)<-rownames(SCcutS[[1]])
      colnames(p.convS)<-colnames(SCcutS[[1]])[3:4]
    }else{
      for(i in 1:nrow(SCcutS[[1]])){
        c(length(which(sapply(SCcutS,function(x) x[i,5])<=0.05))/nsim,
          length(which(sapply(SCcutS,function(x) x[i,6])<=0.05))/nsim)->p.convS[i,]
      }
      rownames(p.convS)<-apply(SCcutS[[1]][,1:2],1,function(x) paste(x[1],x[2],sep="-"))
      colnames(p.convS)<-colnames(SCcutS[[1]])[5:6]
    }

  }else p.convS<-NULL

  list(p.convC,p.convS)->conv.res
  names(conv.res)<-c("clade","state")

  res<-list(mean(real.s),root.conf.int,acefit,conv.res,shift.res,trend.res)
  names(res)<-c("mean.sampling","rootCI","ace.regressions","conv.results","shift.results","trend.results")

  return(res)
}
