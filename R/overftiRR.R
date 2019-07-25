#' @title Testing RRphylo methods overfit
#' @description Testing the robustness of \code{\link{search.trend}} (\cite{Castiglione et al. 2019}) and \code{\link{search.shift}} (\cite{Castiglione et al. 2018}) to sampling effects and phylogenetic uncertainty.
#' @usage overfitRR(RR,y,s=0.25,trend=FALSE,shift.node=NULL,shift.state=NULL,
#' trend.node=NULL,aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,rootV=NULL,nsim=100,clus=.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a named vector of phenotypes.
#' @param s the percentage of tips to be cut off. It is set at 25\% by default.
#' @param trend if \code{TRUE}, the function performs \code{search.trend} on the subsampled tree, searching for trends in either phenoypes or rates through time. If the argument \code{trend.node} is specified (see below), the \code{trend} argument is automatically set to \code{TRUE}.
#' @param shift.node the nodes to be tested in \code{search.shift} under the "clade" condition. Notice that the conditions on \code{search.shift} are exclusive, they cannot be tested together.
#' @param shift.state the named vector of states to be tested in \code{search.shift} under the "sparse" condition. Notice that the conditions on \code{search.shift} are exclusive, they cannot be tested together.
#' @param trend.node the nodes to be tested by \code{search.trend} for trends in phenotypes or rates through time.
#' @param aces if used to produce the \code{RR} object, the vector of those ancestral character values at nodes known in advance must be specified. Names correspond to the nodes in the tree.
#' @param x1 the additional predictor to be specified if the RR object has been created using an additional predictor (i.e. multiple version of \code{RRphylo}). \code{'x1'} vector must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the predictor as well, and taking the vector of ancestral states and tip values to form the \code{x1}.
#' @param aces.x1 a named vector of ancestral character values at nodes for \code{x1}. It must be indicated if the RR object has been created using both \code{aces} and \code{x1}. Names correspond to the nodes in the tree.
#' @param cov if used to produce the \code{RR} object, the covariate must be specified. As in \code{RRphylo}, the covariate vector must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate.
#' @param rootV if used to produce the \code{RR} object, the phenotypic value at the tree root must be specified.
#' @param nsim number of simulations to be performed. It is set at 100 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @return The function returns a 'list' containing:
#' @return \strong{$rootCI} the 95\% confidence interval around the root value.
#' @return \strong{$ace.regressions} the results of linear regression between ancestral state estimates before and after the subsampling.
#' @return \strong{$shift.results} if \code{shift.state} is specified, a list including for each state the p-values at each simulation ($states) and the percentage of simulations producing significant p-value ($p.states). If the argument \code{shift.node} is specified, the list contains for each node the p-values at each simulation ($single clades) and the percentage of simulations producing significant p-value ($p.single.clades), and the same figures by considering all the specified nodes as evolving under a single rate ($all.caldes.togheter and $p.all.clades.together).
#' @return \strong{$trend.results} if \code{trend = TRUE}, a list including slopes and p-values of phenotypes versus age and absolute rates versus age regressions for the entire tree at each simulation ($tree.phenotype and $tree.rates) and the percentage of simulations showing significant p-values for both variables separated by slope sign ($tree.p). If \code{trend.node} is specified, the list also includes the same results at nodes ($node.phenotype, $node.rates and $node.p) and the results for comparison between nodes ($comparison.phenotype, $comparison.rates and $comparison.p).
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @details Methods using a large number of parameters risk being overfit. This usually translates in poor fitting with data and trees other than the those originally used. With \code{RRphylo} methods this risk is usually very low. However, the user can assess how robust the results got by applying \code{search.shift} or \code{search.trend} are by running \code{overfitRR}. With the latter, the original tree and data are subsampled by specifying a \code{s} parameter, that is the proportion of tips to be removed from the tree. Internally, \code{overfitRR} further shuffles the tree by using the function \code{\link{swapONE}}. Thereby, both the potential for overfit and phylogenetic uncertainty are accounted for straight away.
#' @export
#' @importFrom rlist list.append
#' @importFrom ddpcr quiet
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P. (2018). A new method for testing evolutionary rate variation and shifts in phenotypic evolution. \emph{Methods in Ecology and Evolution}, 9: 974-983.doi:10.1111/2041-210X.12954
#' @references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M., Profico, A., Girardi, G., & Raia, P. (2019) Simultaneous detection of macroevolutionary patterns in phenotypic means and rate of change with and within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14: e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' # Extract Pterosaurs tree and data
#'   library(ape)
#'   extract.clade(treedino,748)->treeptero
#'   massdino[match(treeptero$tip.label,names(massdino))]->massptero
#'   massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' \donttest{
#' RRphylo(tree=treedino,y=massdino)->dinoRates
#' RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'
#' # Case 1 search.shift under "clade" condition
#'     search.shift(RR=dinoRates, status.type= "clade",foldername=tempdir())->SSnode
#'     overfitRR(RR=dinoRates,y=massdino,shift.node=rownames(SSnode$single.clades))->overfit.SSnode
#'
#' # Case 2 search.shift under "sparse" condition
#'     search.shift(dinoRates, status.type= "sparse", state=statedino,
#'     foldername=tempdir())->SSstate
#'     overfitRR(RR=dinoRates,y=massdino,shift.state=statedino)->overfit.SSstate
#'
#' # Case 3 search.trend on the entire tree
#'     search.trend(RR=RRptero, y=log(massptero),foldername=tempdir())->STtree
#'     overfitRR(RR=RRptero,y=log(massdino),trend=TRUE)->overfit.STtree
#'
#' # Case 4 search.trend at specified nodes
#'     search.trend(RR=RRptero, y=log(massptero),node=143,foldername=tempdir())->STnode
#'     overfitRR(RR=RRptero,y=log(massdino),trend.node=143)->overfit.STnode
#'
#' # Case 5 overfitRR on multiple RRphylo
#'   data("DataCetaceans")
#'   DataCetaceans$treecet->treecet
#'   DataCetaceans$masscet->masscet
#'   DataCetaceans$brainmasscet->brainmasscet
#'   DataCetaceans$aceMyst->aceMyst
#'
#'   ape::drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#'   masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#'   RRphylo(treecet.multi,masscet.multi)->RRmass.multi
#'   RRmass.multi$aces[,1]->acemass.multi
#'   c(acemass.multi,masscet.multi)->x1.mass
#'
#'   RRphylo(treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#'   search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=0.5,foldername=tempdir())->STcet
#'   overfitRR(RR=RRmulti,y=brainmasscet,trend=TRUE,x1=x1.mass)->overfit.STcet
#'
#' }

overfitRR<-function(RR,y,s=0.25,trend=FALSE,shift.node=NULL,shift.state=NULL,
                    trend.node=NULL,aces=NULL,x1=NULL,aces.x1=NULL,cov=NULL,rootV=NULL,nsim=100,
                    clus=.5)
{
  # require(phytools)
  # require(geiger)
  # require(ddpcr)
  # require(rlist)

  RR$tree->tree
  RR$aces[,1]->y.ace
  tree$node.label<-names(y.ace)
  y[match(tree$tip.label,names(y))]->y

  rootlist<-array()
  acefit<-SScut<-STcut<-list()
  for(k in 1:nsim){
    '%ni%' <- Negate('%in%')


    repeat({
      suppressWarnings(swapONE(tree,0.1,0.1)[[1]])->tree.swap
      sample(y,round((Ntip(tree)*s),0))->offs
      sapply(trend.node,function(x) length(tips(tree,x))-length(which(names(offs)%in%tips(tree,x))))->lt
      sapply(shift.node,function(x) length(tips(tree,x))-length(which(names(offs)%in%tips(tree,x))))->lss
      if(length(which(lt<5))==0&length(which(lss<5))==0) break
    })

    y[match(tree.swap$tip.label,names(y))]->y

    if(s>0){
      tree.swap$edge[tree.swap$edge[,1]==(Ntip(tree.swap)+1),2]->rootdesc
      if(length(which(rootdesc<(Ntip(tree.swap)+1)))>0) tree.swap$tip.label[rootdesc[which(rootdesc<Ntip(tree.swap)+1)]]->saver else saver="xx"
      if(saver%in%names(offs)) offs[-which(names(offs)==saver)]->offs

      y[-which(names(y)%in%names(offs))]->ycut
      drop.tip(tree.swap,which(names(y)%ni%names(ycut)))->treecut
      y.ace[which(names(y.ace)%in%treecut$node.label)]->y.acecut
    }else{
      y->ycut
      tree.swap->treecut
      y.ace->y.acecut
    }

    if(is.null(cov)==FALSE) {
      cov[match(c(tree.swap$node.label,tree.swap$tip.label), names(cov))]->cov
      cov[match(c(names(y.acecut),names(ycut)),names(cov))]->covcut
      names(cov)[1:Nnode(treecut)]<-seq(1:Nnode(treecut))
    }else covcut<-NULL

    if(is.null(x1)==FALSE) {
      x1[match(c(tree.swap$node.label,tree.swap$tip.label), names(x1))]->x1
      x1[match(c(names(y.acecut),names(ycut)),names(x1))]->x1cut
      names(x1cut)[1:Nnode(treecut)]<-seq(1:Nnode(treecut))
    }else x1cut<-NULL

    if(is.null(aces)==FALSE){
      aces->acescut
      drop<-c()
      for(i in 1:length(aces)) {
        if(length(which(tips(tree,names(aces[i]))%in%treecut$tip.label))>1)
          getMRCA(treecut,tips(tree,names(aces[i]))[which(tips(tree,names(aces[i]))%in%treecut$tip.label)])->names(acescut)[i] else
            c(drop,i)->drop
      }
      if(length(drop>0)) acescut[-drop]->acescut

    }else{
      acescut<-NULL
    }

    if(is.null(aces.x1)==FALSE){
      aces.x1->aces.x1cut
      drop<-c()
      for(i in 1:length(aces.x1)) {
        if(length(which(tips(tree,names(aces.x1[i]))%in%treecut$tip.label))>1)
          getMRCA(treecut,tips(tree,names(aces.x1[i]))[which(tips(tree,names(aces.x1[i]))%in%treecut$tip.label)])->names(aces.x1cut)[i] else
            c(drop,i)->drop
      }
      if(length(drop>0)) aces.x1cut[-drop]->aces.x1cut

    }else{
      aces.x1cut<-NULL
    }

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
      shift.state[match(names(ycut),names(shift.state))]->shift.state.cut
      }

    if(is.null(rootV)==FALSE) rootV->rootVcut else rootVcut<-NULL

    RRphylo(treecut,ycut,aces=acescut,x1=x1cut,aces.x1=aces.x1cut,cov=covcut,rootV = rootVcut,clus=clus)->RRcut
    if(trend|is.null(trend.node)==FALSE) quiet(search.trend(RRcut,ycut,x1=x1cut,node=trend.node.cut,foldername=tempdir(),cov=covcut,clus=clus)->stcut->STcut[[k]],all=TRUE)
    if(is.null(shift.node)==FALSE) quiet(search.shift(RRcut,status.type="clade",node=shift.node.cut,foldername=tempdir())->sscut->SScut[[k]],all=TRUE)
    if(is.null(shift.state)==FALSE) quiet(search.shift(RRcut,status.type="sparse",state=shift.state.cut,foldername=tempdir())->sscut->SScut[[k]],all=TRUE)


    RRcut$aces[1,]->rootlist[k]
    # ycut->ylist[[k]]
    summary(lm(y.acecut~RRcut$aces))->acefit[[k]]
  }


  quantile(rootlist,c(0.025,0.975))->CIroot
  data.frame(root=y.ace[1],"CI 2.5"=CIroot[1],"CI 97.5"=CIroot[2])->root.conf.int

  if(is.null(shift.node)==FALSE){
    if(length(SScut[[1]])>2){
      do.call(rbind,lapply(lapply(SScut,"[[",1),function(x) x[,2]))-> p.ran.whole
      apply(p.ran.whole,2,function(x) length(which(x<=0.025|x>=0.975)))/nsim->p.shift.whole
      do.call(rbind,lapply(lapply(SScut,"[[",2),function(x) x[,2]))-> p.ran
      apply(p.ran,2,function(x) length(which(x<=0.025|x>=0.975)))/nsim->p.shift
      names(p.shift)<-colnames(p.ran)<-shift.node
      rownames(p.ran)<-rownames(p.ran.whole)<-seq(1:nsim)
      colnames(p.ran.whole)<-"shift"
      list(p.ran,p.shift,p.ran.whole,p.shift.whole)->shift.res
      names(shift.res)<-c("single.clades","p.single.clades",
                          "all.clades.together","p.all.clades.together")

    }else{
      #do.call(rbind,lapply(lapply(SScut,"[[",2),function(x) x[2]))-> p.ran
      do.call(rbind,lapply(lapply(SScut,"[[",1),function(x) x[2]))-> p.ran
      apply(p.ran,2,function(x) length(which(x<=0.025|x>=0.975)))/nsim->p.shift
      names(p.shift)<-colnames(p.ran)<-shift.node
      rownames(p.ran)<-seq(1:nsim)
      list(p.ran,p.shift)->shift.res
      names(shift.res)<-c("single.clades","p.single.clades")
    }
  }else{
    if(is.null(shift.state)==FALSE){
      p.shift<-array()
      p.ran<-matrix(nrow=nsim,ncol=nrow(SScut[[1]][[1]]))
      for(i in 1:nrow(SScut[[1]][[1]])){
        unlist(lapply(lapply(SScut,"[[",1),function(x) x[i,2]))->pr->p.ran[,i]
        length(which(pr<=0.025|pr>=0.975))/nsim->p.shift[i]
      }
      names(p.shift)<-colnames(p.ran)<-rownames(SScut[[1]][[1]])
      rownames(p.ran)<-seq(1:nsim)
      list(p.ran,p.shift)->shift.res
      names(shift.res)<-c("states","p.states")

    }else shift.res<-NULL
  }

  if(trend|is.null(trend.node)==FALSE){
    do.call(rbind,lapply(STcut,"[[",3))[,c(1,3)]->phen.ran
    do.call(rbind,lapply(STcut,"[[",4))[,c(1,3)]->rat.ran


    rbind(c(length(which(phen.ran$slope>0&phen.ran$p.random<=0.05))/nsim,
      length(which(phen.ran$slope<0&phen.ran$p.random<=0.05))/nsim),
      c(length(which(rat.ran$slope>0&rat.ran$p.random<=0.05))/nsim,
        length(which(rat.ran$slope<0&rat.ran$p.random<=0.05))/nsim))->p.trend

    colnames(p.trend)<-c("p.slope+","p.slope-")
    rownames(p.trend)<-c("phenotype","rates")

    # c(length(which(phen.ran[,2]<=0.05))/nsim,length(which(rat.ran[,2]<=0.05))/nsim)->p.trend
    # names(p.trend)<-c("phenotype","rates")
    rownames(phen.ran)<-rownames(rat.ran)<-seq(1:nsim)
    list(phen.ran,rat.ran,p.trend)->whole.tree.res
    names(whole.tree.res)<-c("tree.phenotype","tree.rates","tree.p")

    if(is.null(trend.node)==FALSE){
      lapply(STcut,"[[",5)->phen.node
      lapply(STcut,"[[",6)->rat.node

      if(length(STcut[[1]])==7){
        if(length(trend.node)>2) {
          phen.comp<-rat.comp<-nod.nam<-list()
          p.comp<-matrix(ncol=2,nrow=nrow(STcut[[1]][[7]][[1]]))
          for(k in 1:nrow(STcut[[1]][[7]][[1]])){
            do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",1),function(w) w[k,]))->pcomp
            pcomp[nsim,1:2]->nod.nam[[k]]
            pcomp[,3:6]->pcomp->phen.comp[[k]]
            do.call(rbind,lapply(lapply(lapply(STcut,"[[",7),"[[",2),function(w) w[k,]))[,3:6]->rcomp
            rcomp->rat.comp[[k]]
            c(length(which(pcomp[,2]<=0.05))/nsim,length(which(rcomp[,2]<=0.05))/nsim)->p.comp[k,]
          }
          colnames(p.comp)<-c("phenotype","rates")
          do.call(rbind, nod.nam)->nod.nam
          t(apply(nod.nam, 1, function(x) unlist(lapply(strsplit(x, "g"), "[[", 2))))->nam.pair
          t(apply(nam.pair,1,function(x) trend.node[match(x,trend.node.cut)]))->nam.pair
          names(phen.comp)<-names(rat.comp)<-rownames(p.comp)<-apply(nam.pair,1, function(x) paste(x[1], x[2], sep="-"))
        }else{
          do.call(rbind,lapply(lapply(STcut,"[[",7),"[[",1))[,3:6]->phen.comp
          do.call(rbind,lapply(lapply(STcut,"[[",7),"[[",2))[,3:6]->rat.comp
          c(length(which(phen.comp[,2]<=0.05))/nsim,length(which(rat.comp[,2]<=0.05))/nsim)->p.comp
          names(p.comp)<-c("phenotype","rates")
        }
      }else{
        phen.comp<-rat.comp<-p.comp<-NULL
      }

      phen.node.ran<-rat.node.ran<-list()
      #p.trend.node<-matrix(ncol=2,nrow=length(trend.node))
      p.phen.node<-matrix(ncol=4,nrow=length(trend.node))
      p.rate.node<-matrix(ncol=4,nrow=length(trend.node))
      for(i in 1:length(trend.node)){
        do.call(rbind,lapply(phen.node,"[[",i))->pnod
        rownames(pnod)<-seq(1:nsim)
        pnod->phen.node.ran[[i]]
        do.call(rbind,lapply(rat.node,"[[",i))->rnod
        rnod->rat.node.ran[[i]]


        c(length(which(pnod$slope>0&pnod$p.slope<=0.05))/nsim,
          length(which(pnod$slope<0&pnod$p.slope<=0.05))/nsim,
          length(which(pnod$emm.difference>0&pnod$p.emm<=0.05))/nsim,
          length(which(pnod$emm.difference<0&pnod$p.emm<=0.05))/nsim)->p.phen.node[i,]

        c(length(which(rnod$emm.difference>0&rnod$p.emm<=0.05))/nsim,
          length(which(rnod$emm.difference<0&rnod$p.emm<=0.05))/nsim,
          length(which(rnod$slope.difference>0&rnod$p.slope<=0.05))/nsim,
          length(which(rnod$slope.difference<0&rnod$p.slope<=0.05))/nsim)->p.rate.node[i,]

      }
      colnames(p.phen.node)<-c("p.slope+","p.slope-","p.emm+","p.emm-")
      colnames(p.rate.node)<-c("p.emm+","p.emm-","p.slope+","p.slope-")
      rownames(p.phen.node)<-rownames(p.rate.node)<-names(phen.node.ran)<-names(rat.node.ran)<-trend.node
      list(p.phen.node,p.rate.node)->p.trend.node
      names(p.trend.node)<-c("phenotype","rates")
      #colnames(p.trend.node)<-c("phenotype","rates")
      #rownames(p.trend.node)<-names(phen.node.ran)<-names(rat.node.ran)<-trend.node
      node.res<-list("node.phenotype"=phen.node.ran,"node.rates"=rat.node.ran,"node.p"=p.trend.node)
      if(length(trend.node)>1) node.res<-list.append(node.res,"comparison.phenotype"=phen.comp,"comparison.rates"=rat.comp,"comparison.p"=p.comp)
      trend.res<-do.call(c,list(whole.tree.res,node.res))
    }else trend.res<-whole.tree.res

  }else trend.res<-NULL

  res<-list(root.conf.int,acefit,shift.res,trend.res)
  names(res)<-c("rootCI","ace.regressions","shift.results","trend.results")

  return(res)
}
