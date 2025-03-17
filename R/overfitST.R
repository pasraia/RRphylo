#' @title Testing search.trend overfit
#' @description Testing the robustness of \code{\link{search.trend}}
#'   (\cite{Castiglione et al. 2019a}) results to sampling effects and
#'   phylogenetic uncertainty.
#' @usage
#' overfitST(RR,y,oveRR,x1=NULL,x1.residuals=FALSE,node=NULL,cov=NULL,clus=0.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a named vector of phenotypes.
#' @param oveRR an object produced by applying \code{\link{overfitRR}} on the
#'   object provided to the function as \code{RR}.
#' @param x1,node,cov arguments as passed to \code{\link{overfitRR}}.
#' @param x1.residuals as passed to \code{\link{search.trend}}. Default is
#'   \code{FALSE}.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{overfitST} set \code{clus} = 0.
#' @return The function returns a 'RRphyloList' object containing:
#' @return \strong{$ST.list} a 'RRphyloList' including the results of each
#'   \code{\link{search.trend}} performed within \code{overfitST}.
#' @return \strong{$trend.results} a list including the percentage of
#'   simulations showing significant p-values for phenotypes versus age and
#'   absolute rates versus age regressions for the entire tree separated by
#'   slope sign ($tree). If one or more nodes are specified within
#'   \code{trend.args}, the list also includes the same results at nodes ($node)
#'   and the results for comparison between nodes ($comparison). For each node
#'   the proportion of tested trees (i.e. where the clade identity was
#'   preserved; always 1 if no \code{phylo.list} is supplied) is also indicated.
#' @return The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Carmela Serio, Giorgia Girardi, Pasquale Raia
#' @export
#' @seealso \href{../doc/overfitRR.html}{\code{overfitST} vignette} ;
#'   \href{../doc/search.trend.html}{\code{search.trend} vignette} ;
#'   \href{../doc/Alternative-trees.html}{\code{Alternative-trees} vignette}
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
#' @examples
#' \dontrun{
#' cc<- 2/parallel::detectCores()
#' library(ape)
#'
#' ## Case 1
#' # load the RRphylo example dataset including Ornithodirans tree and data
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#'
#' # extract Pterosaurs tree and data
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' # perform RRphylo and search.trend on body mass data
#' RRphylo(tree=treeptero,y=log(massptero),clus=cc)->RRptero
#' search.trend(RR=RRptero, y=log(massptero),node=143,clus=cc)->st2
#'
#' ## overfitST routine
#' # generate a list of subsampled and swapped phylogenies setting as node
#' # the clade under testing
#' treeptero.list<-resampleTree(RRptero$tree,s = 0.25,node=143,
#'                              swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#'
#' # test the robustness of search.trend
#' ofRRptero<-overfitRR(RR = RRptero,y=log(massptero),phylo.list=treeptero.list,clus=cc)
#' ofSTptero<-overfitST(RR=RRptero,y=log(massptero),oveRR=ofRRptero,node=143,clus=cc)
#'
#'
#' ## Case 2
#' # load the RRphylo example dataset including Cetaceans tree and data
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#'
#' # cross-reference the phylogenetic tree and body and brain mass data. Remove from
#' # both the tree and vector of body sizes the species whose brain size is missing
#' drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),
#'                                                treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' # peform RRphylo on the variable (body mass) to be used as additional predictor
#' RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#'
#' # create the predictor vector: retrieve the ancestral character estimates
#' # of body size at internal nodes from the RR object ($aces) and collate them
#' # to the vector of species' body sizes to create
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' # peform RRphylo and search.trend on brain mass by using body mass as additional predictor
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc)->STcet
#'
#' ## overfitST routine
#' # generate a list of subsampled and swapped phylogenies to test
#' treecet.list<-resampleTree(RRmulti$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#'
#' # test the robustness of search.trend with and without x1.residuals
#' ofRRcet<-overfitRR(RR = RRmulti,y=brainmasscet,phylo.list=treecet.list,clus=cc,x1 =x1.mass)
#' ofSTcet1<-overfitST(RR=RRmulti,y=brainmasscet,oveRR=ofRRcet,x1 =x1.mass,clus=cc)
#' ofSTcet2<-overfitST(RR=RRmulti,y=brainmasscet,oveRR=ofRRcet,x1 =x1.mass,x1.residuals = TRUE,clus=cc)
#' }

overfitST<-function(RR,y,oveRR,
                    x1=NULL,x1.residuals=FALSE,
                    node=NULL,cov=NULL,
                    clus=0.5)
{
  # require(phytools)
  # require(ddpcr)
  # require(rlist)

  if (!requireNamespace("ddpcr", quietly = TRUE)) {
    stop("Package \"ddpcr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  funcall<-match.call()
  RR$tree->tree
  y <- treedataMatch(tree, y)[[1]]
  oveRR$RR.list->RR.list
  phylo.size<-Ntip(RR.list[[1]]$tree)
  length(RR.list)->nsim

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  STcut<-node.match<-list()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)

    RR.list[[k]]->RRcut
    RRcut$tree->treecut
    y[match(treecut$tip.label,rownames(y)),,drop=FALSE]->ycut

    if(!is.null(cov)){
      treedataMatch(treecut,cov)$y->covcut
      c(RRphylo(treecut,covcut,clus=clus)$aces[,1],covcut)->covcut
    }else covcut<-NULL

    if(!is.null(x1)) {
      as.matrix(x1)->x1
      treedataMatch(treecut,x1)$y->x1cut
      rbind(RRphylo(treecut,x1cut,clus=clus)$aces,x1cut)->x1cut
    }else x1cut<-NULL

    if(!is.null(node)){
      node.cut<-array()
      for(i in 1:length(node)) {
        getMRCA(treecut,tips(tree,node[i])[which(tips(tree,node[i])%in%treecut$tip.label)])->trN
        if(phylo.size==Ntip(tree)){
          length(tips(treecut,trN))/length(tips(tree,node[i]))->sh.tips
          if(sh.tips<=1.1&sh.tips>=0.9) trN->node.cut[i] else NA->node.cut[i]
        } else trN->node.cut[i]
      }
      data.frame(node,node.cut)->node.match[[k]]

      node.cut[which(!is.na(node.cut))]->node.cut
      if(length(node.cut)==0) node.cut<-NULL
    }else node.cut<-NULL


    ddpcr::quiet(search.trend(RRcut,ycut,x1=x1cut,x1.residuals = x1.residuals,node=node.cut,cov=covcut,clus=clus)->stcut->STcut[[k]],all=TRUE)

  }
  close(pb)

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

  if(!is.null(node)){
    mapply(a=node.match,b=lapply(STcut,"[[",4),function(a,b){
      a[match(names(b),a[,2]),1]->names(b)
      b
    },SIMPLIFY = FALSE)->phen.node

    mapply(a=node.match,b=lapply(STcut,"[[",5),function(a,b){
      a[match(names(b),a[,2]),1]->names(b)
      b
    },SIMPLIFY = FALSE)->rat.node

    p.phen.node<-list()
    p.rate.node<-list()
    for(k in 1:length(node)){
      lapply(phen.node,function(j) {
        if(any(as.numeric(names(j))==node[k])) j[[which(as.numeric(names(j))==node[k])]] else NA
      })->pran
      pran[which(!sapply(pran,function(w) all(is.na(w))))]->phen.pran

      lapply(rat.node,function(j) {
        if(any(as.numeric(names(j))==node[k])) j[[which(as.numeric(names(j))==node[k])]] else NA
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
    names(p.phen.node)<-names(p.rate.node)<-node
    list(p.phen.node,p.rate.node)->p.node
    names(p.node)<-c("phenotype","rates")
    node.res<-p.node

    if(length(node)>1){ #### Node comparison ####

      apply(combn(node,2),2,function(j) c(paste(j[1],j[2],sep="-"),paste(j[2],j[1],sep="-")))->tn.pairs

      lapply(STcut,function(j) j$group.comparison)->comptot

      mapply(x=lapply(comptot,"[[",1)[!sapply(comptot,is.null)],
             xx=node.match[!sapply(comptot,is.null)],function(x,xx){
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
                   # data.frame(x[which(revcols==2),c(2,1)],-1*x[which(revcols==2),3],
                   #            x[which(revcols==2),c(4,6,5,7)])->x[which(revcols==2),]
                   x[which(revcols==2),]<-data.frame(x[which(revcols==2),c(2,1)],
                                                     x[which(revcols==2),c(4,3)],
                                                     x[which(revcols==2),5],
                                                     -1*x[which(revcols==2),6],
                                                     x[which(revcols==2),7])
                 }
                 data.frame(x,pair=apply(x[,1:2],1,function(jk) paste(jk,collapse = "-")))
               }
             },SIMPLIFY = FALSE)->pcomptot

      mapply(x=lapply(comptot,"[[",2)[!sapply(comptot,is.null)],
             xx=node.match[!sapply(comptot,is.null)],function(x,xx){
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

    if(length(node)>1) node.res<-list(node=node.res,comparison=p.comp) else node.res<-list(node=node.res)
    trend.res<-do.call(c,list(tree=list(whole.tree.res),node.res))
  }else trend.res<-whole.tree.res

  if(!is.null(STcut)) class(STcut)<-"RRphyloList"

  res<-structure(list(ST.list=STcut,
                      trend.results=trend.res),
                 class = "RRphyloList")
  attr(res,"Call")<-funcall

  res
}
