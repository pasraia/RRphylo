#' @title Testing search.shift overfit
#' @description Testing the robustness of \code{\link{search.shift}}
#'   (\cite{Castiglione et al. 2018}) results to sampling effects and
#'   phylogenetic uncertainty.
#' @usage overfitSS(RR,oveRR,node=NULL,state=NULL)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param oveRR an object produced by applying \code{\link{overfitRR}} on the
#'   object provided to the function as \code{RR}.
#' @param node,state arguments passed to \code{\link{search.shift}}. Arguments
#'   \code{node} and \code{state} can be specified at the same time.
#' @return The function returns a 'RRphyloList' object containing:
#' @return \strong{$SSclade.list} a 'RRphyloList' including the results of each
#'   \code{search.shift - clade condition} performed within \code{overfitSS}.
#' @return \strong{$SSsparse.list} a 'RRphyloList' including the results of each
#'   \code{search.shift - sparse condition} performed within \code{overfitSS}.
#' @return \strong{$shift.results} a list including results for
#'   \code{\link{search.shift}} performed under \code{clade} and \code{sparse}
#'   conditions. If one or more \code{node}s are specified, the
#'   \code{$clade$single.clades} object contains the proportion of simulations
#'   producing significant and positive or significant and negative rate shifts
#'   for each single node, either compared to the rest of the tree
#'   (\code{$singles}) or to the rest of the tree after removing other shifting
#'   clades (\code{$no.others}). The object \code{$clade$all.clades.together}
#'   includes the same proportions obtained by testing all the specified clades
#'   as a whole (i.e. considering them as evolving under a single rate regime).
#'   For each node the proportion of tested trees (i.e. where the clade identity
#'   was preserved) is also indicated. If a \code{state} vector is supplied, the
#'   object \code{$sparse} contains the percentage of simulations producing
#'   significant p-value separated by shift sign (\code{$p.states}).
#' @return The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Carmela Serio, Giorgia Girardi, Pasquale Raia
#' @export
#' @seealso \href{../doc/overfitRR.html}{\code{overfitSS} vignette} ;
#'   \href{../doc/search.shift.html}{\code{search.shift} vignette} ;
#'   \href{../doc/Alternative-trees.html}{\code{Alternative-trees} vignette}
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P. (2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @examples
#' \dontrun{
#' cc<- 2/parallel::detectCores()
#' # load the RRphylo example dataset including Ornithodirans tree and data
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' log(DataOrnithodirans$massdino)->massdino
#' DataOrnithodirans$statedino->statedino
#'
#' # peform RRphylo on Ornithodirans tree and data
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#'
#' # perform search.shift under both "clade" and "sparse" condition
#' search.shift(RR=dinoRates, status.type= "clade")->SSauto
#' search.shift(RR=dinoRates, status.type= "sparse", state=statedino)->SSstate
#'
#' ## overfitSS routine
#' # generate a list of subsampled and swapped phylogenies, setting as categories/node
#' # the state/node under testing
#' treedino.list<-resampleTree(dinoRates$tree,s = 0.25,categories=statedino,
#'                         node=rownames(SSauto$single.clades),swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#'
#' # test the robustness of search.shift
#' ofRRdino<-overfitRR(RR = dinoRates,y=massdino,phylo.list=treedino.list,clus=cc)
#' ofSS<-overfitSS(RR = dinoRates,oveRR = ofRRdino,state=statedino,node=rownames(SSauto$single.clades))
#'
#' }


overfitSS<-function(RR,oveRR,node=NULL,state=NULL)
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
  oveRR$RR.list->RR.list
  phylo.size<-Ntip(RR.list[[1]]$tree)
  length(RR.list)->nsim

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  SScut<-SScutS<-node.match<-list()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)

    RR.list[[k]]->RRcut
    RRcut$tree->treecut

    if(!is.null(node)){
      node.cut<-array()
      for(i in 1:length(node)){
        getMRCA(treecut,tips(tree,node[i])[which(tips(tree,node[i])%in%treecut$tip.label)])->shN
        if(phylo.size==Ntip(tree)){
          length(tips(treecut,shN))/length(tips(tree,node[i]))->sh.tips
          if(sh.tips<=1.1&sh.tips>=0.9) shN->node.cut[i] else NA->node.cut[i]
        } else shN->node.cut[i]
      }
      data.frame(node,node.cut)->node.match[[k]]
      node.cut[which(!is.na(node.cut))]->node.cut
      if(length(node.cut)==0) node.cut<-NULL


      if(!is.null(node.cut))
        ddpcr::quiet(search.shift(RRcut,status.type="clade",node=node.cut)->sscut->SScut[[k]],all=TRUE)
    }
    if(!is.null(state)){
      state[match(treecut$tip.label,names(state))]->state.cut
      ddpcr::quiet(search.shift(RRcut,status.type="sparse",state=state.cut)->sscut->SScutS[[k]],all=TRUE)
    }
  }
  close(pb)

  if(!is.null(node)){
    mapply(a=node.match,b=SScut,function(a,b){
      if(!is.null(b)){
        if(length(b)>1){
          a[match(rownames(b$single.clades$singles),a[,2]),1]->rownames(b$single.clades$singles)
          b$single.clade$singles
        }else{
          a[match(rownames(b$single.clades),a[,2]),1]->rownames(b$single.clades)
          b$single.clade
        }
      }
    },SIMPLIFY = FALSE)->singles

    t(sapply(node,function(k){
      t(sapply(singles,function(j) {
        if(any(as.numeric(rownames(j))==k)) j[which(as.numeric(rownames(j))==k),] else c(NA,NA)
      }))->pran
      pran[which(!is.na(pran[,1])),,drop=FALSE]->pran
      cbind(length(which(pran[,2]>=0.975))/nrow(pran),length(which(pran[,2]<=0.025))/nrow(pran),nrow(pran)/nsim)
    }))->shift.res.clade
    rownames(shift.res.clade)<-node
    colnames(shift.res.clade)<-c("p.shift+","p.shift-","tested.trees")


    mapply(a=node.match,b=SScut,function(a,b){
      if(!is.null(b)){
        if(length(b)>1){
          a[match(rownames(b$single.clades$no.others),a[,2]),1]->rownames(b$single.clades$no.others)
          b$single.clade$no.others
        }else NULL
      }
    },SIMPLIFY = FALSE)->noothers


    t(sapply(node,function(k){
      t(sapply(noothers,function(j) {
        if(any(as.numeric(rownames(j))==k)) j[which(as.numeric(rownames(j))==k),] else c(NA,NA)
      }))->pran
      pran[which(!is.na(pran[,1])),]->pran
      cbind(length(which(pran[,2]>=0.975))/nrow(pran),length(which(pran[,2]<=0.025))/nrow(pran),nrow(pran)/nsim)
    }))->shift.res.noot

    rownames(shift.res.noot)<-node
    colnames(shift.res.noot)<-c("p.shift+","p.shift-","tested.trees")

    lapply(SScut,function(j) j$all.clades)->allcla

    if(!all(sapply(allcla,is.null))){
      as.data.frame(do.call(rbind, allcla))->allcla
      rbind(cbind(length(which(allcla$p.value>=0.975))/nrow(allcla),
                  length(which(allcla$p.value<=0.025))/nrow(allcla),nrow(allcla)/nsim))->shift.res.allcla
      rownames(shift.res.allcla)<-"all.clades"
      colnames(shift.res.allcla)<-c("p.shift+","p.shift-","tested.trees")

      shift.res.clade<-list(all.clades.together=shift.res.allcla,single.clades=list(singles=shift.res.clade,no.others=shift.res.noot))

    }

  }else SScut<-shift.res.clade<-NULL

  if(!is.null(state)){
    p.shift<-matrix(ncol=2,nrow=nrow(SScutS[[1]][[1]]))
    for(i in 1:nrow(SScutS[[1]][[1]])){
      unlist(lapply(lapply(SScutS,"[[",1),function(x) x[i,2]))->pr
      c(length(which(pr>=0.975))/nsim,length(which(pr<=0.025))/nsim)->p.shift[i,]
    }
    rownames(p.shift)<-rownames(SScutS[[1]][[1]])
    colnames(p.shift)<-c("p.shift+","p.shift-")
    p.shift->shift.res.state

  }else SScutS<-shift.res.state<-NULL

  list(shift.res.clade,shift.res.state)->shift.res
  names(shift.res)<-c("clade","sparse")

  if(!is.null(SScut)) class(SScut)<-"RRphyloList"
  if(!is.null(SScutS)) class(SScutS)<-"RRphyloList"

  res<-structure(list(SSclade.list=SScut,
                      SSsparse.list=SScutS,
                      shift.results=shift.res),
                 class = "RRphyloList")
  attr(res,"Call")<-funcall

  res
}
