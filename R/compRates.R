#' @title Comparing average absolute rates between clades
#' @usage compRates(RR,node, nrep = 1000, cov = NULL)
#' @description The function \code{compRates} is an adaptation of
#'   \code{\link{search.shift}} which performs pairwise comparison of average
#'   absolute rates between clades via bootstrap.
#' @param RR an object fitted by the function \code{\link{RRphylo}}.
#' @param node the most recent common ancestors of clades to be tested. The
#'   nodes must be identified on the dicothomized version of the original tree
#'   returned by \code{\link{RRphylo}}. Pairwise comparison between all clades is
#'   performed.
#' @param nrep the number of simulations to be performed for the rate shift
#'   test, by default \code{nrep} is set at 1000.
#' @param cov the covariate vector to be indicated if its effect on rate values
#'   must be accounted for. Contrary to \code{\link{RRphylo}}, \code{cov} needs to be
#'   as long as the number of tips of the tree.
#' @export
#' @return For each node pair, the function returns the average absolute rate
#'   difference (computed as the difference between the average absolute rate
#'   over all branches subtended by the nodes) and related significance level.
#'   Probabilities are derived by contrasting the rate differences to simulated
#'   ones derived by shuffling the rates across the tree branches for a number
#'   of replicates specified by the argument \code{nrep}. Note that the p-values
#'   refer to the number of times the real differences are larger
#'   (p-value>=0.975) or smaller (p-value<=0.025) than the simulated ones,
#'   divided by the number of simulations, hence the test should be considered
#'   as two-tailed.
#'   The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#' @author Silvia Castiglione, Giorgia Girardi
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#' compRates(RR=dinoRates,node=c(696,746))->cr1
#' compRates(RR=dinoRates,node=c(696,746),cov=massdino)->cr2
#'     }

compRates<-function(RR,node,nrep=1000,cov=NULL){
  # require(ape)
  # require(RRphylo)
  # require(phytools)
  funcall <- match.call()
  tree <- RR$tree
  rates <- RR$rates[,,drop=FALSE]
  betas<-RR$multiple.rates[,,drop=FALSE]

  if(!is.null(cov)){
    RRphylo(tree,cov,clus=0)->RRcova
    abs(c(RRcova$aces,cov))->Y
    c(rownames(RRcova$aces),names(cov))->names(Y)
    covRates(Y,betas)->betas

    if(ncol(betas)>1) rates <- as.matrix(apply(betas, 1, function(x) sqrt(sum(x^2)))) else rates<-betas
  }

  combn(node,2)->node.pairs
  diff<-pval<-c()
  pData<-list()
  for(i in 1:ncol(node.pairs)){
    node.pairs[,i]->nns
    nnrates<-lapply(nns,function(k) rates[c(getDescendants(tree,k)[which(getDescendants(tree,k)>Ntip(tree))],tips(tree,k)),])
    lens<-sapply(nnrates,length)
    diff[i]<-realdiff<-diff(sapply(nnrates,function(x) mean(abs(x))))
    pData[[i]]<-randiff<-replicate(nrep,diff(sapply(lens,function(x) mean(abs(sample(unlist(nnrates),x))))))
    pval[i]<-rank(c(realdiff, randiff[-nrep]))[1]/nrep
  }
  data.frame(pair=apply(node.pairs,2,function(x) paste(rev(x),collapse="-")),
             rate.difference=diff,p.value=pval)->res
  plotData<-do.call(cbind,pData)
  colnames(plotData)<-res$pair

  res<-list(results=res,plotData=plotData)
  class(res)<-c("RRphyloList","list")
  attr(res,"hidden")<-"plotData"
  attr(res,"Call")<-funcall

  return(res)
}
