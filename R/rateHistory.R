#' @title Computing historical rates
#' @usage rateHistory(RR)
#' @description The function calculates historical rates for each tip of the
#'   tree. Historical rates represent the sum of rates along subsequent branches
#'   of a lineage. The product of these rates multiplied by the relative branch
#'   lengths represents the phenotypic transformation from one node to the next
#'   one along the path. The function also provides the sum of rate modulus
#'   along lineages. This represents the total amount of evolutionary change per
#'   lineage.
#' @param RR an object fitted by the function \code{\link{RRphylo}}.
#' @export
#' @return The function returns the vector of net historical rates
#'   (\code{$rateHistory$net.rate}) and the sum of rate modulus for each tip
#'   (\code{$rateHistory$norm.rate}), and the matrix of phenotypic changes
#'   from one node to the next along each lineage (\code{phen.path}).
#' @author Pasquale Raia, Silvia Castiglione
#' @examples
#' \dontrun{
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=treecet,y=masscet,clus=cc)->RRcet
#'
#' rateHistory(RR=RRcet)->rh
#'
#'     }

rateHistory<-function(RR){
  RR$tip.path->L->L2
  RR$tree->tree
  RR$multiple.rate->rates

  allrate<-lapply(1:nrow(L), function(i){
    ar<-sapply(1:ncol(rates),function(k)
      rates[match(c(rownames(L)[i],getMommy(tree, rownames(L)[i])),rownames(rates)),k,drop=FALSE])
    rownames(ar)<-c(rownames(L)[i],getMommy(tree, rownames(L)[i]))
    ar
  })
  hrate<-do.call(rbind,lapply(allrate,colSums,na.rm=TRUE))
  norm2vec<-do.call(rbind,lapply(allrate,function(j) apply(j,2,unitV,na.rm=TRUE)))
  rownames(hrate)<-rownames(norm2vec)<-rownames(L)

  LLs<-lapply(1:ncol(rates),function(k){
    lapply(1:length(allrate), function(i){
      L2[i,match(rownames(allrate[[i]]),colnames(L))]<<-allrate[[i]][,k]
     })
      LL<-L2*L
      LL[,1]<-RR$aces[1,k]
      LL
  })

  class(LLs)<-"RRphyloList"

  return(list(rateHistory=list(net.rate=hrate,norm.rate=norm2vec),phen.path=LLs))
}
