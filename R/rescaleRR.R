#' @title Rescaling phylogenetic trees
#' @description The function rescales all branches and leaves of the
#'   phylogenetic tree.
#' @usage
#' rescaleRR(tree,RR=NULL,height=NULL,trend=NULL,delta=NULL,kappa=NULL,lambda=NULL)
#'
#' @param tree the phylogenetic tree to be rescaled.
#' @param RR  is the output of \code{RRphylo} performed on \code{tree}. If this
#'   parameter is indicated, the tree branches are rescaled according to
#'   branch-wise phenotypic evolutionary rates fitted by \code{\link{RRphylo}}.
#'   When a multivariate phenotype is used, rescaling is operated on the norm-2
#'   vector of rates.
#' @param height is the desired height of the tree. If this parameter is
#'   indicated, the tree branches are rescaled to match the total \code{height}.
#' @param trend is a diffusion model with linear trend in rates through time.
#'   The \code{trend} scaling is largely based on package \pkg{geiger}'s
#'   \code{rescale.phylo} function.
#' @param delta if this parameter is indicated, the tree is rescaled according
#'   to Pagel's delta transform (Pagel 1999). Nodes are pushed toward the
#'   present for values of \code{delta} ranging between 0 and 1. The converse
#'   applies for \code{delta} larger than 1. Negative \code{delta} values are
#'   not allowed.
#' @param kappa if this parameter is indicated, the tree is rescaled according
#'   to Pagel's kappa transform (Pagel 1999). At \code{kappa = 1} the tree is
#'   left unmodified. Branches become increasingly closer to 1 as \code{kappa}
#'   approaches 0, making evolution independent from branch lengths. Negative
#'   \code{kappa} values are not allowed.
#' @param lambda if this parameter is indicated, the tree is rescaled according
#'   to Pagel's lambda transform (Pagel 1999). At \code{lambda = 1} the tree is
#'   left unmodified. The tree approaches a star phylogeny as \code{lambda}
#'   approaches zero. \code{lambda} values larger than one are undefined
#'   Negative \code{lambda} values are not allowed.
#' @export
#' @return Rescaled phylogenetic tree.
#' @author Silvia Castiglione, Pasquale Raia
#' @references Castiglione, S., Serio, C., Piccolo, M., Mondanaro, A.,
#'   Melchionna, M., Di Febbraro, M., Sansalone, G., Wroe, S., & Raia, P.
#'   (2020). The influence of domestication, insularity and sociality on the
#'   tempo and mode of brain size evolution in mammals. \emph{Biological Journal
#'   of the Linnean Society},132: 221-231. doi:10.1093/biolinnean/blaa186
#' @references Pagel, M. (1999). Inferring the historical patterns of biological
#'   evolution. \emph{Nature}, 401:877-884.
#' @examples
#' \dontrun{
#' ape::rtree(100)->tree
#' phytools::fastBM(tree)->y
#' max(diag(vcv(tree)))->H
#'
#' RRphylo(tree,y,clus=0)->RR
#' rescaleRR(tree,RR=RR)->treeRR
#'
#' rescaleRR(tree,height=H/3)->tree_height
#'
#' rescaleRR(tree,trend=5)->tree_trend
#'
#' rescaleRR(tree,delta=0.5)->tree_delta05
#' rescaleRR(tree,delta=2)->tree_delta2
#'
#' rescaleRR(tree,kappa=0.5)->tree_kappa
#'
#' rescaleRR(tree,lambda=0.5)->tree_lambda
#' }

rescaleRR<-function(tree,RR=NULL,height=NULL,trend=NULL,delta=NULL,kappa=NULL,lambda=NULL){
  # require(ape)

  formals(rescaleRR)[-1]->allargs

  if(all(sapply(names(allargs),function(x) is.null(get(x)))))
    stop("One of RR, height, trend, delta, kappa, lambda, must be indicated")
  if(sum(sapply(names(allargs),function(x) !is.null(get(x))))>1)
    stop("Only one rescaling type can be indicated")


  tree->tree1
  H<-max(diag(vcv(tree)))

  if(!is.null(kappa)){
    if (kappa<0) stop("kappa must be positive")
    tree$edge.length<-tree$edge.length^kappa
    return(tree)
  }

  if(!is.null(lambda)){
    if (lambda<0) stop("lambda must be positive")
    dN.tips<-dist.nodes(tree)[1:Ntip(tree),(Ntip(tree)+1)]
    tree$edge.length * lambda->bl
    bl[match(names(dN.tips),tree$edge[,2])]<-
      bl[match(names(dN.tips),tree$edge[,2])]+(dN.tips-(dN.tips*lambda))
    tree$edge.length<-bl
    if (any(tree$edge.length<0)) {
      warning("Negative branch lengths generated: please reduce the value of lambda")
    }
  }

  if(!is.null(height)){
    if(!is.null(tree$root.edge)) H+tree$root.edge->H
    tree$edge.length<-(tree$edge.length/H)*height
    if (!is.null(tree$root.edge))
      tree$root.edge <-(tree$root.edge/H)*height
  }

  if(!is.null(trend)|!is.null(delta)){
    H-nodeHeights(tree)->nH
    rownames(nH)<-tree$edge[,2]
    if(!is.null(tree$root.edge)) max(diag(vcv(tree)))+tree$root.edge->H
    data.frame(nH,H-nH[,1],H-nH[,2])->nH
    colnames(nH)<-c("agepar","agedes","head","tail")

    if(!is.null(trend)){
      trend->slope
      nH$br = 1 + nH$head * slope
      nH$er = 1 + nH$tail * slope

      mult<-sapply(1:nrow(nH), function(k) {
        if (rownames(nH)[k] == Ntip(tree) + 1) NA else{
          if (nH$br[k]>0&nH$er[k]>0) mean(c(nH[k,5],nH[k,6])) else
            if (nH$br[k]<0&nH$er[k]<0) 0 else {
              si = -1/slope
              nH$br[k] * (si - nH$head[k])/(2 *(nH$tail[k] - nH$head[k]))
            }
        }
      })
      tree$edge.length<-tree$edge.length*mult
    }

    if(!is.null(delta)){
      if (delta<0) stop("delta must be positive")
      data.frame(nH,bl=tree$edge.length)->nH

      tree$edge.length<-(nH$head+nH$bl)^delta-nH$head^delta
      tree$edge.length<-(tree$edge.length/(H^delta)) * H
    }
  }

  if(!is.null(RR)){
    abs(RR$rates[,1])->rts
    sum(tree$edge.length)->tele
    rts[-1]->rts
    names(rts)[Nnode(tree):length(rts)]<-match(names(rts)[Nnode(tree):length(rts)],tree1$tip.label)
    rts[match(tree$edge[,2],names(rts))]->rts
    tree$edge.length*rts->tree$edge.length
    tele/sum(tree$edge.length)*tree$edge.length->tree$edge.length
  }

  return(tree)
}
