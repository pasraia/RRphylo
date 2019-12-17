#' @title Phylogenetic Generalized Least Square with fossil phylogenies
#' @description The function performs pgls for non-ultrametric trees using either Pagel's lambda transform, Brownian Motion or \code{\link{RRphylo}} rates to change the correlation structure.
#' @usage PGLS_fossil(modform,RR,type=c("Pagel","RRphylo"))
#' @param modform the formula for the regression model
#' @param RR the result of \code{RRphylo} performed on the response variable
#' @param type the correlation structure to be used. If \code{"Pagel"} the function fits Pagel's lambda in the regression for univariate data or uses the tree variance covariance matrix in the multivariate case. If \code{"RRphylo"}, tree branches are rescaled to the absolute branch-wise rate values calculated by \code{RRphylo} to transform the variance-covariance matrix.
#' @importFrom ape corPagel
#' @importFrom nlme varFixed
#' @importFrom geomorph procD.pgls
#' @details With univariate data, the user may want to use either Pagel's lambda or \code{RRphylo} rates to transform the correlation structure. In the former case, the lambda transform is fitted to the data (Revell, 2010). In the latter case, branch lengths are multiplied by absolute rates as computed by \code{RRphylo} to accomodate rate variation across the tree. In the multivariate case, the variance-covariance structure is either left unaltered by specifying \code{"Pagel"} (Adams and Collyer, 2015) or changed according to the norm-2 vector of rates computed for each phenotype by \code{RRphylo}.
#' @export
#' @return Fitted pgls parameters and significance.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references Revell, L.J. (2010). Phylogenetic signal and linear regression on species data. \emph{Methods in Ecology and Evolution}, 1, 319-329. https://doi.org/10.1111/j.2041-210X.2010.00044.x
#' @references Adams, D.C., & Collyer, M. L. (2017). Multivariate phylogenetic comparative methods: evaluations, comparisons, and recommendations. \emph{Systematic Biology}, 67, 14-31. https://doi.org/10.1093/sysbio/syx055
#' @examples
#' library(ape)
#' library(phytools)
#' rtree(100)->tree
#' fastBM(tree)->y
#' fastBM(tree)->x
#' RRphylo(tree,y)->RR
#'
#' PGLS_fossil(y~x,RR=RR,type="Pagel")
#' PGLS_fossil(y~x,RR=RR,type="RRphylo")


PGLS_fossil<-function(modform,RR,type=c("Pagel","RRphylo"))
{
  # require(nlme)
  # require(ape)
  # require(geomorph)
  # require(geiger)

  RR$tree->tree->tree1
  get(as.character(modform)[[2]])->y
  if (class(y) == "data.frame")
    y <- treedata(tree, y, sort = TRUE)[[2]]

  if(type=="Pagel"){
    if(length(y)>Ntip(tree)){
      summary(procD.pgls(modform, Cov = vcv(tree)))->res
    }else{
      co <- corPagel(1, tree)
      v <- diag(vcv(co))
      vf <- varFixed(~ offset(v))
      suppressWarnings(summary(gls(modform, correlation=co, weights=vf))->res)
    }
  }else{
    abs(RR$rates[,1])->rts
    rts[-1]->rts
    names(rts)[Nnode(tree):length(rts)]<-seq(1,Ntip(tree))
    rts[match(tree$edge[,2],names(rts))]->rts
    tree$edge.length*rts->tree1$edge.length

    summary(procD.pgls(modform, Cov = vcv(tree1)))->res
  }
  return(res)
}
