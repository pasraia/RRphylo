#' @title Phylogenetic Generalized Least Square with fossil phylogenies
#' @description The function performs pgls for non-ultrametric trees using
#'   either Pagel's lambda transform, Brownian Motion or \code{\link{RRphylo}}
#'   rates to change the correlation structure.
#' @usage PGLS_fossil(modform,data,tree,RR=NULL)
#' @param modform the formula for the regression model.
#' @param data a list of named vectors including response and predictor
#'   variables as named in \code{modform}.
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param RR the result of \code{RRphylo} performed on the response variable. If
#'   \code{NULL} the function fits Pagel's lambda in the regression for
#'   univariate data or uses the tree variance covariance matrix in the
#'   multivariate case. If \code{RR} is specified, tree branches are rescaled to
#'   the absolute branch-wise rate values calculated by \code{RRphylo} to
#'   transform the variance-covariance matrix.
#' @importFrom ape corPagel
#' @importFrom stats terms
#' @details With univariate data, the user may want to use either Pagel's lambda
#'   or \code{RRphylo} rates to transform the correlation structure. In the
#'   former case, the lambda transform is fitted to the data (Revell, 2010). In
#'   the latter case, branch lengths are multiplied by absolute rates as
#'   computed by \code{RRphylo} to accommodate rate variation across the tree. In
#'   the multivariate case, the variance-covariance structure is either left
#'   unaltered by keeping \code{RR = NULL} (Adams and Collyer, 2015) or changed
#'   according to the norm-2 vector of rates computed for each phenotype by
#'   specifying the \code{RR} object.
#' @export
#' @seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#' @return Fitted pgls parameters and significance.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Revell, L.J. (2010). Phylogenetic signal and linear regression on
#'   species data. \emph{Methods in Ecology and Evolution}, 1, 319-329.
#'   https://doi.org/10.1111/j.2041-210X.2010.00044.x
#' @references Adams, D.C., & Collyer, M. L. (2017). Multivariate phylogenetic
#'   comparative methods: evaluations, comparisons, and recommendations.
#'   \emph{Systematic Biology}, 67, 14-31. https://doi.org/10.1093/sysbio/syx055
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#'
#' rtree(100)->tree
#' fastBM(tree)->resp
#' fastBM(tree,nsim=3)->resp.multi
#' fastBM(tree)->pred1
#' fastBM(tree)->pred2
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree)
#'
#' RRphylo::RRphylo(tree,resp)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree,RR=RR)
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree)
#' cc<- 2/parallel::detectCores()
#' RRphylo::RRphylo(tree,resp.multi,clus=cc)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree,RR=RR)
#' }

PGLS_fossil<-function(modform,data,tree,RR=NULL)
{
  # require(nlme)
  # require(ape)
  # require(geomorph)
  # require(geiger)

  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package \"nlme\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("geomorph", quietly = TRUE)) {
    stop("Package \"geomorph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }

  for(k in 1:length(data)){
    data[[k]]->sam
    if(is.null(nrow(sam))) data[[k]] <- treedata(tree, sam, sort = TRUE)[[2]][,1] else data[[k]] <- treedata(tree, sam, sort = TRUE)[[2]]
  }

  vars <- rownames(attr(terms(modform), "factors"))

  for(i in 1:length(vars)){
    assign(vars[i],data[which(names(data)==vars[i])][[1]],envir = environment())
  }

  data[which(names(data)==vars[1])][[1]]->y

  if(is.null(RR)){
    if(length(y)>Ntip(tree)){
      data->gdf
      class(gdf)<-"geomorph.data.frame"
      cova<- vcv(tree)
      geomorph::procD.pgls(modform,data=gdf, Cov=cova)->res
    }else{
      co <- corPagel(1, tree)
      v <- diag(vcv(co))
      vf <- nlme::varFixed(~ offset(v))
      suppressWarnings(nlme::gls(modform, correlation=co, weights=vf)->res)
    }
  }else{
    tree->tree1
    abs(RR$rates[,1])->rts
    rts[-1]->rts
    names(rts)[Nnode(tree):length(rts)]<-seq(1,Ntip(tree))
    rts[match(tree$edge[,2],names(rts))]->rts
    tree$edge.length*rts->tree1$edge.length

    data->gdf
    class(gdf)<-"geomorph.data.frame"
    cova<- vcv(tree1)
    geomorph::procD.pgls(modform,data=gdf, Cov=cova)->res
  }
  return(res)
}
