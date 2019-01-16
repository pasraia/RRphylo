#' @title Phylogenetic Generalized Least Square with fossil phylogenies
#' @description The function performs pgls for non-ultrametric trees with lambda correlation.
#' @usage PGLS_fossil(tree,x,y)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and fully dichotomous.
#' @param x the predictor variable
#' @param y the response variable
#' @importFrom ape corPagel
#' @importFrom nlme varFixed
#' @export
#' @return Fitted pgls parameters and significance.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' statedino[match(treedino$tip.label,names(statedino))]->statedino
#' massdino[match(treedino$tip.label,names(massdino))]->massdino
#'
#' PGLS_fossil(tree=treedino,x=statedino,y=massdino)
#'


PGLS_fossil<-function(tree=tree,
                      x=x,
                      y=y)
{
  #require(nlme)
  #require(ape)
  co <- corPagel(1, tree)
  v <- diag(vcv(co))
  vf <- varFixed(~ offset(v))
  gls(y ~ x, correlation=co, weights=vf)->res
  return(res)
}
