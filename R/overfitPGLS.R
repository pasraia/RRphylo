#' @title Testing PGLS_fossil overfit
#' @description Testing the robustness of \code{\link{PGLS_fossil}} results to
#'   sampling effects and phylogenetic uncertainty.
#' @usage overfitPGLS(modform,oveRR=NULL,phylo.list=NULL,data=NULL,...)
#' @param modform, data as passed to \code{\link{PGLS_fossil}}.
#' @param oveRR an object produced by applying \code{\link{overfitRR}} to be
#'   provided if \code{\link{PGLS_fossil}} rescaled according to \code{\link{RRphylo}} rates
#'   should be performed.
#' @param phylo.list a list of phylogenetic trees to be
#'   provided if \code{\link{PGLS_fossil}} on unscaled trees should be performed.
#' @param data a data.frame or list including response and predictor variables
#'   as named in \code{modform}. If not found in \code{data}, the variables are
#'   taken from current environment.
#' @param ... further argument passed to \code{\link{PGLS_fossil}}.
#' @return The function returns a list containing two 'RRphyloList' objects
#'   including results of \code{\link{PGLS_fossil}} performed by using the phylogeny as
#'   it is (\code{$tree}) and/or rescaled according to \code{\link{RRphylo}} rates
#'   (\code{$RR}).
#'   The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Carmela Serio, Giorgia Girardi, Pasquale Raia
#' @export
#' @seealso \href{../doc/overfitRR.html}{\code{overfitPGLS} vignette} ;
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
#' library(phytools)
#' library(ape)
#'
#' # generate fictional data to test the function
#' rtree(100)->tree
#' fastBM(tree)->resp
#' fastBM(tree,nsim=3)->resp.multi
#' fastBM(tree)->pred1
#' fastBM(tree)->pred2
#' data.frame(y1=resp,x2=pred1,x1=pred2)->dat
#'
#' # perform RRphylo and PGLS_fossil with univariate/multivariate phenotypic data
#' PGLS_fossil(modform=y1~x1+x2,data=dat,tree=tree)->pgls_noRR
#' RRphylo(tree,resp,clus=cc)->RR
#' PGLS_fossil(modform=resp~pred1+pred2,RR=RR)->pgls_RR
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree)->pgls2_noRR
#' RRphylo(tree,resp.multi,clus=cc)->RR2
#' PGLS_fossil(modform=resp.multi~pred1+pred2,tree=tree,RR=RR2)->pgls2_RR
#'
#' # overfitPGLS routine
#' # generate a list of subsampled and swapped phylogenies to test
#' tree.list<-resampleTree(RR$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#'
#' # test the robustnes of PGLS_fossil with univariate/multivariate phenotypic data
#' ofRR<-overfitRR(RR = RR,y=resp,phylo.list=tree.list,clus=cc)
#' ofPGLS<-overfitPGLS(oveRR = ofRR,phylo.list=tree.list,modform = y1~x1+x2,data=dat)
#'
#' ofRR2<-overfitRR(RR = RR2,y=resp.multi,phylo.list=tree.list,clus=cc)
#' ofPGLS2<-overfitPGLS(oveRR = ofRR2,phylo.list=tree.list,modform = y1~x1+x2,
#'                      data=list(y1=resp.multi,x2=pred1,x1=pred2))
#' }

overfitPGLS<-function(modform,oveRR=NULL,phylo.list=NULL,data=NULL,...)
{
  # require(phytools)
  # require(ddpcr)

  if (!requireNamespace("ddpcr", quietly = TRUE)) {
    stop("Package \"ddpcr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  funcall<-match.call()
  if(is.null(oveRR)&is.null(phylo.list)) stop("One of phylo.list or oveRR must be provided")

  if(!is.null(phylo.list)) length(phylo.list)->nsimP else nsimP<-NULL
  if(!is.null(oveRR)) length(oveRR$RR.list)->nsimRR else nsimRR<-NULL
  if(!is.null(phylo.list)&!is.null(oveRR)&&nsimP!=nsimRR){
    warning("The length of phylo.list and oveRR is different, iterations are performed on the shorter")
    min(c(nsimP,nsimRR))->nsim
  } else min(c(nsimP,nsimRR))->nsim

  data<-model.frame(modform,data=data)
  pgls.args<-list(...)

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  PGLScut<-PGLScutRR<-list()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)
    if(!is.null(phylo.list)){
      phylo.list[[k]]->treecut
      ddpcr::quiet(datacut<-treedataMatch(treecut,data)[[1]])
      ddpcr::quiet(do.call(PGLS_fossil,c(list(modform=modform,data=datacut,tree=treecut),pgls.args))->PGLScut[[k]],all=TRUE)
      attributes(PGLScut[[k]])$Call<-NULL
    }

    if(!is.null(oveRR)){
      oveRR$RR.list[[k]]->RRcut
      ddpcr::quiet(datacut<-treedataMatch(RRcut$tree,data)[[1]])
      ddpcr::quiet(do.call(PGLS_fossil,c(list(modform=modform,data=datacut,RR=RRcut),pgls.args))->PGLScutRR[[k]],all=TRUE)
      attributes(PGLScutRR[[k]])$Call<-NULL
    }
  }
  close(pb)

  if(is.null(phylo.list)) PGLScut<-NULL else class(PGLScut)<-"RRphyloList"
  if(is.null(oveRR)) PGLScutRR<-NULL else class(PGLScutRR)<-"RRphyloList"

  list(PGLScut,PGLScutRR)->pgls.res
  names(pgls.res)<-c("tree","RR")
  attr(pgls.res,"Call")<-funcall

  pgls.res
}
