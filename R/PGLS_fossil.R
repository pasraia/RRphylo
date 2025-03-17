#' @title Phylogenetic Generalized Least Square with phylogenies including
#'   fossils
#' @description The function performs pgls for non-ultrametric trees using a
#'   variety of evolutionary models or \code{\link{RRphylo}} rates to change the
#'   tree correlation structure.
#' @usage PGLS_fossil(modform,data=NULL,tree=NULL,RR=NULL,GItransform=FALSE,
#'   intercept=FALSE,model="BM",method=NULL,permutation="none",...)
#' @param modform the formula for the regression model.
#' @param data a data.frame or list including response and predictor variables
#'   as named in \code{modform}. If not found in \code{data}, the variables are
#'   taken from current environment.
#' @param tree a phylogenetic tree to be indicated for any model except if
#'   \code{\link{RRphylo}} is used to rescale tree branches. The tree needs not to be
#'   ultrametric and fully dichotomous.
#' @param RR the result of \code{\link{RRphylo}} performed on the response variable. If
#'   \code{RR} is specified, tree branches are rescaled to the absolute
#'   branch-wise rate values calculated by \code{\link{RRphylo}} to transform the
#'   variance-covariance matrix.
#' @param GItransform logical indicating whether the PGLS approach as in Garland
#'   and Ives (2000) must be applied.
#' @param intercept under \code{GItransform = TRUE}, indicates whether the
#'   intercept should be included in the model.
#' @param model an evolutionary model as indicated in
#'   \code{\link[phylolm]{phylolm}} (for univariate response variable) or
#'   \code{\link[mvMORPH]{mvgls}} (for multivariate response variable).
#' @param method optional argument to be passed to \code{phylolm} (for
#'   univariate response variable) or \code{mvgls} (for multivariate response
#'   variable). See individual functions for further details.
#' @param permutation passed to \code{\link[mvMORPH]{manova.gls}}
#' @param ... further argument passed to \code{phylolm}, \code{mvgls},
#'   \code{manova.gls}, or \code{\link[stats]{lm}}.
#' @importFrom ape corPagel
#' @importFrom stats terms model.frame update manova
#' @export
#' @seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette};
#' @seealso \code{\link{overfitPGLS}}; \href{../doc/overfit.html#overfitPGLS}{\code{overfitPGLS} vignette}
#' @seealso \code{\link[mvMORPH]{mvgls}} ; \code{\link[mvMORPH]{manova.gls}}
#' @seealso \code{\link[phylolm]{phylolm}}
#' @details The function is meant to work with both univariate or multivariate
#'   data, both low- or high-dimensional. In the first case, \code{PGLS_fossil}
#'   uses \code{phylolm}, returning an object of class "phylolm". In the latter,
#'   regression coefficients are estimated by \code{mvgls}, and statistical
#'   significance is obtained by means of permutations within \code{manova.gls}.
#'   In this case, \code{PGLS_fossil} returns a list including the output of
#'   both analyses. In all cases, for both univariate or multivariate data, if
#'   \code{GItransform = TRUE} the functions returns a standard \code{lm}
#'   output. In the latter case, the output additionally includes the result of
#'   manova applied on the multivariate linear model.
#' @return Fitted pgls parameters and significance. The class of the output
#'   object depends on input data (see details).
#'   The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio, Gabriele Sansalone,
#'   Giorgia Girardi
#' @references Garland, Jr, T., & Ives, A. R. (2000). Using the past to predict
#'   the present: confidence intervals for regression equations in phylogenetic
#'   comparative methods. \emph{The American Naturalist}, 155: 346-364.
#'   doi:10.1086/303327
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#' cc<- 2/parallel::detectCores()
#'
#' rtree(100)->tree
#' fastBM(tree)->resp
#' fastBM(tree,nsim=3)->resp.multi
#' fastBM(tree)->pred1
#' fastBM(tree)->pred2
#'
#' PGLS_fossil(modform=resp~pred1+pred2,tree=tree)->pgls_noRR
#' PGLS_fossil(modform=resp~pred1+pred2,tree=tree,GItransform=TRUE)->GIpgls_noRR
#'
#' RRphylo(tree,resp,clus=cc)->RR
#' PGLS_fossil(modform=resp~pred1+pred2,RR=RR)->pgls_RR
#' PGLS_fossil(modform=resp~pred1+pred2,tree=tree,RR=RR,GItransform=TRUE)->GIpgls_RR
#'
#' # To derive log-likelihood and AIC for outputs of PGLS_fossil applied on univariate
#' # response variables the function AIC can be applied
#' AIC(pgls_noRR)
#' AIC(pgls_RR)
#' AIC(GIpgls_noRR)
#' AIC(GIpgls_RR)
#'
#'
#' PGLS_fossil(modform=resp.multi~pred1+pred2,tree=tree)->pgls2_noRR
#' PGLS_fossil(modform=resp.multi~pred1+pred2,tree=tree,GItransform=TRUE)->GIpgls2_noRR
#'
#' # To evaluate statistical significance of multivariate models, the '$manova'
#' # object must be inspected
#' pgls2_noRR$manova
#' summary(GIpgls2_noRR$manova)
#'
#' RRphylo(tree,resp.multi,clus=cc)->RR2
#' PGLS_fossil(modform=resp.multi~pred1+pred2,RR=RR2)->pgls2_RR
#' PGLS_fossil(modform=resp.multi~pred1+pred2,tree=tree,RR=RR2,GItransform=TRUE)->GIpgls2_RR
#'
#' # To evaluate statistical significance of multivariate models, the '$manova'
#' # object must be inspected
#' pgls2_noRR$manova
#' summary(GIpgls2_noRR$manova)
#' pgls2_RR$manova
#' summary(GIpgls2_RR$manova)
#'
#' logLik(pgls2_noRR$pgls)
#' logLik(pgls2_RR$pgls)
#' }
PGLS_fossil<-function(modform,data=NULL,tree=NULL,RR=NULL,
                      GItransform=FALSE,intercept=FALSE,model="BM",
                      method=NULL,permutation="none",...){
  # require(nlme)
  # require(ape)
  # require(phylolm)

  misspacks<-sapply(c("phylolm","mvMORPH","ddpcr"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  funcall<-match.call()
  if(!is.null(RR)){
    rescaleRR(RR$tree,RR=RR)->tree
    if(model!="BM") warning("The tree is rescaled by using RR rates, 'model' argument is not allowed",immediate. = TRUE)
    model<-"BM"
  }

  mf<-model.frame(modform,data=data)
  if(any(grepl("\\$",colnames(mf)))){
    gsub("\\$","DOLL",colnames(mf))->colnames(mf)
    modform<-eval(parse(text=gsub("\\$","DOLL",deparse(modform))))
  }
  treedataMatch(tree,mf)$y->mf

  if(!GItransform){
    # treemod<-tree
    data<-mf
    if(!is.null(ncol(mf[[1]]))&&ncol(mf[[1]])>1){
      resgls <-  mvMORPH::mvgls(modform,data=data,tree=tree,model=model,method=method,...)
      resman<-mvMORPH::manova.gls(resgls,permutation=permutation,...)
      cat("\n")
      res<-list(pgls=resgls,manova=resman)
    }else{
      res <-phylolm::phylolm(modform,data=data,phy=tree,model=model,method=method,...)
    }
  }else{
    if(model!="BM") warning("Under GItransform, BM model is assumed")
    vcv(tree)->C
    U<-eigen(C)$vectors
    W<-diag(sqrt(eigen(C)$values))
    D<-solve(U%*%W%*%t(U))

    y<-as.matrix(mf[,1,drop=FALSE])
    x<-as.matrix(mf[,-1,drop=FALSE])

    if(intercept){
      x<-cbind("Intercept"=rep(1, nrow(x)), as.matrix(x))
      modform<-update(modform,y~Intercept+.-1)
    } else modform<-update(modform,y~.-1)

    y <- D %*% y
    x <- D %*% x

    rownames(x)<-rownames(y)<-tree$tip.label
    ddpcr::quiet(sapply(1:ncol(x),function(j) assign(colnames(x)[j],x[,j,drop=FALSE],envir = parent.env(environment()))))

    data<-model.frame(modform,data=environment())
    lmmod<-lm(modform,data=data,model=TRUE,method="qr",...)
    if(!is.null(ncol(mf[[1]]))&&ncol(mf[[1]])>1){
      resman<-manova(lmmod)
      res<-list(pgls=lmmod,manova=resman)
    } else res<-lmmod
  }

  attr(res,"Call")<-funcall

  return(res)
}




