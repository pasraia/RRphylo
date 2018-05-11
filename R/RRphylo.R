#' @title Evolutionary rates computation along phylogenies
#' @description The function \code{RRphylo} (\cite{Castiglione et al. 2018}) performs the phylogenetic ridge regression. It takes a tree and a vector of tip data (phenotypes) as entries, calculates the regularization factor, produces the matrices of tip to root (\code{\link{makeL}}), and node to root distances (\code{\link{makeL1}}), the vector of ancestral state estimates, the vector of predicted phenotypes, and the rates vector for all the branches of the tree. For multivariate data, rates are given as both one vector per variable, and as a multivariate vector obtained by computing the Euclidean Norm of individual rate vectors.
#' @usage RRphylo(tree,y,cov=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous.
#' @param y either a single vector variable or a multivariate dataset of class \sQuote{matrix}.
#' @param cov the covariate to be indicated if its effect on the rates must be accounted for. In this case residuals of the covariate versus the rates are used as rates. \code{'cov'} must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate, as in the example below.
#' @export
#' @importFrom ape multi2di Ntip is.binary.tree Nnode
#' @importFrom stats dist lm residuals
#' @importFrom stats4 mle
#' @importFrom geiger treedata tips
#' @return \strong{tree} the tree used by \code{RRphylo}. The fully dichotomous version of the tree argument. For trees with polytomies, the tree is resolved by using \code{multi2di} function in the package \pkg{ape}. If the latter is a dichotomous tree, the two trees will be identical.
#' @return \strong{tip.path} a \eqn{n * m} matrix, where n=number of tips and m=number of branches (i.e. 2*n-1). Each row represents the branch lengths along a root-to-tip path.
#' @return \strong{node.path} a \eqn{n * n} matrix, where n=number of internal branches. Each row represents the branch lengths along a root-to-node path.
#' @return \strong{rates} single rate values computed for each branch of the tree. If \code{y} is a single vector variable, rates are equal to multiple.rates. If \code{y} is a multivariate dataset, rates are computed as the square root of the sum of squares of each row of \code{$multiple.rates}.
#' @return \strong{aces} the phenotypes reconstructed at nodes.
#' @return \strong{predicted.phenotypes} the vector of estimated tip values.
#' @return \strong{multiple.rates} a \eqn{n * m} matrix, where n= number of branches (i.e. n*2-1) and m = number of variables. For each branch, the column entries represent the evolutionary rate.
#' @return \strong{lambda} the regularization factor fitted within \code{RRphylo} by the inner function \code{optL}.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references
#' Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method for testing evolutionary rate variation and shifts in phenotypic evolution. \emph{Methods in Ecology and Evolution}, in press.doi:10.1111/2041-210X.12954
#' @examples
#'  \donttest{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#'
#' # Case 1. "RRphylo" without accounting for the effect of a covariate
#'     RRphylo(tree=treedino,y=massdino)
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#'   # "RRphylo" on the covariate in order to retrieve ancestral state values
#'     RRphylo(tree=treedino,y=massdino)->RRcova
#'     c(RRcova$aces,massdino)->cov.values
#'     c(rownames(RRcova$aces),names(massdino))->names(cov.values)
#'
#'     RRphylo(tree=treedino,y=massdino,cov=cov.values)
#'     }



RRphylo<-function(tree,y,cov=NULL)
  ##### tree does not accept numbers as node labels ####
{
  #require(ape)
  #require(phytools)
  #require(geiger)
  #require(stats4)


  if (is.binary.tree(tree)) t <- tree else t <- multi2di(tree)
  if(class(y)=="data.frame") treedata(tree,y,sort=TRUE)[[2]]->y
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                              2] > Ntip(t))]))
  tippa <- list()
  for (i in 1:length(internals)) {
    tippas <- tips(t, internals[i])
    dato <- data.frame(rep(internals[i], length(tippas)),
                       tippas)
    colnames(dato)[1] <- "node"
    tippa[[i]] <- dato
  }
  Tstr <- do.call(rbind, tippa)
  makeL(t)->L
  makeL1(t)->L1
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                              2] > Ntip(t))]))
  edged <- data.frame(t$edge, t$edge.length)
  optL <- function(lambda) {
    y <- scale(y)
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% as.matrix(y)
    aceRR <- L1 %*% betas[1:Nnode(t), ]
    y.hat <- L %*% betas
    Rvar <- array()
    for (i in 1:Ntip(t)) {
      ace.tip <- betas[match(names(which(L[i, ] != 0)),
                             rownames(betas)), ]
      mat = as.matrix(dist(ace.tip))
      Rvar[i] <- sum(mat[row(mat) == col(mat) + 1])
    }
    abs(1 - (var(Rvar) + (mean(as.matrix(y))/mean(y.hat))))
  }
  h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",upper=10,lower=0.001)
  lambda <- h@coef


  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% t(L)) %*% as.matrix(y)


  aceRR <- L1 %*% betas[1:Nnode(t), ]
  y.hat <- L %*% betas
  rates <- betas
  betasREAL<-betas
  if (length(y) > Ntip(t)) {
    rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
    rates <- as.matrix(rates)
  } else {
    rates <- rates
  }



  if(is.null(cov)){

    rates <- betas
    if (length(y) > Ntip(t)) {
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)
    } else {
      rates <- rates
    }

  }else{
    if (length(y) > Ntip(t)) {
      if(length(which(apply(betas,1,sum)==0))>0){
        which(apply(betas,1,sum)==0)->zeroes
        log(abs(betas))->R
        R[-zeroes,]->R

        abs(cov)->Y
        Y[-zeroes]->Y

        residuals(lm(R~Y))->res
        which(apply(betas,1,sum)!=0)->factOut
        betas[factOut,]<-res
        betas[zeroes,]<-0

      }else {
        log(abs(betas))->R
        abs(cov)->Y
        residuals(lm(R~Y))->res
        as.matrix(res)->betas
      }
      rates<-betas
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)

    }else{

      if(length(which(betas=="0"))>0){
        which(betas=="0")->zeroes
        log(abs(betas))->R
        R[-zeroes]->R

        abs(cov)->Y
        Y[-zeroes]->Y

        residuals(lm(R~Y))->res
        which(betas!="0")->factOut

        betas[factOut]<-res
        betas[zeroes]<-0
      } else {
        log(abs(betas))->R
        abs(cov)->Y
        residuals(lm(R~Y))->res
        as.matrix(res)->betas
      }
      betas->rates
    }
  }

  res <- list(t, L, L1, rates, aceRR,y.hat,betasREAL, lambda)
  names(res) <- c("tree", "tip.path", "node.path", "rates",
                  "aces","predicted.phenotype","multiple.rates", "lambda")
  return(res)
}
