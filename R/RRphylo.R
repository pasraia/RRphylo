#' @title Evolutionary rates computation along phylogenies
#' @description The function \code{RRphylo} (\cite{Castiglione et al. 2018}) performs the phylogenetic ridge regression. It takes a tree and a vector of tip data (phenotypes) as entries, calculates the regularization factor, produces the matrices of tip to root (\code{\link{makeL}}), and node to root distances (\code{\link{makeL1}}), the vector of ancestral state estimates, the vector of predicted phenotypes, and the rates vector for all the branches of the tree. For multivariate data, rates are given as both one vector per variable, and as a multivariate vector obtained by computing the Euclidean Norm of individual rate vectors.
#' @usage RRphylo(tree,y,cov=NULL,rootV=NULL,aces=NULL,clus=0.5)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous.
#' @param y either a single vector variable or a multivariate dataset of class \sQuote{matrix}. In any case, \code{y} must be named.
#' @param cov the covariate to be indicated if its effect on the rates must be accounted for. In this case residuals of the covariate versus the rates are used as rates. \code{'cov'} must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate, as in the example below.
#' @param rootV phenotypic value (values if multivariate) at the tree root. If \code{rootV=NULL} the function 'learns' about the root value from the 10\% tips being closest in time to the tree root, weighted by their temporal distance from the root itself (close tips phenotypes weigh more than more distant tips).
#' @param aces a named vector (or matrix if \code{y} is multivariate) of ancestral character values at nodes. Names correspond to the nodes in the tree.
#' @param clus the proportion of clusters to be used in parallel computing (only if \code{y} is multivariate).
#' @export
#' @importFrom ape multi2di Ntip is.binary.tree Nnode
#' @importFrom stats dist lm residuals weighted.mean
#' @importFrom stats4 mle
#' @importFrom geiger treedata tips
#' @importFrom phytools bind.tip
#' @return \strong{tree} the tree used by \code{RRphylo}. The fully dichotomous version of the tree argument. For trees with polytomies, the tree is resolved by using \code{multi2di} function in the package \pkg{ape}. If the latter is a dichotomous tree, the two trees will be identical.
#' @return \strong{tip.path} a \eqn{n * m} matrix, where n=number of tips and m=number of branches (i.e. 2*n-1). Each row represents the branch lengths along a root-to-tip path.
#' @return \strong{node.path} a \eqn{n * n} matrix, where n=number of internal branches. Each row represents the branch lengths along a root-to-node path.
#' @return \strong{rates} single rate values computed for each branch of the tree. If \code{y} is a single vector variable, rates are equal to multiple.rates. If \code{y} is a multivariate dataset, rates are computed as the square root of the sum of squares of each row of \code{$multiple.rates}.
#' @return \strong{aces} the phenotypes reconstructed at nodes.
#' @return \strong{predicted.phenotypes} the vector of estimated tip values. It is a matrix in the case of multivariate data.
#' @return \strong{multiple.rates} a \eqn{n * m} matrix, where n= number of branches (i.e. n*2-1) and m = number of variables. For each branch, the column entries represent the evolutionary rate.
#' @return \strong{lambda} the regularization factor fitted within \code{RRphylo} by the inner function \code{optL}. With multivariate data, several \code{optL} runs are performed. Hence, the function provides a single lambda for each individual variable.
#' @return \strong{ace values} if \code{aces} are specified, the function returns a dataframe containing the corresponding node number on the \code{RRphylo} tree  for each node , along with estimated values.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references
#' Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method for testing evolutionary rate variation and shifts in phenotypic evolution. \emph{Methods in Ecology and Evolution}, 9: 974-983.doi:10.1111/2041-210X.12954
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
#'
#' # Case 3. "RRphylo" specifying the ancestral states
#'   data("DataCetaceans")
#'   DataCetaceans$treecet->treecet
#'   DataCetaceans$masscet->masscet
#'   DataCetaceans$aceMyst->aceMyst
#'
#'     RRphylo(tree=treecet,y=masscet,aces=aceMyst)
#'
#'     }




RRphylo<-function (tree, y, cov = NULL, rootV = NULL,aces=NULL, clus=0.5)
{
  # require(ape)
  # require(phytools)
  # require(geiger)
  # require(stats4)
  # require(foreach)
  # require(doParallel)
  # require(parallel)

  insert.at <- function(a, pos, ...){
    dots <- list(...)
    stopifnot(length(dots)==length(pos))
    result <- vector("list",2*length(pos)+1)
    result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
    result[c(FALSE,TRUE)] <- dots
    unlist(result)
  }

  optL <- function(lambda) {
    y <- scale(y)
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(y) - rootV)
    aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
    y.hat <- (L %*% betas) + rootV
    Rvar <- array()
    for (i in 1:Ntip(t)) {
      ace.tip <- betas[match(names(which(L[i, ] != 0)),
                             rownames(betas)), ]
      mat = as.matrix(dist(ace.tip))
      Rvar[i] <- sum(mat[row(mat) == col(mat) + 1])
    }
    abs(1 - (var(Rvar) + (mean(as.matrix(y))/mean(y.hat))))
  }

  if (is.binary.tree(tree)) t <- tree else t <- multi2di(tree,random=FALSE)
  if (class(y) == "data.frame")
    y <- treedata(tree, y, sort = TRUE)[[2]]

  if (is.null(rootV)) {
    if (length(y) > Ntip(t)) {
      if(is.binary.tree(tree)==FALSE) u <- data.frame(y, (1/diag(vcv(tree))^2)) else u <- data.frame(y, (1/diag(vcv(t))^2))
      u <- u[order(u[, dim(u)[2]], decreasing = TRUE), ]
      u1 <- u[1:(dim(u)[1] * 0.1), ]
      rootV <- apply(u1[, 1:dim(y)[2]], 2, function(x) weighted.mean(x,
                                                                     u1[, dim(u1)[2]]))
    }else {
      if(is.binary.tree(tree)==FALSE) u <- data.frame(y, (1/diag(vcv(tree))^2)) else u <- data.frame(y, (1/diag(vcv(t))^2))
      u <- u[order(u[, 2], decreasing = TRUE), ]
      u1 <- u[1:(dim(u)[1] * 0.1), ]
      rootV <- weighted.mean(u1[, 1], u1[, 2])
    }
  }else {
    rootV <- rootV
  }


  if(is.null(aces)==FALSE){
    L <- makeL(t)
    aces->aceV
    if(length(y)>Ntip(t)){
      if(is.null(rownames(aceV))) stop("The matrix of aces needs to be named")
      if (is.binary.tree(tree)==FALSE){
        ac<-array()
        for(i in 1:nrow(aceV)){
          getMRCA(t,tips(tree,rownames(aceV)[i]))->ac[i]
        }
        ac->rownames(aceV)
      }
      if(class(aceV)=="data.frame") as.matrix(aceV)->aceV

      aceV->P
      rownames(aceV)->N
      lapply(N,function(x) tips(t,x))->tar.tips
      names(tar.tips)<-N


      t->treeN
      y->ynew
      i=1
      while(i<=length(N)){
        getMRCA(treeN,tar.tips[[i]])->nn
        bind.tip(treeN, tip.label=paste("nod",N[i],sep=""), edge.length=0, where=nn, position=0.001)->treeN
        which(treeN$tip.label==paste("nod",N[i],sep=""))->npos
        if(npos==1) rbind(P[i,],ynew)->ynew
        if(npos==nrow(ynew)+1)  rbind(ynew,P[i,])->ynew else{
          if(npos>1 & npos<(nrow(ynew)+1))  rbind(ynew[1:(npos-1),],unname(P[i,]),ynew[npos:nrow(ynew),])->ynew
        }
        rownames(ynew)[npos]<-paste("nod",N[i],sep="")
        i=i+1
      }
      treeN->t
      ynew->y
    }else{
      if(is.null(names(aceV))) stop("The vector of aces needs to be named")
      if (is.binary.tree(tree)==FALSE){
        ac<-array()
        for(i in 1:length(aceV)){
          getMRCA(t,tips(tree,names(aceV)[i]))->ac[i]
        }
        ac->names(aceV)
      }

      aceV->P
      names(aceV)->N
      lapply(N, function(x) tips(t,x))->tar.tips
      names(tar.tips)<-N

      t->treeN
      y->ynew
      i=1
      while(i<=length(N)){
        getMRCA(treeN,tar.tips[[i]])->nn
        bind.tip(treeN, tip.label=paste("nod",N[i],sep=""), edge.length=0, where=nn, position=0.001)->treeN
        which(treeN$tip.label==paste("nod",N[i],sep=""))->npos
        if(npos==1) c(P[i],ynew)->ynew
        if(npos==length(ynew)+1)  c(ynew,P[i])->ynew else{
          if(npos>1 & npos<(length(ynew)+1))  insert.at(ynew,npos-1,P[i])->ynew
        }
        names(ynew)[npos]<-paste("nod",N[i],sep="")
        i=i+1
      }
      treeN->t
      ynew->y
    }
  }


  L <- makeL(t)
  L1 <- makeL1(t)
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                              2] > Ntip(t))]))
  edged <- data.frame(t$edge, t$edge.length)

  if(length(y)>Ntip(t)){
    dim(y)[2]->k
    y->y.real
    m.aces<-matrix(ncol=k,nrow=Ntip(t)-1)
    m.betas<-matrix(ncol=k,nrow=Ntip(t)+Nnode(t))
    m.yhat<-matrix(ncol=k,nrow=Ntip(t))
    rootV->rv.real


    res<-list()
    cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res <- foreach(i = 1:k,
                   .packages = c("stats4","ape")) %dopar%
                   {

                     gc()

                     rv.real[i]->rootV
                     y.real[,i]->y
                     h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",
                              upper = 10, lower = 0.001)

                     lambda <- h@coef
                     betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                                 t(L)) %*% (as.matrix(y) - rootV)
                     aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
                     y.hat <- (L %*% betas) + rootV
                     list(aceRR, betas, y.hat, lambda)->res[[i]]
                   }
    stopCluster(cl)

    do.call(cbind, lapply(res,"[[",1))->aceRR
    do.call(cbind, lapply(res,"[[",2))->betas
    do.call(cbind, lapply(res,"[[",3))->y.hat
    unname(sapply(res,"[[",4))->lambda
    rv.real->rootV
    y.real->y
    rownames(betas)<-colnames(L)
    rownames(y.hat)<-rownames(y)
    rownames(aceRR)<-colnames(L1)
    colnames(betas)<-colnames(y.hat)<-colnames(aceRR)<-colnames(y)

  }else{

    h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",
             upper = 10, lower = 0.001)
    lambda <- h@coef
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(y) - rootV)
    aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
    y.hat <- (L %*% betas) + rootV

  }


  if(is.null(aces)==FALSE){
    paste("nod",N,sep="")->tip.rem
    nod.rem<-array()
    for(i in 1:length(N)){
      getMRCA(t,c(tip.rem[i],tar.tips[[i]]))->nod.rem[i]
    }
    as.matrix(aceRR[-match(nod.rem,rownames(aceRR)),])->aceRR
    as.matrix(betas[-match(c(nod.rem,tip.rem),rownames(betas)),])->betas
    as.matrix(y.hat[-match(tip.rem,rownames(y.hat)),])->y.hat
    t<-drop.tip(t,tip.rem)
    makeL(t)->L
    makeL1(t)->L1
    rownames(betas)[1:Nnode(t)]<-rownames(aceRR)<-seq((Ntip(t)+1),(Ntip(t)+Nnode(t)),1)
    if(length(y.hat)>Ntip(t)){
      aceRR[match(rownames(aceV),rownames(aceRR)),]->ace.est
      if(nrow(aceV)<2) data.frame(real.node=rownames(aces),RRnode=rownames(aceV),t(as.matrix(ace.est)))->ace.estimates else
        data.frame(real.node=rownames(aces),RRnode=rownames(aceV),ace.estimate=unname(ace.est))->ace.estimates
    }else{
      aceRR[match(names(aceV),rownames(aceRR)),]->ace.est
      data.frame(real.node=names(aces),RRnode=names(ace.est),ace.estimate=unname(ace.est))->ace.estimates
    }
  }

  betasREAL <- betas

  if (is.null(cov)) {
    rates <- betas
    if (length(y) > Ntip(t)) {
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)
    }else {
      rates <- rates
    }
  }else {
    if (length(y) > Ntip(t)) {
      if (length(which(apply(betas, 1, sum) == 0)) > 0) {
        zeroes <- which(apply(betas, 1, sum) == 0)
        R <- log(abs(betas))
        R <- R[-zeroes, ]
        Y <- abs(cov)
        Y <- Y[-zeroes]
        res <- residuals(lm(R ~ Y))
        factOut <- which(apply(betas, 1, sum) != 0)
        betas[factOut, ] <- res
        betas[zeroes, ] <- 0
      }else {
        R <- log(abs(betas))
        Y <- abs(cov)
        res <- residuals(lm(R ~ Y))
        betas <- as.matrix(res)
      }
      rates <- betas
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)
    }else {
      if (length(which(betas == "0")) > 0) {
        zeroes <- which(betas == "0")
        R <- log(abs(betas))
        R <- R[-zeroes]
        Y <- abs(cov)
        Y <- Y[-zeroes]
        res <- residuals(lm(R ~ Y))
        factOut <- which(betas != "0")
        betas[factOut] <- res
        betas[zeroes] <- 0
      }else {
        R <- log(abs(betas))
        Y <- abs(cov)
        res <- residuals(lm(R ~ Y))
        betas <- as.matrix(res)
      }
      rates <- betas
    }
  }

  if(is.null(aces)){
    res <- list(t, L, L1, rates, aceRR, y.hat, betasREAL, lambda)
    names(res) <- c("tree", "tip.path", "node.path", "rates",
                    "aces", "predicted.phenotype", "multiple.rates", "lambda")
  }else{
    res <- list(t, L, L1, rates, aceRR, y.hat, betasREAL, lambda,ace.estimates)
    names(res) <- c("tree", "tip.path", "node.path", "rates",
                    "aces", "predicted.phenotype", "multiple.rates", "lambda","ace values")
  }
  return(res)
}




