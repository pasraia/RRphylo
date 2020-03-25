#' @title Evolutionary rates computation along phylogenies
#' @description The function \code{RRphylo} (\cite{Castiglione et al. 2018}) performs the phylogenetic ridge regression. It takes a tree and a vector of tip data (phenotypes) as entries, calculates the regularization factor, produces the matrices of tip to root (\code{\link{makeL}}), and node to root distances (\code{\link{makeL1}}), the vector of ancestral state estimates, the vector of predicted phenotypes, and the rates vector for all the branches of the tree. For multivariate data, rates are given as both one vector per variable, and as a multivariate vector obtained by computing the Euclidean Norm of individual rate vectors.
#' @usage RRphylo(tree,y,cov=NULL,rootV=NULL,aces=NULL,x1=NULL,aces.x1=NULL,clus=0.5)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous.
#' @param y either a single vector variable or a multivariate dataset of class \sQuote{matrix}. In any case, \code{y} must be named.
#' @param cov the covariate to be indicated if its effect on the rates must be accounted for. In this case residuals of the covariate versus the rates are used as rates. \code{'cov'} must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate, as in the example below.
#' @param rootV phenotypic value (values if multivariate) at the tree root. If \code{rootV=NULL} the function 'learns' about the root value from the 10\% tips being closest in time to the tree root, weighted by their temporal distance from the root itself (close tips phenotypes weigh more than more distant tips).
#' @param aces a named vector (or matrix if \code{y} is multivariate) of ancestral character values at nodes. Names correspond to the nodes in the tree.
#' @param x1 the additional predictor to be indicated to perform the multiple version of \code{RRphylo}. \code{'x1'} vector must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the predictor as well, and taking the vector of ancestral states and tip values to form the \code{x1}.
#' @param aces.x1 a named vector of ancestral character values at nodes for \code{x1}. It must be indicated if both \code{aces} and \code{x1} are specified. Names correspond to the nodes in the tree.
#' @param clus the proportion of clusters to be used in parallel computing (only if \code{y} is multivariate).
#' @export
#' @importFrom ape multi2di Ntip is.binary.tree Nnode dist.nodes drop.tip subtrees nodelabels
#' @importFrom stats dist lm residuals weighted.mean
#' @importFrom stats4 mle
#' @importFrom geiger treedata tips ratematrix
#' @importFrom phytools bind.tip
#' @return \strong{tree} the tree used by \code{RRphylo}. The fully dichotomous version of the tree argument. For trees with polytomies, the tree is resolved by using \code{multi2di} function in the package \pkg{ape}. If the latter is a dichotomous tree, the two trees will be identical.
#' @return \strong{tip.path} a \eqn{n * m} matrix, where n=number of tips and m=number of branches (i.e. 2*n-1). Each row represents the branch lengths along a root-to-tip path.
#' @return \strong{node.path} a \eqn{n * n} matrix, where n=number of internal branches. Each row represents the branch lengths along a root-to-node path.
#' @return \strong{rates} single rate values computed for each branch of the tree. If \code{y} is a single vector variable, rates are equal to multiple.rates. If \code{y} is a multivariate dataset, rates are computed as the square root of the sum of squares of each row of \code{$multiple.rates}.
#' @return \strong{aces} the phenotypes reconstructed at nodes.
#' @return \strong{predicted.phenotypes} the vector of estimated tip values. It is a matrix in the case of multivariate data.
#' @return \strong{multiple.rates} a \eqn{n * m} matrix, where n= number of branches (i.e. n*2-1) and m = number of variables. For each branch, the column entries represent the evolutionary rate.
#' @return \strong{lambda} the regularization factor fitted within \code{RRphylo} by the inner function \code{optL}. With multivariate data, several \code{optL} runs are performed. Hence, the function provides a single lambda for each individual variable.
#' @return \strong{ace.values} if \code{aces} are specified, the function returns a dataframe containing the corresponding node number on the \code{RRphylo} tree  for each node , along with estimated values.
#' @return \strong{x1.rate} if \code{x1} is specified, the function returns the partial regression coefficient for \code{x1}.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method for testing evolutionary rate variation and shifts in phenotypic evolution. \emph{Methods in Ecology and Evolution}, 9: 974-983.doi:10.1111/2041-210X.12954
#' @references Serio, C., Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M., Mondanaro, A., Di Febbraro, M., & Raia, P.(2019). Macroevolution of toothed whales exceptional relative brain size. \emph{Evolutionary Biology}, 46: 332-342. doi:10.1007/s11692-019-09485-7
#' @references  Melchionna, M., Mondanaro, A., Serio, C., Castiglione, S., Di Febbraro, M., Rook, L.,Diniz-Filho,J.A.F., Manzi, G., Profico, A., Sansalone, G., & Raia, P.(2020).Macroevolutionary trends of brain mass in Primates. \emph{Biological Journal of the Linnean Society}, 129: 14-25. doi:10.1093/biolinnean/blz161
#' @examples
#'  \donttest{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#'
#' # Case 1. "RRphylo" without accounting for the effect of a covariate
#' RRphylo(tree=treedino,y=massdino)->RRcova
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#' # "RRphylo" on the covariate in order to retrieve ancestral state values
#' c(RRcova$aces,massdino)->cov.values
#' c(rownames(RRcova$aces),names(massdino))->names(cov.values)
#'
#' RRphylo(tree=treedino,y=massdino,cov=cov.values)->RR
#'
#' # Case 3. "RRphylo" specifying the ancestral states
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#' DataCetaceans$aceMyst->aceMyst
#'
#' RRphylo(tree=treecet,y=masscet,aces=aceMyst)->RR
#'
#' # Case 4. Multiple "RRphylo"
#' library(ape)
#' drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' RRphylo(tree=treecet.multi, y=masscet.multi)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass)->RR
#'     }




RRphylo<-function (tree, y, cov = NULL, rootV = NULL, aces = NULL,x1=NULL,aces.x1=NULL, clus = 0.5)
{
  # library(phytools)
  # library(stats4)
  # library(ape)
  # library(parallel)
  # library(doParallel)

  insert.at <- function(a, pos, ...) {
    dots <- list(...)
    stopifnot(length(dots) == length(pos))
    result <- vector("list", 2 * length(pos) + 1)
    result[c(TRUE, FALSE)] <- split(a, cumsum(seq_along(a) %in%
                                                (pos + 1)))
    result[c(FALSE, TRUE)] <- dots
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
  optLmultiple <- function(lambda){
    scale(c(rootV,y))->yx
    yx[1]->rv
    yx[-1]->y

    scale(y1)->y1
    scale(ace1)->ace1
    ace1[1]->acex

    y-rv->yy
    y1-acex->y1
    cbind(L,y1)->L

    ace1-acex->ace1
    cbind(L1,ace1)->L1


    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(yy))
    y.hat <- (L %*% betas)+rootV
    aceRR <- ((L1 %*% betas[c(1:Nnode(t),length(betas)), ]))+rootV

    abs(1-summary(lm(c(aceRR[1],y.hat)~c(rootV,y)))$coef[2])
  }

  if (is.binary.tree(tree))
    t <- tree else t <- multi2di(tree, random = FALSE)
  # if (inherits(y,"data.frame"))
  #   y <- treedata(tree, y, sort = TRUE)[[2]]
  if(is.null(nrow(y))) y <- treedata(tree, y, sort = TRUE)[[2]][,1] else y <- treedata(tree, y, sort = TRUE)[[2]]

  Loriginal <-L<-makeL(t)
  L1original <-L1<-makeL1(t)

  if (is.null(rootV)) { #### rootV ####
    if (length(y) > Ntip(t)) { #### rootV multi ####
      if (is.binary.tree(tree) == FALSE)
        u <- data.frame(y, (1/diag(vcv(tree))^2))
      else u <- data.frame(y, (1/diag(vcv(t))^2))
      u <- u[order(u[, dim(u)[2]], decreasing = TRUE),
             ]
      u1 <- u[1:(dim(u)[1] * 0.1), ]
      rootV <- apply(u1[, 1:dim(y)[2]], 2, function(x) weighted.mean(x,
                                                                     u1[, dim(u1)[2]]))
    }else { #### rootV uni ####
      if (is.binary.tree(tree) == FALSE)
        u <- data.frame(y, (1/diag(vcv(tree))^2))
      else u <- data.frame(y, (1/diag(vcv(t))^2))
      u <- u[order(u[, 2], decreasing = TRUE), ]
      u1 <- u[1:(dim(u)[1] * 0.1), ]
      rootV <- weighted.mean(u1[, 1], u1[, 2])
    }
  }else {
    if (inherits(rootV,"data.frame")) as.matrix(rootV)->rootV  else  rootV <- rootV

  }

  if(is.null(x1)==FALSE){ #### multiple Ridge Regression ####
    x1[-match(t$tip.label,names(x1))]->ace1
    x1[match(t$tip.label,names(x1))]->y1
  }



  if (is.null(aces) == FALSE) { #### phenotypes at internal nodes ####
    L <- makeL(t)
    aceV <- aces
    if (length(y) > Ntip(t)) { #### aces multi ####
      if (is.null(rownames(aceV)))
        stop("The matrix of aces needs to be named")
      if (is.binary.tree(tree) == FALSE) {
        ac <- array()
        for (i in 1:nrow(aceV)) {
          ac[i] <- getMRCA(t, tips(tree, rownames(aceV)[i]))
        }
        rownames(aceV) <- ac
      }
      if (inherits(aceV,"data.frame"))
        aceV <- as.matrix(aceV)
      P <- aceV
      N <- as.numeric(rownames(aceV))
      tar.tips <- lapply(N, function(x) tips(t, x))
      names(tar.tips) <- N
      treeN <- t
      ynew <- y

      if(is.null(x1)==FALSE){
        Px1<-aces.x1
        y1new <- y1
        ace1new<-ace1
      }

      i = 1
      while (i <= length(N)) {
        nn <- getMRCA(treeN, tar.tips[[i]])
        treeN <- bind.tip(treeN, tip.label = paste("nod",
                                                   N[i], sep = ""), edge.length = 0, where = nn,
                          position = 0.001)
        npos <- which(treeN$tip.label == paste("nod", N[i], sep = ""))
        if (npos == 1) ynew <- rbind(P[i, ], ynew)
        if (npos == nrow(ynew) + 1) ynew <- rbind(ynew, P[i, ]) else {
          if (npos > 1 & npos < (nrow(ynew) + 1)) ynew <- rbind(ynew[1:(npos - 1), ], unname(P[i,]), ynew[npos:nrow(ynew), ])
        }

        if(is.null(x1)==FALSE){
          if (npos == 1) y1new <- c(Px1[i], y1new)

          if (npos == length(y1new) + 1)  y1new <- c(y1new, Px1[i])
          else {
            if (npos > 1 & npos < (length(y1new) + 1)) y1new <- insert.at(y1new, npos - 1, Px1[i])
          }
        }
        rownames(ynew)[npos] <- paste("nod", N[i], sep = "")
        if(is.null(x1)==FALSE) names(y1new)[npos] <- paste("nod", N[i], sep = "")
        i = i + 1
      }

      if(is.null(x1)==FALSE){
        sort(N)->Ns
        Px1[match(Ns,names(Px1))]->Px1

        h=1
        while(h<=length(Ns)){
          np<-which(treeN$tip.label ==
                      paste("nod",Ns[h], sep = ""))
          (getMommy(treeN,np)[1]-(Ntip(treeN)))->npos
          #match(getMommy(treeN,np)[1],names(ace1new))->npos
          if(npos== Nnode(treeN)) c(ace1new,Px1[h])->ace1new else
            insert.at(ace1new,(npos-1),Px1[h])->ace1new
          h=h+1
        }
        names(ace1new)<-seq((Ntip(treeN)+1),(Ntip(treeN)+Nnode(treeN)))
        ace1<-ace1new
        y1<-y1new
      }

      t <- treeN
      y <- ynew
    } else { #### aces uni ####
      if (is.null(names(aceV)))
        stop("The vector of aces needs to be named")
      if (is.binary.tree(tree) == FALSE) {
        ac <- array()
        for (i in 1:length(aceV)) {
          ac[i] <- getMRCA(t, tips(tree, names(aceV)[i]))
        }
        names(aceV) <- ac
      }
      P <- aceV
      N <- as.numeric(names(aceV))
      tar.tips <- lapply(N, function(x) tips(t, x))
      names(tar.tips) <- N
      treeN <- t
      ynew <- y

      if(is.null(x1)==FALSE){
        Px1<-aces.x1
        y1new <- y1
        ace1new<-ace1
      }

      i = 1
      while (i <= length(N)) {
        nn <- getMRCA(treeN, tar.tips[[i]])
        treeN <- bind.tip(treeN,
                          tip.label = paste("nod",N[i], sep = ""),
                          edge.length = 0, where = nn,position = 0.001)
        npos <- which(treeN$tip.label ==
                        paste("nod",N[i], sep = ""))

        if (npos == 1) ynew <- c(P[i], ynew)

        if (npos == length(ynew) + 1) ynew <- c(ynew, P[i]) else {
          if (npos > 1 & npos < (length(ynew) + 1)) ynew <- insert.at(ynew, npos - 1, P[i])
        }

        if(is.null(x1)==FALSE){
          if (npos == 1) y1new <- c(Px1[i], y1new)

          if (npos == length(ynew) + 1)  y1new <- c(y1new, Px1[i])
          else {
            if (npos > 1 & npos < (length(ynew) + 1)) y1new <- insert.at(y1new, npos - 1, Px1[i])
          }
        }
        names(ynew)[npos] <- paste("nod", N[i], sep = "")
        if(is.null(x1)==FALSE) names(y1new)[npos] <- paste("nod", N[i], sep = "")
        i = i + 1
      }

      if(is.null(x1)==FALSE){
        sort(N)->Ns
        Px1[match(Ns,names(Px1))]->Px1

        h=1
        while(h<=length(Ns)){
          np<-which(treeN$tip.label ==
                      paste("nod",Ns[h], sep = ""))
          (getMommy(treeN,np)[1]-(Ntip(treeN)))->npos
          if(npos== Nnode(treeN)) c(ace1new,Px1[h])->ace1new else
            insert.at(ace1new,(npos-1),Px1[h])->ace1new
          h=h+1
        }
        names(ace1new)<-seq((Ntip(treeN)+1),(Ntip(treeN)+Nnode(treeN)))
        ace1<-ace1new
        y1<-y1new
      }

      t <- treeN
      y <- ynew
    }
    L<-makeL(t)
    L1<-makeL1(t)
  }


  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                              2] > Ntip(t))]))
  edged <- data.frame(t$edge, t$edge.length)

  if (length(y) > Ntip(t)) { #### multivariate Ridge Regression ####
    k <- dim(y)[2]
    y.real <- y
    rv.real <- rootV
    res <- list()
    cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res <- foreach(i = 1:k, .packages = c("stats4", "ape")) %dopar% {
      #for(i in 1:k){
      gc()
      rootV <- rv.real[i]
      y <- y.real[, i]
      if(is.null(x1)==FALSE){ #### multiple Ridge Regression ####

        h <- mle(optLmultiple, start = list(lambda = 1), method = "L-BFGS-B",
                 upper = 10, lower = 0.001)
        h@coef->lambda
        cbind(L,y1-ace1[1])->LX
        cbind(L1,ace1-mean(ace1))->LX1

        betas <- (solve(t(LX) %*% LX + lambda * diag(ncol(LX))) %*%
                    t(LX)) %*% (as.matrix(y)-rootV)

        aceRR <- ((LX1 %*% betas[c(1:Nnode(t),length(betas)), ]))+rootV
        y.hat <- (LX %*% betas)+rootV

        betas[length(betas)]->x1.rate
        as.matrix(betas[-length(betas),])->betas
        colnames(betas)<-NULL
        res[[i]] <- list(aceRR, betas, y.hat, lambda,x1.rate)

      }else{

        h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",
                 upper = 10, lower = 0.001)
        lambda <- h@coef
        betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                    t(L)) %*% (as.matrix(y) - rootV)
        aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
        y.hat <- (L %*% betas) + rootV
        res[[i]] <- list(aceRR, betas, y.hat, lambda)
      }
    }
    stopCluster(cl)
    aceRR <- do.call(cbind, lapply(res, "[[", 1))
    betas <- do.call(cbind, lapply(res, "[[", 2))
    y.hat <- do.call(cbind, lapply(res, "[[", 3))
    lambda <- unname(sapply(res, "[[", 4))
    if(is.null(x1)==FALSE) x1.rate <- unname(sapply(res, "[[", 5))
    rootV <- rv.real
    y <- y.real
    rownames(betas) <- colnames(L)
    rownames(y.hat) <- rownames(y)
    rownames(aceRR) <- colnames(L1)
    colnames(betas) <- colnames(y.hat) <- colnames(aceRR) <- colnames(y)

  }else { #### Ridge Regression univariate ####

    if(is.null(x1)==FALSE){ #### multiple Ridge Regression ####

      h <- mle(optLmultiple, start = list(lambda = 1), method = "L-BFGS-B",
               upper = 10, lower = 0.001)
      h@coef->lambda
      cbind(L,y1-ace1[1])->LX
      cbind(L1,ace1-mean(ace1))->LX1

      betas <- (solve(t(LX) %*% LX + lambda * diag(ncol(LX))) %*%
                  t(LX)) %*% (as.matrix(y)-rootV)

      aceRR <- ((LX1 %*% betas[c(1:Nnode(t),length(betas)), ]))+rootV
      y.hat <- (LX %*% betas)+rootV
      betas[length(betas)]->x1.rate
      as.matrix(betas[-length(betas),])->betas
      colnames(betas)<-NULL
    }else{

      h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",
               upper = 10, lower = 0.001)
      lambda <- h@coef
      betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                  t(L)) %*% (as.matrix(y) - rootV)
      aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
      y.hat <- (L %*% betas) + rootV

    }
  }

  if (is.null(aces) == FALSE) {
    tip.rem <- paste("nod", N, sep = "")
    nod.rem <- array()
    for (i in 1:length(N)) {
      nod.rem[i] <- getMRCA(t, c(tip.rem[i], tar.tips[[i]]))
    }
    aceRR <- as.matrix(aceRR[-match(nod.rem, rownames(aceRR)),
                             ])
    betas <- as.matrix(betas[-match(c(nod.rem, tip.rem),
                                    rownames(betas)), ])
    y.hat <- as.matrix(y.hat[-match(tip.rem, rownames(y.hat)),
                             ])
    t <- drop.tip(t, tip.rem)
    rownames(betas)[1:Nnode(t)] <- rownames(aceRR) <- seq((Ntip(t) +
                                                             1), (Ntip(t) + Nnode(t)), 1)
    if (length(y.hat) > Ntip(t)) {
      ace.est <- aceRR[match(rownames(aceV), rownames(aceRR)),
                       ]
      if (nrow(aceV) < 2)
        ace.estimates <- data.frame(real.node = rownames(aces),
                                    RRnode = rownames(aceV), t(as.matrix(ace.est)))
      else ace.estimates <- data.frame(real.node = rownames(aces),
                                       RRnode = rownames(aceV), ace.estimate = unname(ace.est))
    }
    else {
      ace.est <- aceRR[match(names(aceV), rownames(aceRR)),
                       ]
      ace.estimates <- data.frame(real.node = names(aces),
                                  RRnode = names(ace.est), ace.estimate = unname(ace.est))
    }
  }
  betasREAL <- betas

  if (is.null(cov)) {
    rates <- betas
    if (length(y) > Ntip(t)) {
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)
    }
  } else {

    cov[match(rownames(betas),names(cov))]->cov

    if (length(y) > Ntip(t)) { #### Covariate multi ####
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
      } else { #### Covariate uni ####
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
      }
      else {
        R <- log(abs(betas))
        Y <- abs(cov)
        res <- residuals(lm(R ~ Y))
        betas <- as.matrix(res)
      }
      rates <- betas
    }
  }


  if (is.null(aces)) {
    res <- list(t, Loriginal, L1original, rates, aceRR, y.hat, betasREAL,
                lambda)
    names(res) <- c("tree", "tip.path", "node.path", "rates",
                    "aces", "predicted.phenotype", "multiple.rates",
                    "lambda")
  }else {
    res <- list(t, Loriginal, L1original, rates, aceRR, y.hat, betasREAL,
                lambda, ace.estimates)
    names(res) <- c("tree", "tip.path", "node.path", "rates",
                    "aces", "predicted.phenotype", "multiple.rates",
                    "lambda", "ace.values")
  }

  if(is.null(x1)==FALSE) {
    res[[(length(res)+1)]]<-x1.rate
    names(res)[length(res)]<-"x1.rate"
  }

  return(res)
}




