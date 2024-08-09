#' @title Evolutionary rates computation along phylogenies
#' @description The function \code{\link{RRphylo}} (\cite{Castiglione et al.
#'   2018}) performs the phylogenetic ridge regression. It takes a tree and a
#'   vector of tip data (phenotypes) as entries, calculates the regularization
#'   factor, produces the matrices of tip to root (\code{\link{makeL}}), and
#'   node to root distances (\code{\link{makeL1}}), the vector of ancestral
#'   state estimates, the vector of predicted phenotypes, and the rates vector
#'   for all the branches of the tree. For multivariate data, rates are given as
#'   both one vector per variable, and as a multivariate vector obtained by
#'   computing the Euclidean Norm of individual rate vectors.
#' @usage RRphylo(tree,y,cov=NULL,rootV=NULL,aces=NULL,x1=NULL,
#'   aces.x1=NULL,clus=0.5,verbose=FALSE)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous.
#' @param y either a single vector variable or a multivariate dataset. In any
#'   case, \code{y} must be named. In case of categorical variable, this should
#'   be supplied to the function as a numeric vector.
#' @param cov the covariate to be indicated if its effect on the rates must be
#'   accounted for. In this case residuals of the covariate versus the rates are
#'   used as rates. \code{'cov'} must be as long as the number of nodes plus the
#'   number of tips of the tree, which can be obtained by running \code{RRphylo}
#'   on the covariate as well, and taking the vector of ancestral states and tip
#'   values to form the covariate, as in the example below. See
#'   \href{../doc/RRphylo.html#covariate}{RRphylo vignette - covariate} for
#'   details.
#' @param rootV phenotypic value (values if multivariate) at the tree root. If
#'   \code{rootV=NULL} the function 'learns' about the root value from the 10\%
#'   tips being closest in time to the tree root, weighted by their temporal
#'   distance from the root itself (close tips phenotypes weigh more than more
#'   distant tips).
#' @param aces a named vector (or matrix if \code{y} is multivariate) of
#'   ancestral character values at nodes. Names correspond to the nodes in the
#'   tree. See \href{../doc/RRphylo.html#aces}{RRphylo vignette - aces} for
#'   details.
#' @param x1 the additional predictor(s) to be indicated to perform the multiple
#'   version of \code{RRphylo}. \code{'x1'} vector/matrix must be as long as the
#'   number of nodes plus the number of tips of the tree, which can be obtained
#'   by running \code{RRphylo} on the predictors (separately for each predictor)
#'   as well, and taking the vector of ancestral states and tip values to form
#'   the \code{x1}. See \href{../doc/RRphylo.html#predictor}{RRphylo vignette -
#'   predictor} for details.
#' @param aces.x1 a named vector/matrix of ancestral character values at nodes
#'   for \code{x1}. It must be indicated if both \code{aces} and \code{x1} are
#'   specified. Names/rownames correspond to the nodes in the tree.
#' @param clus the proportion of clusters to be used in parallel computing.
#'   Default is 0.5. To run the single-threaded version of \code{RRphylo} set
#'   \code{clus} = 0.
#' @param verbose logical indicating whether a "RRlog.txt" printing progresses
#'   should be stored into the working directory.
#' @export
#' @importFrom ape multi2di Ntip is.binary Nnode dist.nodes drop.tip subtrees
#'   nodelabels
#' @importFrom stats dist lm residuals weighted.mean
#' @importFrom stats4 mle
#' @importFrom phytools bind.tip
#' @return \strong{tree} the tree used by \code{RRphylo}. The fully dichotomous
#'   version of the tree argument. For trees with polytomies, the tree is
#'   resolved by using \code{multi2di} function in the package \pkg{ape}. Note,
#'   tip labels are ordered according to their position in the tree.
#' @return \strong{tip.path} a \eqn{n * m} matrix, where n=number of tips and
#'   m=number of branches (i.e. 2*n-1). Each row represents the branch lengths
#'   along a root-to-tip path.
#' @return \strong{node.path} a \eqn{n * n} matrix, where n=number of internal
#'   branches. Each row represents the branch lengths along a root-to-node path.
#' @return \strong{rates} single rate values computed for each branch of the
#'   tree. If \code{y} is a single vector variable, rates are equal to
#'   multiple.rates. If \code{y} is a multivariate dataset, rates are computed
#'   as the square root of the sum of squares of each row of
#'   \code{$multiple.rates}.
#' @return \strong{aces} the phenotypes reconstructed at nodes.
#' @return \strong{predicted.phenotypes} the vector of estimated tip values. It
#'   is a matrix in the case of multivariate data.
#' @return \strong{multiple.rates} a \eqn{n * m} matrix, where n= number of
#'   branches (i.e. n*2-1) and m = number of variables. For each branch, the
#'   column entries represent the evolutionary rate.
#' @return \strong{lambda} the regularization factor fitted within
#'   \code{RRphylo} by the inner function \code{optL}. With multivariate data,
#'   several \code{optL} runs are performed. Hence, the function provides a
#'   single lambda for each individual variable.
#' @return \strong{ace.values} if \code{aces} are specified, the function
#'   returns a dataframe containing the corresponding node number on the
#'   \code{RRphylo} tree  for each node , along with estimated values.
#' @return \strong{x1.rate} if \code{x1} is specified, the function returns the
#'   partial regression coefficient for \code{x1}.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @references Serio, C., Castiglione, S., Tesone, G., Piccolo, M., Melchionna,
#'   M., Mondanaro, A., Di Febbraro, M., & Raia, P.(2019). Macroevolution of
#'   toothed whales exceptional relative brain size. \emph{Evolutionary
#'   Biology}, 46: 332-342. doi:10.1007/s11692-019-09485-7
#' @references  Melchionna, M., Mondanaro, A., Serio, C., Castiglione, S., Di
#'   Febbraro, M., Rook, L.,Diniz-Filho,J.A.F., Manzi, G., Profico, A.,
#'   Sansalone, G., & Raia, P.(2020).Macroevolutionary trends of brain mass in
#'   Primates. \emph{Biological Journal of the Linnean Society}, 129: 14-25.
#'   doi:10.1093/biolinnean/blz161
#' @references  Castiglione, S., Serio, C., Mondanaro, A., Melchionna, M.,
#'   Carotenuto, F., Di Febbraro, M., Profico, A., Tamagnini, D., & Raia, P.
#'   (2020). Ancestral State Estimation with Phylogenetic Ridge Regression.
#'   \emph{Evolutionary Biology}, 47: 220-232. doi:10.1007/s11692-020-09505-x
#' @references Castiglione, S., Serio, C., Piccolo, M., Mondanaro, A.,
#'   Melchionna, M., Di Febbraro, M., Sansalone, G., Wroe, S., & Raia, P.
#'   (2020). The influence of domestication, insularity and sociality on the
#'   tempo and mode of brain size evolution in mammals. \emph{Biological Journal
#'   of the Linnean Society},132: 221-231. doi:10.1093/biolinnean/blaa186
#' @examples
#'  \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' cc<- 2/parallel::detectCores()
#'
#' # Case 1. "RRphylo" without accounting for the effect of a covariate
#' RRphylo(tree=treedino,y=massdino,clus=cc)->RRcova
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#' # "RRphylo" on the covariate in order to retrieve ancestral state values
#' c(RRcova$aces,massdino)->cov.values
#' c(rownames(RRcova$aces),names(massdino))->names(cov.values)
#'
#' RRphylo(tree=treedino,y=massdino,cov=cov.values,clus=cc)->RR
#'
#' # Case 3. "RRphylo" specifying the ancestral states
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#' DataCetaceans$aceMyst->aceMyst
#'
#' RRphylo(tree=treecet,y=masscet,aces=aceMyst,clus=cc)->RR
#'
#' # Case 4. Multiple "RRphylo"
#' library(ape)
#' drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' RRphylo(tree=treecet.multi, y=masscet.multi,clus=cc)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RR
#'
#' # Case 5. Categorical and multiple "RRphylo" with 2 additional predictors
#' library(phytools)
#'
#' set.seed(1458)
#' rtree(50)->tree
#' fastBM(tree)->y
#' jitter(y)*10->y1
#' rep(1,length(y))->y2
#' y2[sample(1:50,20)]<-2
#' names(y2)<-names(y)
#'
#' c(RRphylo(tree,y1,clus=cc)$aces[,1],y1)->x1
#'
#' RRphylo(tree,y2,clus=cc)->RRcat ### this is the RRphylo on the categorical variable
#' c(RRcat$aces[,1],y2)->x2
#'
#' cbind(c(jitter(mean(y1[tips(tree,83)])),1),
#'       c(jitter(mean(y1[tips(tree,53)])),2))->acex
#' c(jitter(mean(y[tips(tree,83)])),jitter(mean(y[tips(tree,53)])))->acesy
#' names(acesy)<-rownames(acex)<-c(83,53)
#'
#' RRphylo(tree,y,aces=acesy,x1=cbind(x1,x2),aces.x1 = acex,clus=cc)
#'
#'     }




RRphylo<-function (tree, y,
                   cov = NULL,
                   rootV = NULL,
                   aces = NULL,
                   x1=NULL,
                   aces.x1=NULL,
                   clus = 0.5,
                   verbose=FALSE){
  # library(phytools)
  # library(stats4)
  # library(ape)
  # library(parallel)
  # library(doParallel)

  if(verbose){
    sink("RRlog.txt")
    on.exit(sink())
  }

  optL <- function(lambda) {
    y <- scale(y)
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(y) - rootV)
    # aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
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

    apply(y1,2,scale)->y1
    apply(ace1,2,scale)->ace1
    ace1[1,]->acex

    y-rv->yy
    sweep(y1,2,acex)->y1
    cbind(L,y1)->L

    sweep(ace1,2,acex)->ace1
    cbind(L1,ace1)->L1


    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                t(L)) %*% (as.matrix(yy))
    y.hat <- (L %*% betas)+rootV
    aceRR <- ((L1 %*% betas[c(1:Nnode(t),(length(betas)+1-ncol(y1)):length(betas)), ]))+rootV

    abs(1-summary(lm(c(aceRR[1],y.hat)~c(rootV,y)))$coef[2])
  }

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  if (is.binary(tree))
    t <- tree else t <- multi2di(tree, random = FALSE)

  toriginal<-t
  # yoriginal <- y <- treedata(tree, y, sort = TRUE)[[2]]
  yoriginal <- y <- treedataMatch(tree, y)[[1]]

  Loriginal <-L<-makeL(t)
  L1original <-L1<-makeL1(t)

  if(!is.null(x1)){ #### multiple Ridge Regression ####
    if(is.null(nrow(x1))) x1<-as.matrix(x1)
    if(is.null(colnames(x1))) paste("pred",seq(1,ncol(x1)),sep="")->colnames(x1)
    x1[-match(t$tip.label,rownames(x1)),,drop=FALSE]->ace1
    x1[match(t$tip.label,rownames(x1)),,drop=FALSE]->y1
    if(!is.null(aces.x1)){
      as.matrix(aces.x1)->aces.x1
      colnames(x1)->colnames(aces.x1)
    }
  }

  if (!is.null(aces)) { #### phenotypes at internal nodes ####
    #L <- makeL(t)
    aceV <- aces <- as.matrix(aces)
    if (is.null(rownames(aceV)))
      stop("The matrix of aces needs to be named")

    if (!is.binary(tree)){
      sapply(1:nrow(aceV),function(i) getMRCA(t, tips(tree, as.numeric(rownames(aceV)[i]))))->rownames(aceV)
      if(!is.null(aces.x1)) sapply(1:nrow(aces.x1),function(i) getMRCA(t, tips(tree, as.numeric(rownames(aces.x1)[i]))))->rownames(aces.x1)
    }

    t$edge.length[match(rownames(aceV),t$edge[,2])]->aces.bran
    if(any(aces.bran==0)) stop(paste("Error at nodes ",paste(rownames(aceV)[which(aces.bran==0)],collapse = ", "),
                                     ": attempt to set ancestors at nodes with zero-length branches."))

    P <- aceV
    N <- as.numeric(rownames(aceV))
    tar.tips <- lapply(N, function(x) tips(t, x))
    names(tar.tips) <- N
    treeN <- t
    ynew <- y

    if(!is.null(x1)){
      if(is.vector(aces.x1)) t(as.matrix(aces.x1))->Px1 else t(aces.x1)->Px1
      y1new <- y1
      ace1new<-ace1
    }

    ftiplen<-0

    i = 1
    while (i <= length(N)) {
      nn <- getMRCA(treeN, tar.tips[[i]])
      treeN$edge.length[which(treeN$edge[,2]==nn)]->edlen
      if(edlen<=0.001) edlen/10->pp else 0.001->pp
      treeN <- bind.tip(treeN, tip.label = paste("fnd",N[i], sep = ""),
                        edge.length = ftiplen,
                        where = nn,position = pp)
      npos <- which(treeN$tip.label == paste("fnd", N[i], sep = ""))
      if (npos == 1) ynew <- rbind(unname(P[i,,drop=FALSE]), ynew)
      if (npos == nrow(ynew) + 1) ynew <- rbind(ynew, unname(P[i,,drop=FALSE])) else {
        if (npos > 1 & npos < (nrow(ynew) + 1)) ynew <- rbind(ynew[1:(npos - 1), ,drop=FALSE],
                                                              unname(P[i,,drop=FALSE]),
                                                              ynew[npos:nrow(ynew),,drop=FALSE])
      }

      if(!is.null(x1)){
        if (npos == 1) y1new <- rbind(unname(Px1[,i]), y1new)

        if (npos == nrow(y1new) + 1)  y1new <- rbind(y1new, unname(Px1[,i])) else {
          if (npos > 1 & npos < (nrow(y1new) + 1)) y1new <- rbind(y1new[1:(npos - 1), ,drop=FALSE],
                                                                  unname(Px1[,i]),
                                                                  y1new[npos:nrow(y1new),,drop=FALSE])
        }
      }

      rownames(ynew)[npos] <- paste("fnd", N[i], sep = "")
      if(is.null(x1)==FALSE) rownames(y1new)[npos] <- paste("fnd", N[i], sep = "")
      i = i + 1
    }

    if(!is.null(x1)){
      sort(N)->Ns
      Px1[,match(Ns,colnames(Px1)),drop=FALSE]->Px1

      h=1
      while(h<=length(Ns)){
        np<-which(treeN$tip.label ==
                    paste("fnd",Ns[h], sep = ""))
        (getMommy(treeN,np)[1]-(Ntip(treeN)))->npos
        if(npos== Nnode(treeN)) rbind(ace1new,Px1[,h])->ace1new else
          rbind(ace1new[1:(npos - 1), ,drop=FALSE], unname(Px1[,h]), ace1new[npos:nrow(ace1new),,drop=FALSE])->ace1new

        h=h+1
      }
      rownames(ace1new)<-seq((Ntip(treeN)+1),(Ntip(treeN)+Nnode(treeN)))
      ace1<-ace1new
      y1<-y1new
    }

    t <- treeN
    y <- ynew

    L<-makeL(t)
    L1<-makeL1(t)
  }

  internals<-unique(c(t$edge[,1],t$edge[, 2][which(t$edge[,2]>Ntip(t))]))
  edged <- data.frame(t$edge, t$edge.length)

  if(verbose) cat(paste("Initial settings DONE\n"))

  # if(ncol(y)==1) clus<-0
  if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else
    cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
  registerDoParallel(cl)

  {## No missing data
    out.pp<-NULL ### Needed for results

    if (is.null(rootV)) { #### rootV ####
      if (!is.binary(tree)) u <- data.frame(yoriginal, (1/diag(vcv(tree))^2)) else
        u <- data.frame(yoriginal,(1/diag(vcv(toriginal))^2))
      u <- u[order(u[, ncol(u)], decreasing = TRUE),]
      u1 <- u[1:(nrow(u) * 0.1), ,drop=FALSE]
      rootV <- apply(u1[, 1:(ncol(u1)-1),drop=FALSE],2,function(x)
        weighted.mean(x,u1[, dim(u1)[2]]))
    }else  if (inherits(rootV,"data.frame")) as.matrix(rootV)->rootV

    y.real <- y
    rv.real <- rootV
    res <- foreach(i = 1:ncol(y.real), .packages = c("stats4", "ape")) %dopar% {
      if(verbose) sink("RRlog.txt", append=TRUE)
      if(verbose) cat(paste("Variable",i,"- optimization started\n"))
      # for(i in 1:k){
      rootV <- rv.real[i]
      y <- y.real[, i]
      if(!is.null(x1)){ #### multiple Ridge Regression ####
        h <- mle(optLmultiple, start = list(lambda = 1), method = "L-BFGS-B",
                 upper = 10, lower = 0.001)
        if(verbose) cat(paste("Variable",i,"- optimization DONE\n"))

        h@coef->lambda
        cbind(L,sweep(y1,2,ace1[1,]))->LX
        apply(ace1,2,mean)->mean.ace1
        cbind(L1,sweep(ace1,2,mean.ace1))->LX1

        betas <- (solve(t(LX) %*% LX + lambda * diag(ncol(LX))) %*%
                    t(LX)) %*% (as.matrix(y)-rootV)

        aceRR <- ((LX1 %*% betas[c(1:Nnode(t),(length(betas)+1-ncol(y1)):length(betas)), ]))+rootV
        y.hat <- (LX %*% betas)+rootV
        betas[(length(betas)+1-ncol(y1)):length(betas),]->x1.rate
        betas[1:(length(betas)-ncol(y1)),,drop=FALSE]->betas
        colnames(betas)<-NULL
        if(verbose) cat("Variable",i,"- rates and aces estimation DONE\n")

        list(aceRR, betas, y.hat, lambda,x1.rate)

      }else{
        h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",
                 upper = 10, lower = 0.001)
        if(verbose) cat("Variable",i,"- optimization DONE\n")

        lambda <- h@coef
        betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                    t(L)) %*% (as.matrix(y) - rootV)
        aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
        y.hat <- (L %*% betas) + rootV
        if(verbose) cat(paste("Variable",i,"- rates and aces estimation DONE\n"))
        list(aceRR, betas, y.hat, lambda)
      }
    }
    y <- y.real
  }
  stopCluster(cl)

  aceRR <- do.call(cbind, lapply(res, "[[", 1))
  betas <- do.call(cbind, lapply(res, "[[", 2))
  y.hat <- do.call(cbind, lapply(res, "[[", 3))
  lambda <- sapply(res, "[[", 4)
  if(!is.null(x1)){
    x1.rate <- sapply(res,"[[",5)
    if(!is.null(nrow(x1.rate))) colnames(x1.rate)<-colnames(y) else names(x1.rate)<-colnames(y)
  }

  rownames(betas) <- colnames(L)
  rownames(y.hat) <- rownames(y)
  rownames(aceRR) <- colnames(L1)
  colnames(betas) <- colnames(y.hat) <- colnames(aceRR) <- colnames(y)

  if (!is.null(aces)) {
    tip.rem <- paste("fnd", N, sep = "")
    mapply(a=tip.rem,b=tar.tips,function(a,b) getMRCA(t, c(a, b)))->nod.rem
    t <- drop.tip(t, tip.rem)

    if(!is.null(out.pp)){
      out.pp[-match(c(nod.rem, tip.rem), rownames(out.pp)),,drop=FALSE]->out.pp
      rownames(out.pp)[1:Nnode(t)]<-seq((Ntip(t)+1),(Ntip(t)+Nnode(t)),1)
    }

    aceRR <- aceRR[-match(nod.rem, rownames(aceRR)),,drop=FALSE]
    betas <- betas[-match(c(nod.rem, tip.rem),rownames(betas)),,drop=FALSE]
    y.hat <- y.hat[-match(tip.rem, rownames(y.hat)),,drop=FALSE]

    rownames(betas)[1:Nnode(t)] <- rownames(aceRR) <- seq((Ntip(t) +
                                                             1), (Ntip(t) + Nnode(t)), 1)

    ace.est <- aceRR[match(rownames(aceV), rownames(aceRR)),]
    #if (nrow(aceV) < 2)
    ace.estimates <- data.frame(real.node = rownames(aces),
                                RRnode = rownames(aceV), ace.estimate =unname(ace.est))
    # else ace.estimates <- data.frame(real.node = rownames(aces),
    #                                  RRnode = rownames(aceV), ace.estimate = unname(ace.est))

  }
  betasREAL <- betas

  if (!is.null(cov)){
    cov[match(rownames(betas),names(cov))]->cov
    #### Covariate multi ####
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
  }

  if(ncol(y)>1) rates <- as.matrix(apply(betas, 1, function(x) sqrt(sum(x^2)))) else rates<-betas


  res <- list(t, Loriginal, L1original, rates, aceRR, y.hat, betas,
              lambda)
  names(res) <- c("tree", "tip.path", "node.path", "rates",
                  "aces", "predicted.phenotype", "multiple.rates",
                  "lambda")

  if(!is.null(aces)) res <- c(res, ace.values=list(ace.estimates))
  if(!is.null(x1)) res <- c(res, x1.rate=list(x1.rate))
  if(!is.null(out.pp)) res <- c(res, phylopars.res=list(as.matrix(out.pp)))

  return(res)
}
