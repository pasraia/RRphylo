#' @title Matrix of branch lengths along root-to-tip paths
#' @description This function produces a \eqn{n * m} matrix, where n=number of
#'   tips and m=number of branches (i.e. n + number of nodes). Each row
#'   represents the branch lengths aligned along a root-to-tip path.
#' @usage makeL(tree)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @export
#' @importFrom stats median
#' @return The function returns a \eqn{n * m} matrix of branch lengths for all
#'   root-to-tip paths in the tree (one per species).
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' makeL(tree=Tstage)


makeL<-function(tree){

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  internals <- unique(c(tree$edge[, 1], tree$edge[, 2][which(tree$edge[,
                                                                       2] > Ntip(tree))]))


  tippa <- list()
  for (i in 1:length(internals)) {
    tippas <- tips(tree, internals[i])
    dato <- data.frame(rep(internals[i], length(tippas)),
                       tippas)
    colnames(dato)[1] <- "node"
    tippa[[i]] <- dato
  }
  Tstr <- do.call(rbind, tippa)



  L <- matrix(nrow = Ntip(tree), ncol = length(tree$edge.length) +
                1)
  rownames(L) <- tree$tip.label
  colnames(L) <- c(internals, tree$tip.label)
  edged <- data.frame(tree$edge, tree$edge.length)
  order <- edged[, 2][which(edged[, 2] < Ntip(tree) + 1)]
  labs <- tree$tip.label[order]
  edged[, 2][which(edged[, 2] < Ntip(tree) + 1)] <- labs
  tip.path <- list()
  for (i in 1:dim(L)[1]) {
    tip.path[[i]] <- c(Tstr[, 1][which(Tstr[, 2] %in% rownames(L)[i])],
                       which(tree$tip.label == rownames(L)[i]))
  }
  for (j in 1:length(tip.path)) {
    a <- list()
    for (i in 2:length(tip.path[[j]]) - 1) {
      a[[i]] <- tip.path[[j]][c(i, i + 1)]
    }
    b <- do.call(rbind, a)
    b[which(b[, 2] < Ntip(tree) + 1), 2] <- tree$tip.label[b[which(b[,
                                                                     2] < Ntip(tree) + 1), 2]]
    L.match <- b[, 2]
    br.len <- array()
    for (k in 1:dim(b)[1]) {
      br.len[k] <- edged[b[k, 1] == edged[, 1] & b[k, 2] ==
                           edged[, 2], ][, 3]
    }
    d <- data.frame(L.match, br.len)
    L[j, match(d[, 1], colnames(L))] <- d[, 2]
  }
  if (is.null(tree$root.edge) || tree$root.edge==0) L[, 1] <- median(tree$edge.length) else L[, 1] <- tree$root.edge
  L[which(is.na(L))] <- 0
  return(L)
}
