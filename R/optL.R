

optL <- function(lambda) {
  makeL(t)->L
  makeL1(t)->L1
  get("t")->t
  get("y")->y
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
