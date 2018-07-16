#' @title Searching for evolutionary trends in phenotypes and rates
#' @description This function searches for evolutionary trends in phenotypic mean and evolutionary rates .
#' @usage search.trend(RR,y,nsim=100,clus=.5,node=NULL,cov=NULL,foldername,ConfInt=c(FALSE,TRUE))
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y the vector (or matrix if multivariate) of phenotypes.
#' @param node number of nodes to be specifically tested for trends. It is NULL by default. Notice the node number must refer to the dichotomic version of the original tree, as produced by \code{RRphylo}.
#' @param nsim number of simulations to be performed. It is set at 100 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @param foldername the path of the folder where plots are to be found.
#' @param cov the covariate values to be specified if the RR object has been created with covariate. As for \code{RRphylo}, \code{'cov'} must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate (see the example below).
#' @param ConfInt if \code{TRUE}, the function returns 95\% confidence intervals for slopes of phenotype versus age, absolute rates versus age, and relative rates versus age regressions.
#' @return The function returns a ‘list’ object including:
#' @return \strong{$rbt} for each branch of the tree, there are the age of the daughter node/tip (D.age), the age of the parental node (P.age), the \code{RRphylo} rates, and the distance from the tree root (age). If y is multivariate, it also includes the multiple rates for each y component.
#' @return \strong{$pbt} a data frame of phenotypic values and distance from the tree root for each node and tip.
#' @return \strong{$p.trend} results of phenotype versus age regression. It reports a p-value for the regression between the variables (p.real), a p-value for the difference from 0 under standard major axis regression (p.sma0), a p-value computed contrasting the real slope to Brownian motion simulations (p.random), a global p-value corresponding to the least significant between p.real and p.random (p.value; see details for further explanations).
#' @return \strong{$rbt.rateA} results of absolute rate values versus age regression. It reports a p-value for the regression between the variables (p.real), a p-value for the difference from 0 under standard major axis regression (p.sma0), a p-value computed contrasting the real slope to Brownian motion simulations (p.random), and a global p-value (p.value, corresponding to p.random; see details for further explanations).
#' @return \strong{$rbt.rateR} results of rate values versus age regression. It reports a p-value for the regression between the variables (p.real), a p-value for the difference from 0 under standard major axis regression (p.sma0), a p-value computed contrasting the real slope to Brownian motion simulations (p.random), and a global p-value (p.value, corresponding to p.random; see details for further explanations).
#' @return \strong{$ConfInts} the 95\% confidence intervals around the slopes of phenotype versus age ($phenotype), absolute rates versus age ($abs.rate), and relative rates versus age ($rel.rate) regressions.
#' @return If the node argument is specified, the list also includes \strong{$p.trend.nodes}, \strong{$rbt.rateA.nodes}, \strong{$rbt.rateR.nodes}, which return the same results as above for each specified node. Finally, the \strong{$SMA} object contains the comparisons between slopes of regression lines of each pair of nodes, for all the regressions previously performed, under standard major axis regression.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @details The function simultaneously returns the regression of phenotypes and phenotypic evolutionary rates against age (in terms of distance from the tree root). As an output, plots for the absolute and relative rates regressions versus age are saved as a single .pdf file. In the plots, the 95\% confidence intervals of phenotypes and rates simulated under the Brownian motion for each node are plotted as shaded areas. The same is performed for the phenotype versus age regression. Regression slopes are reported for all regressions. To assess significance, slopes are compared to a family of simulated slopes (BMslopes, where the number of simulations is equal to \code{nsim}), generated as to show no phenotypic or rate trend (i.e. Brownian motion), using the \code{\link{setBM}} function. In the simulations, the phenotypic value at the root is the same as calculated for the original data by \code{RRphylo}. In addition, the Brownian rate is calculated for the original data and tree and set for all the Brownian motion simulations. Eventually, a p-value is obtained by comparing the real slope to BMslopes.
#' @importFrom graphics points text title polygon
#' @importFrom stats as.formula coef resid density predict
#' @importFrom binr bins.greedy
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nlme gls varFunc
#' @export
#' @examples
#'  \donttest{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#'
#' # Extract Pterosaurs tree and data
#'   library(ape)
#'   extract.clade(treedino,748)->treeptero
#'   massdino[match(treeptero$tip.label,names(massdino))]->massptero
#'   massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' # Case 1. "RRphylo" whitout accounting for the effect of a covariate
#'   RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'
#'  # Case 1.1. "search.trend" whitout indicating nodes to be tested for trends
#'    search.trend(RR=RRptero, y=log(massptero), nsim=100, clus=0.5,
#'    foldername=tempdir(),cov=NULL,ConfInt=FALSE,node=NULL)
#'
#'  # Case 1.2. "search.trend" indicating nodes to be specifically tested for trends
#'     search.trend(RR=RRptero, y=log(massptero), nsim=100, node=143, clus=0.5,
#'     foldername=tempdir(),cov=NULL,ConfInt=FALSE)
#'
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#'  # "RRphylo" on the covariate in order to retrieve ancestral state values
#'    RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'    c(RRptero$aces,log(massptero))->cov.values
#'    names(cov.values)<-c(rownames(RRptero$aces),names(massptero))
#'    RRphylo(tree=treeptero,y=log(massptero),cov=cov.values)->RRpteroCov
#'
#'  # Case 2.1. "search.trend" whitout indicating nodes to be tested for trends
#'    search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, clus=0.5,
#'    foldername=tempdir(),ConfInt=FALSE,cov=cov.values)
#'
#'  # Case 2.2. "search.trend" indicating nodes to be specifically tested for trends
#'    search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, node=143, clus=0.5,
#'    foldername=tempdir(),ConfInt=FALSE,cov=cov.values)}

search.trend<-function (RR, y, nsim = 100, clus = 0.5, node = NULL, cov = NULL,
                        foldername, ConfInt = c(FALSE, TRUE))
{
  #require(ape)
  #require(phytools)
  #require(geiger)
  #require(stats4)
  #require(foreach)
  #require(doParallel)
  #require(lmtest)
  #require(parallel)
  #require(smatr)
  #require(binr)
  #require(nlme)
  #require(RColorBrewer)

  t <- RR$tree
  if (length(y) > Ntip(t)) {
    if (density(diag(vcv(t)))$bw/max(nodeHeights(t)) < 0.08)
      warning("trend regression test might have low power")
  }else {
    if (density(diag(vcv(t)))$bw/max(nodeHeights(t)) < 0.07)
      warning("trend regression might have low power for a small tree or a modest effect")
  }
  rates <- RR$rates
  betas <- RR$multiple.rates
  aceRR <- RR$aces
  L <- RR$tip.path
  L1 <- RR$node.path
  if (class(y) == "data.frame")
    y <- treedata(t, y, sort = TRUE)[[2]]
  H <- max(nodeHeights(t))
  eds <- t$edge[, 2]
  eds[which(t$edge[, 2] < Ntip(t) + 1)] <- t$tip.label
  eds <- c(Ntip(t) + 1, eds)
  hh <- c(0, nodeHeights(t)[, 2])
  eds <- data.frame(leaf = eds, height = hh)
  if (length(y) > Ntip(t)) {
    rates <- as.data.frame(rates)
    betas <- as.data.frame(betas)
    data <- data.frame(betas = betas[match(eds[, 1], rownames(betas)),
                                     ], rate = rates[match(eds[, 1], rownames(rates)),
                                                     ], age = eds[, 2])
    colnames(data)[1:dim(y)[2]] <- paste("betas", seq(1,
                                                      dim(y)[2], 1), sep = "")
  }else {
    data <- data.frame(rate = rates[match(eds[, 1], rownames(rates)),
                                    ], age = eds[, 2])
    rownames(data) <- rownames(rates)[match(eds[, 1], rownames(rates))]
  }
  B.age <- data.frame(t$edge, nodeHeights(t), t$edge.length)
  B.age <- data.frame(B.age, H - B.age[, 3])
  names(B.age) <- c("parent", "daughter", "rootdist.p", "rootdist.d",
                    "PD.dist", "P.age")
  B.age <- data.frame(B.age, B.age$P.age - B.age$PD.dist)
  colnames(B.age)[7] <- "D.age"
  B.age$D.age <- jitter(B.age$D.age)
  D.death <- findInterval(B.age$D.age, bins.greedy(B.age$D.age,
                                                   nbins = Ntip(t)/3, minpts = 3)$binhi)
  B.age <- data.frame(B.age, D.death)
  b.distrib <- B.age[, c(2, 7, 6, 8)]
  b.distrib[which(b.distrib$daughter < Ntip(t) + 1), 1] <- t$tip.label[b.distrib[which(b.distrib$daughter <
                                                                                         Ntip(t) + 1), 1]]
  rbind(c((Ntip(t)+1),L[1,1],0,0),b.distrib)->b.distrib
  data[,dim(data)[2]]+L[1,1]->data[,dim(data)[2]]
  data <- data.frame(b.distrib, data)
  rbi <- data
  if (length(y) > Ntip(t)) {
    rbi.rate <- rbi[, c(5:(dim(rbi)[2] - 1), dim(rbi)[2])]
    rbi.slopeA <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                             1))
    rbi.slopeR <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                             1))

    REG.betas <- list()
    REGabs.betas <- list()
    for (i in 1:(dim(rbi.rate)[2] - 1)) {
      bet <- rbi.rate[, i]
      age <- rbi.rate[, dim(rbi.rate)[2]]
      aar <- try(gls(bet ~ age, weights = varFunc(~age),
                     control = list(singular.ok = TRUE)))
      if (class(aar) == "try-error")
        aar <- lm(bet ~ age)
      rbi.slopeR[i, ] <- coef(summary(aar))[2, c(1, 4)]
      aaa <- try(gls(abs(bet) ~ age, weights = varFunc(~age),
                     control = list(singular.ok = TRUE)))
      if (class(aaa) == "try-error")
        aaa <- lm(abs(bet) ~ age)
      rbi.slopeA[i, ] <- coef(summary(aaa))[2, c(1, 4)]
      dat <- data.frame(bet, age)
      REG.betas[[i]] <- aar
      REGabs.betas[[i]] <- aaa
    }
    colnames(rbi.slopeR) <- colnames(rbi.slopeA) <- c("slope",
                                                      "p-value")
    rownames(rbi.slopeR) <- rownames(rbi.slopeA) <- names(REG.betas) <- names(REGabs.betas) <- colnames(data)[5:(5 +
                                                                                                                   dim(y)[2])]
  }else {
    rbi.rate <- rbi[, c(5, 6)]
    aar <- try(gls(rate ~ age, data = rbi.rate, weights = varFunc(~age),
                   control = list(singular.ok = TRUE)))
    if (class(aar) == "try-error")
      aar <- lm(rate ~ age, data = rbi.rate)
    rbi.slopeR <- coef(summary(aar))[2, c(1, 4)]
    aaa <- try(gls(abs(rate) ~ age, data = rbi.rate, weights = varFunc(~age),
                   control = list(singular.ok = TRUE)))
    if (class(aaa) == "try-error")
      aaa <- lm(abs(rate) ~ age, data = rbi.rate)
    rbi.slopeA <- coef(summary(aaa))[2, c(1, 4)]
    REG <- aar
    REGabs <- aaa
  }
  nodes <- aceRR[1:Nnode(t), ]
  if (length(y) > Ntip(t)) {
    colnames(nodes)[1:dim(y)[2]] <- paste("y", seq(1, dim(y)[2]),
                                          sep = "")
    P <- rbind(nodes, y)
    PP <- data.frame(P[match(rbi[, 1], rownames(P)), ],
                     rbi$age)
    colnames(PP)[dim(PP)[2]] <- "age"
    trend.reg <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x) summary(lm(x ~
                                                                           PP[, dim(PP)[2]])))
    trend.M <- summary(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                 collapse = "+"), colnames(PP)[length(colnames(PP))],
                                           sep = "~")), data = PP))
    names(trend.reg) <- paste("y", seq(1:dim(y)[2]), sep = "")
    trend.reg[[length(trend.reg) + 1]] <- trend.M
    trend.reg <- lapply(trend.reg, coefficients)
    names(trend.reg)[length(trend.reg)] <- "multiple"
  }else {
    P <- c(nodes, y)
    PP <- data.frame(P[match(rbi[, 1], names(P))], rbi$age)
    colnames(PP) <- c("phenotype", "age")
    trend.reg <- summary(lm(PP))$coeff
  }
  PPtot <- PP
  if (class(node) != "NULL") {
    rbi.sma <- rbi.rate
    PP.sma <- PP
    rbi.sma$group <- rep("NA", dim(rbi.sma)[1])
    rbi.slopeR.sel <- list()
    rbi.slopeA.sel <- list()
    trend.reg.sel <- list()
    REG.betas.y.sel <- list()
    REGabs.betas.y.sel <- list()
    REG.betas.age.sel <- list()
    trend.reg.age.sel <- list()
    trend.reg.y.sel <- list()
    for (j in 1:length(node)) {
      n <- node[j]
      sele <- getDescendants(t, n)
      sele[which(sele < (Ntip(t) + 1))] <- t$tip.label[sele[which(sele <
                                                                    (Ntip(t) + 1))]]
      rbi.sma[match(sele, rownames(rbi.sma)), ]$group <- paste("g",
                                                               n, sep = "")
      rbi.sel <- rbi[match(sele, rownames(rbi)), ]
      if (length(y) > Ntip(t)) {
        rbi.rate <- rbi.sel[, c(5:(dim(rbi.sel)[2] -
                                     1), dim(rbi.sel)[2])]
        rbi.slopeAn <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                  1))
        rbi.slopeRn <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                  1))
        REG.betas.y <- list()
        REG.betas.age <- list()
        REGabs.betas.y <- list()
        for (i in 1:(dim(rbi.rate)[2] - 1)) {
          bet <- rbi.rate[, i]
          age <- rbi.rate[, dim(rbi.rate)[2]]
          aar <- try(gls(bet ~ age, weights = varFunc(~age),
                         control = list(singular.ok = TRUE)))
          if (class(aar) == "try-error")
            aar <- lm(bet ~ age)
          rbi.slopeRn[i, ] <- coef(summary(aar))[2,
                                                 c(1, 4)]
          aaa <- try(gls(abs(bet) ~ age, weights = varFunc(~age),
                         control = list(singular.ok = TRUE)))
          if (class(aaa) == "try-error")
            aaa <- lm(abs(bet) ~ age)
          rbi.slopeAn[i, ] <- coef(summary(aaa))[2,
                                                 c(1, 4)]
          dat <- data.frame(bet, age)
          REG.betas.y[[i]] <- predict(aar)
          REGabs.betas.y[[i]] <- predict(aaa)
          REG.betas.age[[i]] <- age
        }
        colnames(rbi.slopeRn) <- colnames(rbi.slopeAn) <- c("slope",
                                                            "p-value")
        rownames(rbi.slopeRn) <- rownames(rbi.slopeAn) <- names(REG.betas.y) <- names(REG.betas.age) <- names(REGabs.betas.y) <- colnames(data)[5:(5 +
                                                                                                                                                     dim(y)[2])]
      }else {
        rbi.rate <- rbi.sel[, 5:6]
        aar <- try(gls(rate ~ age, data = rbi.rate,
                       weights = varFunc(~age), control = list(singular.ok = TRUE)))
        if (class(aar) == "try-error")
          aar <- lm(rate ~ age, data = rbi.rate)
        rbi.slopeRn <- coef(summary(aar))[2, c(1, 4)]
        aaa <- try(gls(abs(rate) ~ age, data = rbi.rate,
                       weights = varFunc(~age), control = list(singular.ok = TRUE)))
        if (class(aaa) == "try-error")
          aaa <- lm(abs(rate) ~ age, data = rbi.rate)
        rbi.slopeAn <- coef(summary(aaa))[2, c(1, 4)]
        REG.betas.y <- predict(aar)
        REGabs.betas.y <- predict(aaa)
        REG.betas.age <- rbi.rate$age
      }
      rbi.slopeR.sel[[j]] <- rbi.slopeRn
      rbi.slopeA.sel[[j]] <- rbi.slopeAn
      REG.betas.y.sel[[j]] <- REG.betas.y
      REGabs.betas.y.sel[[j]] <- REGabs.betas.y
      REG.betas.age.sel[[j]] <- REG.betas.age
      nodes <- aceRR[1:Nnode(t), ]
      if (length(y) > Ntip(t)) {
        colnames(nodes)[1:dim(y)[2]] <- paste("y", seq(1,
                                                       dim(y)[2]), sep = "")
        P <- rbind(nodes, y)
        PP <- data.frame(P[match(rbi.sel[, 1], rownames(P)),
                           ], rbi.sel$age)
        colnames(PP)[dim(PP)[2]] <- "age"
        trend.regC <- apply(PP[1:(dim(PP)[2] - 1)],
                            2, function(x) summary(lm(x ~ PP[, dim(PP)[2]])))
        trend.reg.y.sel[[j]] <- apply(PP[1:(dim(PP)[2] -
                                              1)], 2, function(x) predict(lm(x ~ PP[, dim(PP)[2]])))
        trend.M <- summary(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                     collapse = "+"), colnames(PP)[length(colnames(PP))],
                                               sep = "~")), data = PP))
        names(trend.regC) <- paste("y", seq(1:dim(y)[2]),
                                   sep = "")
        trend.regC[[length(trend.regC) + 1]] <- trend.M
        trend.regC <- lapply(trend.regC, coefficients)
        names(trend.regC)[length(trend.regC)] <- "multiple"
      }else {
        P <- c(nodes, y)
        PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
                         rbi.sel$age)
        colnames(PP) <- c("phenotype", "age")
        trend.regC <- summary(lm(PP))$coeff
        trend.reg.age.sel[[j]] <- PP$age
        trend.reg.y.sel[[j]] <- predict(lm(PP))
      }
      trend.reg.sel[[j]] <- trend.regC
    }
    names(rbi.slopeR.sel) <- names(rbi.slopeA.sel) <- names(trend.reg.sel) <- node
    rbi.sma$group[which(rbi.sma$group == "NA")] <- "others"
    PP.sma <- cbind(PP.sma, group = rbi.sma[match(rownames(PP.sma),
                                                  rownames(rbi.sma)), ]$group)
    if (length(which(rbi.sma$group == "others")) < 3)
      rbi.sma <- rbi.sma[-which(rbi.sma$group == "others"),
                         ]
    if (length(which(PP.sma$group == "others")) < 3)
      PP.sma <- PP.sma[-which(PP.sma$group == "others"),
                       ]
    rbiRES <- rbi.sma
    if (length(y) > Ntip(t)) {
      sma.resA <- list()
      sma.resR <- list()
      sma.resPP <- list()
      for (i in 1:(dim(rbi.sma)[2] - 2)) {
        group <- rbi.sma[, (dim(rbi.sma)[2])]
        age <- rbi.sma[, (dim(rbi.sma)[2] - 1)]
        bets <- rbi.sma[, i]
        dat <- data.frame(bets, age, group)
        if (length(unique(group)) < 3) {
          SMA <- sma(bets ~ age * group, data = dat,
                     multcomp = FALSE)
          tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                            colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                            SMA$commoncoef$LR, SMA$commoncoef$df, SMA$commoncoef$bs[1,
                                                                                    1], SMA$commoncoef$bs[1, 2])
          colnames(tem) <- c("group_1", "group_2", "Pval",
                             "TestStat", "df", "Slope1", "Slope2")
          sma.resR[[i]] <- tem
          SMA <- sma((abs(bets)) ~ age * group, data = dat,
                     multcomp = FALSE)
          tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                            colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                            SMA$commoncoef$LR, SMA$commoncoef$df, SMA$commoncoef$bs[1,
                                                                                    1], SMA$commoncoef$bs[1, 2])
          colnames(tem) <- c("group_1", "group_2", "Pval",
                             "TestStat", "df", "Slope1", "Slope2")
          sma.resA[[i]] <- tem
        }else {
          sma <- sma(bets ~ age * group, data = dat,
                     multcomp = TRUE)
          sma.resR[[i]] <- sma$multcompresult
          sma((abs(bets)) ~ age * group, data = dat,
              multcomp = TRUE)
          sma.resA[[i]] <- sma$multcompresult
        }
        groupPP <- PP.sma[, (dim(PP.sma)[2])]
        agePP <- PP.sma[, (dim(PP.sma)[2] - 1)]
        if (i <= (dim(PP.sma)[2] - 2)) {
          betPP <- PP.sma[, i]
          datPP <- data.frame(betPP, agePP, groupPP)
          if (length(unique(groupPP)) < 3) {
            SMA <- sma(betPP ~ agePP * groupPP, data = datPP,
                       multcomp = FALSE)
            tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                              colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                              SMA$commoncoef$LR, SMA$commoncoef$df,
                              SMA$commoncoef$bs[1, 1], SMA$commoncoef$bs[1,
                                                                         2])
            colnames(tem) <- c("group_1", "group_2",
                               "Pval", "TestStat", "df", "Slope1", "Slope2")
            sma.resPP[[i]] <- tem
          }else {
            sma <- sma(betPP ~ agePP * groupPP, data = datPP,
                       multcomp = TRUE)
            sma.resPP[[i]] <- sma$multcompresult
          }
        }
      }
    }else {
      if (length(unique(rbi.sma$group)) < 3) {
        SMA <- sma(rate ~ age * group, data = rbi.sma,
                   multcomp = FALSE)
        tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                          colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                          SMA$commoncoef$LR, SMA$commoncoef$df, SMA$commoncoef$bs[1,
                                                                                  1], SMA$commoncoef$bs[1, 2])
        colnames(tem) <- c("group_1", "group_2", "Pval",
                           "TestStat", "df", "Slope1", "Slope2")
        sma.resR <- tem
        SMA <- sma((abs(rate)) ~ age * group, data = rbi.sma,
                   multcomp = FALSE)
        tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                          colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                          SMA$commoncoef$LR, SMA$commoncoef$df, SMA$commoncoef$bs[1,
                                                                                  1], SMA$commoncoef$bs[1, 2])
        colnames(tem) <- c("group_1", "group_2", "Pval",
                           "TestStat", "df", "Slope1", "Slope2")
        sma.resA <- tem
        colnames(PP.sma)[1] <- "y"
        SMA <- sma(y ~ age * group, data = PP.sma, multcomp = FALSE)
        tem <- data.frame(colnames(SMA$commoncoef$bs)[1],
                          colnames(SMA$commoncoef$bs)[2], SMA$commoncoef$p,
                          SMA$commoncoef$LR, SMA$commoncoef$df, SMA$commoncoef$bs[1,
                                                                                  1], SMA$commoncoef$bs[1, 2])
        colnames(tem) <- c("group_1", "group_2", "Pval",
                           "TestStat", "df", "Slope1", "Slope2")
        sma.resPP <- tem
      }else {
        sma <- sma(rate ~ age * group, data = rbi.sma,
                   multcomp = TRUE)
        sma.resR <- sma$multcompresult
        sma <- sma((abs(rate)) ~ age * group, data = rbi.sma,
                   multcomp = TRUE)
        sma.resA <- sma$multcompresult
        colnames(PP.sma)[1] <- "y"
        sma <- sma(y ~ age * group, data = PP.sma, multcomp = TRUE)
        sma.resPP <- sma$multcompresult
      }
    }
    if (class(sma.resA) == "list") {
      names(sma.resR) <- names(sma.resA) <- colnames(rbi.sma)[1:(dim(rbi.sma)[2] -
                                                                   2)]
      names(sma.resPP) <- names(trend.reg)[-length(trend.reg)]
    }
  }
  s2 <- ratematrix(t, y)
  RR$aces[1,]->a
  if (length(y) > Ntip(t)) {
    yy <- list()
    s2 <- diag(s2)
    for (i in 1:dim(y)[2]) {
      yy[[i]] <- setBM(t, s2 = s2[i], a = a[i], nY = nsim, type = "brown")
    }
  }else {
    s2 <- (s2[1, 1])
    yy <- setBM(t, s2 = s2, a = a, nY = nsim, type = "brown")
  }
  res <- list()
  cl <- makeCluster(round((detectCores() * clus), 0))
  registerDoParallel(cl)
  res <- foreach(i = 1:nsim, .packages = c("nlme", "ape",
                                           "geiger", "phytools", "penalized", "doParallel", "lmtest",
                                           "smatr")) %dopar% {
                                             gc()
                                             if (length(y) > Ntip(t)) {
                                               vec <- seq(1:nsim)
                                               y <- lapply(yy, function(x) x[, sample(vec, 1, replace = FALSE)])
                                               y <- do.call(cbind, y)

                                               betas<- matrix(ncol=dim(y)[2], nrow=Ntip(t)+Nnode(t))
                                               aceRR<- matrix(ncol=dim(y)[2], nrow=Nnode(t))
                                               for (s in 1:dim(y)[2]){
                                                 rootV<-a[s]
                                                 lambda <- RR$lambda[s]
                                                 m.betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                                                               t(L)) %*% (as.matrix(y[,s]) - rootV)
                                                 aceRR[,s] <- (L1 %*% m.betas[1:Nnode(t), ]) + rootV
                                                 m.betas->betas[,s]

                                               }
                                               rownames(betas)<-rownames(RR$rates)
                                               rownames(aceRR)<-rownames(RR$aces)
                                             }else {
                                               y <- yy[, i]

                                               rootV<-a
                                               lambda <- RR$lambda
                                               betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                                                           t(L)) %*% (as.matrix(y) - rootV)
                                               aceRR <- (L1 %*% betas[1:Nnode(t), ]) + rootV
                                             }



                                             if (is.null(cov) == FALSE) {
                                               if (length(y) > Ntip(t)) {
                                                 if (length(which(apply(betas, 1, sum) == 0)) >
                                                     0) {
                                                   zeroes <- which(apply(betas, 1, sum) == 0)
                                                   R <- log(abs(betas))
                                                   R <- R[-zeroes, ]
                                                   Y <- abs(cov)
                                                   Y <- Y[-zeroes]
                                                   resi <- residuals(lm(R ~ Y))
                                                   factOut <- which(apply(betas, 1, sum) != 0)
                                                   betas[factOut, ] <- resi
                                                   betas[zeroes, ] <- 0
                                                 }else {
                                                   R <- log(abs(betas))
                                                   Y <- abs(cov)
                                                   resi <- residuals(lm(R ~ Y))
                                                   betas <- as.matrix(resi)
                                                 }
                                               }else {
                                                 if (length(which(betas == "0")) > 0) {
                                                   zeroes <- which(betas == "0")
                                                   R <- log(abs(betas))
                                                   R <- R[-zeroes]
                                                   Y <- abs(cov)
                                                   Y <- Y[-zeroes]
                                                   resi <- residuals(lm(R ~ Y))
                                                   factOut <- which(betas != "0")
                                                   betas[factOut] <- resi
                                                   betas[zeroes] <- 0
                                                 }else {
                                                   R <- log(abs(betas))
                                                   Y <- abs(cov)
                                                   resi <- residuals(lm(R ~ Y))
                                                   betas <- as.matrix(resi)
                                                 }
                                               }
                                             }

                                             nodes <- aceRR[1:Nnode(t), ]
                                             if (length(y) > Ntip(t)) {
                                               colnames(nodes) <- colnames(y)
                                               P <- rbind(nodes, y)
                                               PP <- data.frame(P[match(rbi[, 1], rownames(P)),
                                                                  ], rbi$age)
                                               colnames(PP)[dim(PP)[2]] <- "age"
                                               trend.regR <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x) summary(lm(x ~
                                                                                                                       PP[, dim(PP)[2]])))
                                               trend.M <- summary(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                                                            collapse = "+"), colnames(PP)[length(colnames(PP))],
                                                                                      sep = "~")), data = PP))
                                               names(trend.regR) <- paste("betas", seq(1:dim(y)[2]),
                                                                          sep = "")
                                               trend.regR[[length(trend.regR) + 1]] <- trend.M
                                               trend.regS <- lapply(trend.regR, coefficients)
                                               names(trend.regS)[length(trend.regS)] <- "multiple"

                                             }else {
                                               P <- c(nodes, y)
                                               PP <- data.frame(P[match(rbi[, 1], names(P))], rbi$age)
                                               colnames(PP) <- c("phenotype", "age")
                                               trend.regS <- lm(PP)
                                               trend.regS <- coefficients(trend.regS)
                                             }
                                             PPtot <- PP
                                             eds <- t$edge[, 2]
                                             eds[which(t$edge[, 2] < Ntip(t) + 1)] <- t$tip.label
                                             eds <- c(Ntip(t) + 1, eds)
                                             hh <- c(0, nodeHeights(t)[, 2])
                                             eds <- data.frame(leaf = eds, height = hh)

                                             if (length(y) > Ntip(t)) {
                                               rates <- as.data.frame(apply(betas, 1, function(x) sqrt(sum(x^2))))
                                               data <- data.frame(betas = betas[match(eds[, 1],
                                                                                      rownames(betas)), ], rate = rates[match(eds[,
                                                                                                                                  1], rownames(rates)), ], age = eds[, 2])
                                               colnames(data)[1:dim(y)[2]] <- paste("betas", seq(1,
                                                                                                 dim(y)[2], 1), sep = "")
                                             }else {
                                               rates <- betas
                                               data <- data.frame(rate = rates[match(eds[, 1],
                                                                                     rownames(rates)), ], age = eds[, 2])
                                             }
                                             data[,dim(data)[2]]+L[1,1]->data[,dim(data)[2]]
                                             data <- data.frame(b.distrib, data)
                                             rbi <- data
                                             if (length(y) > Ntip(t)) {
                                               rbi.rate <- rbi[, c(5:(dim(rbi)[2] - 1), dim(rbi)[2])]
                                               rbi.slopeAS <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                                                         1))
                                               rbi.slopeRS <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                                                         1))
                                               for (k in 1:(dim(rbi.rate)[2] - 1)) {
                                                 bet <- rbi.rate[, k]
                                                 age <- rbi.rate[, dim(rbi.rate)[2]]
                                                 betage <- data.frame(bet, age)
                                                 aar <- try(gls(bet ~ age, weights = varFunc(~age),
                                                                control = list(singular.ok = TRUE)))
                                                 if (class(aar) == "try-error")
                                                   aar <- lm(bet ~ age)
                                                 rbi.slopeRS[k, ] <- coef(summary(aar))[2, c(1,
                                                                                             4)]
                                                 aaa <- try(gls(abs(bet) ~ age, weights = varFunc(~age),
                                                                control = list(singular.ok = TRUE)))
                                                 if (class(aaa) == "try-error")
                                                   aaa <- lm(abs(bet) ~ age)
                                                 rbi.slopeAS[k, ] <- coef(summary(aaa))[2, c(1,
                                                                                             4)]
                                                 REGbetas <- aar
                                               }
                                               colnames(rbi.slopeRS) <- colnames(rbi.slopeAS) <- c("slope",
                                                                                                   "p-value")
                                               rownames(rbi.slopeRS) <- rownames(rbi.slopeAS) <- colnames(data)[5:(5 +
                                                                                                                     dim(y)[2])]
                                             }else {
                                               rbi.rate <- rbi[, 5:6]
                                               aar <- try(gls(rate ~ age, data = rbi.rate, weights = varFunc(~age),
                                                              control = list(singular.ok = TRUE)))
                                               if (class(aar) == "try-error")
                                                 aar <- lm(rate ~ age, data = rbi.rate)
                                               rbi.slopeRS <- coef(summary(aar))[2, c(1, 4)]
                                               aaa <- try(gls(abs(rate) ~ age, data = rbi.rate,
                                                              weights = varFunc(~age), control = list(singular.ok = TRUE)))
                                               if (class(aaa) == "try-error")
                                                 aaa <- lm(abs(rate) ~ age, data = rbi.rate)
                                               rbi.slopeAS <- coef(summary(aaa))[2, c(1, 4)]
                                               REG <- aar
                                             }
                                             if (class(node) != "NULL") {
                                               rbi.sma <- rbi.rate
                                               PP.sma <- PP
                                               rbi.sma$group <- rep("NA", dim(rbi.sma)[1])
                                               rbi.slopeRS.sel <- list()
                                               rbi.slopeAS.sel <- list()
                                               trend.reg.SEL <- list()
                                               for (j in 1:length(node)) {
                                                 n <- node[j]
                                                 sele <- getDescendants(t, n)
                                                 sele[which(sele < (Ntip(t) + 1))] <- t$tip.label[sele[which(sele <
                                                                                                               (Ntip(t) + 1))]]
                                                 rbi.sma[match(sele, rownames(rbi.sma)), ]$group <- paste("g",
                                                                                                          n, sep = "")
                                                 rbi.sel <- rbi[match(sele, rownames(rbi)), ]
                                                 if (length(y) > Ntip(t)) {
                                                   rbi.rate <- rbi.sel[, c(5:(dim(rbi.sel)[2] -
                                                                                1), dim(rbi.sel)[2])]
                                                   rbi.slopeAn <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                                                             1))
                                                   rbi.slopeRn <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                                                             1))
                                                   for (i in 1:(dim(rbi.rate)[2] - 1)) {
                                                     bet <- rbi.rate[, i]
                                                     age <- rbi.rate[, dim(rbi.rate)[2]]
                                                     betage <- data.frame(bet, age)
                                                     aar <- try(gls(bet ~ age, weights = varFunc(~age),
                                                                    control = list(singular.ok = TRUE)))
                                                     if (class(aar) == "try-error")
                                                       aar <- lm(bet ~ age)
                                                     rbi.slopeRn[i, ] <- coef(summary(aar))[2,
                                                                                            c(1, 4)]
                                                     aaa <- try(gls(abs(bet) ~ age, weights = varFunc(~age),
                                                                    control = list(singular.ok = TRUE)))
                                                     if (class(aaa) == "try-error")
                                                       aaa <- lm(abs(bet) ~ age)
                                                     rbi.slopeAn[i, ] <- coef(summary(aaa))[2,
                                                                                            c(1, 4)]
                                                     REGbetas <- aar
                                                   }
                                                   colnames(rbi.slopeRn) <- colnames(rbi.slopeAn) <- c("slope",
                                                                                                       "p-value")
                                                   rownames(rbi.slopeRn) <- rownames(rbi.slopeAn) <- colnames(data)[5:(5 +
                                                                                                                         dim(y)[2])]
                                                 }else {
                                                   rbi.rate <- rbi[, 5:6]
                                                   aar <- try(gls(rate ~ age, data = rbi.rate,
                                                                  weights = varFunc(~age), control = list(singular.ok = TRUE)))
                                                   if (class(aar) == "try-error")
                                                     aar <- lm(rate ~ age, data = rbi.rate)
                                                   rbi.slopeRn <- coef(summary(aar))[2, c(1,
                                                                                          4)]
                                                   aaa <- try(gls(abs(rate) ~ age, data = rbi.rate,
                                                                  weights = varFunc(~age), control = list(singular.ok = TRUE)))
                                                   if (class(aaa) == "try-error")
                                                     aaa <- lm(abs(rate) ~ age, data = rbi.rate)
                                                   rbi.slopeAn <- coef(summary(aaa))[2, c(1,
                                                                                          4)]
                                                   REG <- aar
                                                 }
                                                 rbi.slopeRS.sel[[j]] <- rbi.slopeRn
                                                 rbi.slopeAS.sel[[j]] <- rbi.slopeAn
                                                 nodes <- aceRR[1:Nnode(t), ]
                                                 if (length(y) > Ntip(t)) {
                                                   colnames(nodes)[1:dim(y)[2]] <- paste("y",
                                                                                         seq(1, dim(y)[2]), sep = "")
                                                   P <- rbind(nodes, y)
                                                   PP <- data.frame(P[match(rbi.sel[, 1], rownames(P)),
                                                                      ], rbi.sel$age)
                                                   colnames(PP)[dim(PP)[2]] <- "age"
                                                   trend.regC <- apply(PP[1:(dim(PP)[2] - 1)],
                                                                       2, function(x) summary(lm(x ~ PP[, dim(PP)[2]])))
                                                   trend.M <- summary(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                                                                collapse = "+"), colnames(PP)[length(colnames(PP))],
                                                                                          sep = "~")), data = PP))
                                                   names(trend.regC) <- paste("y", seq(1:dim(y)[2]),
                                                                              sep = "")
                                                   trend.regC[[length(trend.regC) + 1]] <- trend.M
                                                   trend.regC <- lapply(trend.regC, coefficients)
                                                   names(trend.regC)[length(trend.regC)] <- "multiple"

                                                 }else {
                                                   P <- c(nodes, y)
                                                   PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
                                                                    rbi.sel$age)
                                                   colnames(PP) <- c("phenotype", "age")
                                                   trend.regC <- lm(PP)
                                                   trend.regC <- coefficients(trend.regC)
                                                 }
                                                 trend.reg.SEL[[j]] <- trend.regC
                                               }
                                               names(rbi.slopeR.sel) <- names(rbi.slopeA.sel) <- names(trend.reg.SEL) <- node
                                               rbi.sma$group[which(rbi.sma$group == "NA")] <- "others"
                                               PP.sma <- cbind(PP.sma, group = rbi.sma[match(rownames(PP.sma),
                                                                                             rownames(rbi.sma)), ]$group)
                                               if (length(which(rbi.sma$group == "others")) < 3)
                                                 rbi.sma <- rbi.sma[-which(rbi.sma$group == "others"),
                                                                    ]
                                               if (length(which(PP.sma$group == "others")) < 3)
                                                 PP.sma <- PP.sma[-which(PP.sma$group == "others"),
                                                                  ]
                                               res[[i]] <- list(trend.regS, rbi.slopeRS, rbi.slopeAS,
                                                                PPtot, rbi, rbi.slopeRS.sel,
                                                                rbi.slopeAS.sel, trend.reg.SEL)
                                             }else {
                                               res[[i]] <- list(trend.regS, rbi.slopeRS, rbi.slopeAS,
                                                                PPtot, rbi)
                                             }
                                           }
  stopCluster(cl)
  if (length(y) > Ntip(t)) {
    p.rbi.slopeR <- array()
    p.rbi.slopeA <- array()
    for (i in 1:(dim(y)[2] + 1)) {
      rbi.slopeRS <- do.call(rbind, lapply(lapply(res,
                                                  "[[", 2), function(x) x[i, 1]))
      rbi.slopeRAS <- do.call(rbind, lapply(lapply(res,
                                                   "[[", 3), function(x) x[i, 1]))

      p.rbi.slopeR[i] <- rank(c(rbi.slopeR[i, 1],
                                rbi.slopeRS[1:(nsim - 1), ]))[1]/nsim

      p.rbi.slopeA[i] <- rank(c(rbi.slopeA[i, 1],
                                rbi.slopeRAS[1:(nsim - 1), ]))[1]/nsim
    }
    p.rbi.slopeA <- unname(p.rbi.slopeA)
    p.rbi.slopeR <- unname(p.rbi.slopeR)
    names(p.rbi.slopeR) <- names(p.rbi.slopeA) <- rownames(rbi.slopeA)
    p.rbi.slopeA <- data.frame(slope = rbi.slopeA[, 1],
                               p.real = rbi.slopeA[, 2], p.random = p.rbi.slopeA)
    p.rbi.slopeR <- data.frame(slope = rbi.slopeR[, 1],
                               p.real = rbi.slopeR[, 2],  p.random = p.rbi.slopeR)
  }else {
    rbi.slopesRS <- do.call(rbind, lapply(res, "[[", 2))[,
                                                         1]
    rbi.slopesRAS <- do.call(rbind, lapply(res, "[[", 3))[,
                                                          1]

    p.rbi.slopeR <- rank(c(rbi.slopeR[1], rbi.slopesRS[1:(nsim -
                                                            1)]))[1]/nsim

    p.rbi.slopeA <- rank(c(rbi.slopeA[1], rbi.slopesRAS[1:(nsim -
                                                             1)]))[1]/nsim
    p.rbi.slopeA <- unname(p.rbi.slopeA)
    p.rbi.slopeR <- unname(p.rbi.slopeR)
    rbi.slopeA <- unname(rbi.slopeA)
    rbi.slopeR <- unname(rbi.slopeR)
    p.rbi.slopeA <- c(slope = rbi.slopeA[1], p.real = rbi.slopeA[2],
                      p.random = p.rbi.slopeA)
    p.rbi.slopeR <- c(slope = rbi.slopeR[1], p.real = rbi.slopeR[2],
                      p.random = p.rbi.slopeR)
  }
  if (class(node) != "NULL") {
    p.slopeA.sel <- list()
    p.slopeR.sel <- list()
    for (k in 1:length(node)) {
      pA <- array()
      pR <- array()
      if (length(y) > Ntip(t)) {
        for (i in 1:(dim(y)[2] + 1)) {
          rbi.slopeRS.sel <- do.call(rbind, lapply(lapply(lapply(res,
                                                                 "[[", 6), "[[", k), function(x) x[i, 1]))
          rbi.slopeRAS.sel <- do.call(rbind, lapply(lapply(lapply(res,
                                                                  "[[", 7), "[[", k), function(x) x[i, 1]))

          pR[i] <- rank(c(rbi.slopeR.sel[[k]][i, 1],
                          rbi.slopeRS.sel[1:(nsim - 1), ]))[1]/nsim

          pA[i] <- rank(c(rbi.slopeA.sel[[k]][i, 1],
                          rbi.slopeRAS.sel[1:(nsim - 1), ]))[1]/nsim
        }
      }else {
        rbi.slopeRS.sel <- do.call(rbind, lapply(lapply(res,
                                                        "[[", 6), "[[", k))[, 1]
        rbi.slopeRAS.sel <- do.call(rbind, lapply(lapply(res,
                                                         "[[", 7), "[[", k))[, 1]

        pR <- rank(c(rbi.slopeR.sel[[k]][1], rbi.slopeRS.sel[1:(nsim -
                                                                  1)]))[1]/nsim

        pA <- rank(c(rbi.slopeA.sel[[k]][1], rbi.slopeRAS.sel[1:(nsim -
                                                                   1)]))[1]/nsim
        pA <- unname(pA)
        pR <- unname(pR)
      }
      p.slopeA.sel[[k]] <- pA
      p.slopeR.sel[[k]] <- pR
    }
    names(p.slopeA.sel) <- names(p.slopeR.sel) <- node
    p.rbi.slopeA.sel <- list()
    p.rbi.slopeR.sel <- list()
    if (length(y) > Ntip(t)) {
      for (p in 1:length(rbi.slopeA.sel)) {
        p.rbi.slopeA.sel[[p]] <- cbind(slope = rbi.slopeA.sel[[p]][,
                                                                   1], p.real = rbi.slopeA.sel[[p]][, 2],
                                       p.random = p.slopeA.sel[[p]])
        p.rbi.slopeR.sel[[p]] <- cbind(slope = rbi.slopeR.sel[[p]][,
                                                                   1], p.real = rbi.slopeR.sel[[p]][, 2],
                                       p.random = p.slopeR.sel[[p]])
      }
    }else {
      for (p in 1:length(rbi.slopeA.sel)) {
        p.rbi.slopeA.sel[[p]] <- c(slope = unname(rbi.slopeA.sel[[p]])[1],
                                   p.real = unname(rbi.slopeA.sel[[p]])[2],
                                   p.random = p.slopeA.sel[[p]])
        p.rbi.slopeR.sel[[p]] <- c(slope = unname(rbi.slopeR.sel[[p]])[1],
                                   p.real = unname(rbi.slopeR.sel[[p]])[2],
                                   p.random = p.slopeR.sel[[p]])
      }
    }
    names(p.rbi.slopeR.sel) <- names(p.rbi.slopeA.sel) <- node
    if (length(y) > Ntip(t)) {
      p.smaR <- list()
      p.smaA <- list()
      p.smaPP <- list()
      for (p in 1:length(sma.resA)) {
        p.smaR[[p]] <- cbind(sma.resR[[p]][, 1:4], sma.resR[[p]][,
                                                                 6:7])
        p.smaA[[p]] <- cbind(sma.resA[[p]][, 1:4], sma.resA[[p]][,
                                                                 6:7])
        if (p <= length(sma.resPP)) {
          p.smaPP[[p]] <- cbind(sma.resPP[[p]][, 1:4],
                                sma.resPP[[p]][, 6:7])
        }
      }
    }else {
      p.smaR <- cbind(sma.resR[, 1:4], sma.resR[, 6:7])
      p.smaA <- cbind(sma.resA[, 1:4], sma.resA[, 6:7])
      p.smaPP <- cbind(sma.resPP[, 1:4], sma.resPP[, 6:7])
    }
    if (class(p.smaR) == "list") {
      names(p.smaR) <- names(p.smaA) <- names(sma.resA)
      names(p.smaPP) <- names(sma.resPP)
    }
  }
  p.trend <- array()
  PP <- PPtot
  if (length(y) > Ntip(t)) {
    trend.slopes <- matrix(ncol = dim(y)[2] + 1, nrow = nsim)
    for (i in 1:(dim(y)[2] + 1)) {
      trend.slopes[, i] <- unlist(lapply(lapply(lapply(res,
                                                       "[[", 1), "[[", i), "[[", 2))

      p.tr <- rank(c(trend.reg[[i]][2, 1], trend.slopes[,
                                                        i][1:(nsim - 1)]))[1]/nsim
      p.trend[i] <- p.tr
    }
    trend.real <- do.call(rbind, lapply(trend.reg, function(x) x[2,
                                                                 c(1, 4)]))
    p.trend <- cbind(trend.real, p.trend)
    colnames(p.trend) <- c("slope", "p.real",
                           "p.random")
    CIphenotype <- list()
    for (i in 1:(dim(y)[2] + 1)) {
      PPci <- apply(do.call(cbind, lapply(lapply(res,
                                                 "[[", 4), function(x) x[, i])), 1, function(u) quantile(u,
                                                                                                         c(0.025, 0.975)))
      colnames(PPci) <- lapply(lapply(res, "[[", 4), rownames)[[1]]
      CIphenotype[[i]] <- t(PPci)
    }
    names(CIphenotype) <- paste("y", seq(1, dim(y)[2]),
                                sep = "")
    if (dim(y)[2] <= 3) {
      pdf(file = paste(foldername, "Phenotypic Trend Test.pdf",
                       sep = "/"))
      par(mar = c(3.5, 3.5, 1, 2))
      par(mfrow = c(dim(y)[2] + 1, 2))
      for (i in 1:(dim(y)[2] + 1)) {
        if (i == length(colnames(PP))) {
          obj <- hist(trend.slopes[, i], xlab = "",
                      ylab = "", main = names(trend.reg[i]), mgp = c(2,
                                                                     0.5, 0))
          title(xlab = "simulated slopes", ylab = "frequency",
                line = 1.5)
          text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
                 0.9, labels = paste("p=", p.trend[i, 3]),
               cex = 1)
          abline(v = trend.reg[[i]][2, 1], lwd = 3,
                 col = "red")
          plot(PP[, length(colnames(PP))], resid(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                                           collapse = "+"), colnames(PP)[length(colnames(PP))],
                                                                     sep = "~")), data = PP)), xlab = "", ylab = "",
               mgp = c(2, 0.5, 0))
          title(xlab = "age", ylab = "residuals of \nmultiple PPT",
                line = 1.5)
          abline(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                           collapse = "+"), colnames(PP)[length(colnames(PP))],
                                     sep = "~")), data = PP), lwd = 4, col = "blue")
        } else {
          PPci <- apply(do.call(cbind, lapply(lapply(res,
                                                     "[[", 4), function(x) x[, i])), 1, function(u) quantile(u,
                                                                                                             c(0.025, 0.975)))
          colnames(PPci) <- lapply(lapply(res, "[[",
                                          4), rownames)[[1]]
          CIphenotype[[i]] <- t(PPci)
          PPci <- cbind(PP[, c(i, dim(PP)[2])], t(PPci[,
                                                       match(rownames(PP), colnames(PPci))]))
          PPci <- PPci[order(PPci[, 2]), ]
          obj <- hist(trend.slopes[, i], xlab = "",
                      ylab = "", main = names(trend.reg[i]), mgp = c(2,
                                                                     0.5, 0))
          title(xlab = "simulated slopes", ylab = "frequency",
                line = 1.5)
          text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
                 0.9, labels = paste("p=", p.trend[i, 3]),
               cex = 1)
          abline(v = trend.reg[[i]][2, 1], lwd = 3,
                 col = "red")
          plot(PP[, c(dim(PP)[2], i)], xlab = "", ylab = "",
               mgp = c(2, 0.5, 0))
          polygon(c(PPci[, 2], rev(PPci[, 2])), c(PPci[,
                                                       3], rev(PPci[, 4])), col = rgb(0.5, 0.5,
                                                                                      0.5, 0.4), border = NA)
          title(xlab = "age", ylab = paste(colnames(PP)[i]),
                line = 1.5)
          points((diag(vcv(t))+PP[1, 2]), y[, i], pch = 21, col = "black",
                 bg = "red")
          if (class(node) != "NULL") {
            for (j in 1:length(node)) {
              cols <- brewer.pal(length(node), "Set2")
              points(trend.reg.age.sel[[j]], trend.reg.y.sel[[j]][,
                                                                  i], lwd = 4, col = cols[j], type = "l")
            }
            abline(lm(PP[, i] ~ age, data = PP), lwd = 3,
                   col = "blue", lty = 2)
            if (i == 1)
              legend(min(PP$age), max(PP[, i]), legend = node,
                     fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                      0, 0, 0), border = NA, x.intersp = 0.25)
          }else {
            abline(lm(PP[, i] ~ age, data = PP), lwd = 4,
                   col = "blue")
          }
        }
      }
    }else {
      pdf(file = paste(foldername, "Phenotypic Trend Test.pdf",
                       sep = "/"))
      par(mar = c(3.5, 3.5, 1, 2))
      par(mfrow = c(2, 1))
      i <- length(colnames(PP))
      obj <- hist(trend.slopes[, i], xlab = "", ylab = "",
                  main = "Phenotypic Trend Test", mgp = c(2, 0.5,
                                                          0))
      title(xlab = "simulated slopes", ylab = "frequency",
            line = 1.5)
      text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
             0.9, labels = paste("p=", p.trend[i, 3]), cex = 1)
      abline(v = trend.reg[[i]][2, 1], lwd = 3, col = "red")
      plot(PP[, length(colnames(PP))], resid(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                                                       collapse = "+"), colnames(PP)[length(colnames(PP))],
                                                                 sep = "~")), data = PP)), xlab = "", ylab = "",
           mgp = c(2, 0.5, 0))
      title(xlab = "age", ylab = "residuals of \nmultiple PPT",
            line = 1.5)
      abline(lm(as.formula(paste(paste(colnames(PP)[-length(colnames(PP))],
                                       collapse = "+"), colnames(PP)[length(colnames(PP))],
                                 sep = "~")), data = PP), lwd = 4, col = "blue")
    }
    dev.off()
  }else {
    trend.slopes <- do.call(rbind, lapply(res, "[[", 1))[,
                                                         2]

    p.trend <- rank(c(trend.reg[2, 1], trend.slopes[1:(nsim -
                                                         1)]))[1]/nsim
    p.trend <- c(trend.reg[2, c(1, 4)],  p.trend)
    names(p.trend) <- c("slope", "p.real", "p.random")
    pdf(file = paste(foldername, "Phenotypic Trend Test.pdf",
                     sep = "/"))
    par(mar = c(3.5, 3.5, 1, 2))
    par(mfrow = c(2, 1))
    obj <- hist(trend.slopes, xlab = "simulated slopes",
                ylab = "frequency", main = "Phenotypic Trend Test",
                mgp = c(2, 0.5, 0))
    text(quantile(trend.slopes, 0.01), max(obj$counts) *
           0.8, labels = paste("p=", p.trend[3]), cex = 1)
    abline(v = trend.reg[2, 1], lwd = 3, col = "red")
    plot(PP[, c(2, 1)], mgp = c(2, 0.5, 0))
    PPci <- apply(do.call(cbind, lapply(lapply(res, "[[",
                                               4), function(x) x[, 1])), 1, function(u) quantile(u,
                                                                                                 c(0.025, 0.975)))
    colnames(PPci) <- lapply(lapply(res, "[[", 4), rownames)[[1]]
    CIphenotype <- t(PPci)
    PPci <- cbind(PP, t(PPci[, match(rownames(PP), colnames(PPci))]))
    PPci <- PPci[order(PPci[, 2]), ]
    polygon(c(PPci[, 2], rev(PPci[, 2])), c(PPci[, 3], rev(PPci[,
                                                                4])), col = rgb(0.5, 0.5, 0.5, 0.4), border = NA)
    points((diag(vcv(t))+PP[1, 2]), y, pch = 21, col = "black", bg = "red")
    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        cols <- brewer.pal(length(node), "Set2")
        points(trend.reg.age.sel[[j]], trend.reg.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(lm(PP), lwd = 3, col = "blue", lty = 2)
      legend(min(PP[, 2]), max(PP[, 1]), legend = node,
             fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                              0, 0, 0), border = NA, x.intersp = 0.25)
    }else {
      abline(lm(PP), lwd = 4, col = "blue")
    }
    dev.off()
  }
  if (class(node) != "NULL") {
    p.trend.sel <- list()
    for (u in 1:length(node)) {
      if (length(y) > Ntip(t)) {
        p.sele <- list()
        for (i in 1:(dim(y)[2] + 1)) {
          slopeR <- unlist(lapply(lapply(lapply(lapply(res,
                                                       "[[", 8), "[[", u), "[[", i), function(x) x[2,
                                                                                                   1]))

          p.sel <- rank(c(trend.reg.sel[[u]][[i]][2,
                                                  1], slopeR[1:(nsim - 1)]))[1]/nsim
          p.sele[[i]] <- c(slope = trend.reg.sel[[u]][[i]][2,
                                                           1], p.real = trend.reg.sel[[u]][[i]][2,
                                                                                                4], p.random = p.sel)
        }
        names(p.sele) <- names(trend.reg.sel[[u]])
        p.selt <- do.call(rbind, p.sele)
        p.trend.sel[[u]] <- p.selt
      } else {
        slopeR <- unlist(lapply(lapply(lapply(res, "[[",
                                              8), "[[", u), function(x) x[2]))

        p.sel <- rank(c(trend.reg.sel[[u]][2, 1],
                        slopeR[1:(nsim - 1)]))[1]/nsim
        p.trend.sel[[u]] <- c(slope = trend.reg.sel[[u]][2,
                                                         1], p.real = trend.reg.sel[[u]][2, 4],
                              p.random = p.sel)
      }
    }
    names(p.trend.sel) <- names(trend.reg.sel)
  }
  rbi <- rbi[, -4]
  if (length(y) > Ntip(t)) {
    A <- rbi[, 4:dim(rbi)[2]]
    CIrelative <- CIabsolute <- list()
    if (dim(y)[2] <= 3) {
      pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                       sep = "/"))
      par(mfrow = c(dim(y)[2] + 1, 2))
      par(mar = c(3.5, 3.5, 1, 1))
      for (i in 1:(dim(y)[2] + 1)) {
        if (i == dim(y)[2] + 1)
          ynam <- "rate" else ynam <- paste("betas", i, sep = "")
          bet <- A[, i]
          age <- A[, dim(A)[2]]
          RBTRci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                              "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                                function(k) k[, i])), 1, function(u) quantile(u,
                                                                                              c(0.025, 0.975)))
          RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                              "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                                function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                                   c(0.025, 0.975)))
          colnames(RBTAci) <- colnames(RBTRci) <- lapply(lapply(res,
                                                                "[[", 5), rownames)[[1]]
          CIabsolute[[i]] <- t(RBTAci)
          CIrelative[[i]] <- t(RBTRci)
          RBTRci <- cbind(A[, c(i, dim(A)[2])], t(RBTRci[,
                                                         match(rownames(A), colnames(RBTRci))]))
          RBTAci <- cbind(A[, c(i, dim(A)[2])], t(RBTAci[,
                                                         match(rownames(A), colnames(RBTAci))]))
          RBTRci <- RBTRci[order(RBTRci[, 2]), ]
          RBTAci <- RBTAci[order(RBTAci[, 2]), ]
          plot(abs(bet) ~ age, main = "absolute rate",
               ylab = ynam, mgp = c(2, 0.5, 0))
          polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                             3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5,
                                                                                              0.5, 0.4), border = NA)
          if (class(node) != "NULL") {
            for (j in 1:length(node)) {
              cols <- brewer.pal(length(node), "Set2")
              points(REG.betas.age.sel[[j]][[i]], REGabs.betas.y.sel[[j]][[i]],
                     lwd = 4, col = cols[j], type = "l")
            }
            abline(REGabs.betas[[i]], lwd = 3, col = "red",
                   lty = 2)
            if (i == 1)
              legend(min(age), max(abs(bet)), legend = node,
                     fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                      0, 0, 0), border = NA, x.intersp = 0.25)
          } else {
            abline(REGabs.betas[[i]], lwd = 4, col = "red")
          }
          plot(bet ~ age, main = "rate", ylab = ynam,
               mgp = c(2, 0.5, 0))
          polygon(c(RBTRci[, 2], rev(RBTRci[, 2])), c(RBTRci[,
                                                             3], rev(RBTRci[, 4])), col = rgb(0.5, 0.5,
                                                                                              0.5, 0.4), border = NA)
          if (class(node) != "NULL") {
            for (j in 1:length(node)) {
              points(REG.betas.age.sel[[j]][[i]], REG.betas.y.sel[[j]][[i]],
                     lwd = 4, col = cols[j], type = "l")
            }
            abline(REG.betas[[i]], lwd = 3, col = "blue",
                   lty = 2)
          } else {
            abline(REG.betas[[i]], lwd = 4, col = "blue")
          }
      }
    } else {
      for (i in 1:(dim(y)[2] + 1)) {
        RBTRci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                            "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                              function(k) k[, i])), 1, function(u) quantile(u,
                                                                                            c(0.025, 0.975)))
        RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                            "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                              function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                                 c(0.025, 0.975)))
        colnames(RBTAci) <- colnames(RBTRci) <- lapply(lapply(res,
                                                              "[[", 5), rownames)[[1]]
        CIabsolute[[i]] <- t(RBTAci)
        CIrelative[[i]] <- t(RBTRci)
      }
      pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                       sep = "/"))
      par(mfrow = c(2, 1))
      bet <- A[, (dim(y) + 1)[2]]
      age <- A[, dim(A)[2]]
      RBTRci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                          "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                            function(k) k[, i])), 1, function(u) quantile(u,
                                                                                          c(0.025, 0.975)))
      RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                          "[[", 5), function(x) x[, c(5:dim(x)[2])]),
                                            function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                               c(0.025, 0.975)))
      colnames(RBTAci) <- colnames(RBTRci) <- lapply(lapply(res,
                                                            "[[", 5), rownames)[[1]]
      RBTRci <- cbind(A[, c(i, dim(A)[2])], t(RBTRci[,
                                                     match(rownames(A), colnames(RBTRci))]))
      RBTAci <- cbind(A[, c(i, dim(A)[2])], t(RBTAci[,
                                                     match(rownames(A), colnames(RBTAci))]))
      RBTRci <- RBTRci[order(RBTRci[, 2]), ]
      RBTAci <- RBTAci[order(RBTAci[, 2]), ]
      plot(abs(bet) ~ age, ylab = "absolute rate", mgp = c(2,
                                                           0.5, 0))
      polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                         3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                          0.4), border = NA)
      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          cols <- brewer.pal(length(node), "Set2")
          points(REG.betas.age.sel[[j]][[(dim(y) + 1)[2]]],
                 REGabs.betas.y.sel[[j]][[(dim(y) + 1)[2]]],
                 lwd = 4, col = cols[j], type = "l")
        }
        abline(REGabs.betas[[(dim(y) + 1)[2]]], lwd = 3,
               col = "red", lty = 2)
        legend(min(age), max(abs(bet)), legend = node,
               fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                0, 0, 0), border = NA, x.intersp = 0.25)
      }else {
        abline(REGabs.betas[[(dim(y) + 1)[2]]], lwd = 4,
               col = "red")
      }
      plot(bet ~ age, ylab = "rate", mgp = c(2, 0.5, 0))
      polygon(c(RBTRci[, 2], rev(RBTRci[, 2])), c(RBTRci[,
                                                         3], rev(RBTRci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                          0.4), border = NA)
      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          points(REG.betas.age.sel[[j]][[(dim(y) + 1)[2]]],
                 REG.betas.y.sel[[j]][[(dim(y) + 1)[2]]],
                 lwd = 4, col = cols[j], type = "l")
        }
        abline(REG.betas[[(dim(y) + 1)[2]]], lwd = 3,
               col = "blue", lty = 2)
      }else {
        abline(REG.betas[[(dim(y) + 1)[2]]], lwd = 4,
               col = "blue")
      }
    }
    names(CIabsolute) <- names(CIrelative) <- c(paste("betas",
                                                      seq(1, dim(y)[2]), sep = ""), "rate")
    dev.off()
  }else {
    A <- rbi
    AA <- A[, 4:5]
    RBTRci <- apply(do.call(cbind, lapply(lapply(res, "[[",
                                                 5), function(x) x[, 5])), 1, function(u) quantile(u,
                                                                                                   c(0.025, 0.975)))
    RBTAci <- apply(do.call(cbind, lapply(lapply(res, "[[",
                                                 5), function(x) abs(x[, 5]))), 1, function(u) quantile(u,
                                                                                                        c(0.025, 0.975)))
    colnames(RBTAci) <- colnames(RBTRci) <- lapply(lapply(res,
                                                          "[[", 5), rownames)[[1]]
    CIrelative <- t(RBTRci)
    CIabsolute <- t(RBTAci)
    RBTRci <- cbind(AA, t(RBTRci[, match(rownames(AA), colnames(RBTRci))]))
    RBTAci <- cbind(AA, t(RBTAci[, match(rownames(AA), colnames(RBTAci))]))
    RBTRci <- RBTRci[order(RBTRci[, 2]), ]
    RBTAci <- RBTAci[order(RBTAci[, 2]), ]
    pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                     sep = "/"))
    par(mfrow = c(2, 1))
    par(mar = c(3, 3.5, 1, 1.5))
    plot(abs(rate) ~ age, data = AA, ylab = "absolute rate",
         mgp = c(2, 0.5, 0))
    polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                       3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                        0.4), border = NA)
    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        cols <- brewer.pal(length(node), "Set2")
        points(REG.betas.age.sel[[j]], REGabs.betas.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(REGabs, col = "red3", lwd = 3, lty = 2)
      legend(min(AA$age), max(abs(AA$rate)), legend = node,
             fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                              0, 0, 0), border = NA, x.intersp = 0.25)
    }else {
      abline(REGabs, col = "red3", lwd = 4)
    }
    plot(rate ~ age, data = AA, ylab = "rate", mgp = c(2,
                                                       0.5, 0))
    polygon(c(RBTRci[, 2], rev(RBTRci[, 2])), c(RBTRci[,
                                                       3], rev(RBTRci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                        0.4), border = NA)
    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        points(REG.betas.age.sel[[j]], REG.betas.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(REG, lwd = 3, col = "blue", lty = 2)
    } else {
      abline(REG, lwd = 4, col = "blue")
    }
    dev.off()
  }
  colnames(rbi)[1] <- "branch"

  if (length(y) > Ntip(t)) {
    p.trend <- cbind(p.trend, p.value = rep(NA, dim(p.trend)[1]))
    p.rbi.slopeA <- cbind(p.rbi.slopeA, p.value = rep(NA,
                                                      dim(p.rbi.slopeA)[1]))
    p.rbi.slopeR <- cbind(p.rbi.slopeR, p.value = rep(NA,
                                                      dim(p.rbi.slopeR)[1]))
    for (i in 1:(dim(y)[2] + 1)) {
      if (p.trend[i, 3] > 0.5)
        p.trend[i, 3] <- 1 - p.trend[i, 3]
      if (p.rbi.slopeA[i, 3] > 0.5)
        p.rbi.slopeA[i, 3] <- 1 - p.rbi.slopeA[i, 3]
      if (p.rbi.slopeR[i, 3] > 0.5)
        p.rbi.slopeR[i, 3] <- 1 - p.rbi.slopeR[i, 3]
      p.trend[i, 4] <- p.trend[i, 3]
      p.rbi.slopeA[i, 4] <- p.rbi.slopeA[i, 3]
      p.rbi.slopeR[i, 4] <- p.rbi.slopeR[i, 3]
      if (p.trend[i, 4] <= 0.05) {
        if (i == dim(y)[2] + 1) {
          if (p.trend[i, 1] > 0)
            print("There is a trend for increase in phenotype y.multi")else print("There is a trend for increase in phenotype y.multi",
                                                                                  i, sep = "")
        } else {
          if (p.trend[i, 1] > 0)
            print(paste("There is a trend for increase in phenotype",
                        colnames(y)[i])) else print(paste("There is a trend for decrease in phenotype",
                                                          colnames(y)[i]))
        }
      }
      if (p.rbi.slopeA[i, 4] <= 0.05) {
        if (i == dim(y)[2] + 1) {
          if (p.rbi.slopeA[i, 1] > 0)
            print("There is a trend for increase in absolute evolutionary rates for variable y.multi") else print("There is a trend for decrease in absolute evolutionary rates for variable y")
        }else {
          if (p.rbi.slopeA[i, 1] > 0)
            print(paste("There is a trend for increase in absolute evolutionary rates for variable",
                        colnames(y)[i])) else print(paste("There is a trend for decrease in absolute evolutionary rates for variable",
                                                          colnames(y)[i]))
        }
      }
      if (p.rbi.slopeR[i, 4] <= 0.05) {
        if (i == dim(y)[2] + 1) {
          if (p.rbi.slopeR[i, 1] > 0)
            print("There is a trend for increase in relative evolutionary rates for variable y.multi")else print("There is a trend for decrease in relative evolutionary rates for variable y")
        }else {
          if (p.rbi.slopeR[i, 1] > 0)
            print(paste("There is a trend for increase in relative evolutionary rates for variable",
                        colnames(y)[i]))else print(paste("There is a trend for decrease in relative evolutionary rates for variable",
                                                         colnames(y)[i]))
        }
      }
    }
  }else {
    if (p.trend[3] > 0.5)
      p.trend[3] <- 1 - p.trend[3]
    if (p.rbi.slopeA[3] > 0.5)
      p.rbi.slopeA[3] <- 1 - p.rbi.slopeA[3]
    if (p.rbi.slopeR[3] > 0.5)
      p.rbi.slopeR[3] <- 1 - p.rbi.slopeR[3]
    p.trend[4] <- p.trend[3]
    p.rbi.slopeA[4] <- p.rbi.slopeA[3]
    p.rbi.slopeR[4] <- p.rbi.slopeR[3]
    names(p.rbi.slopeA)[4] <- names(p.rbi.slopeR)[4] <- names(p.trend)[4] <- "p.value"
    if (p.trend[4] <= 0.05) {
      if (p.trend[1] > 0)
        print("There is a trend for increase in phenotype")else print("There is a trend for decrease in phenotype")
    }
    if (p.rbi.slopeA[4] <= 0.05) {
      if (p.rbi.slopeA[1] > 0)
        print("There is a trend for increase in absolute evolutionary rates")else print("There is a trend for decrease in absolute evolutionary rates")
    }
    if (p.rbi.slopeR[4] <= 0.05) {
      if (p.rbi.slopeR[1] > 0)
        print("There is a trend for increase in relative evolutionary rates")else print("There is a trend for decrease in relative evolutionary rates")
    }
  }
  if (class(node) != "NULL") {
    if (length(y) > Ntip(t)) {
      for (j in 1:length(p.trend.sel)) {
        p.trend.sel[[j]] <- cbind(p.trend.sel[[j]],
                                  p.value = rep(NA, dim(p.trend.sel[[j]])[1]))
        p.rbi.slopeA.sel[[j]] <- cbind(p.rbi.slopeA.sel[[j]],
                                       p.value = rep(NA, dim(p.rbi.slopeA.sel[[j]])[1]))
        p.rbi.slopeR.sel[[j]] <- cbind(p.rbi.slopeR.sel[[j]],
                                       p.value = rep(NA, dim(p.rbi.slopeR.sel[[j]])[1]))
        for (i in 1:(dim(y)[2] + 1)) {
          if (p.trend.sel[[j]][i, 3] > 0.5)
            p.trend.sel[[j]][i, 3] <- 1 - p.trend.sel[[j]][i,
                                                           3]
          if (p.rbi.slopeA.sel[[j]][i, 3] > 0.5)
            p.rbi.slopeA.sel[[j]][i, 3] <- 1 - p.rbi.slopeA.sel[[j]][i,
                                                                     3]
          if (p.rbi.slopeR.sel[[j]][i, 3] > 0.5)
            p.rbi.slopeR.sel[[j]][i, 3] <- 1 - p.rbi.slopeR.sel[[j]][i,
                                                                     3]
          p.trend.sel[[j]][i, 4] <- p.trend.sel[[j]][i,3]
          p.rbi.slopeA.sel[[j]][i, 4] <- p.rbi.slopeA.sel[[j]][i,3]
          p.rbi.slopeR.sel[[j]][i, 4] <- p.rbi.slopeR.sel[[j]][i,3]
        }
      }
      for (j in 1:length(p.smaPP)) {
        colnames(p.smaPP[[j]])[3] <- "p.value.sma"
        colnames(p.smaA[[j]])[3] <- "p.value.sma"
        colnames(p.smaR[[j]])[3] <- "p.value.sma"
        for (i in 1:dim(p.smaPP[[j]])[1]) {
          if (as.character(p.smaPP[[j]][i, 2]) == "others") {
            if (p.smaPP[[j]][i, 3] <= 0.05)
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is different from regression through others for variable",
                          colnames(y)[j]))
            if (p.smaA[[j]][i, 3] <= 0.05)
              print(paste("Absolute evolutionary rates regression through node",
                          lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                   1]), "g"), "[[", 2), "is different from regression through others for variable",
                          colnames(y)[j]))
            if (p.smaR[[j]][i, 3] <= 0.05)
              print(paste("Relative evolutionary rates regression through node",
                          lapply(strsplit(as.character(p.smaR[[j]][i,
                                                                   1]), "g"), "[[", 2), "is different from regression through others for variable",
                          colnames(y)[j]))
          }else {
            if (p.smaPP[[j]][i, 3] <= 0.05)
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is different from regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    2]), "g"), "[[", 2), "for variable",
                          colnames(y)[j]))
            if (p.smaA[[j]][i, 3] <= 0.05)
              print(paste("Absolute evolutionary rates regression through node",
                          lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                   1]), "g"), "[[", 2), "is different from regression through node",
                          lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                   2]), "g"), "[[", 2), "for variable",
                          colnames(y)[j]))
            if (p.smaR[[j]][i, 3] <= 0.05)
              print(paste("Relative evolutionary rates regression through node",
                          lapply(strsplit(as.character(p.smaR[[j]][i,
                                                                   1]), "g"), "[[", 2), "is different from regression through node",
                          lapply(strsplit(as.character(p.smaR[[j]][i,
                                                                   2]), "g"), "[[", 2), "for variable",
                          colnames(y)[j]))
          }
        }
        if (as.character(p.smaPP[[j]][i, 2]) == "others") {
          if (p.smaA[[length(p.smaA)]][i, 3] <= 0.05)
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[[length(p.smaA)]][i,
                                                                              1]), "g"), "[[", 2), "is different from regression through others for variable y.multi"))
          if (p.smaR[[length(p.smaA)]][i, 3] <= 0.05)
            print(paste("Relative evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaR[[length(p.smaA)]][i,
                                                                              1]), "g"), "[[", 2), "is different from regression through others for variable y.multi"))
        }else {
          if (p.smaA[[length(p.smaA)]][i, 3] <= 0.05)
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[[length(p.smaA)]][i,
                                                                              1]), "g"), "[[", 2), "is different from regression through node",
                        lapply(strsplit(as.character(p.smaA[[length(p.smaA)]][i,
                                                                              2]), "g"), "[[", 2), "for variable y.multi"))
          if (p.smaR[[length(p.smaA)]][i, 3] <= 0.05)
            print(paste("Relative evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaR[[length(p.smaA)]][i,
                                                                              1]), "g"), "[[", 2), "is different from regression through node",
                        lapply(strsplit(as.character(p.smaR[[length(p.smaA)]][i,
                                                                              2]), "g"), "[[", 2), "for variable y.multi"))
        }
      }
      colnames(p.smaA[[length(p.smaA)]])[3] <- "p.value.sma"
      colnames(p.smaR[[length(p.smaR)]])[3] <- "p.value.sma"
    }else {
      for (j in 1:length(p.trend.sel)) {
        if (p.trend.sel[[j]][3] > 0.5)
          p.trend.sel[[j]][3] <- 1 - p.trend.sel[[j]][3]
        if (p.rbi.slopeA.sel[[j]][3] > 0.5)
          p.rbi.slopeA.sel[[j]][3] <- 1 - p.rbi.slopeA.sel[[j]][3]
        if (p.rbi.slopeR.sel[[j]][3] > 0.5)
          p.rbi.slopeR.sel[[j]][3] <- 1 - p.rbi.slopeR.sel[[j]][3]
        p.trend.sel[[j]][4] <- p.trend.sel[[j]][3]
        p.rbi.slopeA.sel[[j]][4] <- p.rbi.slopeA.sel[[j]][3]
        p.rbi.slopeR.sel[[j]][4] <- p.rbi.slopeR.sel[[j]][3]
        names(p.rbi.slopeA.sel[[j]])[4] <- names(p.rbi.slopeR.sel[[j]])[4] <- names(p.trend.sel[[j]])[4] <- "p.value"
      }
      for (i in 1:dim(p.smaPP)[1]) {
        if (as.character(p.smaPP[i, 2]) == "others") {
          if (p.smaPP[i, 3] <= 0.05)
            print(paste("Phenotypic regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             1]), "g"), "[[", 2), "is different from regression through others"))
          if (p.smaA[i, 3] <= 0.05)
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            1]), "g"), "[[", 2), "is different from regression through others"))
          if (p.smaR[i, 3] <= 0.05)
            print(paste("Relative evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaR[i,
                                                            1]), "g"), "[[", 2), "is different from regression through others"))
        }else {
          if (p.smaPP[i, 3] <= 0.05)
            print(paste("Phenotypic regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             1]), "g"), "[[", 2), "is different from regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             2]), "g"), "[[", 2)))
          if (p.smaA[i, 3] <= 0.05)
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            1]), "g"), "[[", 2), "is different from regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            2]), "g"), "[[", 2)))
          if (p.smaR[i, 3] <= 0.05)
            print(paste("Relative evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaR[i,
                                                            1]), "g"), "[[", 2), "is different from regression through node",
                        lapply(strsplit(as.character(p.smaR[i,
                                                            2]), "g"), "[[", 2)))
        }
      }
    }
    SMA.res <- list(p.smaPP, p.smaA, p.smaR)
    CInts <- list(CIphenotype, CIabsolute, CIrelative)
    names(SMA.res) <- c("phenotype", "abs.rate", "rel.rate")
    names(CInts) <- c("phenotype", "abs.rate", "rel.rate")
    if (ConfInt == TRUE) {
      res <- list(rbiRES, PPtot, p.trend, p.rbi.slopeA,
                  p.rbi.slopeR, p.trend.sel, p.rbi.slopeA.sel,
                  p.rbi.slopeR.sel, SMA.res, CInts)
      names(res) <- c("rbt", "pbt", "p.trend", "rbt.rateA",
                      "rbt.rateR", "p.trend.nodes", "rbt.rateA.nodes",
                      "rbt.rateR.nodes", "SMA", "ConfInts")
    }else {
      res <- list(rbiRES, PPtot, p.trend, p.rbi.slopeA,
                  p.rbi.slopeR, p.trend.sel, p.rbi.slopeA.sel,
                  p.rbi.slopeR.sel, SMA.res)
      names(res) <- c("rbt", "pbt", "p.trend", "rbt.rateA",
                      "rbt.rateR", "p.trend.nodes", "rbt.rateA.nodes",
                      "rbt.rateR.nodes", "SMA")
    }
  }else {
    if (ConfInt == TRUE) {
      CInts <- list(CIphenotype, CIabsolute, CIrelative)
      names(CInts) <- c("phenotype", "abs.rate", "rel.rate")
      res <- list(rbi, PPtot, p.trend, p.rbi.slopeA, p.rbi.slopeR,
                  CInts)
      names(res) <- c("rbt", "pbt", "p.trend", "rbt.rateA",
                      "rbt.rateR", "ConfInts")
    }else {
      res <- list(rbi, PPtot, p.trend, p.rbi.slopeA, p.rbi.slopeR)
      names(res) <- c("rbt", "pbt", "p.trend", "rbt.rateA",
                      "rbt.rateR")
    }
  }
  return(res)
}
