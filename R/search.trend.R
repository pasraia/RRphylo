#' @title Searching for evolutionary trends in phenotypes and rates
#' @description This function searches for evolutionary trends in the phenotypic mean and the evolutionary rates for the entire tree and individual clades.
#' @usage search.trend(RR,y,nsim=100,clus=.5,node=NULL,cov=NULL,foldername,ConfInt=FALSE)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y the named vector (or matrix if multivariate) of phenotypes.
#' @param node the node number of individual clades to be specifically tested and contrasted to each other. It is \code{NULL} by default. Notice the node number must refer to the dichotomic version of the original tree, as produced by \code{RRphylo}.
#' @param nsim number of simulations to be performed. It is set at 100 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @param foldername the path of the folder where plots are to be found.
#' @param cov the covariate values to be specified if the RR object has been created using a  covariate for rates calculation.  As for \code{RRphylo}, \code{'cov'} must be as long as the number of nodes plus the number of tips of the tree, which can be obtained by running \code{RRphylo} on the covariate as well, and taking the vector of ancestral states and tip values to form the covariate (see the example below).
#' @param ConfInt if \code{TRUE}, the function returns 95\% confidence intervals around phenotypes and rates produced according to the Brownian motion model of evolution. It is \code{FALSE} by default.
#' @return The function returns a ‘list’ object including:
#' @return \strong{$rbt} for each branch of the tree, there are the age of the daughter node/tip (D.age), the age of the parental node (P.age), the \code{RRphylo} rates, and the distance from the tree root (age). If y is multivariate, it also includes the multiple rates for each y vector.
#' @return \strong{$pbt} a data frame of phenotypic values and their distance from the tree root for each node (i.e. ancestral states) and tip of the tree.
#' @return \strong{$phenotypic.regression} results of phenotype versus age regression. It reports a p-value for the regression slope between the variables (p.real) a p-value computed contrasting the real slope to Brownian motion simulations (p.random). Only the latter should be inspected to assess significance.
#' @return \strong{$rate.regression} results of the rates (absolute values) versus age regression. It reports a p-value for the regression between the variables (p.real), a p-value computed contrasting the real slope to Brownian motion simulations (p.random). Only the latter should be inspected to assess significance.
#' @return \strong{$ConfInts} the 95\% confidence intervals around phenotypes and rates produced according to the Brownian motion model of evolution.
#' @return If specified, individual nodes are tested as the whole tree, the results are summarized in the objects \strong{$node.phenotypic.regression} and \strong{$node.rate.regression}. If more than one node is specified, the object  \strong{$group.comparison} reports the results obtained by comparing individual clades to each other.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @details The function simultaneously returns the regression of phenotypes and phenotypic evolutionary rates against age (in terms of distance from the tree root) tested against Brownian motion simulations to assess signficance. It stores the rates (absolute values) versus age regression and the phenotype versus age regression plots as .pdf files. In the plots, the 95\% confidence intervals of phenotypes and rates simulated under the Brownian motion for each node are plotted as shaded areas. Regression lines are printed for all regressions. To assess significance, slopes are compared to a family of simulated slopes (BMslopes, where the number of simulations is equal to \code{nsim}), generated under the Brownian motion, using the \code{fastBM} function in the package \pkg{phytools}. Individual nodes are compared to the rest of the tree in different ways depending on whether phenotypes or rates versus age regressions are tested. With the former, the regression slopes for individual clades and the slope difference between clades is contrasted slopes obtained through Brownian motion simulations. For the latter, regression models are tested and contrasted to each other referring to estimated marginal means, by using the \code{emmeans} function in the package \pkg{emmeans}.
#' @importFrom graphics points text title polygon pairs
#' @importFrom stats as.formula coef resid density predict
#' @importFrom binr bins.greedy
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nlme gls varFunc
#' @importFrom outliers outlier
#' @importFrom car outlierTest
#' @importFrom emmeans emmeans
#' @export
#' @references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M., Profico, A., Girardi, G., & Raia, P. (2019) Simultaneous detection of macroevolutionary patterns in phenotypic means and rate of change with and within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14: e0210101. https://doi.org/10.1371/journal.pone.0210101
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
                        foldername, ConfInt = FALSE)
{
  #require(ape)
  #require(phytools)
  #require(geiger)
  #require(stats4)
  #require(foreach)
  #require(doParallel)
  #require(lmtest)
  #require(parallel)
  #require(binr)
  #require(nlme)
  #require(RColorBrewer)
  #require(emmeans)
  #require(outliers)
  #require(car)

  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  t <- RR$tree
  if(min(diag(vcv(t)))/max(diag(vcv(t)))>=0.9) stop("not enough fossil information")
  rates <- RR$rates
  betas <- RR$multiple.rates
  aceRR <- RR$aces
  L <- RR$tip.path
  L1 <- RR$node.path
  if (length(y) > Ntip(t)&is.null(rownames(y))) stop("The matrix of phenotypes needs to be named")
  if (length(y) == Ntip(t)&is.null(names(y))) stop("The vector of phenotypes needs to be named")

  if (class(y) == "data.frame")
    y <- treedata(t, y, sort = TRUE)[[2]]
  H <- max(nodeHeights(t))
  eds <- t$edge[, 2]
  eds[which(t$edge[, 2] < Ntip(t) + 1)] <- t$tip.label
  eds <- c(Ntip(t) + 1, eds)
  hh <- c(0.0001, nodeHeights(t)[, 2])
  eds <- data.frame(leaf = eds, height = hh)
  if (length(y) > Ntip(t)) {
    y.multi <- (L %*% rates)
    as.matrix(as.data.frame(y.multi[match(rownames(y),rownames(y.multi)),]))->y.multi
    aceRR.multi <- (L1 %*% rates[1:Nnode(t), ])
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
  data[match(b.distrib[,1],rownames(data)),]->data
  data <- data.frame(b.distrib, data)
  rbi <- data


  if (length(y) > Ntip(t)) {##### Rate Trend Real Multi #####
    rbi.rate <- rbi[, c(5:(dim(rbi)[2] - 1), dim(rbi)[2])]
    rbi.slopeA <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                             1))
    REGabs.betas <- list()
    e1<-array()
    for (i in 1:(dim(rbi.rate)[2] - 1)) {
      bet <- rbi.rate[, i]
      age <- rbi.rate[, dim(rbi.rate)[2]]
      names(bet)<-names(age)<-rownames(rbi.rate)
      as.matrix(as.data.frame(abs(bet)))->rts
      rts->rtsA
      log(range01(rts))->rts
      age->ageC

      sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->e1[i]

      c(which(rts=="-Inf"),which(age=="-Inf"))->outs
      if(length(outs)>0)
      {
        as.matrix(as.data.frame(rts[-outs,]))->rts
        age[-outs]->age
      }
      lm(rts~age)->bb

      residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
      if((Ntip(t)+1)%in%names(resout)){
        as.matrix(as.data.frame(rts[-1,]))->rts
        age[-1]->age
      }

      lm(rts~age)->regr.1
      rbi.slopeA[i, ] <- coef(summary(regr.1))[2, c(1, 4)]
      REGabs.betas[[i]] <- regr.1
    }
    colnames(rbi.slopeA) <- c("slope","p-value")
    rownames(rbi.slopeA) <- names(REGabs.betas)<-colnames(data)[5:(5 +dim(y)[2])]
  }else { #### Rate Trend Real Uni #####
    rbi.rate <- rbi[, c(5, 6)]
    bet <- rbi.rate[, 1]
    age <- rbi.rate[, 2]
    names(bet)<-names(age)<-rownames(rbi.rate)
    as.matrix(as.data.frame(abs(bet)))->rts
    rts->rtsA
    log(range01(rts))->rts
    age->ageC

    sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->e1

    c(which(rts=="-Inf"),which(age=="-Inf"))->outs
    if(length(outs)>0)
    {
      as.matrix(as.data.frame(rts[-outs,]))->rts
      age[-outs]->age
    }
    lm(rts~age)->bb

    residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
    if((Ntip(t)+1)%in%names(resout)){
      as.matrix(as.data.frame(rts[-1,]))->rts
      age[-1]->age
    }

    lm(rts~age)->regr.1
    rbi.slopeA<- coef(summary(regr.1))[2, c(1, 4)]
    REGabs <- regr.1
  }

  nodes <- aceRR[1:Nnode(t), ]

  if (length(y) > Ntip(t)) {##### Phenotypic Trend Real Multi #####
    nodes <- cbind(aceRR[1:Nnode(t), ],aceRR.multi[1:Nnode(t), ])
    colnames(nodes)<- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
    P <- rbind(nodes, cbind(y,y.multi))
    PP <- data.frame(P[match(rbi[, 1], rownames(P)), ],
                     rbi$age)
    colnames(PP)[dim(PP)[2]] <- "age"
    trendR <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x) lm(range01(x) ~ PP[, dim(PP)[2]]))
    trend.reg <- lapply(trendR, function(x) coefficients(summary(x)))
    names(trend.reg) <- c(paste("y", seq(1:dim(y)[2]), sep = ""),"multiple")
    dev<-array()
    for(i in 1:length(trend.reg)){
      if(trend.reg[[i]][2,1]<0) (min(predict(trendR[[i]]))-mean(range01(PP[,i])))/sd(range01(PP[,i]))->dev[i] else (max(predict(trendR[[i]]))-mean(range01(PP[,i])))/sd(range01(PP[,i]))->dev[i]
    }
    names(dev)<-names(trend.reg)

  }else {##### Phenotypic Trend Real Uni #####
    P <- c(nodes, y)
    PP <- data.frame(P[match(rbi[, 1], names(P))], rbi$age)
    colnames(PP) <- c("phenotype", "age")
    lm(range01(PP[,1])~PP[,2])->trendR
    summary(trendR)$coef->trend.reg
    if(trend.reg[2,1]<0) (min(predict(trendR))-mean(range01(PP[,1])))/sd(range01(PP[,1]))->dev else (max(predict(trendR))-mean(range01(PP[,1])))/sd(range01(PP[,1]))->dev

  }
  PPtot <- PP

  #### Nodes Real ####
  if (class(node) != "NULL") {
    rbi.sma <- rbi.rate
    PP.sma <- PP
    rbi.sma$group <- rep("NA", dim(rbi.sma)[1])
    trend.reg.sel <- list()
    REGabs.betas.y.sel <- list()
    REG.betas.age.sel <- list()
    trend.reg.age.sel <- list()
    trend.reg.y.sel <- list()
    for (j in 1:length(node)) {
      n <- node[j]
      sele <- getDescendants(t, n)
      sele[which(sele < (Ntip(t) + 1))] <- t$tip.label[sele[which(sele <
                                                                    (Ntip(t) + 1))]]
      rbi.sma[match(c(n,sele), rownames(rbi.sma)), ]$group <- paste("g",
                                                                    n, sep = "")
      rbi.sel <- rbi[match(c(n,sele), rownames(rbi)), ]

      if (length(y) > Ntip(t)) { #### Rate Trend Real Node Multi ####
        rbi.rate <- rbi.sel[, c(5:(dim(rbi.sel)[2] -
                                     1), dim(rbi.sel)[2])]
        rbi.slopeAn <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                  1))
        REG.betas.age <- list()
        REGabs.betas.y <- list()
        for (i in 1:(dim(rbi.rate)[2] - 1)) {
          bet <- rbi.rate[, i]
          age <- rbi.rate[, dim(rbi.rate)[2]]
          names(bet)<-names(age)<-rownames(rbi.rate)
          as.matrix(as.data.frame(abs(bet)))->rts
          rts->rtsA
          REG.betas.age[[i]]<-age
          REGabs.betas.y[[i]]<-predict(lm(rtsA~age))
          log(range01(rts))->rts
          age->ageC

          c(which(rts=="-Inf"),which(age=="-Inf"))->outs
          if(length(outs)>0)
          {
            as.matrix(as.data.frame(rts[-outs,]))->rts
            age[-outs]->age
          }
          lm(rts~age)->bb

          residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
          if((Ntip(t)+1)%in%names(resout)){
            as.matrix(as.data.frame(rts[-1,]))->rts
            age[-1]->age
          }

          lm(rts~age)->regr.1
          rbi.slopeAn[i, ] <- coef(summary(regr.1))[2, c(1, 4)]


        }

        colnames(rbi.slopeAn) <- c("slope","p-value")
        rownames(rbi.slopeAn) <- names(REG.betas.age) <- names(REGabs.betas.y) <- colnames(data)[5:(5 +dim(y)[2])]
      }else {#### Rate Trend Real Node Uni ####

        rbi.rate <- rbi.sel[, 5:6]

        bet <- rbi.rate[, 1]
        age <- rbi.rate[, 2]
        names(bet)<-names(age)<-rownames(rbi.rate)
        as.matrix(as.data.frame(abs(bet)))->rts
        rts->rtsA
        log(range01(rts))->rts
        age->ageC

        REG.betas.age<-age
        REGabs.betas.y<-predict(lm(rtsA~age))

        c(which(rts=="-Inf"),which(age=="-Inf"))->outs
        if(length(outs)>0)
        {
          as.matrix(as.data.frame(rts[-outs,]))->rts
          age[-outs]->age
        }
        lm(rts~age)->bb

        residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
        if((Ntip(t)+1)%in%names(resout)){
          as.matrix(as.data.frame(rts[-1,]))->rts
          age[-1]->age
        }

        lm(rts~age)->regr.1
        rbi.slopeAn<- coef(summary(regr.1))[2, c(1, 4)]

      }

      REGabs.betas.y.sel[[j]] <- REGabs.betas.y
      REG.betas.age.sel[[j]] <- REG.betas.age

      nodes <- aceRR[1:Nnode(t), ]
      if (length(y) > Ntip(t)) {##### Phenotypic Trend Real Node Multi #####
        nodes <- cbind(aceRR[1:Nnode(t), ],aceRR.multi[1:Nnode(t), ])
        colnames(nodes)<- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
        P <- rbind(nodes, cbind(y,y.multi))
        PP <- data.frame(P[match(rbi.sel[, 1], rownames(P)),
                           ], rbi.sel$age)
        colnames(PP)[dim(PP)[2]] <- "age"
        trend.regC <- apply(PP[1:(dim(PP)[2] - 1)],
                            2, function(x) summary(lm(range01(x) ~ PP[, dim(PP)[2]])))
        trend.regC <- lapply(trend.regC, coefficients)
        names(trend.regC) <- c(paste("y", seq(1:dim(y)[2]),sep = ""),"multiple")
        trend.reg.y.sel[[j]] <- apply(PP[1:(dim(PP)[2] -
                                              1)], 2, function(x) predict(lm(x ~ PP[, dim(PP)[2]])))

        trend.reg.age.sel[[j]] <- PP$age

      }else {##### Phenotypic Trend Real Node Multi #####
        P <- c(nodes, y)
        PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
                         rbi.sel$age)
        colnames(PP) <- c("phenotype", "age")
        summary(lm(range01(PP[,1])~PP[,2]))$coef->trend.regC
        trend.reg.age.sel[[j]] <- PP$age
        trend.reg.y.sel[[j]] <- predict(lm(PP))
      }
      trend.reg.sel[[j]] <- trend.regC
    }

    names(trend.reg.sel)<-node

    rbi.sma$group[which(rbi.sma$group == "NA")] <- "others"
    PP.sma <- cbind(PP.sma, group = rbi.sma[match(rownames(PP.sma),
                                                  rownames(rbi.sma)), ]$group)
    if (length(which(rbi.sma$group == "others")) < 3)
      rbi.sma <- rbi.sma[-which(rbi.sma$group == "others"),]

    PP.sma <- PP.sma[-which(PP.sma$group == "others"),]

    rbiRES <- rbi.sma
    if (length(y) > Ntip(t)) { #### Node Comparison Multi ####
      sma.resA <- list()
      sma.resPP <- list()

      group <- rbi.sma[, (dim(rbi.sma)[2])]
      age <- rbi.sma[, (dim(rbi.sma)[2] - 1)]
      groupPP <- PP.sma[, (dim(PP.sma)[2])]
      agePP <- PP.sma[, (dim(PP.sma)[2] - 1)]
      for (i in 1:(dim(rbi.sma)[2] - 2)) {
        bets <- rbi.sma[, i]
        dat <- data.frame(bets, age, group)
        suppressMessages(as.data.frame(pairs(emmeans(lm(abs(bets)~age*group,data=rbi.sma),specs="group")))->smaA)
        data.frame(do.call(rbind,strsplit(as.character(smaA[,1])," - ")),smaA[,-c(1,3,4,5)])->sma.resA[[i]]

        colnames(sma.resA[[i]])[1:2]<-c("group_1","group_2")

        if(length(node)>1){
          sapply(lapply(trend.reg.sel,"[[",i),"[[",2)->slope.tot
          names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")

          combn(sort(unique(as.character(groupPP))),2)->pair
          slope.diff<-array()
          for(jj in 1:dim(pair)[2]){
            slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

          }
          data.frame(t(pair),slope.diff)->sma.resPP[[i]]
          colnames(sma.resPP[[i]])<-colnames(sma.resA[[i]])[1:3]

        }else{
          sma.resPP<-NULL
        }
      }


      lapply(sma.resA,function(x) subset(x,x$group_2=="others"))->n.ot
      lapply(sma.resA,function(x) subset(x,x$group_2!="others"))->sma.resA

      sapply(strsplit(as.character(n.ot[[1]][,1]),"g"),"[[",2)->grn
      lapply(n.ot,function(x) x[match(node,grn),])->n.ot
      rbi.slopeA.sel<-list()
      for(i in 1:length(node)){
        do.call(rbind,lapply(n.ot,function(x) x[i,3:4]))->rbi.slopeA.sel[[i]]
        rownames(rbi.slopeA.sel[[i]])<-rownames(rbi.slopeA)
      }



    }else {#### Node Comparison Uni ####

      colnames(PP.sma)[1] <- "y"
      suppressMessages(sma.resA<-as.data.frame(pairs(emmeans(lm(abs(rate)~age*group,data=rbi.sma),specs="group"))))
      data.frame(do.call(rbind,strsplit(as.character(sma.resA[,1])," - ")),sma.resA[,-c(1,3,4,5)])->sma.resA
      colnames(sma.resA)[1:2]<-c("group_1","group_2")

      subset(sma.resA,sma.resA$group_2=="others")->n.ot
      subset(sma.resA,sma.resA$group_2!="others")->sma.resA

      sapply(strsplit(as.character(n.ot[,1]),"g"),"[[",2)->grn
      n.ot[match(node,grn),]->n.ot
      rownames(n.ot)<-NULL
      rbi.slopeA.sel<-list()
      for(i in 1:nrow(n.ot)){
        n.ot[i,c(3,4)]->rbi.slopeA.sel[[i]]
        rownames(rbi.slopeA.sel[[i]])<-NULL
      }


      if(dim(sma.resA)[1]>0){
        PP.sma[,3]<-as.character(PP.sma[,3])
        sapply(trend.reg.sel,"[[",2)->slope.tot
        names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")

        combn(sort(unique(PP.sma[,3])),2)->pair
        slope.diff<-array()
        for(jj in 1:dim(pair)[2]){
          slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

        }
        data.frame(t(pair),slope.diff)->sma.resPP
        colnames(sma.resPP)<-colnames(sma.resA)[1:3]

      }else{
        sma.resPP<-NULL
        sma.resA<-NULL
      }
    }
    names(rbi.slopeA.sel)<-node
    if (class(sma.resA) == "list") {
      if(is.null(sma.resPP)==FALSE){
        names(sma.resA) <- colnames(rbi.sma)[1:(dim(rbi.sma)[2] - 2)]
        names(sma.resPP) <- names(trend.reg)
      }
    }
  }


  #### Random ####
  RR$aces[1,]->a

  if (length(y) > Ntip(t)) {
    yyD <- list()
    yyT <- list()
    for (i in 1:dim(y)[2]) {
      yyD[[i]] <- replicate(nsim,fastBM(t, sig2 = 1, a = a[i], bounds=c(min(y[,i]),max(y[,i]))))
      yyT[[i]] <- replicate(nsim,fastBM(t, sig2 = 1, a = mean(y[,i]), bounds=c(min(y[,i]),max(y[,i]))))
    }
  }else {
    yyD<-replicate(nsim,fastBM(t,sig2=1,a=a,bounds=c(min(y),max(y))))
    yyT<-replicate(nsim,fastBM(t,sig2=1,a=mean(y),bounds=c(min(y),max(y))))
  }


  res <- list()
  cl <- makeCluster(round((detectCores() * clus), 0))
  registerDoParallel(cl)
  res <- foreach(i = 1:nsim, .packages = c("car","outliers","nlme", "ape",
                                           "geiger", "phytools", "penalized", "doParallel", "lmtest"
  )) %dopar% {
    gc()
    if (length(y) > Ntip(t)) {
      vec <- seq(1:nsim)
      yD <- lapply(yyD, function(x) x[, sample(vec, 1, replace = FALSE)])
      yD <- do.call(cbind, yD)

      yT <- lapply(yyT, function(x) x[, sample(vec, 1, replace = FALSE)])
      yT <- do.call(cbind, yT)

      betasT<- matrix(ncol=dim(yT)[2], nrow=Ntip(t)+Nnode(t))
      aceRRT<- matrix(ncol=dim(yT)[2], nrow=Nnode(t))
      betasD<- matrix(ncol=dim(yD)[2], nrow=Ntip(t)+Nnode(t))
      aceRRD<- matrix(ncol=dim(yD)[2], nrow=Nnode(t))
      for (s in 1:dim(yT)[2]){
        rootV<-a[s]
        lambda <- RR$lambda[s]
        m.betasT <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                       t(L)) %*% (as.matrix(yT[,s]) - rootV)
        aceRRT[,s] <- (L1 %*% m.betasT[1:Nnode(t), ]) + rootV
        m.betasT->betasT[,s]

        m.betasD <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                       t(L)) %*% (as.matrix(yD[,s]) - rootV)
        aceRRD[,s] <- (L1 %*% m.betasD[1:Nnode(t), ]) + rootV
        m.betasD->betasD[,s]


      }
      rownames(betasT)<-rownames(betasD)<-rownames(RR$rates)
      rownames(aceRRT)<-rownames(aceRRD)<-rownames(RR$aces)

      ratesT<-as.matrix(as.data.frame(apply(betasT, 1, function(x) sqrt(sum(x^2)))))
      aceRRTmulti <- (L1 %*% ratesT[1:Nnode(t), ])
      yTmulti <- (L %*% ratesT)
      yTmulti<-as.matrix(as.data.frame(yTmulti[match(rownames(yT),rownames(yTmulti)),]))

      ratesD<-as.matrix(as.data.frame(apply(betasD, 1, function(x) sqrt(sum(x^2)))))
      aceRRDmulti <- (L1 %*% ratesD[1:Nnode(t), ])
      yDmulti <- (L %*% ratesD)+aceRRDmulti[1,]
      yDmulti<-as.matrix(as.data.frame(yDmulti[match(rownames(yD),rownames(yDmulti)),]))

    }else {
      yT <- yyT[, i]
      yD <- yyD[, i]

      rootV<-a
      lambda <- RR$lambda
      betasT <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                   t(L)) %*% (as.matrix(yT) - rootV)
      aceRRT <- (L1 %*% betasT[1:Nnode(t), ]) + rootV

      betasD <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                   t(L)) %*% (as.matrix(yD) - rootV)
      aceRRD <- (L1 %*% betasD[1:Nnode(t), ]) + rootV
    }


    #### Covariate ####
    if (is.null(cov) == FALSE) {
      if (length(yT) > Ntip(t)) {
        if (length(which(apply(betasT, 1, sum) == 0)) >
            0) {
          zeroes <- which(apply(betasT, 1, sum) == 0)
          R <- log(abs(betasT))
          R <- R[-zeroes, ]
          Y <- abs(cov)
          Y <- Y[-zeroes]
          resi <- residuals(lm(R ~ Y))
          factOut <- which(apply(betasT, 1, sum) != 0)
          betasT[factOut, ] <- resi
          betasT[zeroes, ] <- 0
        }else {
          R <- log(abs(betasT))
          Y <- abs(cov)
          resi <- residuals(lm(R ~ Y))
          betasT <- as.matrix(resi)
        }
      }else {
        if (length(which(betasT == "0")) > 0) {
          zeroes <- which(betasT == "0")
          R <- log(abs(betasT))
          R <- R[-zeroes]
          Y <- abs(cov)
          Y <- Y[-zeroes]
          resi <- residuals(lm(R ~ Y))
          factOut <- which(betasT != "0")
          betasT[factOut] <- resi
          betasT[zeroes] <- 0
        }else {
          R <- log(abs(betasT))
          Y <- abs(cov)
          resi <- residuals(lm(R ~ Y))
          betasT <- as.matrix(resi)
        }
      }
    }
    #### End of Covariate ####

    eds <- t$edge[, 2]
    eds[which(t$edge[, 2] < Ntip(t) + 1)] <- t$tip.label
    eds <- c(Ntip(t) + 1, eds)
    hh <- c(0.0001, nodeHeights(t)[, 2])
    eds <- data.frame(leaf = eds, height = hh)

    if (length(y) > Ntip(t)) {
      data <- data.frame(betas = betasT[match(eds[, 1],
                                              rownames(betasT)), ], rate = ratesT[match(eds[,
                                                                                            1], rownames(ratesT)), ], age = eds[, 2])
      colnames(data)[1:dim(yT)[2]] <- paste("betas", seq(1,
                                                         dim(yT)[2], 1), sep = "")
    }else {
      rates <- betasT
      data <- data.frame(rate = rates[match(eds[, 1],
                                            rownames(rates)), ], age = eds[, 2])
    }
    data[,dim(data)[2]]+L[1,1]->data[,dim(data)[2]]
    data[match(b.distrib[,1],rownames(data)),]->data
    data <- data.frame(b.distrib, data)
    rbi <- data
    if (length(yT) > Ntip(t)) {#### Rate Trend Random Multi ####
      rbi.rate <- rbi[, c(5:(dim(rbi)[2] - 1), dim(rbi)[2])]
      rbi.slopeAS <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                1))
      ee<-array()
      for (k in 1:(dim(rbi.rate)[2] - 1)) {
        bet <- rbi.rate[, k]
        age <- rbi.rate[, dim(rbi.rate)[2]]
        names(bet)<-names(age)<-rownames(rbi.rate)

        as.matrix(as.data.frame(abs(bet)))->rts
        rts->rtsA
        log(range01(rts))->rts
        age->ageC

        sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->ee[k]

        c(which(rts=="-Inf"),which(age=="-Inf"))->outs
        if(length(outs)>0)
        {
          as.matrix(as.data.frame(rts[-outs,]))->rts
          age[-outs]->age
        }
        lm(rts~age)->bb

        residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
        if((Ntip(t)+1)%in%names(resout)){
          as.matrix(as.data.frame(rts[-1,]))->rts
          age[-1]->age
        }

        lm(rts~age)->regr.1
        rbi.slopeAS[k, ] <- coef(summary(regr.1))[2, c(1, 4)]

      }
      colnames(rbi.slopeAS) <- c("slope","p-value")
      rownames(rbi.slopeAS) <- colnames(data)[5:(5 +dim(y)[2])]
    }else {#### Rate Trend Random Uni ####
      rbi.rate <- rbi[, 5:6]
      bet <- rbi.rate[, 1]
      age <- rbi.rate[, 2]
      names(bet)<-names(age)<-rownames(rbi.rate)
      as.matrix(as.data.frame(abs(bet)))->rts
      rts->rtsA
      log(range01(rts))->rts
      age->ageC

      sd(range01(rtsA)[ageC<0.5*max(ageC)])/sd(range01(rtsA)[ageC>0.5*max(ageC)])->ee

      c(which(rts=="-Inf"),which(age=="-Inf"))->outs
      if(length(outs)>0)
      {
        as.matrix(as.data.frame(rts[-outs,]))->rts
        age[-outs]->age
      }
      lm(rts~age)->bb

      residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
      if((Ntip(t)+1)%in%names(resout)){
        as.matrix(as.data.frame(rts[-1,]))->rts
        age[-1]->age
      }

      lm(rts~age)->regr.1
      rbi.slopeAS <- coef(summary(regr.1))[2, c(1, 4)]

    }
    rbi.sma <- rbi.rate

    nodesD <- aceRRD[1:Nnode(t), ]
    if (length(yD) > Ntip(t)) { #### Phenotypic Trend Random Multi ####
      nodesD<-cbind(nodesD,aceRRDmulti[1:Nnode(t), ])
      colnames(nodesD) <- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
      P <- rbind(nodesD, cbind(yD,yDmulti))
      PP <- data.frame(P[match(rbi[, 1], rownames(P)),
                         ], rbi$age)
      colnames(PP)[dim(PP)[2]] <- "age"
      trend.regR <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x)
        summary(lm(range01(x) ~PP[, dim(PP)[2]])))
      trend.regS <- lapply(trend.regR, coefficients)
      names(trend.regR) <- c(paste("betas", seq(1:dim(y)[2]),sep = ""),"multiple")
      cbind(yT,yTmulti)->yT
      cbind(yD,yDmulti)->yD

    }else { #### Phenotypic Trend Random Uni ####
      P <- c(nodesD, yD)
      PP <- data.frame(P[match(rbi[, 1], names(P))], rbi$age)
      colnames(PP) <- c("phenotype", "age")
      trend.regS <- summary(lm(range01(PP[,1])~PP[,2]))$coef
    }

    PPtot <- PP

    if (class(node) != "NULL") {
      PP.sma <- PP
      PP.sma$group <- rep("NA", dim(PP.sma)[1])
      rbi.slopeAS.sel <- list()
      trend.reg.SEL <- list()
      for (j in 1:length(node)) {
        n <- node[j]
        sele <- getDescendants(t, n)
        sele[which(sele < (Ntip(t) + 1))] <- t$tip.label[sele[which(sele <
                                                                      (Ntip(t) + 1))]]
        PP.sma[match(c(n,sele), rownames(PP.sma)), ]$group <- paste("g",
                                                                    n, sep = "")
        rbi.sel <- rbi[match(c(n,sele), rownames(rbi)), ]

        nodesD <- aceRRD[1:Nnode(t), ]
        if (length(y) > Ntip(t)) {#### Phenotypic Trend Random Node Multi ####
          nodesD<-cbind(nodesD,aceRRDmulti[1:Nnode(t),])
          colnames(nodesD)<- c(paste("y",seq(1, dim(y)[2]), sep = ""),"multiple")
          P <- rbind(nodesD, yD)
          PP <- data.frame(P[match(rbi.sel[, 1], rownames(P)),
                             ], rbi.sel$age)
          colnames(PP)[dim(PP)[2]] <- "age"
          trend.regC <- apply(PP[1:(dim(PP)[2] - 1)],
                              2, function(x) summary(lm(range01(x) ~ PP[, dim(PP)[2]])))
          trend.regC <- lapply(trend.regC, coefficients)
          names(trend.regC) <- c(paste("y", seq(1:dim(y)[2]),sep = ""),"multiple")

        }else {
          P <- c(nodesD, yD)
          PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
                           rbi.sel$age)
          colnames(PP) <- c("phenotype", "age")
          trend.regC <- summary(lm(range01(PP[,1])~PP[,2]))$coef

        }
        trend.reg.SEL[[j]] <- trend.regC
      }
      names(trend.reg.SEL) <- node


      if(is.null(sma.resPP)==FALSE){
        PP.sma <- PP.sma[-which(PP.sma$group == "NA"),]

        if (length(y) > Ntip(t)) {#### Node Comparison Random Multi ####
          groupPP <- PP.sma[, (dim(PP.sma)[2])]
          agePP <- PP.sma[, (dim(PP.sma)[2] - 1)]

          sma.resPP.R<-list()
          for (k in 1:(dim(rbi.rate)[2] - 1)) {

            sapply(lapply(trend.reg.SEL,"[[",k),"[[",2)->slope.tot
            names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")


            combn(sort(unique(as.character(groupPP))),2)->pair
            slope.diff<-array()
            for(jj in 1:dim(pair)[2]){
              slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

            }
            data.frame(t(pair),slope.diff)->sma.resPP.R[[k]]
            colnames(sma.resPP.R[[k]])<-c("group_1","group_2","estimate")

          }
          names(sma.resPP.R)<-colnames(PP.sma)[1:(dim(PP.sma)[2] - 2)]
        }else{#### Node Comparison Random Uni ####

          PP.sma[,3]<-as.character(PP.sma[,3])
          sapply(trend.reg.SEL,"[[",2)->slope.tot
          names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")

          combn(sort(unique(PP.sma[,3])),2)->pair
          slope.diff<-array()
          for(jj in 1:dim(pair)[2]){
            slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

          }
          data.frame(t(pair),slope.diff)->sma.resPP.R
          colnames(sma.resPP.R)<-c("group_1","group_2","estimate")

        }
      }else{
        sma.resPP.R<-NULL
      }
      res[[i]] <- list(trend.regS, rbi.slopeAS,
                       PPtot, rbi,yT,yD,ee,
                       rbi.slopeAS.sel, trend.reg.SEL,sma.resPP.R)
    }else {
      res[[i]] <- list(trend.regS, rbi.slopeAS,
                       PPtot, rbi,yT,yD,ee)
    }
  }
  stopCluster(cl)
  #### End of Random ####

  if (length(y) > Ntip(t)) {#### p Rate Trend Multi####
    p.rbi.slopeA <- array()
    spread<-array()
    for (i in 1:(dim(y)[2] + 1)) {
      cbind(y,y.multi)->yTot
      rbi.slopeRAS <- do.call(rbind, lapply(lapply(res,"[[", 2), function(x) x[i, 1]))
      yRT<- lapply(lapply(res,"[[",5), function(x) x[, i])

      nsim-length(which(rbi.slopeRAS<rbi.slopeA[i, 1]))->pp

      sapply(lapply(res,"[[",7),"[[",i)->ee

      if(pp>(nsim*.5) & length(which(ee>e1[i]))>=(nsim*.95) & rbi.slopeA[i, 1]>0){
        nsim-pp->pp
      } else {
        if(pp<(nsim*.5) & length(which(ee>e1[i]))<=(nsim*.05) & rbi.slopeA[i, 1]<0)
          nsim-pp->pp
      }
      pp/nsim->p.rbi.slopeA[i]

      (diff(range(yTot[,i]))/diff(range(yTot[,i][which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))*
        mean(sapply(yRT,function(x)(diff(range(x))/diff(range(x[which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))))->spread[i]

    }
    names(spread)<-names(p.rbi.slopeA) <- rownames(rbi.slopeA)
    p.rbi.slopeA <- data.frame(slope = rbi.slopeA[, 1],
                               p.real = rbi.slopeA[, 2], p.random = p.rbi.slopeA,spread=spread)
  }else {#### p Rate Trend Uni ####
    rbi.slopesRAS <- do.call(rbind, lapply(res, "[[", 2))[,
                                                          1]
    yRT<-lapply(res, "[[", 5)

    nsim-length(which(rbi.slopesRAS<rbi.slopeA[1]))->pp

    sapply(res,"[[",7)->ee

    if(pp>(nsim*.5) & length(which(ee>e1))>=(nsim*.95) & rbi.slopeA[1]>0){
      nsim-pp->pp
    } else {
      if(pp<(nsim*.5) & length(which(ee>e1))<=(nsim*.05) & rbi.slopeA[1]<0)
        nsim-pp->pp
    }
    pp/nsim->p.rbi.slopeA

    (diff(range(y))/diff(range(y[which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))*
      mean(sapply(yRT,function(x)(diff(range(x))/diff(range(x[which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))))->spread

    rbi.slopeA <- unname(rbi.slopeA)
    p.rbi.slopeA <- c(slope = rbi.slopeA[1], p.real = rbi.slopeA[2],
                      p.random = p.rbi.slopeA,spread=spread)
  }


  if (class(node) != "NULL") {
    p.rbi.slopeA.sel<-rbi.slopeA.sel
    p.smaA <- sma.resA

    if(is.null(sma.resPP)==FALSE){
      if (length(y) > Ntip(t)) {#### p Phenotypic Slope Comparison Multi  ####
        p.smaPP<-list()
        for (k in 1:(dim(y)[2] + 1)) {
          p.smaR<-array()
          for(i in 1:nrow(sma.resPP[[k]])){
            rank(c(sma.resPP[[k]][i,3],sapply(lapply(lapply(res,"[[",10),"[[",k),function(x) x[i,3])[1:(nsim-1)]))[1]/nsim->p.smaR[i]
          }
          data.frame(sma.resPP[[k]],p.value=1-p.smaR)->p.smaPP[[k]]
        }
        names(p.smaPP)<-names(trend.reg)
      }else{
        p.smaR<-array()
        for(i in 1:nrow(sma.resPP)){
          rank(c(sma.resPP[i,3],sapply(lapply(res,"[[",10),function(x) x[i,3])[1:(nsim-1)]))[1]/nsim->p.smaR[i]
        }
        p.smaPP <-data.frame(sma.resPP,p.value=1-p.smaR)
      }
    }else{
      p.smaPP<-NULL
    }
  }

  p.trend <- array()
  PP <- PPtot
  if (length(y) > Ntip(t)) { #### p Phenotypic Trend Multi and pdf ####
    trend.slopes <- matrix(ncol = dim(y)[2] + 1, nrow = nsim)
    for (i in 1:(dim(y)[2] + 1)) {
      trend.slopes[, i] <- unlist(lapply(lapply(lapply(res,
                                                       "[[", 1), "[[", i), function(x) x[2,1]))

      p.trend[i] <- (nsim-length(which(trend.slopes[, i]<trend.reg[[i]][2,1])))/nsim
    }
    trend.real <- do.call(rbind, lapply(trend.reg, function(x) x[2,
                                                                 c(1, 4)]))
    p.trend <- cbind(trend.real, p.trend,dev)
    colnames(p.trend) <- c("slope", "p.real",
                           "p.random","dev")

    CIphenotype <- list()
    if (dim(y)[2] <= 3) {
      pdf(file = paste(foldername, "Phenotypic Trend Test.pdf",
                       sep = "/"))
      par(mar = c(3.5, 3.5, 1, 2))
      par(mfrow = c(dim(y)[2] + 1, 2))
      for (i in 1:(dim(y)[2] + 1)) {
        PPci <- apply(do.call(cbind, lapply(lapply(res,
                                                   "[[", 3), function(x) x[, i])), 1, function(u) quantile(u,
                                                                                                           c(0.025, 0.975)))
        colnames(PPci) <- lapply(lapply(res, "[[",
                                        3), rownames)[[1]]
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
        points((diag(vcv(t))+PP[1, dim(PP)[2]]), yTot[, i], pch = 21, col = "black",
               bg = "red")
        if (class(node) != "NULL") {
          for (j in 1:length(node)) {
            cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
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
    }else {
      pdf(file = paste(foldername, "Phenotypic Trend Test.pdf",
                       sep = "/"))
      par(mar = c(3.5, 3.5, 1, 2))
      par(mfrow = c(2, 1))
      i <- dim(yTot)[2]
      obj <- hist(trend.slopes[, i], xlab = "", ylab = "",
                  main = "Phenotypic Trend Test", mgp = c(2, 0.5,
                                                          0))
      title(xlab = "simulated slopes", ylab = "frequency",
            line = 1.5)
      text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
             0.9, labels = paste("p=", p.trend[i, 3]), cex = 1)
      abline(v = trend.reg[[i]][2, 1], lwd = 3, col = "red")
      plot(PP[, c(dim(PP)[2], i)], xlab = "", ylab = "",mgp = c(2, 0.5, 0))
      points((diag(vcv(t))+PP[1, dim(PP)[2]]), yTot[, i], pch = 21, col = "black",
             bg = "red")
      title(xlab = "age", ylab = "y.multi",
            line = 1.5)
      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
          points(trend.reg.age.sel[[j]], trend.reg.y.sel[[j]][,
                                                              i], lwd = 4, col = cols[j], type = "l")
        }
        abline(lm(PP[, i] ~ age, data = PP), lwd = 3,
               col = "blue", lty = 2)
        legend(min(PP$age), max(PP[, i]), legend = node,
               fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                0, 0, 0), border = NA, x.intersp = 0.25)
      }else {
        abline(lm(PP[, i] ~ age, data = PP), lwd = 4,
               col = "blue")
      }
    }
    dev.off()
  }else {#### p Phenotypic Trend Uni and pdf ####
    trend.slopes <- do.call(rbind, lapply(lapply(res, "[[", 1),function(x) x[2,1]))

    p.trend <- (nsim-length(which(trend.slopes<trend.reg[2,1])))/nsim

    p.trend <- c(trend.reg[2, c(1, 4)],  p.trend, dev)
    names(p.trend) <- c("slope", "p.real", "p.random","dev")
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
                                               3), function(x) x[, 1])), 1, function(u) quantile(u,
                                                                                                 c(0.025, 0.975)))
    colnames(PPci) <- lapply(lapply(res, "[[", 3), rownames)[[1]]
    CIphenotype <- t(PPci)
    PPci <- cbind(PP, t(PPci[, match(rownames(PP), colnames(PPci))]))
    PPci <- PPci[order(PPci[, 2]), ]
    polygon(c(PPci[, 2], rev(PPci[, 2])), c(PPci[, 3], rev(PPci[,
                                                                4])), col = rgb(0.5, 0.5, 0.5, 0.4), border = NA)
    points((diag(vcv(t))+PP[1, 2]), y, pch = 21, col = "black", bg = "red")
    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
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
      if (length(y) > Ntip(t)) {#### p Phenotypic Trend Node Multi ####
        p.sele <- list()
        for (i in 1:(dim(y)[2] + 1)) {
          slopeR <- unlist(lapply(lapply(lapply(lapply(res,
                                                       "[[", 9), "[[", u), "[[", i), function(x) x[2,
                                                                                                   1]))

          p.sel <- (nsim-length(which(slopeR<trend.reg.sel[[u]][[i]][2,1])))/nsim
          p.sele[[i]] <- c(estimate = trend.reg.sel[[u]][[i]][2,
                                                              1],p.value = p.sel)
        }
        names(p.sele) <- names(trend.reg.sel[[u]])
        p.selt <- do.call(rbind, p.sele)
        p.trend.sel[[u]] <- p.selt
      } else {
        slopeR <- unlist(lapply(lapply(lapply(res, "[[",
                                              9), "[[", u), function(x) x[2,1]))

        p.sel <- (nsim-length(which(slopeR<trend.reg.sel[[u]][2,1])))/nsim
        p.trend.sel[[u]] <- c(estimate = trend.reg.sel[[u]][2,1],p.value = p.sel)
      }
    }
    names(p.trend.sel) <- names(trend.reg.sel)
  }

  rbi <- rbi[, -4]
  if (length(y) > Ntip(t)) { #### Rate pdf Multi ####
    A <- rbi[, 4:dim(rbi)[2]]
    CIabsolute <- list()
    if (dim(y)[2] <= 3) {
      pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                       sep = "/"))
      par(mfrow = c(dim(y)[2] + 1, 1))
      par(mar = c(3.5, 3.5, 1, 1))
      for (i in 1:(dim(y)[2] + 1)) {
        if (i == dim(y)[2] + 1) ynam <- "rate" else ynam <- paste("betas", i, sep = "")
        bet <- A[, i]
        age <- A[, dim(A)[2]]
        names(bet)<-names(age)<-rownames(A)
        RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                            "[[", 4), function(x) x[, c(5:dim(x)[2])]),
                                              function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                                 c(0.025, 0.975)))
        colnames(RBTAci) <-  lapply(lapply(res,"[[", 4), rownames)[[1]]
        CIabsolute[[i]] <- t(RBTAci)
        RBTAci <- cbind(A[, c(i, dim(A)[2])], t(RBTAci[,
                                                       match(rownames(A), colnames(RBTAci))]))
        RBTAci <- RBTAci[order(RBTAci[, 2]), ]
        plot(abs(bet) ~ age, main = " ",
             ylab = ynam, mgp = c(2, 0.5, 0))
        polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                           3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5,
                                                                                            0.5, 0.4), border = NA)
        points(age[which(names(age)%in%t$tip.label)], abs(bet[which(names(bet)%in%t$tip.label)]), pch = 21, col = "black",
               bg = "red")

        if (class(node) != "NULL") {
          for (j in 1:length(node)) {
            cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
            points(REG.betas.age.sel[[j]][[i]], REGabs.betas.y.sel[[j]][[i]],
                   lwd = 4, col = cols[j], type = "l")
          }
          abline(lm(abs(bet) ~ age), lwd = 3, col = "blue",
                 lty = 2)
          if (i == 1)
            legend(min(age), max(abs(bet)), legend = node,
                   fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                    0, 0, 0), border = NA, x.intersp = 0.25)
        } else {
          abline(lm(abs(bet) ~ age), lwd = 4, col = "blue")
        }

      }
    } else {
      for (i in 1:(dim(y)[2] + 1)) {
        RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                            "[[", 4), function(x) x[, c(5:dim(x)[2])]),
                                              function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                                 c(0.025, 0.975)))
        colnames(RBTAci) <- lapply(lapply(res,"[[", 4), rownames)[[1]]
        CIabsolute[[i]] <- t(RBTAci)
      }
      pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                       sep = "/"))
      bet <- A[, (dim(y) + 1)[2]]
      age <- A[, dim(A)[2]]
      names(bet)<-names(age)<-rownames(A)
      RBTAci <- apply(do.call(cbind, lapply(lapply(lapply(res,
                                                          "[[", 4), function(x) x[, c(5:dim(x)[2])]),
                                            function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                               c(0.025, 0.975)))
      colnames(RBTAci) <- lapply(lapply(res,"[[", 4), rownames)[[1]]
      RBTAci <- cbind(A[, c(i, dim(A)[2])], t(RBTAci[,
                                                     match(rownames(A), colnames(RBTAci))]))
      RBTAci <- RBTAci[order(RBTAci[, 2]), ]
      plot(abs(bet) ~ age, ylab = "absolute rate", mgp = c(2,
                                                           0.5, 0))
      polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                         3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                          0.4), border = NA)

      points(age[which(names(age)%in%t$tip.label)], abs(bet[which(names(bet)%in%t$tip.label)]), pch = 21, col = "black",
             bg = "red")
      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
          points(REG.betas.age.sel[[j]][[(dim(y) + 1)[2]]],
                 REGabs.betas.y.sel[[j]][[(dim(y) + 1)[2]]],
                 lwd = 4, col = cols[j], type = "l")
        }
        abline(lm(abs(bet)~age), lwd = 3,
               col = "blue", lty = 2)
        legend(min(age), max(abs(bet)), legend = node,
               fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                0, 0, 0), border = NA, x.intersp = 0.25)
      }else {
        abline(lm(abs(bet)~age), lwd = 4,
               col = "blue")
      }
    }
    names(CIabsolute) <- c(paste("betas",seq(1, dim(y)[2]), sep = ""), "rate")
    dev.off()
  }else {#### Rate pdf Uni ####
    A <- rbi
    AA <- A[, 4:5]
    RBTAci <- apply(do.call(cbind, lapply(lapply(res, "[[",
                                                 4), function(x) abs(x[, 5]))), 1, function(u) quantile(u,
                                                                                                        c(0.025, 0.975)))
    colnames(RBTAci) <- lapply(lapply(res,"[[", 4), rownames)[[1]]
    CIabsolute <- t(RBTAci)
    RBTAci <- cbind(AA, t(RBTAci[, match(rownames(AA), colnames(RBTAci))]))
    RBTAci <- RBTAci[order(RBTAci[, 2]), ]
    pdf(file = paste(foldername, "Evolutionary Rate Trend Test.pdf",
                     sep = "/"))
    plot(abs(rate) ~ age, data = AA, ylab = "absolute rate",
         mgp = c(2, 0.5, 0))
    polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                       3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                        0.4), border = NA)

    points(AA$age[which(rownames(AA)%in%t$tip.label)], abs(AA$rate[which(rownames(AA)%in%t$tip.label)]), pch = 21, col = "black",
           bg = "red")

    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        cols <- suppressWarnings(brewer.pal(length(node), "Set2"))
        points(REG.betas.age.sel[[j]], REGabs.betas.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(lm(abs(rate) ~ age, data = AA), col = "blue", lwd = 3, lty = 2)
      legend(min(AA$age), max(abs(AA$rate)), legend = node,
             fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                              0, 0, 0), border = NA, x.intersp = 0.25)
    }else {
      abline(lm(abs(rate) ~ age, data = AA), col = "blue", lwd = 4)
    }
    dev.off()
  }
  colnames(rbi)[1] <- "branch"

  if (length(y) > Ntip(t)) {
    for (i in 1:(dim(y)[2] + 1)) {
      if (i == dim(y)[2] + 1){
        if (p.trend[i, 3] <= 0.05) print("There is a trend for increase in phenotype y.multi")
        if (p.trend[i, 3] >= 0.95) print("There is a trend for decrease in phenotype y.multi")
        if (p.rbi.slopeA[i, 3] <= 0.05) print("There is a trend for increase in absolute evolutionary rates for variable y.multi")
        if (p.rbi.slopeA[i, 3] >= 0.95) print("There is a trend for decrease in absolute evolutionary rates for variable y.multi")
      }else{
        if (p.trend[i, 3] <= 0.05) print(paste("There is a trend for increase in phenotype y",i, sep = ""))
        if (p.trend[i, 3] >= 0.95) print(paste("There is a trend for decrease in phenotype y",i, sep = ""))
        if (p.rbi.slopeA[i, 3] <= 0.05) print(paste("There is a trend for increase in absolute evolutionary rates for variable y",i, sep = ""))
        if (p.rbi.slopeA[i, 3] >= 0.95) print(paste("There is a trend for decrease in absolute evolutionary rates for variable y",i, sep = ""))
      }
    }
  }else {
    if (p.trend[3] <= 0.05) print("There is a trend for increase in phenotype")
    if (p.trend[3] >= 0.95) print("There is a trend for decrease in phenotype")
    if (p.rbi.slopeA[3] <= 0.05) print("There is a trend for increase in absolute evolutionary rates")
    if (p.rbi.slopeA[3] >= 0.95) print("There is a trend for decrease in absolute evolutionary rates")
  }

  if (class(node) != "NULL") {
    if (length(y) > Ntip(t)) {
      for (i in 1:length(p.rbi.slopeA.sel)) {
        for(k in 1:dim(p.rbi.slopeA.sel[[i]])[1]){
          if(k==dim(p.rbi.slopeA.sel[[i]])[1]) ".multi"-> vv else k->vv
          if(p.rbi.slopeA.sel[[i]][k,2]<=.05&p.rbi.slopeA.sel[[i]][k,1]>0)
            print(paste("There is a trend for increase in absolute evolutionary rates for variable",paste("y",vv,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))
          if(p.rbi.slopeA.sel[[i]][k,2]<=.05&p.rbi.slopeA.sel[[i]][k,1]<0)
            print(paste("There is a trend for decrease in absolute evolutionary rates for variable",paste("y",vv,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))
          if(p.trend.sel[[i]][k,2]<=.05)
            print(paste("There is a trend for increase in phenotype for variable",paste("y",vv,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))
          if(p.trend.sel[[i]][k,2]>=.95)
            print(paste("There is a trend for decrease in phenotype for variable",paste("y",vv,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))



        }
      }


      if(is.null(p.smaPP)==FALSE){
        for (j in 1:length(p.smaPP)) {
          if(j==length(p.smaPP)) ".multi"->vv else j->vv
          for (i in 1:dim(p.smaPP[[j]])[1]) {
            if (p.smaPP[[j]][i, 4] >= 0.95){
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is lower than regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    2]), "g"), "[[", 2), "for variable",
                          paste("y",vv,sep="")))
            }
            if (p.smaPP[[j]][i, 4] <= 0.05){
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is higher than regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    2]), "g"), "[[", 2), "for variable",
                          paste("y",vv,sep="")))

            }
            if (p.smaA[[j]][i, 4] <= 0.05){
              if(p.smaA[[j]][i, 3] >0){
                print(paste("Absolute evolutionary rates regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     1]), "g"), "[[", 2), "is higher than regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     2]), "g"), "[[", 2), "for variable",
                            paste("y",vv,sep="")))
              }else{
                print(paste("Absolute evolutionary rates regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     1]), "g"), "[[", 2), "is lower than regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     2]), "g"), "[[", 2), "for variable",
                            paste("y",vv,sep="")))
              }
            }
          }
        }
      }
    }else {
      for (i in 1:length(p.rbi.slopeA.sel)) {
        if(p.rbi.slopeA.sel[[i]][,2]<=.05&p.rbi.slopeA.sel[[i]][,1]>0)
          print(paste("There is a trend for increase in absolute evolutionary rates through node",names(p.rbi.slopeA.sel)[i]))
        if(p.rbi.slopeA.sel[[i]][,2]<=.05&p.rbi.slopeA.sel[[i]][,1]<0)
          print(paste("There is a trend for decrease in absolute evolutionary rates through node",names(p.rbi.slopeA.sel)[i]))
        if(p.trend.sel[[i]][2]<=.05)
          print(paste("There is a trend for increase in phenotype through node",names(p.rbi.slopeA.sel)[i]))
        if(p.trend.sel[[i]][2]>=.95)
          print(paste("There is a trend for decrease in phenotype through node",names(p.rbi.slopeA.sel)[i]))
      }

      if(is.null(p.smaPP)==FALSE){
        for(i in 1:dim(p.smaPP)[1]){
          if (p.smaPP[i, 4] >= 0.95){
            print(paste("Phenotypic regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             1]), "g"), "[[", 2), "is lower than regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             2]), "g"), "[[", 2)))
          }

          if (p.smaPP[i, 4] <= 0.05){
            print(paste("Phenotypic regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             1]), "g"), "[[", 2), "is higher than regression through node",
                        lapply(strsplit(as.character(p.smaPP[i,
                                                             2]), "g"), "[[", 2)))

          }
        }

        if (p.smaA[i, 4] <= 0.05){
          if (p.smaA[i, 3]>0){
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            1]), "g"), "[[", 2), "is higher than regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            2]), "g"), "[[", 2)))

          }else{
            print(paste("Absolute evolutionary rates regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            1]), "g"), "[[", 2), "is lower than regression through node",
                        lapply(strsplit(as.character(p.smaA[i,
                                                            2]), "g"), "[[", 2)))

          }
        }
      }

    }

    if(is.null(p.smaPP)){
      SMA.res<-NULL
    }else{
      SMA.res <- list(p.smaPP, p.smaA)
      names(SMA.res) <- c("phenotype", "rate")
    }
    CInts <- list(CIphenotype, CIabsolute)
    names(CInts) <- c("phenotype", "rate")
    if (ConfInt == TRUE) {
      res <- list(rbiRES, PPtot, p.trend, p.rbi.slopeA,
                  p.trend.sel, p.rbi.slopeA.sel,
                  SMA.res, CInts)
      names(res) <- c("rbt", "pbt", "phenotypic.regression", "rate.regression",
                      "node.phenotypic.regression", "node.rate.regression",
                      "group.comparison", "ConfInts")
    }else {
      res <- list(rbiRES, PPtot, p.trend, p.rbi.slopeA,
                  p.trend.sel, p.rbi.slopeA.sel,
                  SMA.res)
      names(res) <- c("rbt", "pbt", "phenotypic.regression", "rate.regression",
                      "node.phenotypic.regression", "node.rate.regression",
                      "group.comparison")
    }
  }else {
    if (ConfInt == TRUE) {
      CInts <- list(CIphenotype, CIabsolute)
      names(CInts) <- c("phenotype", "abs.rate")
      res <- list(rbi, PPtot, p.trend, p.rbi.slopeA,
                  CInts)
      names(res) <- c("rbt", "pbt", "phenotypic.regression", "rate.regression",
                      "ConfInts")
    }else {
      res <- list(rbi, PPtot, p.trend, p.rbi.slopeA)
      names(res) <- c("rbt", "pbt", "phenotypic.regression", "rate.regression")
    }
  }

  if(length(which(sapply(res,function(x) is.null(x))))>0)
    res<-res[-which(sapply(res,function(x) is.null(x)))]
  return(res)
}
