#'@title Searching for evolutionary trends in phenotypes and rates
#'@description This function searches for evolutionary trends in the phenotypic
#'  mean and the evolutionary rates for the entire tree and individual clades.
#'@usage search.trend(RR,y,x1=NULL,x1.residuals = FALSE,
#'  node=NULL,cov=NULL,nsim=100,clus=0.5)
#'@param RR an object produced by \code{\link{RRphylo}}.
#'@param y the named vector (or matrix if multivariate) of phenotypes.
#'@param x1 the additional predictor to be specified if the RR object has been
#'  created using an additional predictor (i.e. multiple version of
#'  \code{\link{RRphylo}}). \code{'x1'} vector must be as long as the number of nodes
#'  plus the number of tips of the tree, which can be obtained by running
#'  \code{\link{RRphylo}} on the predictor as well, and taking the vector of ancestral
#'  states and tip values to form the \code{x1}. Note: only one predictor at
#'  once can be specified.
#'@param x1.residuals logical specifying whether the residuals of regression
#'  between \code{y} and \code{x1} should be inspected for a phenotypic trend
#'  (see details and examples below). Default is \code{FALSE}.
#'@param node the node number of individual clades to be specifically tested and
#'  contrasted to each other. It is \code{NULL} by default. Notice the node
#'  number must refer to the dichotomic version of the original tree, as
#'  produced by \code{\link{RRphylo}}.
#'@param cov the covariate values to be specified if the RR object has been
#'  created using a  covariate for rates calculation.  As for \code{\link{RRphylo}},
#'  \code{'cov'} must be as long as the number of nodes plus the number of tips
#'  of the tree, which can be obtained by running \code{\link{RRphylo}} on the
#'  covariate as well, and taking the vector of ancestral states and tip values
#'  to form the covariate (see the example below).
#'@param nsim number of simulations to be performed. It is set at 100 by
#'  default.
#'@param clus the proportion of clusters to be used in parallel computing. To
#'  run the single-threaded version of \code{search.trend} set \code{clus} = 0.
#'@return The function returns a list object containing:
#'@return \strong{$trends.data} a 'RRphyloList' object including:
#'  \enumerate{\item{\code{$phenotypeVStime}}: a data frame of phenotypic values
#'  (or \code{y} versus \code{x1} regression residuals if
#'  \code{x1.residuals=TRUE}) and their distance from the tree root for each
#'  node (i.e. ancestral states) and tip of the tree.
#'  \item{\code{$absrateVStime}}: a data frame of \code{\link{RRphylo}} rates and the
#'  distance from the tree root (age). If y is multivariate, it also includes
#'  the multiple rates for each y vector. If \code{node} is specified, each
#'  branch is classified as belonging to an indicated clade.
#'  \item{\code{$rescaledrateVStime}}: a data frame of rescaled \code{\link{RRphylo}}
#'  rates and the distance from the tree root (age). If y is multivariate, it
#'  also includes the multiple rates for each y vector. If \code{node} is
#'  specified, each branch is classified as belonging to an indicated clade. NAs
#'  correspond either to very small values or to outliers which are excluded
#'  from the analysis.}
#'@return \strong{$phenotypic.regression} results of phenotype (\code{y} versus
#'  \code{x1} regression residuals) versus age regression. It reports a p-value
#'  for the regression slope between the variables (p.real), a p-value computed
#'  contrasting the real slope to Brownian motion simulations (p.random), and a
#'  parameter indicating the deviation of the phenotypic mean from the root
#'  value in terms of the number of standard deviations of the trait
#'  distribution (dev). dev is 0 under Brownian Motion. Only p.random should be
#'  inspected to assess significance.
#'@return \strong{$rate.regression} results of the rates (rescaled absolute
#'  values) versus age regression. It reports a p-value for the regression
#'  between the variables (p.real), a p-value computed contrasting the real
#'  slope to Brownian motion simulations (p.random), and a parameter indicating
#'  the ratio between the range of phenotypic values and the range of such
#'  values halfway along the tree height, divided to the same figure under
#'  Brownian motion (spread). spread is 1 under Brownian Motion. Only p.random
#'  should be inspected to assess significance.
#'@return \strong{$ConfInts} a 'RRphyloList' object including the 95\%
#'  confidence intervals around regression slopes of phenotypes and rates (both
#'  rescaled and unscaled absolute rates) produced according to the Brownian
#'  motion model of evolution.
#'@return If specified, individual nodes are tested as the whole tree, the
#'  results are summarized in the objects:
#'@return \strong{$node.phenotypic.regression} results of phenotype (or \code{y}
#'  versus \code{x1} regression residuals) versus age regression through node.
#'  It reports the slope for the regression between the variables at node
#'  (slope), a p-value computed contrasting the real slope to Brownian motion
#'  simulations (p.random), the difference between estimated marginal means
#'  predictions for the group and for the rest of the tree (emm.difference), and
#'  a p-value for the emm.difference (p.emm).
#'@return \strong{$node.rate.regression} results of the rates (absolute values)
#'  versus age regression through node. It reports the difference between
#'  estimated marginal means predictions for the group and for the rest of the
#'  tree (emm.difference), a p-value for the emm.difference (p.emm), the
#'  regression slopes for the group (slope.node) and for the rest of the tree
#'  (slope.others), and a p-value for the difference between such slopes
#'  (p.slope).
#'@return If more than one node is specified, the object
#'  \strong{$group.comparison} reports the same results as
#'  $node.phenotypic.regression and $node.rate.regression obtained by comparing
#'  individual clades to each other.
#'@return The output always has an attribute "Call" which returns an unevaluated call to the function.
#'@author Silvia Castiglione, Carmela Serio, Pasquale Raia, Alessandro
#'  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'  Carotenuto
#'@details The function simultaneously returns the regression of phenotypes and
#'  phenotypic evolutionary rates against age tested against Brownian motion
#'  simulations to assess significance. To this aim rates are rescaled in the
#'  0-1 range and then logged. To assess significance, slopes are compared to a
#'  family of simulated slopes (BMslopes, where the number of simulations is
#'  equal to \code{nsim}), generated under the Brownian motion, using the
#'  \code{fastBM} function in the package \pkg{phytools}. Individual nodes are
#'  compared to the rest of the tree in different ways depending on whether
#'  phenotypes or rates (always unscaled in this case) versus age regressions
#'  are tested. With the former, the regression slopes for individual clades and
#'  the slope difference between clades is contrasted to slopes obtained through
#'  Brownian motion simulations. For the latter, regression models are tested
#'  and contrasted to each other referring to estimated marginal means, by using
#'  the \code{emmeans} function in the package \pkg{emmeans}.
#'
#'  The \href{../doc/RRphylo.html#predictor}{multiple regression version of
#'  RRphylo} allows to incorporate the effect of an additional predictor in the
#'  computation of evolutionary rates without altering the ancestral character
#'  estimation. Thus, when a multiple \code{\link{RRphylo}} output is fed to
#'  \code{search.trend}, the predictor effect is accounted for on the absolute
#'  evolutionary rates, but not on the phenotype. However, in some situations
#'  the user might want to factor out the predictor effect on phenotypes as
#'  well. Under the latter circumstance, by setting the argument
#'  \code{x1.residuals = TRUE}, the residuals of the response to predictor
#'  regression are used as to represent the phenotype.
#'@importFrom stats as.formula coef resid density predict cor
#'@importFrom phytools nodeHeights
#'@importFrom parallel makeCluster detectCores stopCluster clusterEvalQ clusterExport parLapply
#'@importFrom doParallel registerDoParallel
#'@importFrom foreach foreach %dopar%
#'@importFrom utils combn
#'@importFrom emmeans emmeans emtrends
#'@export
#'@seealso \href{../doc/search.trend.html}{\code{search.trend} vignette}
#'@seealso \code{\link{overfitST}}; \href{../doc/overfit.html#overfitST}{\code{overfitST} vignette}
#'@seealso \code{\link{plotTrend}}; \href{../doc/Plotting-tools.html#plotTrend}{\code{plotTrend} vignette}
#'@references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M.,
#'  Profico, A., Girardi, G., & Raia, P. (2019) Simultaneous detection of
#'  macroevolutionary patterns in phenotypic means and rate of change with and
#'  within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14:
#'  e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @examples
#'  \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' cc<- 2/parallel::detectCores()
#'
#' # Extract Pterosaurs tree and data
#' library(ape)
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' # Case 1. "RRphylo" whitout accounting for the effect of a covariate
#' RRphylo(tree=treeptero,y=log(massptero),clus=cc)->RRptero
#'
#' # Case 1.1. "search.trend" whitout indicating nodes to be tested for trends
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, clus=cc,cov=NULL,node=NULL)->st1
#'
#' # Case 1.2. "search.trend" indicating nodes to be specifically tested for trends
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, node=143, clus=cc,cov=NULL)->st2
#'
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#' # "RRphylo" on the covariate in order to retrieve ancestral state values
#' c(RRptero$aces,log(massptero))->cov.values
#' names(cov.values)<-c(rownames(RRptero$aces),names(massptero))
#' RRphylo(tree=treeptero,y=log(massptero),cov=cov.values,clus=cc)->RRpteroCov
#'
#' # Case 2.1. "search.trend" whitout indicating nodes to be tested for trends
#' search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, clus=cc,cov=cov.values)->st3
#'
#' # Case 2.2. "search.trend" indicating nodes to be specifically tested for trends
#' search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, node=143, clus=cc,cov=cov.values)->st4
#'
#'
#' # Case 3. "search.trend" on multiple "RRphylo"
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#' DataCetaceans$aceMyst->aceMyst
#'
#' drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#'
#' # incorporating the effect of body size at inspecting trends in absolute evolutionary rates
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc)->STcet
#'
#' # incorporating the effect of body size at inspecting trends in both absolute evolutionary
#' # rates and phenotypic values (by using brain versus body mass regression residuals)
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,clus=cc)->st5
#'    }


search.trend<-function (RR,y,
                        x1=NULL,x1.residuals=FALSE,
                        node = NULL, cov = NULL,
                        nsim = 100, clus = 0.5)
{
  # require(ape)
  # require(phytools)
  # require(stats4)
  # require(foreach)
  # require(doParallel)
  # require(parallel)
  # require(nlme)
  # require(emmeans)
  # require(car)

  misspacks<-sapply(c("car","nlme"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  funcall <- match.call()

  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

  t <- RR$tree
  if(min(diag(vcv(t)))/max(diag(vcv(t)))>=0.9) stop("not enough fossil information")
  rates <- RR$rates
  betas <- RR$multiple.rates
  aceRR <- RR$aces
  L <- RR$tip.path
  L1 <- RR$node.path
  if (!is.null(ncol(y))&is.null(rownames(y))) stop("The matrix of phenotypes needs to be named")
  if (is.null(ncol(y))&is.null(names(y))) stop("The vector of phenotypes needs to be named")
  ynam<-deparse(substitute(y))
  # if(is.null(nrow(y))) y <- treedata(t, y, sort = TRUE)[[2]][,1] else y <- treedata(t, y, sort = TRUE)[[2]]
  y <- treedataMatch(t, y)[[1]]

  H <- max(nodeHeights(t))
  eds <- data.frame(leaf = c(Ntip(t)+1,t$edge[, 2]),height=c(0,nodeHeights(t)[, 2]))
  eds[which(eds[,1]<=Ntip(t)),1] <- t$tip.label[eds[which(eds[,1]<=Ntip(t)),1]]

  if (ncol(y) > 1) {
    y.multi <- (L %*% rates)
    y.multi[match(rownames(y),rownames(y.multi)),,drop=FALSE]->y.multi
    aceRR.multi <- (L1 %*% rates[1:Nnode(t), ])
    nodes <- cbind(aceRR,aceRR.multi)
    if(is.null(colnames(y))) colnames(nodes)<- c(paste("y", seq(1,ncol(y)),sep = ""),"y.multi") else
      colnames(nodes)<-c(colnames(y),"y.multi")
    P <- rbind(nodes, cbind(y,y.multi))
    if(isTRUE(x1.residuals)) apply(P,2,function(x) residuals(lm(x~x1)))->P

    rates <- as.data.frame(rates)
    betas <- as.data.frame(betas)

    rate.data <- data.frame(betas = betas[match(eds[, 1], rownames(betas)),],
                            rate = rates[match(eds[, 1], rownames(rates)),],
                            age = eds[, 2])
    if(is.null(colnames(y))) colnames(rate.data)[1:ncol(y)] <- paste("betas", seq(1,ncol(y)), sep = "") else
      colnames(rate.data)[1:ncol(y)]<-colnames(y)

  }else{
    P <- rbind(aceRR, y)
    if(isTRUE(x1.residuals)) as.matrix(residuals(lm(P~x1)))->P
    # colnames(P)<-"y"
    colnames(P)<-ynam

    rate.data <- data.frame(rate = rates[match(eds[, 1], rownames(rates)),],
                            age = eds[, 2])
    rownames(rate.data) <- rownames(rates)[match(eds[, 1], rownames(rates))]

  }

  rate.data[,ncol(rate.data)]+L[1,1]->rate.data[,ncol(rate.data)]
  phen.data <- data.frame(P[match(rownames(rate.data), rownames(P)),,drop=FALSE],age=rate.data$age,check.names = FALSE)
  rate.dataRES<-rate.data
  phen.dataRES <- phen.data
  {##### Rate Trend Real Multi #####
    rate.reg <- scalrat.data<-rate.coef <-list()
    e1<-array()
    for (i in 1:(ncol(rate.data)-1)) {
      bet <- rate.data[, i,drop=FALSE]
      age <- rate.data[, ncol(rate.data),drop=FALSE]
      abs(bet)->rts->rtsA
      log(range01(rts))->rts

      c(which(rts=="-Inf"),which(age=="-Inf"))->outs
      if(length(outs)>0){
        rts[-outs,,drop=FALSE]->rts
        age[-outs,,drop=FALSE]->age
      }
      age->ageC

      sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->e1[i]

      if(!is.null(x1)){
        rts[-1,,drop=FALSE]->rts
        age[-1,,drop=FALSE]->age

        car::outlierTest(lm(as.matrix(rts)~as.matrix(age)))->ouT
        if(any(ouT$bonf.p<=0.05)){
          rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
          age[-match(names(which(ouT$bonf.p<=0.05)),rownames(age)),,drop=FALSE]->age
        }
      }else{
        lm(as.matrix(rts)~as.matrix(age))->bb

        residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
        if((Ntip(t)+1)%in%names(resout)){
          rts[-1,,drop=FALSE]->rts
          age[-1,,drop=FALSE]->age
        }
      }
      lm(as.matrix(rts)~as.matrix(age))->regr.1
      rate.coef[[i]] <- coef(summary(regr.1))[2, c(1, 4)]
      rate.reg[[i]] <- regr.1

      scalrat.data[[i]]<-rts[match(rownames(rate.data),rownames(rts)),,drop=FALSE]
      rownames(scalrat.data[[i]])<-rownames(rate.data)
    }
    do.call(rbind,rate.coef)->rate.coef
    colnames(rate.coef) <- c("slope","p-value")
    do.call(cbind,scalrat.data)->scalrat.data
    colnames(scalrat.data)<-colnames(rate.data)[1:(ncol(rate.data)-1)]
    data.frame(scalrat.data,age=rate.data$age)->scalrat.data

    rownames(rate.coef) <- names(rate.reg)<-colnames(rate.data)[1:(ncol(rate.data)-1)]

  }

  {##### Phenotypic Trend Real Multi #####
    phen.reg <-apply(phen.data[,1:(ncol(phen.data)- 1),drop=FALSE], 2, function(x) lm(range01(x) ~ phen.data[,ncol(phen.data)]))
    phen.coef <- lapply(phen.reg, function(x) coefficients(summary(x))[2,c(1,4)])
    # if(ncol(y)>1){
    #   if(is.null(colnames(y)))
    #     names(phen.coef) <- c(paste("y", seq(1:ncol(y)), sep = ""),"y.multiple") else
    #       names(phen.coef) <- c(colnames(y),"y.multiple")
    # }else names(phen.coef)<-"y"

    names(phen.coef) <- colnames(P)

    sapply(1:length(phen.coef),function(i){
      if(phen.coef[[i]][1]<0)
        (min(predict(phen.reg[[i]]))-mean(range01(phen.data[,i])))/sd(range01(phen.data[,i]))
      else (max(predict(phen.reg[[i]]))-mean(range01(phen.data[,i])))/sd(range01(phen.data[,i]))
    })->dev
    names(dev)<-names(phen.coef)
    do.call(rbind,phen.coef)->phen.coef
  }


  #### Nodes Real ####
  if (!is.null(node)) {
    rate.dataN <- rate.data
    phen.dataN <- phen.data
    rate.dataN$group <- rep("others", nrow(rate.dataN))
    phen.reg.sel<-phen.reg.age.sel<-phen.reg.y.sel<-list()
    rate.reg.y.sel<-rate.reg.age.sel<-list()
    for (j in 1:length(node)) {
      n <- node[j]
      sele <- getDescendants(t, n)
      sele[which(sele <= Ntip(t))]<-t$tip.label[sele[which(sele<=Ntip(t))]]
      rate.dataN[match(c(n,sele), rownames(rate.dataN)), ]$group <- paste("g",
                                                                          n, sep = "")
      rep("others",nrow(rate.dataN))->gg
      gg[match(c(n,sele), rownames(rate.dataN))]<-paste("g",n, sep = "")
      data.frame(rate.dataN,gg)->rate.dataN
      rate.data.sel <- rate.data[match(c(n,sele), rownames(rate.data)),,drop=FALSE]

      {#### Rate Trend Real Node Multi ####
        rate.reg.age <- rate.data.sel[,ncol(rate.data.sel),drop=FALSE]

        lapply(1:(ncol(rate.data.sel)-1),function(i){
          bet <- rate.data.sel[,i,drop=FALSE]
          age <- rate.data.sel[,ncol(rate.data.sel),drop=FALSE]
          abs(bet)->rts
          predict(lm(as.matrix(rts)~as.matrix(age)))->pp
          names(pp)<-rownames(rts)
          pp
        })->rate.reg.y

        # names(rate.reg.y) <- colnames(rate.data)[1:(ncol(rate.data)-1)]
        names(rate.reg.y) <- rownames(rate.coef)
      }
      rate.reg.y.sel[[j]] <- do.call(cbind,rate.reg.y)
      rate.reg.age.sel[[j]] <- rate.reg.age

      {##### Phenotypic Trend Real Node Multi #####
        PPsel <- phen.data[match(rownames(rate.data.sel), rownames(phen.data)),,drop=FALSE]
        phen.regN <- apply(PPsel[1:(ncol(PPsel)-1)],2, function(x)
          coefficients(summary(lm(range01(x) ~ PPsel[, ncol(PPsel)])))[2,c(1,4)])

        # if(ncol(y)>1){
        #   if(is.null(colnames(y)))
        #     colnames(phen.regN) <- c(paste("y", seq(1,dim(y)[2]), sep = ""),"y.multiple") else
        #       colnames(phen.regN) <- c(colnames(y),"y.multiple")
        # }else colnames(phen.regN)<-"y"

        phen.reg.y.sel[[j]] <- apply(PPsel[1:(ncol(PPsel)-1)], 2,function(x)
          predict(lm(x ~ PPsel[, ncol(PPsel)])))
        phen.reg.age.sel[[j]] <- PPsel$age
      }
      phen.reg.sel[[j]] <- phen.regN
    }

    names(phen.reg.sel)<-node

    #colnames(rate.dataN)[4:ncol(rate.dataN)]<-paste("group",seq(1,(ncol(rate.dataN)-3)),sep="")
    phen.dataN <- cbind(phen.dataN, rate.dataN[,(ncol(rate.data)+1):ncol(rate.dataN)])
    if (length(which(rate.dataN$group == "others")) < 3)
      rate.dataN <- rate.dataN[which(rate.dataN$group != "others"),]


    rate.dataRES<-rate.dataN[,1:(ncol(rate.dataN)-length(node))]

    { #### Node Comparison ####
      rate.NvsN<-rate.NvsO<-list()
      phen.NvsN<-phen.NvsN.emm<-phen.NvsO<-list()

      groupPP <- phen.dataN[which(phen.dataN$group != "others"), ]$group
      agePP <- phen.dataN$age

      for (i in 1:(ncol(rate.data)-1)) {
        bets <- rate.dataN[, i,drop=FALSE]
        yPP<-phen.dataN[,i,drop=FALSE]

        rate.comp<-list()
        phen.comp<-list()
        for (w in 1:(length(node))){
          data.frame(rate=unname(bets),age=rate.dataN$age,group=rate.dataN[,(w+ncol(rate.data)+1)])->emdat
          suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(rate)~age+group,data=emdat),specs="group"))))
          # mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(rate)~age*group,data=emdat),specs="group",var="age")))
          # mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
          # data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->rate.comp[[w]]

          emmeans::emtrends(lm(abs(rate)~age*group,data=emdat),specs="group",var="age")->mtrends.test
          mtrends<-as.data.frame(pairs(mtrends.test))
          mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
          data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],
                     t(as.data.frame(mtrends.test)[,2,drop=FALSE]),mtrends[,6])->rate.comp[[w]]

          data.frame(yPP,age=phen.dataN$age,group=phen.dataN[,(w+ncol(rate.data)+1)])->PPemdat
          if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
            data.frame(PPemdat,x1=x1[match(rownames(PPemdat),rownames(x1)),])->PPemdat
            suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP[,1])~age+x1+group,data=PPemdat),specs="group"))))
          }else{
            suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP[,1])~age+group,data=PPemdat),specs="group"))))
          }
          data.frame(do.call(rbind,strsplit(as.character(PPmeans[,1])," - ")),PPmeans[,-c(1,3,4,5)])->phen.comp[[w]]

        }
        do.call(rbind,rate.comp)->rate.NvsO[[i]]
        # colnames(rate.NvsO[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
        colnames(rate.NvsO[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.node","slope.others","p.slope")

        do.call(rbind,phen.comp)->phen.comp
        colnames(phen.comp)<-c("group_1","group_2","mean","p.mean")
        sapply(strsplit(as.character(phen.comp[,1]),"g"),"[[",2)->PPnam
        phen.comp[,3:4]->phen.comp
        rownames(phen.comp)<-PPnam
        phen.comp[match(node,rownames(phen.comp)),]->phen.comp
        phen.comp->phen.NvsO[[i]]

        dat <- data.frame(bets=unname(bets), age=rate.dataN$age, group=rate.dataN$group)
        suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(bets)~age+group,data=dat),specs="group"))))
        # mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(bets)~age*group,data=dat),specs="group",var="age")))
        # mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
        # data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->rate.NvsN[[i]]
        # colnames(rate.NvsN[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
        # rate.NvsN[[i]][-which(rate.NvsN[[i]]$group_2=="others"),]->rate.NvsN[[i]]

        emmeans::emtrends(lm(abs(bets)~age*group,data=dat),specs="group",var="age")->mtrends.test
        as.data.frame(mtrends.test)[,1:2]->mtrends.slope
        mtrends<-as.data.frame(pairs(mtrends.test))
        mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
        data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)])->mtrends.dat
        data.frame(mtrends.dat,mtrends.slope[match(mtrends.dat[,1],mtrends.slope[,1]),2],
                   mtrends.slope[match(mtrends.dat[,2],mtrends.slope[,1]),2],mtrends[,6])->rate.NvsN[[i]]
        colnames(rate.NvsN[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.group1","slope.group2","p.slope")
        rate.NvsN[[i]][which(rate.NvsN[[i]]$group_2!="others"),]->rate.NvsN[[i]]

        dat <- data.frame(yPP=unname(yPP), age=phen.dataN$age, group=phen.dataN$group)
        if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
          data.frame(dat,x1=x1[match(rownames(dat),rownames(x1)),])->dat
          suppressMessages(PPpairs<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+x1+group,data=dat),specs="group"))))
        }else{
          suppressMessages(PPpairs<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+group,data=dat),specs="group"))))
        }

        data.frame(do.call(rbind,strsplit(as.character(PPpairs[,1])," - ")),
                   PPpairs[,-c(1,3,4,5)])->Pcomp.emm
        colnames(Pcomp.emm)<-c("group_1","group_2","mean","p.mean")
        Pcomp.emm[which(Pcomp.emm$group_2!="others"),]->Pcomp.emm
        Pcomp.emm->phen.NvsN.emm[[i]]

        if(length(node)>1){
          sapply(phen.reg.sel,function(w) w[1,i])->slope.tot
          names(slope.tot)<-paste("g",names(phen.reg.sel),sep="")

          combn(sort(unique(as.character(groupPP))),2)->pair

          # sapply(1:ncol(pair), function(jj){
          #   slope.tot[match(pair[1,jj],names(slope.tot))]-
          #     slope.tot[match(pair[2,jj],names(slope.tot))]
          # })->slope.diff
          # data.frame(t(pair),slope.diff)->phen.NvsN[[i]]
          # colnames(phen.NvsN[[i]])<-c(colnames(rate.NvsN[[i]])[1:2],"estimate")

          do.call(rbind,lapply(1:ncol(pair), function(jj){
            data.frame(slope.tot[match(pair[1,jj],names(slope.tot))],
                       slope.tot[match(pair[2,jj],names(slope.tot))])
          }))->slopes

          data.frame(t(pair),slopes)->phen.NvsN[[i]]
          colnames(phen.NvsN[[i]])<-c(colnames(rate.NvsN[[i]])[1:2],"slope.group_1","slope.group_2")

        }else{
          phen.NvsN<-NULL
          phen.NvsN.emm<-NULL
        }
      }

      if(!is.null(phen.NvsN.emm)) names(phen.NvsN.emm)<-colnames(phen.dataN)[1:ncol(rate.data)-1]

      sapply(strsplit(as.character(rate.NvsO[[1]][,1]),"g"),"[[",2)->grn
      lapply(rate.NvsO,function(x) x[match(node,grn),])->rate.NvsO

      lapply(1:length(node),function(u){
        # do.call(rbind,lapply(rate.NvsO,function(x) x[u,3:6]))->rcs
        do.call(rbind,lapply(rate.NvsO,function(x) x[u,3:7]))->rcs
        rownames(rcs)<-rownames(rate.coef)
        rcs
      })->rate.coef.sel

    }

    names(rate.coef.sel)<-node
    names(rate.NvsO) <- colnames(rate.dataN)[1:(ncol(rate.data)-1)]

    if(!is.null(phen.NvsN)){
      names(rate.NvsN)<-names(rate.NvsO)
      names(phen.NvsN.emm) <- names(phen.NvsN) <- names(phen.reg)
    }
  }


  #### Random ####
  RR$aces[1,,drop=FALSE]->a
  if(!is.null(x1)) x1[1:Ntip(t)]->y1
  {
    yyD <- list()
    yyT <- list()

    if(is.null(x1)){
      for (i in 1:ncol(y)) {
        yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = a[i], bounds=c(min(y[,i]),max(y[,i])))))
        yyT[[i]] <- replicate(nsim,fastBM(t, sig2 = 1, a = mean(y[,i]), bounds=c(min(y[,i]),max(y[,i]))))
      }
    }else{
      if(isFALSE(x1.residuals)){
        for (i in 1:ncol(y))
          yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = a[i], bounds=c(min(y[,i]),max(y[,i])))))
      }else{
        phen.data[which(rownames(phen.data)==(Ntip(t)+1)),]->ares
        phen.data[which(rownames(phen.data)%in%t$tip.label),]->yres
        for (i in 1:ncol(y))
          yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = ares[1,i], bounds=c(min(yres[,i]),max(yres[,i])))))
      }

      matrix(ncol=nsim,nrow=Ntip(t))->yy1T
      matrix(ncol=nsim,nrow=Nnode(t))->ace1T
      for(k in 1:nsim){
        phenb<-list()
        for(i in 1:ncol(y)) fastBM(t, sig2 = 1, a = mean(y[,i]), bounds=c(min(y[,i]),max(y[,i])),internal = TRUE)->phenb[[i]]
        fastBM(t,a=mean(y1),bounds=c(min(y1),max(y1)),internal = TRUE)->phen1b
        do.call(cbind,phenb)->phenb

        cor(cbind(y,y1))->cr
        cbind(phenb,phen1b)->xx1
        cr->m
        t(chol(m))->U
        U%*%t(xx1)->xx3
        xx3[-nrow(xx3),1:Ntip(t),drop=FALSE]->yyT[[k]]
        if(ncol(yyT[[k]])>nrow(yyT[[k]])) t(yyT[[k]])->yyT[[k]]
        xx3[nrow(xx3),1:Ntip(t)]->yy1T[,k]
        xx3[nrow(xx3),(Ntip(t)+1):(Ntip(t)+Nnode(t))]->ace1T[,k]
      }
      rownames(yy1T)<-rownames(y)
      rownames(ace1T)<-rownames(aceRR)
    }
  }
  #### End of Random ####
  core.chunk<-  expression({
    {
      yD <- lapply(yyD, function(x) x[,ii])
      yD <- do.call(cbind, yD)

      if(is.null(x1)){
        yT <- lapply(yyT, function(x) x[,ii])
        yT <- do.call(cbind, yT)
      }else{
        yT <-yyT[[ii]]
        y1T<-yy1T[,ii]
        aceT<-ace1T[,ii]
        cbind(L,y1T-aceT[1])->LX
        cbind(L1,aceT-mean(aceT))->L1X
      }

      betasT<- matrix(ncol=ncol(yT), nrow=Ntip(t)+Nnode(t))
      aceRRT<- matrix(ncol=ncol(yT), nrow=Nnode(t))
      betasD<- matrix(ncol=ncol(yD), nrow=Ntip(t)+Nnode(t))
      aceRRD<- matrix(ncol=ncol(yD), nrow=Nnode(t))
      for (s in 1:ncol(yT)){
        rootV<-a[s]
        # lambda <- RR$lambda[[s]]@coef
        lambda <- RR$lambda[[s]]$par
        if(!is.null(x1)){

          m.betasT <- (solve(t(LX) %*% LX + lambda * diag(ncol(LX))) %*%
                         t(LX)) %*% (as.matrix(yT[,s]) - rootV)
          aceRRT[,s] <- (L1X %*% m.betasT[c(1:Nnode(t),length(m.betasT)), ]) + rootV
          as.matrix(m.betasT[-length(m.betasT),])->betasT[,s]

        }else{
          m.betasT <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                         t(L)) %*% (as.matrix(yT[,s]) - rootV)

          aceRRT[,s] <- (L1 %*% m.betasT[1:Nnode(t), ]) + rootV
          m.betasT->betasT[,s]
        }

        if(isTRUE(x1.residuals)) rootVD<-mean(range(yres[,s])) else rootVD<-mean(range(y[,s])) ### CONTROLLARE CON x1
        m.betasD <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                       t(L)) %*% (as.matrix(yD[,s]) - rootVD)
        aceRRD[,s] <- (L1 %*% m.betasD[1:Nnode(t), ]) + rootVD
        m.betasD->betasD[,s]


      }
      rownames(betasT)<-rownames(betasD)<-rownames(RR$rates)
      rownames(aceRRT)<-rownames(aceRRD)<-rownames(RR$aces)

    }

    betasT->betasTreal
    #### Covariate ####
    if (!is.null(cov)) {
      cov[match(rownames(betas),names(cov))]->cov
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
      ratesT<-as.matrix(as.data.frame(apply(betasT, 1, function(x) sqrt(sum(x^2)))))
    }
    #### End of Covariate ####

    if (ncol(y)>1) {
      betasTreal->betasT
      ratesT<-as.matrix(as.data.frame(apply(betasT, 1, function(x) sqrt(sum(x^2)))))
      aceRRTmulti <- (L1 %*% ratesT[1:Nnode(t), ])
      yTmulti <- (L %*% ratesT)
      rate.data <- data.frame(betas = betasT[match(eds[, 1],rownames(betasT)), ],
                              rate = ratesT[match(eds[,1], rownames(ratesT)), ],
                              age = eds[, 2])
      colnames(rate.data)[1:ncol(yT)] <- paste("betas", seq(1,ncol(yT), 1), sep = "")

      ratesD<-as.matrix(as.data.frame(apply(betasD, 1, function(x) sqrt(sum(x^2)))))
      aceRRDmulti <- (L1 %*% ratesD[1:Nnode(t), ])
      yDmulti <- (L %*% ratesD)+aceRRDmulti[1,]

      nodesD<-cbind(aceRRD,aceRRDmulti)
      colnames(nodesD) <- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
      P <- rbind(nodesD, cbind(yD,yDmulti))
      if(isTRUE(x1.residuals)) apply(P,2,function(x) residuals(lm(x~x1)))->P

    }else {
      rate.data <- data.frame(rate = betasT[match(eds[, 1],rownames(betasT)), ],
                              age = eds[, 2])

      P <- rbind(aceRRD, yD)
      colnames(P)<-"y"
      if(isTRUE(x1.residuals)) as.matrix(residuals(lm(P~x1)))->P
    }
    rate.data[,ncol(rate.data)]+L[1,1]->rate.data[,ncol(rate.data)]
    phen.data <- data.frame(P[match(rownames(rate.data), rownames(P)), ],age=rate.data$age)

    {##### Rate Trend Random Multi #####
      rate.reg <- scalrat.data<-rate.coef <-list()
      ee<-array()
      for (i in 1:(ncol(rate.data)-1)) {
        bet <- rate.data[, i,drop=FALSE]
        age <- rate.data[, ncol(rate.data),drop=FALSE]
        abs(bet)->rts->rtsA
        log(range01(rts))->rts

        c(which(rts=="-Inf"),which(age=="-Inf"))->outs
        if(length(outs)>0){
          rts[-outs,,drop=FALSE]->rts
          age[-outs,,drop=FALSE]->age
        }
        age->ageC

        sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->ee[i]

        if(!is.null(x1)){
          rts[-1,,drop=FALSE]->rts
          age[-1,,drop=FALSE]->age

          car::outlierTest(lm(as.matrix(rts)~as.matrix(age)))->ouT
          if(any(ouT$bonf.p<=0.05)){
            rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
            age[-match(names(which(ouT$bonf.p<=0.05)),rownames(age)),,drop=FALSE]->age
          }
        }else{
          lm(as.matrix(rts)~as.matrix(age))->bb

          residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
          if((Ntip(t)+1)%in%names(resout)){
            rts[-1,,drop=FALSE]->rts
            age[-1,,drop=FALSE]->age
          }
        }
        lm(as.matrix(rts)~as.matrix(age))->regr.1
        rate.coef[[i]] <- coef(summary(regr.1))[2, c(1, 4)]
        rate.reg[[i]] <- regr.1

        scalrat.data[[i]]<-rts[match(rownames(rate.data),rownames(rts)),,drop=FALSE]
      }
      do.call(rbind,rate.coef)->rate.coef
      colnames(rate.coef) <- c("slope","p-value")
      do.call(cbind,scalrat.data)->scalrat.data
      colnames(scalrat.data)<-colnames(rate.data)[1:(ncol(rate.data)-1)]
      data.frame(scalrat.data,age=rate.data$age)->scalrat.data
      rownames(scalrat.data)<-rownames(rate.data)

      if(is.null(colnames(y))|ncol(y)==1) rownames(rate.coef) <- names(rate.reg)<-colnames(rate.data)[1:(ncol(rate.data)-1)] else
        rownames(rate.coef) <- names(rate.reg)<-c(colnames(y),"y.multiple")
    }


    { #### Phenotypic Trend Random Multi ####
      phen.reg <- apply(phen.data[1:(dim(phen.data)[2] - 1)], 2, function(x)
        summary(lm(range01(x) ~phen.data[, dim(phen.data)[2]])))
      phen.coef <- lapply(phen.reg, coefficients)
      if(ncol(y)>1){
        if(is.null(colnames(y)))
          names(phen.coef) <- c(paste("y", seq(1:ncol(y)), sep = ""),"y.multiple") else
            names(phen.coef) <- c(colnames(y),"y.multiple")
          cbind(yT,yTmulti)->yT
          cbind(yD,yDmulti)->yD
      }else names(phen.coef)<-"y"

    }


    #### Nodes Random ####
    if (!is.null(node)) {
      phen.dataN <- phen.data
      phen.dataN$group <- rep("NA", nrow(phen.dataN))
      phen.reg.sel <- list()
      for (j in 1:length(node)) {
        n <- node[j]
        sele <- getDescendants(t, n)
        sele[which(sele <= Ntip(t))]<-t$tip.label[sele[which(sele<=Ntip(t))]]
        phen.dataN[match(c(n,sele), rownames(phen.dataN)), ]$group <- paste("g",
                                                                            n, sep = "")
        {#### Phenotypic Trend Random Node Multi ####
          PPsel <- phen.data[match(c(n,sele), rownames(phen.data)),,drop=FALSE]
          phen.regN <- apply(PPsel[1:(ncol(PPsel)- 1)],
                             2, function(x) coefficients(summary(lm(range01(x) ~ PPsel[, ncol(PPsel)])))[2,c(1,4)])
        }
        phen.reg.sel[[j]] <- phen.regN
      }
      names(phen.reg.sel) <- node


      if(!is.null(phen.NvsN)){
        phen.dataN <- phen.dataN[-which(phen.dataN$group == "NA"),]

        {#### Node Comparison Random ####
          groupPP <- phen.dataN[, (ncol(phen.dataN))]
          agePP <- phen.dataN[, (ncol(phen.dataN)-1)]

          phen.NvsN.ran<-list()
          for (k in 1:(ncol(rate.data)-1)) {

            sapply(phen.reg.sel,function(w) w[1,k])->slope.tot
            names(slope.tot)<-paste("g",names(phen.reg.sel),sep="")


            combn(sort(unique(as.character(groupPP))),2)->pair
            sapply(1:ncol(pair), function(jj){
              slope.tot[match(pair[1,jj],names(slope.tot))]-
                slope.tot[match(pair[2,jj],names(slope.tot))]
            })->slope.diff

            data.frame(t(pair),slope.diff)->phen.NvsN.ran[[k]]
            colnames(phen.NvsN.ran[[k]])<-c("group_1","group_2","estimate")

          }
          names(phen.NvsN.ran)<-colnames(phen.dataN)[1:(dim(phen.dataN)[2] - 2)]
        }
      }else{
        phen.NvsN.ran<-NULL
      }
      list(phen.coef, rate.coef,
           phen.data, rate.data,yT,yD,ee,scalrat.data,
           phen.reg.sel,phen.NvsN.ran)
    }else {
      list(phen.coef, rate.coef,
           phen.data, rate.data,yT,yD,ee,scalrat.data)
    }
  })

  cldef<-({
    'cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
     clusterEvalQ(cl, {library(ape)\nlibrary(phytools)\nlibrary(nlme)\nlibrary(car)})
     clusterExport(cl=cl,varlist = list("cov"),envir = environment())
    '
  })


  ii=NULL
  res=NULL
  if(round((detectCores() * clus), 0)>1)
    eval(parse(text=paste0(cldef,'\nres<-parLapply(cl=cl,1:nsim,function(ii)',
                           core.chunk,")\n stopCluster(cl)"))) else
                             eval(parse(text=paste0('res<-lapply(1:nsim,function(ii)',core.chunk,")")))

  {#### p Rate Trend Multi####
    p.rate<-spread<-array()
    rate.coefR<-list()
    scalrat.CI <- CIabsolute <- list()
    if(ncol(y)>1) cbind(y,y.multi)->yTot else y->yTot
    for (i in 1:(ncol(rate.data)-1)){
      scalrat.CI[[i]]<-rate.coefR <- do.call(rbind, lapply(lapply(res,"[[", 2), function(x) x[i, 1]))
      # quantile(rate.coefR,c(0.025,0.975))
      yRT<- lapply(lapply(res,"[[",5), function(x) x[, i])

      rank(c(rate.coef[i, 1],rate.coefR[-1]))[1]/nsim->pp->p.rate[i]
      sapply(lapply(res,"[[",7),"[[",i)->ee

      if(i<=ncol(y)){
        ifelse(rate.coef[i, 1]>0,{
          if(pp<=0.025) dir.rate<-"increase slower" else
            if(pp>=0.975) dir.rate<-"increase faster" else
              dir.rate<-"nosig"
        },{
          if(pp>=0.975) dir.rate<-"decrease slower" else
            if(pp<=0.025) dir.rate<-"decrease faster" else
              dir.rate<-"nosig"
        })

        if(dir.rate!="nosig"){
          if(ncol(y)>1) yprint<-paste(" for",colnames(phen.data)[i]) else yprint<-NULL
          print(paste("Absolute evolutionary rates ",
                      dir.rate," than BM expectation",yprint,sep=""))
        }

      }

      (diff(range(yTot[,i]))/diff(range(yTot[,i][which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))*
        mean(sapply(yRT,function(x)(diff(range(x))/diff(range(x[which(diag(vcv(t))<diff(range(diag(vcv(t)))))])))))->spread[i]

      as.matrix(sapply(lapply(res,"[[", 4),function(k) unname(coef(lm(abs(k[,i])~k$age))[2])))->CIabsolute[[i]]

    }
    names(spread)<-names(p.rate) <- rownames(rate.coef)
    p.rate <- data.frame(slope = rate.coef[, 1],p.real = rate.coef[, 2],
                         p.random = p.rate,spread=spread)
  }

  { #### p Phenotypic Trend Multi and pdf ####
    CIphenotype<-list()
    p.phen <- array()
    for (i in 1:(ncol(phen.data)-1)) {
      CIphenotype[[i]]<-sapply(lapply(lapply(res,"[[", 1), "[[", i), function(x) x[2,1])
      p.phen[i] <- rank(c(phen.coef[i,1],CIphenotype[[i]][-1]))[1]/nsim

      if(i<=ncol(y)){
        ifelse(phen.coef[i,1]>0,{
          if(p.phen[i]<=0.025) dir.phen<-"increases less" else
            if(p.phen[i]>=0.975) dir.phen<-"increases more" else
              dir.phen<-"nosig"
        },{
          if(p.phen[i]<=0.025) dir.phen<-"decreases more" else
            if(p.phen[i]>=0.975) dir.phen<-"decreases less"else
              dir.phen<-"nosig"
        })

        if(dir.phen!="nosig"){
          if(ncol(y)>1) yprint<-paste(" for",colnames(phen.data)[i]) else yprint<-NULL
          print(paste("Phenotypic regression slope ",dir.phen ," than BM expectation",yprint,sep=""))
        }
      }
    }

    p.phen <- cbind(phen.coef, p.phen,dev)
    colnames(p.phen) <- c("slope", "p.real",
                          "p.random","dev")

  }

  if (!is.null(node)) {
    p.phen.sel <- list()
    for (u in 1:length(node)) {
      {#### p Phenotypic Trend Node Multi ####
        p.sele <- list()
        for (i in 1:(ncol(phen.data)-1)) {
          slopeR <- sapply(lapply(lapply(res,"[[", 9), "[[", u),function(a) a[1,i])
          p.sel <- rank(c(phen.reg.sel[[u]][1,i],slopeR[-1]))[1]/nsim

          p.sele[[i]] <- c(slope = phen.reg.sel[[u]][1,i],
                           p.slope = p.sel,
                           emm.difference=phen.NvsO[[i]][match(names(phen.reg.sel)[u],rownames(phen.NvsO[[i]])),1],
                           p.emm=phen.NvsO[[i]][match(names(phen.reg.sel)[u],rownames(phen.NvsO[[i]])),2])

        }
        names(p.sele) <- colnames(phen.reg.sel[[u]])
        p.phen.sel[[u]]  <- do.call(rbind, p.sele)
      }
    }
    names(p.phen.sel) <- names(phen.reg.sel)


    p.rate.sel<-rate.coef.sel
    p.rate.comp <- rate.NvsN
    if(!is.null(phen.NvsN)){
      {#### p Phenotypic Slope Comparison Multi  ####
        p.phen.comp<-list()
        for (k in 1:(ncol(phen.data)-1)) {
          p.nodeR<-mean.diff<-array()
          for(i in 1:nrow(phen.NvsN[[k]])){
            rank(c(phen.NvsN[[k]][i,3]-phen.NvsN[[k]][i,4],sapply(lapply(lapply(res,"[[",10),"[[",k),function(x) x[i,3])[-1]))[1]/nsim->p.nodeR[i]
            if((phen.NvsN[[k]][i,3]-phen.NvsN[[k]][i,4])>0&p.nodeR[i]<=0.05|(phen.NvsN[[k]][i,3]-phen.NvsN[[k]][i,4])<0&p.nodeR[i]>=0.95) p.nodeR[i]<-0.5
            if(as.character(interaction(phen.NvsN[[k]][i,c(1,2)]))==as.character(interaction(phen.NvsN.emm[[k]][i,c(1,2)])))
              phen.NvsN.emm[[k]][i,3]->mean.diff[i] else
                -1*phen.NvsN.emm[[k]][i,3]->mean.diff[i]
          }
          # data.frame(phen.NvsN[[k]][,c(1,2)],slope.difference=phen.NvsN[[k]][,3],p.slope=1-p.nodeR,
          #            emm.difference=mean.diff,p.emm=phen.NvsN.emm[[k]][,4])->p.phen.comp[[k]]
          data.frame(phen.NvsN[[k]][,c(1,2)],slope.group_1=phen.NvsN[[k]][,3],slope.group_2=phen.NvsN[[k]][,4],p.slope=p.nodeR,
                     emm.difference=mean.diff,p.emm=phen.NvsN.emm[[k]][,4])->p.phen.comp[[k]]
        }
        names(p.phen.comp)<-names(phen.reg)
      }
    }else{
      p.phen.comp<-NULL
    }

    data.frame(scalrat.data,group=rate.dataRES[match(rownames(scalrat.data),rownames(rate.dataRES)),]$group)->scalrat.data
  }

  CInts <- list(CIphenotype, CIabsolute,scalrat.CI)
  names(CInts) <- c("phenotype", "absolute_rate","rescaled_rate")
  class(CInts)<-"RRphyloList"

  list(phen.dataRES,rate.dataRES,scalrat.data)->tabs
  names(tabs) <- c("phenotypeVStime", "rateVStime","rescaledrateVStime")
  class(tabs)<-"RRphyloList"

  if (!is.null(node)) {
    {
      for (i in 1:length(p.rate.sel)) {
        for(k in 1:ncol(y)){
          ifelse(p.phen.sel[[i]][k,1]>0,{
            if(p.phen.sel[[i]][k,2]<=0.025) phenN.dir<-" increases less" else
              if(p.phen.sel[[i]][k,2]>=0.975) phenN.dir<-" increases more" else
                phenN.dir<-"nosig"
          },{
            if(p.phen.sel[[i]][k,2]<=0.025) phenN.dir<-" decreases more" else
              if(p.phen.sel[[i]][k,2]>=0.975) phenN.dir<-" decreases less" else
                phenN.dir<-"nosig"
          })

          if(phenN.dir!="nosig"){
            if(ncol(y)>1) yprint<-paste(" for",rownames(p.phen.sel[[i]])[k]) else yprint<-NULL
            print(paste("Phenotypic regression slope through node ",names(p.rate.sel)[i],
                        phenN.dir," than BM expectation",yprint,sep=""))
          }


          if(p.rate.sel[[i]][k,3]>p.rate.sel[[i]][k,4]) rateN.dir<-"higher" else rateN.dir<-"lower"
          if(p.rate.sel[[i]][k,5]<=0.05){
            if(ncol(y)>1) yprint<-paste(" for",rownames(p.phen.sel[[i]])[k]) else yprint<-NULL
            print(paste("Absolute evolutionary rates regression slope",yprint,
                        " through node ",names(p.rate.sel)[i],
                        " is ",rateN.dir," than for the rest of the tree",sep=""))
          }
        }
      }

      if(!is.null(p.phen.comp)){
        for (j in 1:ncol(y)) {
          for (i in 1:nrow(p.phen.comp[[j]])) {
            if ((p.phen.comp[[j]][i,3]-p.phen.comp[[j]][i,4])>0&p.phen.comp[[j]][i,5]>=0.95) phen.comp.dir<-"higher" else
              if ((p.phen.comp[[j]][i,3]-p.phen.comp[[j]][i,4])<0&p.phen.comp[[j]][i,5]<=0.05) phen.comp.dir<-"lower" else
                phen.comp.dir<-NULL
              if(!is.null(phen.comp.dir)){
                gsub("g","",p.phen.comp[[j]][i,1:2])->nprint
                if(ncol(y)>1) yprint<-paste(" for",names(p.phen.comp)[j]) else yprint<-NULL
                print(paste("Phenotypic regression through node ",nprint[1],
                            " is ",phen.comp.dir," than regression through node ",nprint[2],
                            yprint,sep=""))
              }

              if (p.rate.comp[[j]][i, 7] <= 0.05){
                if (p.rate.comp[[j]][i, 5]<p.rate.comp[[j]][i, 6]) rate.comp.dir<-"lower" else rate.comp.dir<-"higher"
                gsub("g","",p.rate.comp[[j]][i,1:2])->nprint
                if(ncol(y)>1) yprint<-paste(" for",names(p.rate.comp)[j]) else yprint<-NULL
                print(paste("Absolute evolutionary rates regression slope through node ",nprint[1],
                            " is ",rate.comp.dir,"than regression through node ",nprint[2],yprint,sep=""))
              }
          }
        }
      }
    }

    if(is.null(p.phen.comp)){
      node.comp.res<-NULL
    }else{
      if(length(p.phen.comp)==1) p.phen.comp[[1]]->p.phen.comp
      if(length(p.rate.comp)==1) p.rate.comp[[1]]->p.rate.comp
      node.comp.res <- list(p.phen.comp, p.rate.comp)
      names(node.comp.res) <- c("phenotype", "rate")
    }

    res <- list(tabs, p.phen, p.rate,
                p.phen.sel, p.rate.sel,
                node.comp.res, CInts)
    names(res) <- c("trend.data", "phenotypic.regression", "rate.regression",
                    "node.phenotypic.regression", "node.rate.regression",
                    "group.comparison", "ConfInts")



  }else {
    res <- list(tabs, p.phen, p.rate,
                CInts)
    names(res) <- c("trend.data", "phenotypic.regression", "rate.regression",
                    "ConfInts")
  }

  if(length(which(sapply(res,function(x) is.null(x))))>0)
    res<-res[-which(sapply(res,function(x) is.null(x)))]

  attr(res,"Call")<-funcall
  return(res)
}
