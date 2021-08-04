#'@title Searching for evolutionary trends in phenotypes and rates
#'@description This function searches for evolutionary trends in the phenotypic
#'  mean and the evolutionary rates for the entire tree and individual clades.
#'@usage search.trend(RR,y,x1=NULL,x1.residuals = FALSE,
#'  node=NULL,cov=NULL,nsim=100,clus=0.5,ConfInt=FALSE,foldername=NULL,filename)
#'@param RR an object produced by \code{\link{RRphylo}}.
#'@param y the named vector (or matrix if multivariate) of phenotypes.
#'@param x1 the additional predictor to be specified if the RR object has been
#'  created using an additional predictor (i.e. multiple version of
#'  \code{RRphylo}). \code{'x1'} vector must be as long as the number of nodes
#'  plus the number of tips of the tree, which can be obtained by running
#'  \code{RRphylo} on the predictor as well, and taking the vector of ancestral
#'  states and tip values to form the \code{x1}. Note: only one predictor at
#'  once can be specified.
#'@param x1.residuals logical specifying whether the residuals of regression
#'  between \code{y} and \code{x1} should be inspected for a phenotypic trend
#'  (see details and examples below). Default is \code{FALSE}.
#'@param node the node number of individual clades to be specifically tested and
#'  contrasted to each other. It is \code{NULL} by default. Notice the node
#'  number must refer to the dichotomic version of the original tree, as
#'  produced by \code{RRphylo}.
#'@param cov the covariate values to be specified if the RR object has been
#'  created using a  covariate for rates calculation.  As for \code{RRphylo},
#'  \code{'cov'} must be as long as the number of nodes plus the number of tips
#'  of the tree, which can be obtained by running \code{RRphylo} on the
#'  covariate as well, and taking the vector of ancestral states and tip values
#'  to form the covariate (see the example below).
#'@param nsim number of simulations to be performed. It is set at 100 by
#'  default.
#'@param clus the proportion of clusters to be used in parallel computing. To
#'  run the single-threaded version of \code{search.trend} set \code{clus} = 0.
#'@param ConfInt if \code{TRUE}, the function returns 95\% confidence intervals
#'  around phenotypes and rates produced according to the Brownian motion model
#'  of evolution. It is \code{FALSE} by default.
#'@param foldername has been deprecated; please see the argument \code{filename}
#'  instead.
#'@param filename a character indicating the name of the pdf file and the path
#'  where it is to be saved. If no path is indicated the file is stored in the
#'  working directory
#'@return The function returns a list object containing:
#'@return \strong{$trends.data} a 'RRphyloList' object including:
#'  \enumerate{\item{\code{$phenotypeVStime}}: a data frame of phenotypic values (or
#'  \code{y} versus \code{x1} regression residuals if \code{x1.residuals=TRUE})
#'  and their distance from the tree root for each node (i.e. ancestral states)
#'  and tip of the tree. \item{\code{$absrateVStime}}: a data frame of
#'  \code{RRphylo} rates and the distance from the tree root (age). If y is
#'  multivariate, it also includes the multiple rates for each y vector. If
#'  \code{node} is specified, each branch is classified as belonging to an
#'  indicated clade. \item{\code{$rescaledrateVStime}}: a data frame of rescaled
#'  \code{RRphylo} rates and the distance from the tree root (age). If y is
#'  multivariate, it also includes the multiple rates for each y vector. If
#'  \code{node} is specified, each branch is classified as belonging to an
#'  indicated clade. NAs correspond either to very small values or to outliers
#'  which are excluded from the analysis.}
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
#'  confidence intervals around phenotypes and rates (both rescaled and unscaled
#'  absolute rates) produced according to the Brownian motion model of
#'  evolution.
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
#'  difference between regression slopes for the group and for the rest of the
#'  tree (slope.difference), and a p-value for the slope.difference (p.slope).
#'@return If more than one node is specified, the object
#'  \strong{$group.comparison} reports the same results as
#'  $node.phenotypic.regression and $node.rate.regression obtained by comparing
#'  individual clades to each other.
#'@author Silvia Castiglione, Carmela Serio, Pasquale Raia, Alessandro
#'  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'  Carotenuto
#'@details The function simultaneously returns the regression of phenotypes and
#'  phenotypic evolutionary rates against age tested against Brownian motion
#'  simulations to assess significance. To this aim rates are rescaled in the
#'  0-1 range and then logged. The function stores the rates (both rescaled and
#'  unscaled absolute values) versus age regression and the phenotype versus age
#'  regression plots as .pdf files. In the plots, the 95\% confidence intervals
#'  of phenotypes and rates simulated under the Brownian motion for each node
#'  are plotted as shaded areas. Regression lines are printed for all
#'  regressions. To assess significance, slopes are compared to a family of
#'  simulated slopes (BMslopes, where the number of simulations is equal to
#'  \code{nsim}), generated under the Brownian motion, using the \code{fastBM}
#'  function in the package \pkg{phytools}. Individual nodes are compared to the
#'  rest of the tree in different ways depending on whether phenotypes or rates
#'  (always unscaled in this case) versus age regressions are tested. With the
#'  former, the regression slopes for individual clades and the slope difference
#'  between clades is contrasted to slopes obtained through Brownian motion
#'  simulations. For the latter, regression models are tested and contrasted to
#'  each other referring to estimated marginal means, by using the
#'  \code{emmeans} function in the package \pkg{emmeans}.
#'
#'  The \href{../doc/RRphylo.html#predictor}{multiple regression version of
#'  RRphylo} allows to incorporate the effect of an additional predictor in the
#'  computation of evolutionary rates without altering the ancestral character
#'  estimation. Thus, when a multiple \code{RRphylo} output is fed to
#'  \code{search.trend}, the predictor effect is accounted for on the absolute
#'  evolutionary rates, but not on the phenotype. However, in some situations
#'  the user might want to factor out the predictor effect on phenotypes
#'  as well. Under the latter circumstance, by setting the argument
#'  \code{x1.residuals = TRUE}, the residuals of the response to predictor
#'  regression are used as to represent the phenotype.
#'@importFrom graphics points text title polygon pairs plot
#'@importFrom stats as.formula coef resid density predict cor
#'@importFrom phytools nodeHeights
#'@importFrom parallel makeCluster detectCores stopCluster
#'@importFrom doParallel registerDoParallel
#'@importFrom foreach foreach %dopar%
#'@importFrom grDevices pdf dev.off
#'@importFrom utils combn
#'@importFrom emmeans emmeans emtrends
#'@export
#'@seealso \href{../doc/search.trend.html}{\code{search.trend} vignette}
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
#' RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'
#' # Case 1.1. "search.trend" whitout indicating nodes to be tested for trends
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, clus=cc,
#'              filename=paste(tempdir(), "ST", sep="/"),cov=NULL,ConfInt=FALSE,node=NULL)
#'
#' # Case 1.2. "search.trend" indicating nodes to be specifically tested for trends
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, node=143, clus=cc,
#'              filename=paste(tempdir(), "STnode", sep="/"),cov=NULL,ConfInt=FALSE)
#'
#'
#' # Case 2. "RRphylo" accounting for the effect of a covariate
#' # "RRphylo" on the covariate in order to retrieve ancestral state values
#' RRphylo(tree=treeptero,y=log(massptero))->RRptero
#' c(RRptero$aces,log(massptero))->cov.values
#' names(cov.values)<-c(rownames(RRptero$aces),names(massptero))
#' RRphylo(tree=treeptero,y=log(massptero),cov=cov.values)->RRpteroCov
#'
#' # Case 2.1. "search.trend" whitout indicating nodes to be tested for trends
#' search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, clus=cc,
#'              filename=paste(tempdir(), "ST_cov", sep="/"),ConfInt=FALSE,cov=cov.values)
#'
#' # Case 2.2. "search.trend" indicating nodes to be specifically tested for trends
#' search.trend(RR=RRpteroCov, y=log(massptero), nsim=100, node=143, clus=cc,
#'              filename=paste(tempdir(), "STnode_cov", sep="/"),ConfInt=FALSE,cov=cov.values)
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
#' RRphylo(tree=treecet.multi,y=masscet.multi)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#'
#' # incorporating the effect of body size at inspecting trends in absolute evolutionary rates
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc,
#'              filename=paste(tempdir(), "STmulti_rate", sep="/"))
#'
#' # incorporating the effect of body size at inspecting trends in both absolute evolutionary
#' # rates and phenotypic values (by using brain versus body mass regression residuals)
#' search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,clus=cc,
#'              filename=paste(tempdir(), "STmulti", sep="/"))
#'    }


search.trend<-function (RR,y,
                        x1=NULL,x1.residuals=FALSE,
                        node = NULL, cov = NULL,
                        nsim = 100, clus = 0.5, ConfInt = FALSE,
                        foldername=NULL,filename)
{
  # require(ape)
  # require(phytools)
  # require(geiger)
  # require(stats4)
  # require(foreach)
  # require(doParallel)
  # require(parallel)
  # require(nlme)
  # require(RColorBrewer)
  # require(emmeans)
  # require(outliers)
  # require(car)

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package \"car\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!missing(foldername)){
    stop("argument foldername is deprecated; please use filename instead.",
         call. = FALSE)
  }

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

  if(is.null(nrow(y))) y <- treedata(t, y, sort = TRUE)[[2]][,1] else y <- treedata(t, y, sort = TRUE)[[2]]

  H <- max(nodeHeights(t))
  eds <- t$edge[, 2]
  eds[which(t$edge[, 2] < Ntip(t) + 1)] <- t$tip.label
  eds <- c(Ntip(t) + 1, eds)
  hh <- c(0, nodeHeights(t)[, 2])
  eds <- data.frame(leaf = eds, height = hh)
  if (length(y) > Ntip(t)) {
    y.multi <- (L %*% rates)
    as.matrix(as.data.frame(y.multi[match(rownames(y),rownames(y.multi)),]))->y.multi
    aceRR.multi <- (L1 %*% rates[1:Nnode(t), ])
    rates <- as.data.frame(rates)
    betas <- as.data.frame(betas)

    rbi.rate <- data.frame(betas = betas[match(eds[, 1], rownames(betas)),
    ], rate = rates[match(eds[, 1], rownames(rates)),
    ], age = eds[, 2])
    colnames(rbi.rate)[1:dim(y)[2]] <- paste("betas", seq(1,
                                                          dim(y)[2], 1), sep = "")

  }else {
    rbi.rate <- data.frame(rate = rates[match(eds[, 1], rownames(rates)),
    ], age = eds[, 2])
    rownames(rbi.rate) <- rownames(rates)[match(eds[, 1], rownames(rates))]

  }

  rbi.rate[,dim(rbi.rate)[2]]+L[1,1]->rbi.rate[,dim(rbi.rate)[2]]

  if (length(y) > Ntip(t)) {##### Rate Trend Real Multi #####
    rbi.slopeA <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                             1))
    REGabs.betas <- list()
    e1<-array()
    scalrat.data<-matrix(ncol=(ncol(y)+1),nrow=nrow(rbi.rate))
    for (i in 1:(dim(rbi.rate)[2] - 1)) {
      bet <- rbi.rate[, i]
      age <- rbi.rate[, dim(rbi.rate)[2]]
      names(bet)<-names(age)<-rownames(rbi.rate)
      as.matrix(as.data.frame(abs(bet)))->rts
      rts->rtsA
      log(range01(rts))->rts

      c(which(rts=="-Inf"),which(age=="-Inf"))->outs
      if(length(outs)>0)
      {
        as.matrix(as.data.frame(rts[-outs,]))->rts
        age[-outs]->age
      }
      age->ageC

      sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->e1[i]

      if(!is.null(x1)){
        rts[-1,,drop=FALSE]->rts
        age[-1]->age

        car::outlierTest(lm(rts~age))->ouT
        if(length(which(ouT$bonf.p<=0.05))>0){
          rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
          age[-match(names(which(ouT$bonf.p<=0.05)),names(age))]->age
        }

      }else{
        lm(rts~age)->bb

        residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
        if((Ntip(t)+1)%in%names(resout)){
          as.matrix(as.data.frame(rts[-1,]))->rts
          age[-1]->age
        }
      }
      lm(rts~age)->regr.1
      rbi.slopeA[i, ] <- coef(summary(regr.1))[2, c(1, 4)]
      REGabs.betas[[i]] <- regr.1

      ##############
      scalrat.data[,i]<-rts[match(rownames(rbi.rate),rownames(rts)),]

    }
    colnames(scalrat.data)<-colnames(rbi.rate)[1:(ncol(y)+1)]
    rownames(scalrat.data)<-rownames(rbi.rate)
    data.frame(as.data.frame(scalrat.data),age=rbi.rate$age)->scalrat.data

    colnames(rbi.slopeA) <- c("slope","p-value")
    if(is.null(colnames(y))) rownames(rbi.slopeA) <- names(REGabs.betas)<-colnames(rbi.rate)[1:(ncol(y)+1)] else
      rownames(rbi.slopeA) <- names(REGabs.betas)<-c(colnames(y),"multiple")
  }else { #### Rate Trend Real Uni #####
    bet <- rbi.rate[, 1]
    age <- rbi.rate[, 2]
    names(bet)<-names(age)<-rownames(rbi.rate)
    as.matrix(as.data.frame(abs(bet)))->rts
    rts->rtsA
    log(range01(rts))->rts

    c(which(rts=="-Inf"),which(age=="-Inf"))->outs
    if(length(outs)>0)
    {
      as.matrix(as.data.frame(rts[-outs,]))->rts
      age[-outs]->age
    }
    age->ageC

    sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->e1

    if(!is.null(x1)){
      rts[-1,,drop=FALSE]->rts

      age[-1]->age

      car::outlierTest(lm(rts~age))->ouT
      if(length(which(ouT$bonf.p<=0.05))>0){
        rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
        age[-match(names(which(ouT$bonf.p<=0.05)),names(age))]->age
      }
    }else{

      lm(rts~age)->bb

      residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
      if((Ntip(t)+1)%in%names(resout)){
        as.matrix(as.data.frame(rts[-1,]))->rts
        age[-1]->age
      }
    }
    lm(rts~age)->regr.1
    rbi.slopeA<- coef(summary(regr.1))[2, c(1, 4)]
    REGabs <- regr.1

    data.frame(rate=rts[match(rownames(rbi.rate),rownames(rts)),1],age=age[match(rownames(rbi.rate),names(age))],row.names =rownames(rbi.rate))->scalrat.data
    rownames(scalrat.data)<-rownames(rbi.rate)
  }

  nodes <- aceRR[1:Nnode(t), ]
  rbiRES<-rbi.rate

  if (length(y) > Ntip(t)) {##### Phenotypic Trend Real Multi #####
    nodes <- cbind(aceRR[1:Nnode(t), ],aceRR.multi[1:Nnode(t), ])
    colnames(nodes)<- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
    P <- rbind(nodes, cbind(y,y.multi))
    if(isTRUE(x1.residuals)) apply(P,2,function(x) residuals(lm(x~x1)))->P
    #PP <- data.frame(P[match(rbi[, 1], rownames(P)), ],rbi$age)
    # PP <- data.frame(P[match(rownames(data), rownames(P)), ],data$age)
    PP <- data.frame(P[match(rownames(rbi.rate), rownames(P)), ],rbi.rate$age)
    colnames(PP)[dim(PP)[2]] <- "age"
    trendR <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x) lm(range01(x) ~ PP[, dim(PP)[2]]))
    trend.reg <- lapply(trendR, function(x) coefficients(summary(x)))
    if(is.null(colnames(y))) names(trend.reg) <- c(paste("y", seq(1:dim(y)[2]), sep = ""),"multiple") else
      names(trend.reg) <- c(colnames(y),"multiple")
    dev<-array()
    for(i in 1:length(trend.reg)){
      if(trend.reg[[i]][2,1]<0) (min(predict(trendR[[i]]))-mean(range01(PP[,i])))/sd(range01(PP[,i]))->dev[i] else (max(predict(trendR[[i]]))-mean(range01(PP[,i])))/sd(range01(PP[,i]))->dev[i]
    }
    names(dev)<-names(trend.reg)

  }else {##### Phenotypic Trend Real Uni #####
    P <- c(nodes, y)
    if(isTRUE(x1.residuals)) P<-residuals(lm(P~x1))
    # PP <- data.frame(P[match(rbi[, 1], names(P))], rbi$age)
    # PP <- data.frame(P[match(rownames(data), names(P))], data$age)
    PP <- data.frame(P[match(rownames(rbi.rate), names(P))], rbi.rate$age)
    colnames(PP) <- c("phenotype", "age")
    lm(range01(PP[,1])~PP[,2])->trendR
    summary(trendR)$coef->trend.reg
    if(trend.reg[2,1]<0) (min(predict(trendR))-mean(range01(PP[,1])))/sd(range01(PP[,1]))->dev else (max(predict(trendR))-mean(range01(PP[,1])))/sd(range01(PP[,1]))->dev

  }
  PPtot <- PP

  #### Nodes Real ####
  if (!is.null(node)) {
    rbi.sma <- rbi.rate
    PP.sma <- PP
    rbi.sma$group <- rep("others", dim(rbi.sma)[1])
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
      rep("others",nrow(rbi.sma))->gg
      gg[match(c(n,sele), rownames(rbi.sma))]<-paste("g",n, sep = "")
      data.frame(rbi.sma,gg)->rbi.sma
      # rbi.sel <- rbi[match(c(n,sele), rownames(rbi)), ]
      # rbi.sel <- data[match(c(n,sele), rownames(data)), ]
      rbi.sel <- rbi.rate[match(c(n,sele), rownames(rbi.rate)), ]

      if (length(y) > Ntip(t)) { #### Rate Trend Real Node Multi ####
        # rbi.rate <- rbi.sel[, c(5:(dim(rbi.sel)[2] -
        #                              1), dim(rbi.sel)[2])]
        ## rbi.rate <- rbi.sel
        REG.betas.age <- list()
        REGabs.betas.y <- list()
        # for (i in 1:(dim(rbi.rate)[2] - 1)) {
        #   bet <- rbi.rate[, i]
        #   age <- rbi.rate[, dim(rbi.rate)[2]]
        #   names(bet)<-names(age)<-rownames(rbi.rate)
        for (i in 1:(dim(rbi.sel)[2] - 1)) {
          bet <- rbi.sel[, i]
          age <- rbi.sel[, dim(rbi.sel)[2]]
          names(bet)<-names(age)<-rownames(rbi.sel)
          as.matrix(as.data.frame(abs(bet)))->rts
          REG.betas.age[[i]]<-age
          REGabs.betas.y[[i]]<-predict(lm(rts~age))
        }
        #names(REG.betas.age) <- names(REGabs.betas.y) <- colnames(data)[5:(5 +dim(y)[2])]
        #names(REG.betas.age) <- names(REGabs.betas.y) <- colnames(data)[1:(ncol(y)+1)]
        names(REG.betas.age) <- names(REGabs.betas.y) <- colnames(rbi.rate)[1:(ncol(y)+1)]
      }else {#### Rate Trend Real Node Uni ####
        # rbi.rate <- rbi.sel[, 5:6]
        ## rbi.rate <- rbi.sel
        # bet <- rbi.rate[, 1]
        # age <- rbi.rate[, 2]
        # names(bet)<-names(age)<-rownames(rbi.rate)
        bet <- rbi.sel[, 1]
        age <- rbi.sel[, 2]
        names(bet)<-names(age)<-rownames(rbi.sel)
        as.matrix(as.data.frame(abs(bet)))->rts
        REG.betas.age<-age
        REGabs.betas.y<-predict(lm(rts~age))

      }

      REGabs.betas.y.sel[[j]] <- REGabs.betas.y
      REG.betas.age.sel[[j]] <- REG.betas.age

      # nodes <- aceRR[1:Nnode(t), ]
      if (length(y) > Ntip(t)) {##### Phenotypic Trend Real Node Multi #####
        PPsel <- PP[match(rownames(rbi.sel), rownames(PP)),]
        # colnames(PP)[dim(PP)[2]] <- "age"
        trend.regC <- apply(PPsel[1:(dim(PPsel)[2] - 1)],
                            2, function(x) summary(lm(range01(x) ~ PPsel[, dim(PPsel)[2]])))
        trend.regC <- lapply(trend.regC, coefficients)
        #names(trend.regC) <- c(paste("y", seq(1:dim(y)[2]),sep = ""),"multiple")
        if(is.null(colnames(y))) names(trend.regC) <- c(paste("y", seq(1:dim(y)[2]), sep = ""),"multiple") else
          names(trend.regC) <- c(colnames(y),"multiple")
        trend.reg.y.sel[[j]] <- apply(PPsel[1:(dim(PPsel)[2] -
                                                 1)], 2, function(x) predict(lm(x ~ PPsel[, dim(PPsel)[2]])))
        trend.reg.age.sel[[j]] <- PPsel$age

      }else {##### Phenotypic Trend Real Node Multi #####
        # P <- c(nodes, y)
        # PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
        #                  rbi.sel$age)
        # colnames(PP) <- c("phenotype", "age")
        # summary(lm(range01(PP[,1])~PP[,2]))$coef->trend.regC
        # trend.reg.age.sel[[j]] <- PP$age
        # trend.reg.y.sel[[j]] <- predict(lm(PP))

        PPsel <- PP[match(rownames(rbi.sel), rownames(PP)),]
        # colnames(PPsel) <- c("phenotype", "age")
        summary(lm(range01(PPsel[,1])~PPsel[,2]))$coef->trend.regC
        trend.reg.age.sel[[j]] <- PPsel$age
        trend.reg.y.sel[[j]] <- predict(lm(PPsel))
      }
      trend.reg.sel[[j]] <- trend.regC
    }

    names(trend.reg.sel)<-node

    #colnames(rbi.sma)[4:ncol(rbi.sma)]<-paste("group",seq(1,(ncol(rbi.sma)-3)),sep="")
    if(length(y) > Ntip(t))  PP.sma <- cbind(PP.sma, rbi.sma[,(ncol(y)+3):ncol(rbi.sma)]) else
      PP.sma <- cbind(PP.sma, rbi.sma[,3:ncol(rbi.sma)])
    if (length(which(rbi.sma$group == "others")) < 3)
      rbi.sma <- rbi.sma[-which(rbi.sma$group == "others"),]


    rbiRES<-rbi.sma[,1:(ncol(rbi.sma)-length(node))]

    if (length(y) > Ntip(t)) { #### Node Comparison Multi ####
      sma.resA <- list()
      sma.resPP <- list()
      sma.resPPemm<-list()
      PPmeans.multi<-list()
      n.ot<-list()

      # group <- rbi.sma$group
      # age <- rbi.sma$age
      groupPP <- PP.sma[-which(PP.sma$group == "others"), ]$group
      agePP <- PP.sma$age


      for (i in 1:(dim(y)[2]+1)) {
        bets <- rbi.sma[, i]
        yPP<-PP.sma[,i]

        names(bets)<-rownames(rbi.sma)
        names(yPP)<-rownames(PP.sma)

        nn.ot<-list()
        PPn.ot<-list()
        for (w in 1:(length(node))){
          data.frame(rate=bets,age=rbi.sma$age,group=rbi.sma[,(w+ncol(y)+3)])->emdat
          suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(rate)~age+group,data=emdat),specs="group"))))
          mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(rate)~age*group,data=emdat),specs="group",var="age")))
          mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
          data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->nn.ot[[w]]

          data.frame(yPP,age=PP.sma$age,group=PP.sma[,(w+ncol(y)+3)])->PPemdat
          if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
            data.frame(PPemdat,x1=x1[match(rownames(PPemdat),names(x1))])->PPemdat
            suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+x1+group,data=PPemdat),specs="group"))))
          }else{
            suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+group,data=PPemdat),specs="group"))))
          }
          data.frame(do.call(rbind,strsplit(as.character(PPmeans[,1])," - ")),PPmeans[,-c(1,3,4,5)])->PPn.ot[[w]]

        }
        do.call(rbind,nn.ot)->n.ot[[i]]
        do.call(rbind,PPn.ot)->PPmeans
        colnames(n.ot[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
        colnames(PPmeans)<-c("group_1","group_2","mean","p.mean")

        dat <- data.frame(bets, age=rbi.sma$age, group=rbi.sma$group)
        suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(bets)~age+group,data=dat),specs="group"))))
        mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(bets)~age*group,data=dat),specs="group",var="age")))
        mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
        data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->sma.resA[[i]]
        colnames(sma.resA[[i]])<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
        sma.resA[[i]][-which(sma.resA[[i]]$group_2=="others"),]->sma.resA[[i]]

        dat <- data.frame(yPP, age=PP.sma$age, group=PP.sma$group)
        if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
          data.frame(dat,x1=x1[match(rownames(dat),names(x1))])->dat
          suppressMessages(PPpairs<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+x1+group,data=dat),specs="group"))))
        }else{
          suppressMessages(PPpairs<-as.data.frame(pairs(emmeans::emmeans(lm(range01(yPP)~age+group,data=dat),specs="group"))))
        }

        data.frame(do.call(rbind,strsplit(as.character(PPpairs[,1])," - ")),PPpairs[,-c(1,3,4,5)])->sma.resPPemt
        colnames(sma.resPPemt)<-c("group_1","group_2","mean","p.mean")
        #subset(sma.resPPemt,sma.resPPemt$group_2=="others")->PPmeans
        #subset(sma.resPPemt,sma.resPPemt$group_2!="others")->sma.resPPemt
        sma.resPPemt[-which(sma.resPPemt$group_2=="others"),]->sma.resPPemt
        sma.resPPemt->sma.resPPemm[[i]]


        sapply(strsplit(as.character(PPmeans[,1]),"g"),"[[",2)->PPnam
        PPmeans[,3:4]->PPmeans
        rownames(PPmeans)<-PPnam
        PPmeans[match(node,rownames(PPmeans)),]->PPmeans
        PPmeans->PPmeans.multi[[i]]


        if(length(node)>1){
          sapply(lapply(trend.reg.sel,"[[",i),"[[",2)->slope.tot
          names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")

          combn(sort(unique(as.character(groupPP))),2)->pair
          slope.diff<-array()
          for(jj in 1:dim(pair)[2]){
            slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

          }
          data.frame(t(pair),slope.diff)->sma.resPP[[i]]
          colnames(sma.resPP[[i]])<-c(colnames(sma.resA[[i]])[1:2],"estimate")

        }else{
          sma.resPP<-NULL
          sma.resPPemm<-NULL
        }
      }

      names(PPmeans.multi)<-colnames(PP.sma)[1:(dim(y)[2]+1)]

      # lapply(sma.resA,function(x) subset(x,x$group_2=="others"))->n.ot
      # lapply(sma.resA,function(x) subset(x,x$group_2!="others"))->sma.resA

      sapply(strsplit(as.character(n.ot[[1]][,1]),"g"),"[[",2)->grn
      lapply(n.ot,function(x) x[match(node,grn),])->n.ot
      rbi.slopeA.sel<-list()
      for(i in 1:length(node)){
        do.call(rbind,lapply(n.ot,function(x) x[i,3:6]))->rbi.slopeA.sel[[i]]
        rownames(rbi.slopeA.sel[[i]])<-rownames(rbi.slopeA)
      }

    }else {#### Node Comparison Uni ####
      colnames(PP.sma)[1] <- "y"
      n.ot<-list()
      PPn.ot<-list()
      for (w in 1:(length(node))){
        data.frame(rate=rbi.sma$rate,age=rbi.sma$age,group=rbi.sma[,(w+3)])->emdat
        suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(rate)~age+group,data=emdat),specs="group"))))
        mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(rate)~age*group,data=emdat),specs="group",var="age")))
        mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
        data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->n.ot[[w]]

        data.frame(y=PP.sma$y,age=PP.sma$age,group=PP.sma[,(w+3)])->PPemdat
        if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
          #data.frame(PPemdat,x1=x1[match(rownames(PPemdat),names(x1))])->PPemdat
          data.frame(PPemdat,x1=x1[match(rownames(PP.sma),names(x1))])->PPemdat
          suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(y)~age+x1+group,data=PPemdat),specs="group"))))
        }else{
          suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(y)~age+group,data=PPemdat),specs="group"))))
        }
        data.frame(do.call(rbind,strsplit(as.character(PPmeans[,1])," - ")),PPmeans[,-c(1,3,4,5)])->PPn.ot[[w]]
      }

      do.call(rbind,n.ot)->n.ot
      do.call(rbind,PPn.ot)->PPn.ot
      colnames(n.ot)<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
      colnames(PPn.ot)<-c("group_1","group_2","mean","p.mean")

      suppressMessages(mmeans<-as.data.frame(pairs(emmeans::emmeans(lm(abs(rate)~age+group,data=rbi.sma),specs="group"))))
      mtrends<-as.data.frame(pairs(emmeans::emtrends(lm(abs(rate)~age*group,data=rbi.sma),specs="group",var="age")))
      mtrends[match(mmeans[,1],mtrends[,1]),]->mtrends
      data.frame(do.call(rbind,strsplit(as.character(mmeans[,1])," - ")),mmeans[,-c(1,3,4,5)],mtrends[,c(2,6)])->sma.resA

      colnames(sma.resA)<-c("group_1","group_2","emm.difference","p.emm","slope.difference","p.slope")
      sma.resA[-which(sma.resA$group_2=="others"),]->sma.resA

      sapply(strsplit(as.character(n.ot[,1]),"g"),"[[",2)->grn
      n.ot[match(node,grn),]->n.ot
      rownames(n.ot)<-NULL
      rbi.slopeA.sel<-list()
      for(i in 1:nrow(n.ot)){
        n.ot[i,c(3,4,5,6)]->rbi.slopeA.sel[[i]]
        rownames(rbi.slopeA.sel[[i]])<-NULL
      }


      if((!is.null(x1))&isFALSE(x1.residuals)){ #### emmeans multiple ####
        data.frame(PP.sma,x1=x1[match(rownames(PP.sma),names(x1))])->PP.sma
        suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(y)~age+x1+group,data=PP.sma),specs="group"))))
      }else{
        suppressMessages(PPmeans<-as.data.frame(pairs(emmeans::emmeans(lm(range01(y)~age+group,data=PP.sma),specs="group"))))
      }

      data.frame(do.call(rbind,strsplit(as.character(PPmeans[,1])," - ")),PPmeans[,-c(1,3,4,5)])->sma.resPPemm
      colnames(sma.resPPemm)<-c("group_1","group_2","mean","p.mean")

      PPn.ot->PPmeans
      sma.resPPemm[-which(sma.resPPemm$group_2=="others"),]->sma.resPPemm
      sapply(strsplit(as.character(PPmeans[,1]),"g"),"[[",2)->PPnam
      PPmeans[,c(3,4)]->PPmeans
      rownames(PPmeans)<-PPnam
      PPmeans[match(node,rownames(PPmeans)),]->PPmeans


      if(dim(sma.resA)[1]>0){
        PP.sma[,3]<-as.character(PP.sma[,3])
        sapply(trend.reg.sel,"[[",2)->slope.tot
        names(slope.tot)<-paste("g",names(trend.reg.sel),sep="")

        combn(sort(unique(PP.sma[-which(PP.sma$group == "others"),3])),2)->pair
        slope.diff<-array()
        for(jj in 1:dim(pair)[2]){
          slope.tot[match(pair[1,jj],names(slope.tot))]-slope.tot[match(pair[2,jj],names(slope.tot))]->slope.diff[jj]

        }
        data.frame(t(pair),slope.diff)->sma.resPP
        colnames(sma.resPP)<-c(colnames(sma.resA)[1:2],"estimate")

      }else{
        sma.resPP<-NULL
        sma.resPPemm<-NULL
        sma.resA<-NULL
      }
    }
    names(rbi.slopeA.sel)<-node
    if (class(sma.resA) == "list") {
      if(is.null(sma.resPP)==FALSE){
        if(is.null(colnames(y))) names(sma.resA) <- colnames(rbi.sma)[1:(dim(y)[2]+1)] else
          names(sma.resA) <-c(colnames(y),"multiple")

        names(sma.resPP) <- names(sma.resPPemm) <- names(trend.reg)
      }
    }
  }


  #### Random ####
  RR$aces[1,]->a
  if(!is.null(x1)) x1[1:Ntip(t)]->y1
  if (length(y) > Ntip(t)) {
    yyD <- list()
    yyT <- list()

    if(is.null(x1)){
      for (i in 1:dim(y)[2]) {
        yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = a[i], bounds=c(min(y[,i]),max(y[,i])))))
        yyT[[i]] <- replicate(nsim,fastBM(t, sig2 = 1, a = mean(y[,i]), bounds=c(min(y[,i]),max(y[,i]))))
      }
    }else{
      if(isFALSE(x1.residuals)){
        for (i in 1:dim(y)[2])
          yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = a[i], bounds=c(min(y[,i]),max(y[,i])))))
      }else{
        PP[which(rownames(PP)==(Ntip(t)+1)),]->ares
        PP[which(rownames(PP)%in%t$tip.label),]->yres
        for (i in 1:dim(y)[2])
          yyD[[i]] <- suppressWarnings(replicate(nsim,fastBM(t, sig2 = 1, a = ares[1,i], bounds=c(min(yres[,i]),max(yres[,i])))))
      }

      matrix(ncol=nsim,nrow=Ntip(t))->yy1T
      matrix(ncol=nsim,nrow=Nnode(t))->ace1T
      for(k in 1:nsim){
        phenb<-list()
        for(i in 1:dim(y)[2]) fastBM(t, sig2 = 1, a = mean(y[,i]), bounds=c(min(y[,i]),max(y[,i])),internal = TRUE)->phenb[[i]]
        fastBM(t,a=mean(y1),bounds=c(min(y1),max(y1)),internal = TRUE)->phen1b
        do.call(cbind,phenb)->phenb

        cor(cbind(y,y1))->cr
        cbind(phenb,phen1b)->xx1
        cr->m
        #cor(xx1)->m
        #m[1:length(cr),ncol(m)]<-cr
        #m[nrow(m),1:length(cr)]<-cr
        t(chol(m))->U
        U%*%t(xx1)->xx3
        t(xx3[-nrow(xx3),1:Ntip(t)])->yyT[[k]]
        xx3[nrow(xx3),1:Ntip(t)]->yy1T[,k]
        xx3[nrow(xx3),(Ntip(t)+1):(Ntip(t)+Nnode(t))]->ace1T[,k]
      }
      rownames(yy1T)<-rownames(y)
      rownames(ace1T)<-rownames(aceRR)
    }
  }else {
    if(is.null(x1)){
      suppressWarnings(yyD<-replicate(nsim,fastBM(t,sig2=1,a=a,bounds=c(min(y),max(y)))))
      yyT<-replicate(nsim,fastBM(t,sig2=1,a=mean(y),bounds=c(min(y),max(y))))
    } else {
      if(isFALSE(x1.residuals)) suppressWarnings(yyD<-replicate(nsim,fastBM(t,sig2=1,a=a,bounds=c(min(y),max(y))))) else{
        PP[which(rownames(PP)%in%t$tip.label),1]->yres
        PP[which(rownames(PP)==(Ntip(t)+1)),1]->ares
        suppressWarnings(yyD<-replicate(nsim,fastBM(t,sig2=1,a=ares,bounds=c(min(yres),max(yres)))))
      }

      matrix(ncol=nsim,nrow=Ntip(t))->yyT
      matrix(ncol=nsim,nrow=Ntip(t))->yy1T
      matrix(ncol=nsim,nrow=Nnode(t))->ace1T

      for(i in 1:nsim){
        fastBM(t,a=mean(y),bounds=c(min(y),max(y)),internal = TRUE)->phenb
        fastBM(t,a=mean(y1),bounds=c(min(y1),max(y1)),internal = TRUE)->phen1b


        cor(y,y1)->cr
        rbind(phenb,phen1b)->xx1
        cor(t(xx1))->m
        m[-diag(m)]<-cr
        diag(m)<-rep(1,length(diag(m)))
        t(chol(m))->U
        U%*%xx1->xx3

        xx3[1,1:Ntip(t)]->yyT[,i]
        #xx3[1,(Ntip(t)+1):(Ntip(t)+Nnode(t))]->ace
        xx3[2,1:Ntip(t)]->yy1T[,i]
        xx3[2,(Ntip(t)+1):(Ntip(t)+Nnode(t))]->ace1T[,i]
      }
      rownames(yyT)<-rownames(yy1T)<-names(y)
      rownames(ace1T)<-rownames(aceRR)
    }
  }

  res <- list()
  if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else
    cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
  registerDoParallel(cl)
  res <- foreach(i = 1:nsim, .packages = c("car","nlme", "ape",
                                           "geiger", "phytools", "doParallel"
  )) %dopar% {
    gc()

    if (length(y) > Ntip(t)) {
      vec <- seq(1:nsim)
      yD <- lapply(yyD, function(x) x[, sample(vec, 1, replace = FALSE)])
      yD <- do.call(cbind, yD)

      if(is.null(x1)){
        yT <- lapply(yyT, function(x) x[, sample(vec, 1, replace = FALSE)])
        yT <- do.call(cbind, yT)
      }else{
        yT <-yyT[[i]]
        y1T<-yy1T[,i]
        aceT<-ace1T[,i]
        cbind(L,y1T-aceT[1])->LX
        cbind(L1,aceT-mean(aceT))->L1X
      }

      betasT<- matrix(ncol=dim(yT)[2], nrow=Ntip(t)+Nnode(t))
      aceRRT<- matrix(ncol=dim(yT)[2], nrow=Nnode(t))
      betasD<- matrix(ncol=dim(yD)[2], nrow=Ntip(t)+Nnode(t))
      aceRRD<- matrix(ncol=dim(yD)[2], nrow=Nnode(t))
      for (s in 1:dim(yT)[2]){
        rootV<-a[s]
        lambda <- RR$lambda[s]

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


        if(isTRUE(x1.residuals)) rootVD<-ares[,s] else rootVD<-a[s]
        m.betasD <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                       t(L)) %*% (as.matrix(yD[,s]) - rootVD)
        aceRRD[,s] <- (L1 %*% m.betasD[1:Nnode(t), ]) + rootVD
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
      if(!is.null(x1)) {
        y1T <- yy1T[, i]
        aceT<-ace1T[,i]

        cbind(L,y1T-aceT[1])->LX
        cbind(L1,aceT-mean(aceT))->L1X

        betasT <- (solve(t(LX) %*% LX + lambda * diag(ncol(LX))) %*%
                     t(LX)) %*% (as.matrix(yT) - rootV)
        aceRRT <- (L1X %*% betasT[c(1:Nnode(t),length(betasT)), ]) + rootV
        as.matrix(betasT[-length(betasT),])->betasT

      }else{

        betasT <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                     t(L)) %*% (as.matrix(yT) - rootV)
        aceRRT <- (L1 %*% betasT[1:Nnode(t), ]) + rootV
      }

      if(isTRUE(x1.residuals)) rootVD<-ares else rootVD<-a
      betasD <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%
                   t(L)) %*% (as.matrix(yD) - rootVD)
      aceRRD <- (L1 %*% betasD[1:Nnode(t), ]) + rootVD
    }

    betasT->betasTreal
    #### Covariate ####
    if (!is.null(cov)) {
      cov[match(rownames(betas),names(cov))]->cov

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
        ratesT<-as.matrix(as.data.frame(apply(betasT, 1, function(x) sqrt(sum(x^2)))))
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
    hh <- c(0, nodeHeights(t)[, 2])
    eds <- data.frame(leaf = eds, height = hh)

    if (length(y) > Ntip(t)) {
      betasTreal->betasT
      # data <- data.frame(betas = betasT[match(eds[, 1],
      #                                         rownames(betasT)), ], rate = ratesT[match(eds[,
      #                                                                                       1], rownames(ratesT)), ], age = eds[, 2])
      # colnames(data)[1:dim(yT)[2]] <- paste("betas", seq(1,
      #                                                    dim(yT)[2], 1), sep = "")

      rbi.rate <- data.frame(betas = betasT[match(eds[, 1],
                                                  rownames(betasT)), ], rate = ratesT[match(eds[,
                                                                                                1], rownames(ratesT)), ], age = eds[, 2])
      colnames(rbi.rate)[1:dim(yT)[2]] <- paste("betas", seq(1,
                                                             dim(yT)[2], 1), sep = "")
    }else {
      rates <- betasT
      # data <- data.frame(rate = rates[match(eds[, 1],
      #                                       rownames(rates)), ], age = eds[, 2])
      rbi.rate <- data.frame(rate = rates[match(eds[, 1],
                                                rownames(rates)), ], age = eds[, 2])
    }
    #data[,dim(data)[2]]+L[1,1]->data[,dim(data)[2]]
    rbi.rate[,dim(rbi.rate)[2]]+L[1,1]->rbi.rate[,dim(rbi.rate)[2]]

    if (length(yT) > Ntip(t)) {#### Rate Trend Random Multi ####
      rbi.slopeAS <- matrix(ncol = 2, nrow = (dim(rbi.rate)[2] -
                                                1))
      ee<-array()
      scalrat.data<-matrix(ncol=(ncol(y)+1),nrow=nrow(rbi.rate))
      for (k in 1:(dim(rbi.rate)[2] - 1)) {
        bet <- rbi.rate[, k]
        age <- rbi.rate[, dim(rbi.rate)[2]]
        names(bet)<-names(age)<-rownames(rbi.rate)

        as.matrix(as.data.frame(abs(bet)))->rts
        rts->rtsA
        log(range01(rts))->rts

        c(which(rts=="-Inf"),which(age=="-Inf"))->outs
        if(length(outs)>0)
        {
          as.matrix(as.data.frame(rts[-outs,]))->rts
          age[-outs]->age
        }
        age->ageC

        sd(range01(rtsA[ageC<0.5*max(ageC)]))/sd(range01(rtsA)[ageC>0.5*max(ageC)])->ee[k]

        if(!is.null(x1)){
          rts[-1,,drop=FALSE]->rts
          age[-1]->age

          car::outlierTest(lm(rts~age))->ouT
          if(length(which(ouT$bonf.p<=0.05))>0){
            rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
            age[-match(names(which(ouT$bonf.p<=0.05)),names(age))]->age
          }

        }else{
          lm(rts~age)->bb

          residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
          if((Ntip(t)+1)%in%names(resout)){
            as.matrix(as.data.frame(rts[-1,]))->rts
            age[-1]->age
          }
        }
        lm(rts~age)->regr.1
        rbi.slopeAS[k, ] <- coef(summary(regr.1))[2, c(1, 4)]
        scalrat.data[,k]<-rts[match(rownames(rbi.rate),rownames(rts)),1]

      }
      colnames(scalrat.data)<-colnames(rbi.rate)[1:(ncol(y)+1)]
      rownames(scalrat.data)<-rownames(rbi.rate)
      data.frame(as.data.frame(scalrat.data),age=rbi.rate$age)->scalrat.data


      colnames(rbi.slopeAS) <- c("slope","p-value")
      # rownames(rbi.slopeAS) <- colnames(data)[5:(5 +dim(y)[2])]
      # rownames(rbi.slopeAS) <- colnames(data)[1:(ncol(y)+1)]
      rownames(rbi.slopeAS) <- colnames(rbi.rate)[1:(ncol(y)+1)]
    }else {#### Rate Trend Random Uni ####
      # rbi.rate <- rbi[, 5:6]
      bet <- rbi.rate[, 1]
      age <- rbi.rate[, 2]
      names(bet)<-names(age)<-rownames(rbi.rate)
      as.matrix(as.data.frame(abs(bet)))->rts
      rts->rtsA
      log(range01(rts))->rts
      c(which(rts=="-Inf"),which(age=="-Inf"))->outs
      if(length(outs)>0)
      {
        as.matrix(as.data.frame(rts[-outs,]))->rts
        age[-outs]->age
      }
      age->ageC


      sd(range01(rtsA)[ageC<0.5*max(ageC)])/sd(range01(rtsA)[ageC>0.5*max(ageC)])->ee

      if(!is.null(x1)){
        rts[-1,,drop=FALSE]->rts
        age[-1]->age

        car::outlierTest(lm(rts~age))->ouT
        if(length(which(ouT$bonf.p<=0.05))>0){
          rts[-match(names(which(ouT$bonf.p<=0.05)),rownames(rts)),,drop=FALSE]->rts
          age[-match(names(which(ouT$bonf.p<=0.05)),names(age))]->age
        }

      }else{
        lm(rts~age)->bb

        residuals(bb)[order(residuals(bb),decreasing=TRUE)][1:(Ntip(t)/15)]->resout
        if((Ntip(t)+1)%in%names(resout)){
          as.matrix(as.data.frame(rts[-1,]))->rts
          age[-1]->age
        }
      }
      lm(rts~age)->regr.1
      rbi.slopeAS <- coef(summary(regr.1))[2, c(1, 4)]

      data.frame(rate=rts[match(rownames(rbi.rate),rownames(rts)),1],age=age[match(rownames(rbi.rate),names(age))],row.names =rownames(rbi.rate))->scalrat.data
      rownames(scalrat.data)<-rownames(rbi.rate)
    }

    nodesD <- aceRRD[1:Nnode(t), ]
    if (length(yD) > Ntip(t)) { #### Phenotypic Trend Random Multi ####
      nodesD<-cbind(nodesD,aceRRDmulti[1:Nnode(t), ])
      colnames(nodesD) <- c(paste("y", seq(1, dim(y)[2]),sep = ""),"y.multi")
      P <- rbind(nodesD, cbind(yD,yDmulti))
      # PP <- data.frame(P[match(rbi[, 1], rownames(P)),
      #                    ], rbi$age)
      if(isTRUE(x1.residuals)) apply(P,2,function(x) residuals(lm(x~x1)))->P
      # PP <- data.frame(P[match(rownames(data), rownames(P)), ],data$age)
      PP <- data.frame(P[match(rownames(rbi.rate), rownames(P)), ],rbi.rate$age)
      colnames(PP)[dim(PP)[2]] <- "age"
      trend.regR <- apply(PP[1:(dim(PP)[2] - 1)], 2, function(x)
        summary(lm(range01(x) ~PP[, dim(PP)[2]])))
      trend.regS <- lapply(trend.regR, coefficients)
      names(trend.regR) <- c(paste("betas", seq(1:dim(y)[2]),sep = ""),"multiple")
      cbind(yT,yTmulti)->yT
      cbind(yD,yDmulti)->yD

    }else { #### Phenotypic Trend Random Uni ####
      P <- c(nodesD, yD)
      if(isTRUE(x1.residuals)) P<-residuals(lm(P~x1))
      PP <- data.frame(P[match(rownames(rbi.rate), names(P))], rbi.rate$age)
      colnames(PP) <- c("phenotype", "age")
      trend.regS <- summary(lm(range01(PP[,1])~PP[,2]))$coef
    }

    PPtot <- PP

    #### Nodes Random ####
    if (!is.null(node)) {
      # rbi.sma <- rbi.rate
      PP.sma <- PP
      PP.sma$group <- rep("NA", dim(PP.sma)[1])
      trend.reg.SEL <- list()
      for (j in 1:length(node)) {
        n <- node[j]
        sele <- getDescendants(t, n)
        sele[which(sele < (Ntip(t) + 1))] <- t$tip.label[sele[which(sele <
                                                                      (Ntip(t) + 1))]]
        PP.sma[match(c(n,sele), rownames(PP.sma)), ]$group <- paste("g",
                                                                    n, sep = "")
        # rbi.sel <- rbi[match(c(n,sele), rownames(rbi)), ]

        # nodesD <- aceRRD[1:Nnode(t), ]
        if (length(y) > Ntip(t)) {#### Phenotypic Trend Random Node Multi ####
          # nodesD<-cbind(nodesD,aceRRDmulti[1:Nnode(t),])
          # colnames(nodesD)<- c(paste("y",seq(1, dim(y)[2]), sep = ""),"multiple")
          # P <- rbind(nodesD, yD)
          # PP <- data.frame(P[match(rbi.sel[, 1], rownames(P)),
          #                    ], rbi.sel$age)
          # colnames(PP)[dim(PP)[2]] <- "age"
          # trend.regC <- apply(PP[1:(dim(PP)[2] - 1)],
          #                     2, function(x) summary(lm(range01(x) ~ PP[, dim(PP)[2]])))

          # PPsel <- PP[match(rownames(data), rownames(PP)),]
          PPsel <- PP[match(rownames(rbi.rate), rownames(PP)),]
          # colnames(PP)[dim(PP)[2]] <- "age"
          trend.regC <- apply(PPsel[1:(dim(PPsel)[2] - 1)],
                              2, function(x) summary(lm(range01(x) ~ PPsel[, dim(PPsel)[2]])))
          trend.regC <- lapply(trend.regC, coefficients)
          names(trend.regC) <- c(paste("y", seq(1:dim(y)[2]),sep = ""),"multiple")

        }else {
          # P <- c(nodesD, yD)
          # PP <- data.frame(P[match(rbi.sel[, 1], names(P))],
          #                  rbi.sel$age)
          # colnames(PP) <- c("phenotype", "age")
          # trend.regC <- summary(lm(range01(PP[,1])~PP[,2]))$coef

          # PPsel <- PP[match(rownames(data), rownames(PP)),]
          PPsel <- PP[match(rownames(rbi.rate), rownames(PP)),]
          # colnames(PPsel) <- c("phenotype", "age")
          trend.regC <- summary(lm(range01(PPsel[,1])~PPsel[,2]))$coef

        }
        trend.reg.SEL[[j]] <- trend.regC
      }
      names(trend.reg.SEL) <- node


      if(!is.null(sma.resPP)){
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
                       PPtot, rbi.rate,yT,yD,ee,scalrat.data,
                       trend.reg.SEL,sma.resPP.R)
    }else {
      res[[i]] <- list(trend.regS, rbi.slopeAS,
                       PPtot, rbi.rate,yT,yD,ee,scalrat.data)
    }
  }
  stopCluster(cl)
  #### End of Random ####

  if (length(y) > Ntip(t)) {#### p Rate Trend Multi####
    p.rbi.slopeA <- array()
    spread<-array()
    cbind(y,y.multi)->yTot
    for (i in 1:(dim(y)[2] + 1)) {
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


    scalrat.CI <- CIabsolute <- list()
    for (i in 1:(dim(y)[2] + 1)) {
      scalrat.ci <- apply(do.call(cbind, lapply(lapply(res,"[[", 8),function(k) k[, i])), 1,
                          function(u){
                            if(all(is.na(u)))  rep(NA,2) else quantile(na.omit(u),c(0.025, 0.975))
                          })
      colnames(scalrat.ci) <- lapply(lapply(res,"[[", 8), rownames)[[1]]
      scalrat.CI[[i]] <- t(scalrat.ci)

      RBTAci <- apply(do.call(cbind, lapply(lapply(res,
                                                   "[[", 4),
                                            function(k) abs(k[, i]))), 1, function(u) quantile(u,
                                                                                               c(0.025, 0.975)))
      colnames(RBTAci) <- lapply(lapply(res,"[[", 4), rownames)[[1]]
      CIabsolute[[i]] <- t(RBTAci)
    }
    names(CIabsolute) <- names(scalrat.CI) <- c(paste("betas",seq(1, dim(y)[2]), sep = ""), "rate")

    scalrat.A <- scalrat.data
    A <- rbi.rate
    if(is.null(x1)){
      scalrat.age <- max(na.omit(scalrat.data[, ncol(scalrat.data)]))-scalrat.A[, ncol(scalrat.A)]
      age <- max(rbi.rate[, dim(rbi.rate)[2]])-A[, dim(A)[2]]
      names(scalrat.age)<-names(age)<-rownames(A)
    }else{
      scalrat.age <- max(scalrat.data[-which(rownames(scalrat.data)==(Ntip(t)+1)), ncol(scalrat.data)])-
        scalrat.A[-which(rownames(scalrat.A)==(Ntip(t)+1)), ncol(scalrat.A)]
      age <- max(rbi.rate[-which(rownames(rbi.rate)==(Ntip(t)+1)), ncol(rbi.rate)])-A[-which(rownames(A)==(Ntip(t)+1)), ncol(A)]
      names(scalrat.age)<-names(age)<-rownames(A)[-which(rownames(A)==(Ntip(t)+1))]
    }

    # pdf(file = paste(filename, "Rates.pdf",
    #                  sep = "-"))

    if (dim(y)[2] <= 2) {
      pdf(file = paste(filename, "Rates.pdf",sep = "-"))
      par(mfrow = c(dim(y)[2] + 1, 2),mar = c(2.5, 3, 2, 1))
      for (i in 1:(dim(y)[2] + 1)) {
        if (i == dim(y)[2] + 1) ynam <- "rate" else ynam <- paste("betas", i, sep = "")
        if(is.null(x1)) {
          scalrat.bet <- scalrat.A[, i]
          bet <- A[, i]
          names(bet)<-names(scalrat.bet)<-rownames(A)
          CIabsolute[[i]]->RBTAci
          scalrat.CI[[i]]->scalrat.ci
        }else{
          bet<-A[-which(rownames(A)==(Ntip(t)+1)),i]
          scalrat.bet<-scalrat.A[-which(rownames(scalrat.A)==(Ntip(t)+1)),i]
          names(bet)<-names(scalrat.bet)<-rownames(A)[-which(rownames(A)==(Ntip(t)+1))]
          CIabsolute[[i]][-which(rownames(CIabsolute[[i]])==(Ntip(t)+1)),]->RBTAci
          scalrat.CI[[i]][-which(rownames(scalrat.CI[[i]])==(Ntip(t)+1)),]->scalrat.ci
        }

        apply(scalrat.ci,2,function(j) all(is.na(j)))->NAci
        if(any(NAci)) scalrat.ci[,which(NAci)]<-scalrat.AA[match(names(which(NAci)),rownames(scalrat.AA)),1]
        cbind(scalrat.bet,scalrat.age,scalrat.ci)->scalrat.ci
        scalrat.ci <- scalrat.ci[order(scalrat.ci[, 2]), ]
        scalrat.ci[which(!is.na(scalrat.ci[, 1])), ]->scalrat.ci

        if(i==1) maintitle<-"rescaled rate" else maintitle<-" "
        plot(scalrat.bet ~ scalrat.age,data=scalrat.ci, main = maintitle,xlab="age", cex.axis=0.8,mgp = c(1.2, 0.4, 0),
             ylab = ynam,xlim=c(max(na.omit(scalrat.age)),min(na.omit(scalrat.age))))

        polygon(c(scalrat.ci[, 2], rev(scalrat.ci[, 2])),
                c(scalrat.ci[,3], rev(scalrat.ci[, 4])), col = rgb(0.5, 0.5, 0.5,0.4), border = NA)

        points(scalrat.age[which(rownames(scalrat.A)%in%t$tip.label)], scalrat.bet[which(rownames(scalrat.A)%in%t$tip.label)],
               pch = 21, col = "black",bg = "red")

        abline(lm(scalrat.bet ~ scalrat.age), lwd = 4, col = "blue")


        cbind(bet,age,RBTAci)->RBTAci
        RBTAci <- RBTAci[order(RBTAci[, 2]), ]
        if(i==1) maintitle<-"absolute rate" else maintitle<-" "
        plot(abs(bet) ~ age, main = maintitle, cex.axis=0.8,mgp = c(1.2, 0.4, 0),
             ylab = ynam,xlim=c(max(age),min(age)))
        polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                           3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5,
                                                                                            0.5, 0.4), border = NA)
        points(age[which(names(age)%in%t$tip.label)], abs(bet[which(names(bet)%in%t$tip.label)]), pch = 21, col = "black",
               bg = "red")

        if (class(node) != "NULL") {
          for (j in 1:length(node)) {
            cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
            points(max(rbi.rate$age)-REG.betas.age.sel[[j]][[i]], REGabs.betas.y.sel[[j]][[i]],
                   lwd = 4, col = cols[j], type = "l")
          }
          abline(lm(abs(bet) ~ age), lwd = 3, col = "blue",
                 lty = 2)
          if (i == 1)
            legend(max(age), max(abs(bet)), legend = node,
                   fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                    0, 0, 0), border = NA, x.intersp = 0.25)
        } else {
          abline(lm(abs(bet) ~ age), lwd = 4, col = "blue")
        }

      }
    } else {

      pdf(file = paste(filename, "Rates.pdf",sep = "-"),width=4.7,height=7)
      i=ncol(A)-1
      par(mfrow = c(2, 1),mar = c(2.5, 3, 2, 1))
      if(is.null(x1)) {
        scalrat.bet <- scalrat.A[, i]
        bet <- A[, i]
        names(bet)<-names(scalrat.bet)<-rownames(A)
        CIabsolute[[i]]->RBTAci
        scalrat.CI[[i]]->scalrat.ci
      }else{
        bet<-A[-which(rownames(A)==(Ntip(t)+1)),i]
        scalrat.bet<-scalrat.A[-which(rownames(scalrat.A)==(Ntip(t)+1)),i]
        names(bet)<-names(scalrat.bet)<-rownames(A)[-which(rownames(A)==(Ntip(t)+1))]
        CIabsolute[[i]][-which(rownames(CIabsolute[[i]])==(Ntip(t)+1)),]->RBTAci
        scalrat.CI[[i]][-which(rownames(scalrat.CI[[i]])==(Ntip(t)+1)),]->scalrat.ci
      }

      apply(scalrat.ci,2,function(j) all(is.na(j)))->NAci
      if(any(NAci)) scalrat.ci[,which(NAci)]<-scalrat.AA[match(names(which(NAci)),rownames(scalrat.AA)),1]
      cbind(scalrat.bet,scalrat.age,scalrat.ci)->scalrat.ci
      scalrat.ci <- scalrat.ci[order(scalrat.ci[, 2]), ]
      scalrat.ci[which(!is.na(scalrat.ci[, 1])), ]->scalrat.ci

      plot(scalrat.bet ~ scalrat.age,data=scalrat.ci, xlab="age",cex.axis=0.8,
           ylab = "rescaled rate", mgp = c(1.2, 0.4, 0),xlim=c(max(na.omit(scalrat.age)),min(na.omit(scalrat.age))))

      polygon(c(scalrat.ci[, 2], rev(scalrat.ci[, 2])),
              c(scalrat.ci[,3], rev(scalrat.ci[, 4])), col = rgb(0.5, 0.5, 0.5,0.4), border = NA)

      points(scalrat.age[which(rownames(scalrat.A)%in%t$tip.label)], scalrat.bet[which(rownames(scalrat.A)%in%t$tip.label)],
             pch = 21, col = "black",bg = "red")

      abline(lm(scalrat.bet ~ scalrat.age), lwd = 4, col = "blue")


      cbind(bet,age,RBTAci)->RBTAci
      RBTAci <- RBTAci[order(RBTAci[, 2]), ]
      par(mar = c(3, 3, 1, 1))
      plot(abs(bet) ~ age,cex.axis=0.8,
           ylab = "absolute rate", mgp = c(1.2, 0.4, 0),xlim=c(max(age),min(age)))
      polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                         3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5,
                                                                                          0.5, 0.4), border = NA)
      points(age[which(names(age)%in%t$tip.label)], abs(bet[which(names(bet)%in%t$tip.label)]), pch = 21, col = "black",
             bg = "red")

      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
          points(max(age)-REG.betas.age.sel[[j]][[i]], REGabs.betas.y.sel[[j]][[i]],
                 lwd = 4, col = cols[j], type = "l")
        }
        abline(lm(abs(bet) ~ age), lwd = 3, col = "blue",
               lty = 2)
        if (i == 1)
          legend(max(age), max(abs(bet)), legend = node,
                 fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                  0, 0, 0), border = NA, x.intersp = 0.25)
      } else {
        abline(lm(abs(bet) ~ age), lwd = 4, col = "blue")
      }

    }
    dev.off()
  }else {#### p Rate Trend Uni and pdf ####
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


    pdf(file = paste(filename, "Rates.pdf",
                     sep = "-"),width=4.7,height=7)
    par(mfrow = c(2, 1),mar = c(2.5, 3, 2, 1))
    scalrat.AA<-scalrat.data
    max(na.omit(scalrat.AA$age))-scalrat.AA$age->scalrat.AA$age
    scalrat.ci <- apply(do.call(cbind,lapply(lapply(res, "[[",8),function(x) x[, 1])), 1,
                        function(u){
                          if(all(is.na(u)))  rep(NA,2) else quantile(na.omit(u),c(0.025, 0.975))
                        })
    colnames(scalrat.ci) <- lapply(lapply(res,"[[", 8), rownames)[[1]]
    scalrat.CI <- t(scalrat.ci)

    apply(scalrat.ci,2,function(j) all(is.na(j)))->NAci
    if(any(NAci)) scalrat.ci[,which(NAci)]<-scalrat.AA[match(names(which(NAci)),rownames(scalrat.AA)),1]
    scalrat.ci <- cbind(scalrat.AA, t(scalrat.ci[, match(rownames(scalrat.AA), colnames(scalrat.ci))]))
    scalrat.ci <- scalrat.ci[order(scalrat.ci[, 2]), ]

    if(!is.null(x1)){
      scalrat.AA[-which(rownames(scalrat.AA)==(Ntip(t)+1)),]->scalrat.AA
      scalrat.ci[-which(rownames(scalrat.ci)==(Ntip(t)+1)),]->scalrat.ci
    }

    scalrat.ci[which(!is.na(scalrat.ci[, 1])), ]->scalrat.ci
    plot(rate ~ age, data = scalrat.ci, ylab = "rescaled rate",cex.axis=0.8,mgp = c(1.2, 0.4, 0),
         xlim=c(max(na.omit(scalrat.AA$age)),min(na.omit(scalrat.AA$age))))

    polygon(c(scalrat.ci[, 2], rev(scalrat.ci[, 2])),
            c(scalrat.ci[,3], rev(scalrat.ci[, 4])), col = rgb(0.5, 0.5, 0.5,0.4), border = NA)

    points(scalrat.AA$age[which(rownames(scalrat.AA)%in%t$tip.label)], scalrat.AA$rate[which(rownames(scalrat.AA)%in%t$tip.label)],
           pch = 21, col = "black",bg = "red")

    abline(lm(rate ~ age, data = scalrat.AA), col = "blue", lwd = 3)

    AA<-rbi.rate
    max(rbi.rate$age)-AA$age->AA$age
    RBTAci <- apply(do.call(cbind,lapply(lapply(res, "[[",4),function(x) abs(x[, 1]))), 1,
                    function(u) quantile(u,c(0.025, 0.975)))
    colnames(RBTAci) <- lapply(lapply(res,"[[", 4), rownames)[[1]]
    CIabsolute <- t(RBTAci)
    RBTAci <- cbind(AA, t(RBTAci[, match(rownames(AA), colnames(RBTAci))]))
    RBTAci <- RBTAci[order(RBTAci[, 2]), ]
    if(!is.null(x1)){
      AA[-which(rownames(AA)==(Ntip(t)+1)),]->AA
      RBTAci[-which(rownames(RBTAci)==(Ntip(t)+1)),]->RBTAci
    }

    par(mar = c(3, 3, 1, 1))
    plot(abs(rate) ~ age, data = AA, ylab = "absolute rate",
         mgp = c(2, 0.5, 0),xlim=c(max(AA$age),min(AA$age)),cex.axis=0.8,mgp = c(1.2, 0.4, 0))
    polygon(c(RBTAci[, 2], rev(RBTAci[, 2])), c(RBTAci[,
                                                       3], rev(RBTAci[, 4])), col = rgb(0.5, 0.5, 0.5,
                                                                                        0.4), border = NA)

    points(AA$age[which(rownames(AA)%in%t$tip.label)], abs(AA$rate[which(rownames(AA)%in%t$tip.label)]), pch = 21, col = "black",
           bg = "red")

    if (class(node) != "NULL") {
      for (j in 1:length(node)) {
        cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
        points(max(rbi.rate$age)-REG.betas.age.sel[[j]], REGabs.betas.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(lm(abs(rate) ~ age, data = AA), col = "blue", lwd = 3, lty = 2)
      legend(max(AA$age), max(abs(AA$rate)), legend = node,
             fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                              0, 0, 0), border = NA, x.intersp = 0.25)
    }else {
      abline(lm(abs(rate) ~ age, data = AA), col = "blue", lwd = 4)
    }
    dev.off()

  }



  p.trend <- array()
  if (length(y) > Ntip(t)) { #### p Phenotypic Trend Multi and pdf ####
    trend.slopes <- matrix(ncol = dim(y)[2] + 1, nrow = nsim)
    CIphenotype <- list()
    for (i in 1:(dim(y)[2] + 1)) {
      trend.slopes[, i] <- unlist(lapply(lapply(lapply(res,
                                                       "[[", 1), "[[", i), function(x) x[2,1]))

      p.trend[i] <- (nsim-length(which(trend.slopes[, i]<trend.reg[[i]][2,1])))/nsim

      trend.real <- do.call(rbind, lapply(trend.reg, function(x) x[2,
                                                                   c(1, 4)]))


      PPci <- apply(do.call(cbind, lapply(lapply(res,
                                                 "[[", 3), function(x) x[, i])), 1, function(u) quantile(u,
                                                                                                         c(0.025, 0.975)))
      colnames(PPci) <- lapply(lapply(res, "[[",
                                      3), rownames)[[1]]
      CIphenotype[[i]] <- t(PPci)
    }
    p.trend <- cbind(trend.real, p.trend,dev)
    colnames(p.trend) <- c("slope", "p.real",
                           "p.random","dev")

    PP[,dim(PP)[2]]<-max(PPtot[,dim(PPtot)[2]])-PP[,dim(PP)[2]]
    if (dim(y)[2] <= 2) {
      pdf(file = paste(filename, "Phenotypes.pdf",sep = "-"))
      par(mar = c(2.5, 3, 2, 1),mfrow = c(dim(y)[2] + 1, 2))
      for (i in 1:(dim(y)[2] + 1)) {
        PPci <- cbind(PP[, c(i, dim(PP)[2])], CIphenotype[[i]])
        PPci <- PPci[order(PPci[, 2]), ]
        obj <- hist(trend.slopes[, i], xlab = "",
                    ylab = "", main = names(trend.reg[i]), cex.axis=0.8,mgp = c(1.2, 0.4, 0))
        title(xlab = "simulated slopes", ylab = "frequency",
              line = 1.5)
        if(p.trend[i,3]>0.5) 1-p.trend[i,3]->pp else p.trend[i,3]->pp
        text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
               0.9, labels = paste("p=", pp),
             cex = 1)
        abline(v = trend.reg[[i]][2, 1], lwd = 3,
               col = "red")
        plot(PP[, c(dim(PP)[2], i)], xlab = "", ylab = "",
             cex.axis=0.8,mgp = c(1.2, 0.4, 0),xlim=c(max(PP[,dim(PP)[2]]),min(PP[,dim(PP)[2]])))
        polygon(c(PPci[,2], rev(PPci[, 2])), c(PPci[,3], rev(PPci[, 4])), col = rgb(0.5, 0.5,0.5, 0.4), border = NA)
        title(xlab = "age", ylab = paste(colnames(PP)[i]),
              line = 1.5)
        # points((max(diag(vcv(t)))-diag(vcv(t))), yTot[, i],xlim=c(max(PP[,dim(PP)[2]]),min(PP[,dim(PP)[2]])), pch = 21, col = "black",
        #        bg = "red")
        points((max(diag(vcv(t)))-diag(vcv(t))), PP[match(t$tip.label,rownames(PP)),i],xlim=c(max(PP[,dim(PP)[2]]),min(PP[,dim(PP)[2]])), pch = 21, col = "black",
               bg = "red")
        if (class(node) != "NULL") {
          for (j in 1:length(node)) {
            cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
            points(max(PPtot[,dim(PPtot)[2]])-trend.reg.age.sel[[j]], trend.reg.y.sel[[j]][,
                                                                                           i], lwd = 4, col = cols[j], type = "l")
          }
          abline(lm(PP[, i] ~ age, data = PP), lwd = 3,
                 col = "blue", lty = 2)
          if (i == 1)
            legend(max(PP$age),max(PP[, i]), legend = node,
                   fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                                    0, 0, 0), border = NA, x.intersp = 0.25)
        }else {
          abline(lm(PP[, i] ~ age, data = PP), lwd = 4,
                 col = "blue")
        }

      }
    }else {
      pdf(file = paste(filename, "Phenotypes.pdf",sep = "-"),width=4.7,height=7)
      par(mfrow = c(2, 1),mar = c(2.5, 3, 2, 1))
      i <- dim(yTot)[2]
      obj <- hist(trend.slopes[, i], xlab = "", ylab = "",
                  main = "Phenotypic Trend Test", cex.axis=0.8,mgp = c(1.2, 0.4, 0))
      title(xlab = "simulated slopes", ylab = "frequency",
            line = 1.5)
      if(p.trend[i,3]>0.5) 1-p.trend[i,3]->pp else p.trend[i,3]->pp
      text(quantile(trend.slopes[, i], 0.01), max(obj$counts) *
             0.9, labels = paste("p=", pp), cex = 1)
      abline(v = trend.reg[[i]][2, 1], lwd = 3, col = "red")
      par(mar = c(3, 3, 1, 1))
      plot(PP[, c(dim(PP)[2], i)], xlab = "", ylab = "",
           cex.axis=0.8,mgp = c(1.2, 0.4, 0),xlim=c(max(PP[,dim(PP)[2]]),min(PP[,dim(PP)[2]])))
      # points((max(diag(vcv(t)))-diag(vcv(t))), yTot[, i], pch = 21, col = "black",
      #        bg = "red")
      points((max(diag(vcv(t)))-diag(vcv(t))),PP[match(t$tip.label,rownames(PP)),i], pch = 21, col = "black",
             bg = "red")

      title(xlab = "age", ylab = "y.multi",
            line = 1.5)
      if (class(node) != "NULL") {
        for (j in 1:length(node)) {
          cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
          points(max(PPtot[,dim(PPtot)[2]])-trend.reg.age.sel[[j]], trend.reg.y.sel[[j]][,
                                                                                         i], lwd = 4, col = cols[j], type = "l")
        }
        abline(lm(PP[, i] ~ age, data = PP), lwd = 3,
               col = "blue", lty = 2)
        legend(max(PP$age), max(PP[, i]), legend = node,
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
    pdf(file = paste(filename, "Phenotypes.pdf",sep = "-"),width=4.7,height=7)
    par(mfrow = c(2, 1),mar = c(2.5, 3, 2, 1))
    obj <- hist(trend.slopes, xlab = "simulated slopes",
                ylab = "frequency", main = "Phenotypic Trend Test",
                cex.axis=0.8,mgp = c(1.2, 0.4, 0))
    text(quantile(trend.slopes, 0.01), max(obj$counts) *
           0.8, labels = paste("p=", p.trend[3]), cex = 1)
    abline(v = trend.reg[2, 1], lwd = 3, col = "red")
    max(PPtot[,2])-PP[,2]->PP[,2]
    par(mar = c(3, 3, 1, 1))
    plot(PP[, c(2, 1)],cex.axis=0.8,mgp = c(1.2, 0.4, 0),xlim=c(max(PP[,2]),min(PP[,2])))
    PPci <- apply(do.call(cbind, lapply(lapply(res, "[[",3),function(x) x[, 1])), 1, function(u) quantile(u,c(0.025, 0.975)))
    colnames(PPci) <- lapply(lapply(res, "[[", 3), rownames)[[1]]
    CIphenotype <- t(PPci)
    PPci <- cbind(PP, t(PPci[, match(rownames(PP), colnames(PPci))]))
    PPci <- PPci[order(PPci[, 2]), ]
    polygon(c(PPci[, 2], rev(PPci[, 2])), c(PPci[, 3], rev(PPci[,
                                                                4])), col = rgb(0.5, 0.5, 0.5, 0.4), border = NA)
    # points((max(diag(vcv(t)))-diag(vcv(t))), y, pch = 21, col = "black", bg = "red")
    points((max(diag(vcv(t)))-diag(vcv(t))), PP[match(t$tip.label,rownames(PP)),1], pch = 21, col = "black", bg = "red")
    if (!is.null(node)) {
      for (j in 1:length(node)) {
        cols <- suppressWarnings(RColorBrewer::brewer.pal(length(node), "Set2"))
        points(max(PPtot[,2])-trend.reg.age.sel[[j]], trend.reg.y.sel[[j]],
               lwd = 4, col = cols[j], type = "l")
      }
      abline(lm(PP), lwd = 3, col = "blue", lty = 2)
      legend(max(PP[, 2]), max(PP[, 1]), legend = node,
             fill = cols, bg = rgb(0, 0, 0, 0), box.col = rgb(0,
                                                              0, 0, 0), border = NA, x.intersp = 0.25)
    }else {
      abline(lm(PP), lwd = 4, col = "blue")
    }
    dev.off()
  }

  if (!is.null(node)) {
    p.trend.sel <- list()
    for (u in 1:length(node)) {
      if (length(y) > Ntip(t)) {#### p Phenotypic Trend Node Multi ####
        p.sele <- list()
        for (i in 1:(dim(y)[2] + 1)) {
          slopeR <- unlist(lapply(lapply(lapply(lapply(res,
                                                       "[[", 9), "[[", u), "[[", i), function(x) x[2,
                                                                                                   1]))

          p.sel <- (nsim-length(which(slopeR<trend.reg.sel[[u]][[i]][2,1])))/nsim

          p.sele[[i]] <- c(slope = trend.reg.sel[[u]][[i]][2,1],
                           p.slope = p.sel,emm.difference=PPmeans.multi[[i]][match(names(trend.reg.sel)[u],rownames(PPmeans.multi[[i]])),1],
                           p.emm=PPmeans.multi[[i]][match(names(trend.reg.sel)[u],rownames(PPmeans.multi[[i]])),2])

        }
        names(p.sele) <- names(trend.reg.sel[[u]])
        p.selt <- do.call(rbind, p.sele)
        p.trend.sel[[u]] <- p.selt
      } else {#### p Phenotypic Trend Node Uni ####
        slopeR <- unlist(lapply(lapply(lapply(res, "[[",
                                              9), "[[", u), function(x) x[2,1]))

        p.sel <- (nsim-length(which(slopeR<trend.reg.sel[[u]][2,1])))/nsim
        p.trend.sel[[u]] <- c(slope = trend.reg.sel[[u]][2,1],p.slope = p.sel,
                              emm.difference=unname(PPmeans[match(names(trend.reg.sel)[u],rownames(PPmeans)),1]),
                              p.emm=unname(PPmeans[match(names(trend.reg.sel)[u],rownames(PPmeans)),2]))
      }
    }
    names(p.trend.sel) <- names(trend.reg.sel)


    p.rbi.slopeA.sel<-rbi.slopeA.sel
    p.smaA <- sma.resA
    if(!is.null(sma.resPP)){
      if (length(y) > Ntip(t)) {#### p Phenotypic Slope Comparison Multi  ####
        p.smaPP<-list()
        for (k in 1:(dim(y)[2] + 1)) {
          p.smaR<-array()
          mean.diff<-array()
          for(i in 1:nrow(sma.resPP[[k]])){
            rank(c(sma.resPP[[k]][i,3],sapply(lapply(lapply(res,"[[",10),"[[",k),function(x) x[i,3])[1:(nsim-1)]))[1]/nsim->p.smaR[i]
            if(as.character(interaction(sma.resPP[[k]][i,c(1,2)]))==as.character(interaction(sma.resPPemm[[k]][i,c(1,2)]))) sma.resPPemm[[k]][i,3]->mean.diff[i] else
              -1*sma.resPPemm[[k]][i,3]->mean.diff[i]
          }
          data.frame(sma.resPP[[k]][,c(1,2)],slope.difference=sma.resPP[[k]][,3],p.slope=1-p.smaR,emm.difference=mean.diff,p.emm=sma.resPPemm[[k]][,4])->p.smaPP[[k]]
        }
        names(p.smaPP)<-names(trend.reg)
      }else{#### p Phenotypic Slope Comparison Uni  ####
        p.smaR<-array()
        mean.diff<-array()
        for(i in 1:nrow(sma.resPP)){
          rank(c(sma.resPP[i,3],sapply(lapply(res,"[[",10),function(x) x[i,3])[1:(nsim-1)]))[1]/nsim->p.smaR[i]
          if(as.character(interaction(sma.resPP[i,c(1,2)]))==as.character(interaction(sma.resPPemm[i,c(1,2)]))) sma.resPPemm[i,3]->mean.diff[i] else
            -1*sma.resPPemm[i,3]->mean.diff[i]
        }
        p.smaPP <-data.frame(sma.resPP[,1:2],slope.difference=sma.resPP[,3],p.slope=1-p.smaR,emm.difference=mean.diff,p.emm=sma.resPPemm[,4])
      }
    }else{
      p.smaPP<-NULL
    }

    data.frame(scalrat.data,group=rbiRES[match(rownames(scalrat.data),rownames(rbiRES)),]$group)->scalrat.data
  }


  if (length(y) > Ntip(t)) {
    for (i in 1:dim(y)[2]) {
      if (p.trend[i, 3] <= 0.05) print(paste("There is a trend for increase in phenotype y",i, sep = ""))
      if (p.trend[i, 3] >= 0.95) print(paste("There is a trend for decrease in phenotype y",i, sep = ""))
      if (p.rbi.slopeA[i, 1]>0 & p.rbi.slopeA[i, 3] <= 0.05) print(paste("Absolute evolutionary rates increase faster than BM for variable y",i, sep = ""))
      if (p.rbi.slopeA[i, 1]<0 & p.rbi.slopeA[i, 3] <= 0.05) print(paste("Absolute evolutionary rates decrease slower than BM for variable y",i, sep = ""))
      if (p.rbi.slopeA[i, 1]>0 & p.rbi.slopeA[i, 3] >= 0.95) print(paste("Absolute evolutionary rates increase slower than BM for variable y",i, sep = ""))
      if (p.rbi.slopeA[i, 1]<0 & p.rbi.slopeA[i, 3] >= 0.95) print(paste("Absolute evolutionary rates decrease faster than BM for variable y",i, sep = ""))

    }
    for(i in 1:nrow(p.trend)){
      if(p.trend[i,3]>=0.5) 1-p.trend[i,3]->p.trend[i,3]
      if(p.rbi.slopeA[i,3]>=0.5) 1-p.rbi.slopeA[i,3]->p.rbi.slopeA[i,3]
      #p.rbi.slopeA[i,3]*2->p.rbi.slopeA[i,3]
    }
  }else {
    if (p.trend[3] <= 0.05) print("There is a trend for increase in phenotype")
    if (p.trend[3] >= 0.95) print("There is a trend for decrease in phenotype")
    if (p.rbi.slopeA[1]>0 & p.rbi.slopeA[3] <= 0.05) print("Absolute evolutionary rates increase faster than BM")
    if (p.rbi.slopeA[1]<0 & p.rbi.slopeA[3] <= 0.05) print("Absolute evolutionary rates decrease slower than BM")
    if (p.rbi.slopeA[1]>0 & p.rbi.slopeA[3] >= 0.95) print("Absolute evolutionary rates increase slower than BM")
    if (p.rbi.slopeA[1]<0 & p.rbi.slopeA[3] >= 0.95) print("Absolute evolutionary rates decrease faster than BM")

    if(p.trend[3] >= 0.5) 1-p.trend[3]->p.trend[3]
    if(p.rbi.slopeA[3] >= 0.5) 1-p.rbi.slopeA[3]->p.rbi.slopeA[3]
    #p.rbi.slopeA[3]*2->p.rbi.slopeA[3]
  }

  CInts <- list(CIphenotype, CIabsolute,scalrat.CI)
  names(CInts) <- c("phenotype", "absolute rate","rescaled rate")
  class(CInts)<-"RRphyloList"

  list(PPtot,rbiRES,scalrat.data)->tabs
  names(tabs) <- c("phenotypeVStime", "absrateVStime","rescaledrateVStime")
  class(tabs)<-"RRphyloList"

  if (!is.null(node)) {
    if (length(y) > Ntip(t)) {
      for (i in 1:length(p.rbi.slopeA.sel)) {
        for(k in 1:dim(p.rbi.slopeA.sel[[i]])[1]){
          if(k<dim(p.rbi.slopeA.sel[[i]])[1]){
            # if(k==dim(p.rbi.slopeA.sel[[i]])[1]) ".multi"-> vv else k->vv
            if(p.trend.sel[[i]][k,2]<=.05)
              print(paste("There is a trend for increase in phenotype for variable",paste("y",k,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))
            if(p.trend.sel[[i]][k,2]>=.95)
              print(paste("There is a trend for decrease in phenotype for variable",paste("y",k,sep=""),"through node",names(p.rbi.slopeA.sel)[i]))
            if(p.rbi.slopeA.sel[[i]][k,4]<=.05&p.rbi.slopeA.sel[[i]][k,3]>0)
              print(paste("Absolute evolutionary rates regression slope for variable",paste("y",k,sep=""),"through node",names(p.rbi.slopeA.sel)[i],"is higher than for the rest of the tree"))
            if(p.rbi.slopeA.sel[[i]][k,4]<=.05&p.rbi.slopeA.sel[[i]][k,3]<0)
              print(paste("Absolute evolutionary rates regression slope for variable",paste("y",k,sep=""),"through node",names(p.rbi.slopeA.sel)[i],"is lower than for the rest of the tree"))
          }
          if(p.trend.sel[[i]][k,2]>=.5) 1-p.trend.sel[[i]][k,2]->p.trend.sel[[i]][k,2]
        }
      }


      if(is.null(p.smaPP)==FALSE){
        for (j in 1:(length(p.smaPP)-1)) {
          #if(j==length(p.smaPP)) ".multi"->vv else j->vv
          for (i in 1:dim(p.smaPP[[j]])[1]) {
            if (p.smaPP[[j]][i, 4] >= 0.95){
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is lower than regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    2]), "g"), "[[", 2), "for variable",
                          paste("y",j,sep="")))
            }
            if (p.smaPP[[j]][i, 4] <= 0.05){
              print(paste("Phenotypic regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    1]), "g"), "[[", 2), "is higher than regression through node",
                          lapply(strsplit(as.character(p.smaPP[[j]][i,
                                                                    2]), "g"), "[[", 2), "for variable",
                          paste("y",j,sep="")))

            }

            if (p.smaPP[[j]][i, 4] >= 0.5) 1-p.smaPP[[j]][i, 4]->p.smaPP[[j]][i, 4]

            if (p.smaA[[j]][i, 6] <= 0.05){
              if(p.smaA[[j]][i, 5] >0){
                print(paste("Absolute evolutionary rates regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     1]), "g"), "[[", 2), "is higher than regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     2]), "g"), "[[", 2), "for variable",
                            paste("y",j,sep="")))
              }else{
                print(paste("Absolute evolutionary rates regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     1]), "g"), "[[", 2), "is lower than regression through node",
                            lapply(strsplit(as.character(p.smaA[[j]][i,
                                                                     2]), "g"), "[[", 2), "for variable",
                            paste("y",j,sep="")))
              }
            }
          }
        }
      }
    }else {
      for (i in 1:length(p.rbi.slopeA.sel)) {
        if(p.trend.sel[[i]][2]<=.05)
          print(paste("There is a trend for increase in phenotype through node",names(p.rbi.slopeA.sel)[i]))
        if(p.trend.sel[[i]][2]>=.95)
          print(paste("There is a trend for decrease in phenotype through node",names(p.rbi.slopeA.sel)[i]))
        if(p.rbi.slopeA.sel[[i]][,2]<=.05&p.rbi.slopeA.sel[[i]][,1]>0)
          print(paste("Absolute evolutionary rates regression through node",names(p.rbi.slopeA.sel)[i],"is higher than for the rest of the tree"))
        if(p.rbi.slopeA.sel[[i]][,2]<=.05&p.rbi.slopeA.sel[[i]][,1]<0)
          print(paste("Absolute evolutionary rates regression through node",names(p.rbi.slopeA.sel)[i],"is lower than for the rest of the tree"))
        if(p.trend.sel[[i]][2]>=.5) 1-p.trend.sel[[i]][2]->p.trend.sel[[i]][2]
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
          if (p.smaPP[i, 4] >= 0.5) 1-p.smaPP[i, 4]->p.smaPP[i, 4]
        }

        if (p.smaA[i, 6] <= 0.05){
          if (p.smaA[i, 5]>0){
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

    if (isTRUE(ConfInt)) {
      res <- list(tabs, p.trend, p.rbi.slopeA,
                  p.trend.sel, p.rbi.slopeA.sel,
                  SMA.res, CInts)
      names(res) <- c("trend.data", "phenotypic.regression", "rate.regression",
                      "node.phenotypic.regression", "node.rate.regression",
                      "group.comparison", "ConfInts")
    }else {
      res <- list(tabs, p.trend, p.rbi.slopeA,
                  p.trend.sel, p.rbi.slopeA.sel,
                  SMA.res)
      names(res) <- c("trend.data", "phenotypic.regression", "rate.regression",
                      "node.phenotypic.regression", "node.rate.regression",
                      "group.comparison")
    }

  }else {
    if (ConfInt == TRUE) {
      res <- list(tabs, p.trend, p.rbi.slopeA,
                  CInts)
      names(res) <- c("trend.data", "phenotypic.regression", "rate.regression",
                      "ConfInts")
    }else {
      res <- list(tabs, p.trend, p.rbi.slopeA)
      names(res) <- c("trend.data", "phenotypic.regression", "rate.regression")
    }
  }

  if(length(which(sapply(res,function(x) is.null(x))))>0)
    res<-res[-which(sapply(res,function(x) is.null(x)))]
  return(res)
}
