#' @title Rate By Clade
#' @description The \code{'RBC'} function (\cite{Piras et al. 2018}) provides a test for Brownian rate variation via the non-censored approach (\cite{O'Meara et al. 2006}).
#' @usage RBC(RR,y,n.shift=c("clust","fix"),NS=3,clus=.5,f=NULL,foldername)
#' @param RR an object fitted by the function \code{\link{RRphylo}}.
#' @param y either a single vector variable or a multivariate dataset of class \sQuote{matrix}.
#' @param n.shift it specifies whether the function automatically searches for shifts (\code{"clust"}) or the number of shifts to search for has to be set by the user (\code{"fix"}).
#' @param NS the number of shifts to be specified for \code{n.shift = "fix"}. It is set at 3 by default.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @param f the size of the smallest clade to be tested in \code{RBC}. By default, nodes subtending to at least one tenth of the tree tips are tested.
#' @param foldername the path of the folder where plots are to be found.
#' @export
#' @importFrom ape Ntip dist.nodes drop.tip subtrees nodelabels
#' @importFrom stats pchisq
#' @importFrom geiger tips ratematrix
#' @importFrom phytools nodeHeights make.simmap plotSimmap brownieREML
#' @importFrom pvclust pvclust pvpick
#' @importFrom utils combn
#' @importFrom mvMORPH mvBM LRT
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @return \strong{rbc} brownian rate values for all possible clades identified (i.e. subtrees having number of tips > \code{f}).
#' @return \strong{best.model} best Brownian rate variation ("non-censored") model selected by means of likelihood ratio test.
#' @return \strong{alternative.models} alternative Brownian rate variation ("non-censored") model selected by means of LRT.
#' @return \strong{n.shifts} the number of shifts for the \code{RBC} analysis.
#' @details The user indicates the minimum size of the clades to be tested for rate variation must be (by setting the argument \code{'f'}, which is by default set at one-tenth of the tree size). Individual nodes are arranged according to their rates (i.e. in descending Brownian rate value). Then, the user is left with two different options to locate the number of potential shifts. First (\code{n.shift = "fix"}), it is possible to specify the number \emph{n} (=\code{'NS'}) of shifts to be searched for all combinations of the \emph{n} nodes with the \emph{n} largest sigma2 (Brownian rate) value, with size 1 to \emph{n}. For instance, with \code{'NS'} = 3 \code{RRphylo} will search through all the eight possible combinations of the 3 nodes with the largest sigma2 values (three combinations with one shift only, one for each node; three combinations of two shifts at two different nodes; and a single combination including all the three shifts for all \code{'NS'}=3 nodes, plus Brownian motion, which means no shift applied). Alternatively (\code{n.shift = "clust"}), all nodes subtending a number of tips as large as \code{f} are partitioned in groups according to their patristic distance, and the number of distinct groups with potential shifts is established via bootstrapped cluster analysis of the internodes distances, by using \pkg{pvclust} package in R (\cite{Suzuki and Shimodaira 2015}).  This way the number of potential shifts are located in topologically distinct parts of the tree. The resulting number of groups \emph{k} is thus taken to be equivalent to the number of shifts to be searched, by examining all possible combinations of the \emph{k} nodes with the largest simga2 values. Of course, it is still possible (and in fact tested) that more than one shift falls in the same region of the tree. Once potential shifts are located, their combinations represent different rate variation models, which are compared to each other (and to a single rate, Brownian motion model) by means of restricted maximum likelihood fitted with the function \code{brownieREML} in \pkg{phytools} (\cite{Revell 2012}), or \code{mvBM} in \pkg{mvMORPH} (\cite{Clavel et al. 2015}) in the multivariate case. The likelihoods of individual models are contrasted to each other to find the best model by means of likelihood ratio test.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @references
#' Piras, P., Silvestro, D., Carotenuto, F., Castiglione, S., Kotsakis, A., Maiorino, L., Melchionna, M., Mondanaro, A., Sansalone, G., Serio, C., Vero, V.A., & Raia, P. (2018). Evolution of the sabertooth mandible: A deadly ecomorphological specialization. \emph{Palaeogeography, Palaeoclimatology, Palaeoecology}, in press
#' @references
#' O'Meara, B. C., Ané, C., Sanderson, M. J., & Wainwright, P. C. (2006). Testing for different rates of continuous trait evolution using likelihood. \emph{Evolution}, 60: 922-933.
#' @references
#' Suzuki, R., & Shimodaira, H. (2015). pvclust: Hierarchical Clustering with P-Values via Multiscale Bootstrap Resampling. R package version 2.0-0. https://CRAN.R-project.org/package=pvclust
#' @references
#' Revell, L. J. (2012). phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods in Ecology and Evolution}, 3: 217-223.doi:10.1111/j.2041-210X.2011.00169.x
#' @references
#' Clavel, J., Escarguel, G., & Merceron, G.(2015). mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. \emph{Methods in Ecology and Evolution}, 6: 1311-1319. doi: 10.1111/2041-210X.12420
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#'
#'  \donttest{
#'
#'     RRphylo(tree=treedino,y=massdino)->RR
#'
#' # Case 1. RBC fixing the number of shifts at 2
#'     RBC(RR=RR,y=massdino,n.shift="fix",NS=2,foldername=tempdir())
#'
#' # Case 2. RBC automatically searching for shifts
#'     RBC(RR=RR,y=massdino,n.shift="clust",foldername=tempdir())
#'
#'     }


RBC<-function(RR,y,n.shift=c("clust","fix"),NS=3,clus=.5,f=NULL,foldername)
{
  #require(ape)
  #require(phytools)
  #require(geiger)
  #require(foreach)
  #require(doParallel)
  #require(binr)
  #require(pvclust)
  #require(mvMORPH)
  #require(parallel)


  RR$tree->t
  if(is.null(f)) f<-Ntip(t)/10
  RR$rates->rates
  RR$tip.path->L

  if(class(y)=="data.frame") treedata(t,y,sort=TRUE)[[2]]->y

  if (length(y) > Ntip(t)) {
    multi.rate <-rates[,1]
    names(multi.rate)<-rownames(rates)
    multi.rate.n <- multi.rate[-which(names(multi.rate) %in% t$tip.label)]
    n.tip <- list()
    for (i in 1:length(multi.rate.n)) {
      n.tip[[i]] <- tips(t, as.numeric(names(multi.rate.n)[i]))
    }
    names(n.tip) <- names(multi.rate.n)
    multi.rate.sel <- multi.rate.n[which(names(multi.rate.n) %in% names(which(lapply(n.tip, length) > round(f))))]
    multi.rate.sel <- multi.rate.sel[-which(names(multi.rate.sel) ==
                                              Ntip(t) + 1)]
    multi.rate.sel <- sort(multi.rate.sel, decreasing = TRUE)
    match.arg(n.shift)
    if (n.shift == "clust") {
      d <- dist.nodes(t)[names(multi.rate.sel), names(multi.rate.sel)]
      d2 <- pvclust(d, parallel = TRUE)
      n.shift <- length(pvpick(d2, alpha = 0.99)$edges) -
        1
      if (n.shift < 2)
        n.shift <- 1
    } else {
      NS <- NS
      n.shift = NS
    }
    YT<-array()
    clade.rateACE <- list()
    for (i in 1:length(multi.rate.sel)) {
      y.tips <- tips(t, as.numeric(names(multi.rate.sel)[i]))
      length(y.tips)->YT[i]
      y.tippa <- y[which(rownames(y) %in% y.tips),
                   ]
      t.tippa <- drop.tip(t, rownames(y)[-which(rownames(y) %in%
                                                  y.tips) ])
      clade.rateACE[[i]] <- diag(ratematrix(t.tippa,y.tippa))

    }
    if(length(which(YT==2))>0) clade.rateACE[-which(YT==2)]->clade.rateACE
    if(length(which(YT==2))>0) multi.rate.sel[-which(YT==2)]->multi.rate.sel


    clade.rateACE <- do.call(rbind, clade.rateACE)
    rownames(clade.rateACE) <- names(multi.rate.sel)
    rbc <- clade.rateACE
    rbc[order(apply(rbc,1,mean),decreasing=TRUE),]->rbc
    rate.sel <- rbc[1:n.shift,]
    if(n.shift==1){
      sum(rate.sel)->rate.sel
      names(rate.sel)<-rownames(rbc)[1]
    } else {

      apply(rate.sel,1,sum)->rate.sel
    }


    st <- rep("A", Ntip(t))
    names(st) <- t$tip.label
    M.co <- list()
    for (i in 1:length(rate.sel)) {
      co <- combn(names(rate.sel), i)
      M <- list()
      for (j in 1:dim(co)[2]) {
        nodes <- as.numeric(co[, j])
        st <- rep("A", Ntip(t))
        names(st) <- t$tip.label
        data.frame(nodeHeights(t)[match(nodes,t$edge[,2]),2],nodes)->ordering
        names(ordering)<-c("H","node")
        ordering[order(ordering$H),][,2]->nodes

        for (l in 1:length(nodes)) {
          st1 <- st
          st1[match(tips(t, nodes[l]), names(st1))] <- letters[l +
                                                                 1]
          st <- st1
        }
        tsy <- make.simmap(t, st1)
        M[[j]] <- mvBM(tsy, y, model = "BMM", optimization = "L-BFGS-B",
                       method = "rpf")
      }
      if (max(apply(co, 2, length)) == 1) {
        names(M) <- co
      } else {
        names(M) <- apply(co, 2, function(x) paste(x,
                                                   collapse = "_"))
      }
      if (i == 1) M.co <- M else M.co <- append(M.co, M)
    }
    mod.liks <- as.data.frame(unlist(lapply(M.co, "[[",
                                            1)))
    mod.aicc <- as.data.frame(unlist(lapply(M.co, "[[",
                                            3)))
    AICcs <- data.frame(model = rownames(mod.aicc), AICc = mod.aicc[,
                                                                    1])
    bestAICc <- AICcs[order(AICcs[, 2]), ][1:n.shift,
                                           ]
    if (nrow(bestAICc) == 1) {
      selected.models <- bestAICc
      best.models <- "no.alternative"
    } else {
      if (nrow(bestAICc) == 2) {
        iss <- LRT(M.co[[match(bestAICc[1, 1], names(M.co))]], M.co[[match(bestAICc[2, 1], names(M.co))]])

        comp <- iss$p
      } else {
        compare <- combn(bestAICc[, 1], 2)
        compare <- compare[, which(compare[1, ] ==
                                     bestAICc[1, 1])]
        comp <- array()
        for (i in 1:ncol(compare)) {
          iss <- LRT(M.co[[match(compare[1, i], names(M.co))]],
                     M.co[[match(compare[2, i], names(M.co))]])
          comp[i] <- iss$p
        }
      }
      best <- data.frame(bestAICc, mod.likelihoods = mod.liks[match(bestAICc[,1], rownames(mod.liks)),],lrt.test = c(1,comp))


      selected.models <- best[which(best[, 4] < 0.05),
                              ]
      best.models <- best[which(best[, 4] > 0.05),
                          ]
      if (dim(selected.models)[1] == 0)
        selected.models <- best[1, ]
      if (dim(best.models)[1] == 0)
        best.models <- "no.alternative"
    }
    rownames(selected.models) <- selected.models[, 1]
  } else {
    rates[,1]->rat
    names(rat)<-rownames(rates)
    rat->rates
    y.hat <- L %*% rates
    ic <- ratematrix(t, y.hat)
    internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                                2] > Ntip(t))]))
    clade.rateACE <- array()
    for (i in 1:length(internals)) {
      tl <- length(tips(t, internals[i]))
      if (tl < 10) {
        clade.rateACE[i] <- 0
      } else {
        y.clade <- y[names(y) %in% tips(t, internals[i])]
        t.clade <- subtrees(t)[[i]]
        clade.rateACE[i] <- ratematrix(t.clade, y.clade)
      }
    }
    mR <- data.frame(internals, clade.rateACE)
    mR <- mR[-which(mR[, 2] == 0), ]
    mR[, 2] <- as.numeric(as.character(mR[, 2]))
    mR <- mR[order(abs(mR[, 2]), decreasing = TRUE),]
    CS <- array()
    for (i in 1:dim(mR)[1]) {
      CS[i] <- length(tips(t, mR[i, 1]))
    }
    if (length(which(CS > round(f))) == dim(mR)[1]) sel <- mR else sel <- mR[which(CS > round(f)), ]
    rbc <- mR
    shift.sel <- sel
    match.arg(n.shift)
    if (n.shift == "clust") {
      d <- dist.nodes(t)[shift.sel[, 1], shift.sel[,1]]
      if(length(d)==1){
        n.shift <-1
      }else{
        d2 <- pvclust(d, parallel = TRUE)
        n.shift <- length(pvpick(d2, alpha = 0.99)$edges) - 1
      }
      if (n.shift < 2) n.shift <- 1
    } else {
      NS <- NS
      n.shift = NS
    }

    if ((Ntip(t) + 1) %in% shift.sel$int) {
      shift.sel <- shift.sel[-match((Ntip(t) + 1), shift.sel$int), ]
      data.frame(nodeHeights(t)[match(shift.sel[,1],t$edge[,2]),][,2],shift.sel[,1])->ordering
      names(ordering)<-c("H","node")
      ordering[order(ordering$H),][,2]->ordered
      shift.sel <- shift.sel[match(ordered, shift.sel[,1]), ]
    } else {
      data.frame(nodeHeights(t)[match(shift.sel[,1],t$edge[,2]),][,2],shift.sel[,1])->ordering
      names(ordering)<-c("H","node")
      ordering[order(ordering$H),][,2]->ordered
      shift.sel <- shift.sel[match(ordered, shift.sel[,1]), ]

    }
    if(dim(shift.sel)[1]==0)
    {
      rbc <- "no.rbc"
      selected.models <- cbind("brown","no.rbc.estimates")
      rownames(selected.models)<-"brown_brown"
      best.models <- cbind("brown","no.rbc.estimates")
      rownames(best.models)<-"brown_brown"

      warning("no potential shift located")
    } else {
      shi <- seq(1, n.shift)
      cl <- makeCluster(round((detectCores() * clus), 0))
      registerDoParallel(cl)
      gc()
      shi.LL <- list()
      for (k in 1:length(shi)) {
        combinations <- combn(shift.sel$int[1:(2 * n.shift)], shi[k])

        LL.matrix <- list()
        LL.matrix <- foreach(j = 1:dim(combinations)[2]) %dopar%
        {
          node <- combinations[, j]
          hh <- phytools::nodeHeights(t)[which(t$edge[,
                                                      2] %in% node)]
          names(hh) <- node
          node <- node[order(hh, decreasing = TRUE)]
          rs <- rep("r1", ape::Ntip(t))
          reg <- list()
          for (i in 1:length(node)) {
            reg[[i]] <- geiger::tips(t, node[i])
          }
          for (i in 1:length(node)) {
            rs[which(names(y) %in% reg[[i]])] <- paste("r",
                                                       node[i], sep = "")
          }
          names(rs) <- t$tip.label
          tryCatch(phy.simm <- phytools::make.simmap(t,
                                                     rs), error = function(e) NULL)
          tryCatch(t(phytools::brownieREML(phy.simm,
                                           y)[c(4, 3, 7, 8)]), error = function(e) NULL)
        }
        names(LL.matrix) <- apply(combinations, 2, function(x) paste(x,
                                                                     collapse = "_"))
        shi.LL[[k]] <- LL.matrix
      }
      stopCluster(cl)
      shift.likelihoods <- list()
      for (i in 1:length(shi.LL)) {
        shi.LL[[i]] <- Filter(Negate(is.null), shi.LL[[i]])
        shi.ll <- do.call(rbind, shi.LL[[i]])
        rownames(shi.ll) <- names(shi.LL[[i]])
        shift.likelihoods[[i]] <- shi.ll
      }
      shift.likelihoods <- do.call(rbind, shift.likelihoods)
      if (dim(shift.likelihoods)[2] == 1) {
        best.models <- "REML failed"
        selected.models <- "REML failed"
      } else {
        if (length(which(shift.likelihoods[, 3] == "NULL")) > 0)
          shift.likelihoods <- shift.likelihoods[-which(shift.likelihoods[,3] == "NULL"), ]
        else shift.likelihoods <- shift.likelihoods
        shift.likelihoods <- shift.likelihoods[order(as.numeric(as.character(shift.likelihoods[,3])), decreasing = TRUE),]
        shift.likelihoods <- as.data.frame(shift.likelihoods)
        restrict.s <- shift.likelihoods[shift.likelihoods[,3] > shift.likelihoods[[1, 3]] - 4, ]
        "brown" <- restrict.s[1, c(2, 1)]
        brown <- data.frame(brown)
        names(brown) <- colnames(restrict.s)[3:4]
        rownames(brown) <- "brown"
        restrict.s <- rbind(brown, restrict.s[, 3:4])
        restrict.s[, 1] <- as.numeric(as.character(restrict.s[,1]))
        restrict.s[, 2] <- as.numeric(as.character(restrict.s[,2]))
        n.par <- array()
        for (i in 1:dim(restrict.s)[1]) {
          if (rownames(restrict.s)[i] == "brown") {
            n.par[i] <- 2
          } else {
            n.par[i] <- length(strsplit(rownames(restrict.s[i,]), split = "_")[[1]]) + 2
          }
        }
        restrict.s <- data.frame(restrict.s, n.par)
        best <- restrict.s[which.max(restrict.s[, 1]),
                           ]
        res.s <- restrict.s[-which.max(restrict.s[, 1]),
                            ]
        pX2 <- array()
        for (i in 1:nrow(res.s)) {
          pX2[i] <- 1 - pchisq((best[1, 1] - res.s[i,
                                                   1]) * 2, abs(best[1, 3] - res.s[i, 3]))
        }
        res.s <- data.frame(res.s, pX2)
        res.s <- rbind(data.frame(best, pX2 = 0), res.s)
        best.models <- rbind(res.s[1, ], res.s[res.s$p >
                                                 0.05, ])
        selected.models <- best.models[which(best.models$n.par ==
                                               min(best.models$n.par)), ]
        LogL.brownian <- res.s[which(rownames(res.s) ==
                                       "brown"), ][, 1]
        selected.models <- data.frame(selected.models,
                                      rep(LogL.brownian, dim(selected.models)[1]))
        colnames(selected.models)[5] <- "LogL.brownian"
      }
    }
  }
  node <- list()
  for (j in 1:nrow(selected.models)) {
    node[[j]] <- strsplit(rownames(selected.models)[j], "_")[[1]]
  }
  for (j in 1:length(node)) {
    if (length(y) > Ntip(t)) nod <- node[[1]] else nod <- node[[j]]
    if ("brown" %in% nod) {
      pdf(paste(paste(foldername,rownames(selected.models)[j],sep="/"),".pdf"))
      plot(t)
      nodelabels(bg = "white", cex = 0.5)
      dev.off()
    } else {
      rs <- rep("r1", Ntip(t))
      reg <- list()
      for (i in 1:length(nod)) {
        reg[[i]] <- tips(t, as.numeric(as.character(nod[[i]])))
        if (length(y) > Ntip(t))
          rs[which(rownames(y) %in% reg[[i]])] <- paste("r",
                                                        nod[i], sep = "")
        else rs[which(names(y) %in% reg[[i]])] <- paste("r",
                                                        nod[i], sep = "")
        names(rs) <- t$tip.label
      }
      phy.simm <- make.simmap(t, rs)
      if (length(y) > Ntip(t)) {
        pdf(paste(paste(foldername,as.character(selected.models[1, 1]),sep="/"),
                  ".pdf"))
        plotSimmap(phy.simm)
        dev.off()
      } else {
        sim <- brownieREML(phy.simm, y)
        ss <- sprintf("%4.3f", sim$sig2.multiple)
        pdf(paste(paste(foldername,rownames(selected.models)[j],sep="/"), ".pdf"))
        plotSimmap(phy.simm)
        nodelabels(node = c(Ntip(t) + 1, as.numeric(unlist(nod))),
                   text = paste(ss), frame = "none", cex = 1.5)
        dev.off()
      }
    }
  }
  res <- list(rbc, selected.models, best.models,
              n.shift)
  names(res) <- c("rbc", "best.model", "alternative.models",
                  "n.shifts")
  return(res)
}
