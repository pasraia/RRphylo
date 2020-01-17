#' @title Producing simulated phenotypes with trends
#'
#' @description The function \code{setBM} is wrapper around \pkg{phytools} \code{fastBM} function, which generates BM simulated phenotypes with or without a trend.
#' @usage setBM(tree, nY = 1, s2 = 1, a = 0, type = c("", "brown","trend", "drift"),
#' trend.type = c("linear", "stepwise"),tr = 10, t.shift = 0.5, es=2, ds=1)
#' @param tree a phylogenetic tree.
#' @param nY the number of phenotypes to simulate.
#' @param s2 value of the Browian rate to use in the simulations.
#' @param a the phenotype at the tree root.
#' @param type the type of phenotype to simulate. With the option \code{"brown"} the phenotype will have no trend in the phenotypic mean or in the rate of evolution. A variation in the phenotypic mean over time (a phenotypic trend) is obtained by selecting the option \code{"drift"}. A trend in the rate of evolution produces an increased variance in the residuals over time. This is obtained by specifying the option \code{"trend"}.
#' @param trend.type two kinds of heteroscedastic residuals are generated under the \code{"trend"} type. The option \code{"linear"} produces an exponential linear increase (or decrease) in heteroscedasticity, whereas the \code{"stepwise"} option produces an increase (or decrease) after a specified point in time.
#' @param tr the intensity of the trend with the \code{"stepwise"} option is controlled by the \code{tr} argument. The scalar \code{tr} is the multiplier of the branches extending after the shift point as indicated by \code{t.shift}.
#' @param t.shift the relative time distance from the tree root where the stepwise change in the rate of evolution is indicated to apply.
#' @param es when \code{trend.type="linear"}, \code{es} is a scalar representing the exponent at which the evolutionary time (i.e. distance from the root) scales to change to phenotypic variance. With \code{es} = 1 the phenotypic rate will be trendless, with \code{es} < 1 the variance of the phenotypes will decrease exponentially towards the present and the other way around with \code{es} > 1.
#' @param ds a scalar indicating the change in phenotypic mean in the unit time, in \code{type="drift"} case. With \code{ds} = 0 the phenotype will be trendless, with \code{ds} < 0 the phenotypic mean will decrease exponentially towards the present and the other way around with \code{ds} > 0.
#' @details Note that \code{setBM} differs from \code{fastBM} in that the produced phenotypes are checked for the existence of a temporal trend in the phenotype. The user may specify whether she wants trendless data (option \code{"brown"}), phenotypes trending in time (option \code{"drift"}), or phenotypes whose variance increases/decreases exponentially over time, consistently with the existence of a trend in the rate of evolution (option \code{"trend"}). In the latter case, the user may indicate the intensity of the trend (by applying different values of \code{es}), and whether it should occur after a given proportion of the tree height (hence a given point back in time, specified by the argument \code{t.shift}). Trees in \code{setBM} are treated as non ultrametric. If an ultrametric tree is fed to the function, \code{setBM} alters slightly the leaf lengths multiplying randomly half of the leaves by 1 * 10e-3,in order to make it non-ultrametric.
#' @return Either an object of class \code{'array'} containing a single phenotype or an object of class \code{'matrix'} of \emph{n} phenotypes as columns, where \emph{n} is indicated as \code{nY} = \emph{n}.
#' @export
#' @importFrom phytools fastBM
#' @importFrom ape vcv
#' @importFrom lmtest bptest
#' @importFrom stats coefficients
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#'  \donttest{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#'
#' setBM(tree=treedino, nY= 1, type="brown")
#' setBM(tree=treedino, nY= 1, type="drift", ds=2)
#' setBM(tree=treedino, nY= 1, type="trend", trend.type="linear", es=2)}



setBM<-function (tree, nY = 1, s2 = 1, a = 0,
                 type = c("", "brown","trend", "drift"),
                 trend.type = c("linear", "stepwise"),
                 tr = 10, t.shift = 0.5, es = 2, ds = 1)
{
  if (type == "")
    stop("argument 'type' must be defined")
  switch(type, brown = {
    i = 1
    times <- diag(vcv(tree))
    if (var(times) < 1e-04) {
      tree <- makeFossil(tree, ex = 1.001)
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      if (coefficients(summary(lm(y ~ times)))[8] > 0.05) {
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  }, trend = {
    i = 1
    times <- diag(vcv(tree))
    if (var(times) < 1e-04) {
      tree <- makeFossil(tree, ex = 1.001)
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      if (coefficients(summary(lm(y ~ times)))[8] > 0.05) {
        match.arg(trend.type)
        if (trend.type == "linear") {
          y <- (diag(vcv(tree))^es)/(diag(vcv(tree))) * y
        } else {
          y[match(names(which(times > t.shift * max(nodeHeights(tree)))),
                  names(y))] <- y[match(names(which(times >
                                                      t.shift * max(nodeHeights(tree)))), names(y))] *
            tr
        }
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  }, drift = {
    i = 1
    times <- diag(vcv(tree))
    if (var(times) < 1e-04) {
      tree <- makeFossil(tree, ex = 1.001)
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      if (coefficients(summary(lm(y ~ times)))[8] > 0.05) {
        y <- y + diag(vcv(tree)) * ds
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  })
  return(yy)
}
