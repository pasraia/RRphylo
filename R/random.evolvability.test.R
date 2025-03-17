#' @title Randomization test for phylogenetic structuring in evolvability
#' @description The function is a wrapper around the function
#'   \code{\link[evolqg]{MeanMatrixStatistics}} from the package \pkg{evolqg} (\cite{Melo et
#'   al. 2015}). It estimates ancestral character at internal nodes either
#'   according to Brownian Motion or by means of \code{\link{RRphylo}} (see the
#'   argument \code{node.estimation}), then performs \code{MeanMatrixStatistics}
#'   to calculate: Mean Squared Correlation, ICV, Autonomy,
#'   ConditionalEvolvability, Constraints, Evolvability, Flexibility,
#'   Pc1Percent, and Respondability. To assess the importance of phylogenetic
#'   structuring (signal) on Respondability Evolvability, and Flexibility. The
#'   function performs a randomization test by randomly shuffling the species on
#'   tree and replicating the analyses \code{nsim} times. A p-value is computed
#'   by contrasting the real metrics to the ones derived by randomization.
#' @usage random.evolvability.test(tree,data,node.estimation=c("RR","BM"),
#'   aces=NULL,iterations=1000,nsim=100,clus=0.5)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous.
#' @param data a matrix or data.frame of phenotypic data having species as rownames
#' @param node.estimation specify the method to compute ancestral character at
#'   nodes. It can be one of \code{"RR"}, to compute ancestral states by mean of
#'   \code{RRphylo}, or \code{"BM"}, to use \pkg{phytools}' function \code{\link[phytools]{fastAnc}}
#'   (\cite{Paradis & Schliep 2019}) to estimate ancestral characters at nodes
#'   according to Brownian Motion.
#' @param aces a named matrix of known ancestral character values at nodes.
#'   Names correspond to the nodes in the tree.
#' @param iterations the iterations argument to be indicated in
#'   \code{MeanMatrixStatistics}
#' @param nsim the number of simulations to be performed for the randomization
#'   test, by default \code{nrep} is set at 100.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{random.evolvability.test} set \code{clus} = 0.
#' @export
#' @seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#' @importFrom stats cov
#' @importFrom phytools fastAnc
#' @return The function returns a list object including (\code{$means}) the mean
#'   values for all statistics as produced by \code{MeanMatrixStatistics} and
#'   (\code{$means}) the significance levels for Respondability, Evolvability,
#'   and Flexibility.
#'   The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @references Melo, D., Garcia, G., Hubbe, A., Assis, A. P., & Marroig, G.
#'   (2015). EvolQG-An R package for evolutionary quantitative genetics.
#'   \emph{F1000Research}, 4.
#' @references Revell, L. J. (2012) phytools: An R package for phylogenetic
#'   comparative biology (and other things).
#'   \emph{Methods in Ecology and Evolution}, 3, 217-223.
#' @author Silvia Castiglione, Gabriele Sansalone, Pasquale Raia
#' @examples
#'  \dontrun{
#'  library(ape)
#'  library(phytools)
#'
#'  rtree(30)->phy
#'  fastBM(phy,nsim=4)->phen
#'
#'  random.evolvability.test(tree=phy,data=phen,node.estimation="RR")->rEvTest
#'
#'     }

random.evolvability.test<-function(tree,data,node.estimation=c("RR","BM"),aces=NULL,iterations=1000,nsim=100,clus=0.5){
  # require(doParallel)
  # require(evolqg)
  # require(ape)
  # require(phytools)
  # require(RRphylo)
  # require(ddpcr)

  misspacks<-sapply(c("evolqg","ddpcr"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  funcall <- match.call()
  if(!is.binary(tree)){
    tree->mutree
    multi2di(tree,random=FALSE)->tree
    tree$edge[which(tree$edge.length==0),2]->len0
    lapply(names(which(table(mutree$edge[,1])>2)),function(x) tips(mutree,as.numeric(x)))->specs
    lapply(specs,function(w) getMRCA(tree,w))->newn
    unlist(lapply(newn,function(w)getDescendants(tree,w)))->des
    if(any(len0%in%des)) len0[which(len0%in%des)]->remn
  } else remn<-NULL

  # treedata(tree,data)[[2]]->data
  treedataMatch(tree,data)[[1]]->data

  if(node.estimation=="RR") rec<-RRphylo(tree,data,aces=aces,clus=clus)$aces else{
    rec<-NULL
    for(i in 1:ncol(data)){
      if(!is.null(aces)){
        aces[,i]->ancs
        names(ancs)<-rownames(aces)
      }else ancs<-NULL
      reci<-fastAnc(tree,data[,i],anc.states=ancs,CI=T)$ace
      rec<-cbind(rec,reci)
    }
  }
  if(!is.null(remn)) rec<-rec[-which(rownames(rec)%in%remn),]
  cv<-cov(t(rec))

  if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
  registerDoParallel(cl)
  ddpcr::quiet(xx<-evolqg::MeanMatrixStatistics(cv,iterations=iterations,full.results = TRUE,parallel = TRUE))
  xx$mean->means

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  random.means<-list()
  j=1
  repeat({
    setTxtProgressBar(pb,j)
    data->dataR
    rownames(dataR)<-sample(rownames(data))
    # dataR<-treedata(tree,dataR)[[2]]
    dataR<-treedataMatch(tree,dataR)[[1]]

    if(node.estimation=="RR") rec<-RRphylo(tree,dataR,aces=aces,clus=clus)$aces else{
      rec<-NULL
      for(i in 1:ncol(dataR)){
        if(!is.null(aces)){
          aces[,i]->ancs
          names(ancs)<-rownames(aces)
        }else ancs<-NULL
        reci<-fastAnc(tree,dataR[,i],anc.states=ancs,CI=T)$ace
        rec<-cbind(rec,reci)
      }

    }
    if(!is.null(remn)) rec<-rec[-which(rownames(rec)%in%remn),]
    cv<-cov(t(rec))

    if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1, setup_strategy = "sequential") else cl <- makeCluster(round((detectCores() * clus), 0), setup_strategy = "sequential")
    registerDoParallel(cl)
    ddpcr::quiet(try(xx<-evolqg::MeanMatrixStatistics(cv,iterations=iterations,full.results = T,parallel = TRUE),silent = TRUE)->trytest)
    if(inherits(trytest,"try-error")) j=j else{
      xx$mean->random.means[[j]]
      if(j==(nsim-1)) break else j=j+1
    }
  })

  do.call(cbind,random.means)->random.means

  cbind(means[c("respondability","evolvability","flexibility")],
        random.means[c("respondability","evolvability","flexibility"),])->totmeans
  apply(totmeans,1,function(w) rank(w)[1]/nsim)->p
  res<-list(means=means,p.value=p)
  attr(res,"Call")<-funcall

  return(res)
}

