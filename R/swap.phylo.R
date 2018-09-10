#' @title Test the effect of phylogenetic uncertainty on rate shifts found at a particular node
#' @description The function uses a number of alternative phylogenies with altered (as compared to the reference tree) topology and branch lengths tests whether the tips descending from the specified node (\code{'node'}) have statistically different rates from the rest of the tree. A phenotypic vector \code{'y'} must be supplied. Eventually, the effect of a covariate could be included.
#' @usage swap.phylo(tree,si=0.5,si2=0.5,node,y,rts,nrep=100,cov=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or fully dichotomous.
#' @param si the proportion of tips whose topologic arrangement will be swapped.
#' @param si2 the proportion of nodes whose age will be changed.
#' @param node the focal node to be tested.
#' @param y the phenotype under testing.
#' @param rts the rates found by \code{\link{RRphylo}} on the original tree.
#' @param nrep the number of simulated trees to be produced.
#' @param cov the covariate to be indicated if its effect on rate values must be accounted for. Contrary to \code{RRphylo}, \code{'cov'} needs to be as long as the number of tips of the tree.
#' @importFrom graphics hist par
#' @importFrom grDevices rgb
#' @importFrom stats t.test
#' @details \code{swap.phylo} changes the tree topology and branch lengths to a level specified by the user. Up to half of the tips, and half of the branch lengths can be changed randomly. The function provides a 'swapped' tree, yet, importantly, once a shift in the rate of evolution has been found by \code{\link{RRphylo}}, this function can be used to test whether the shift depends on the tree topology and branch lengths. It runs \code{RRphylo} on swapped trees (default is 100) and then calculates the absolute rate difference between all the branches of the shifted node and the rest of the tree. A t-test is eventually performed to assess significance.
#' @export
#' @return The function returns a 'list' object containing:
#' @return \strong{$p.swap} the probability that the rates at \code{'node'} are different from rates at the rest of the tree.
#' @return \strong{$rates} the distribution of rates per branch as calculated by \code{RRphylo} on 'swapped' phylogenies.
#' @examples
#'   data("DataApes")
#'   DataApes$PCstage->PCstage
#'   DataApes$Tstage->Tstage
#'   DataApes$CentroidSize->CS
#'
#' \donttest{
#' # Case 1. swap.phylo without accounting for the effect of a covariate
#'   RRphylo(tree=Tstage,y=PCstage)->RR
#'   RR$rates->rr
#'   swap.phylo(Tstage,node=61,y=PCstage,rts=rr)

#' # Case 2. swap.phylo accounting for the effect of a covariate
#'   RRphylo(tree=Tstage,y=CS)->RRcova
#'   c(RRcova$aces,CS)->cov.values
#'   c(rownames(RRcova$aces),names(CS))->names(cov.values)
#'   RRphylo(tree=Tstage,y=PCstage,cov=cov.values)->RR
#'   RR$rates->rr
#'   swap.phylo(Tstage,node=61,y=PCstage,rts=rr,cov=CS)
#' }




swap.phylo<-
  function(tree,si=0.5,si2=0.5,node,y,rts,nrep=100,cov=NULL)


  {


    maxN <- function(x, N=2){
      len <- length(x)
      if(N>len){
        warning('N greater than length(x).  Setting N=length(x)')
        N <- length(x)
      }
      sort(x,partial=len-N+1)[len-N+1]
    }

    tree->tree1
    if(class(y)=="data.frame") treedata(tree,y,sort=TRUE)[[2]]->y
    Rrts<-list()
    if(is.null(cov))
    {
      replicate(nrep,RRphylo(swapONE(tree,si,si2)[[1]],y)$rates)->Rrts
    }else{
      RRphylo(tree,cov)->RRcov
      c(RRcov$aces,cov)->cova
      names(cova)<-c(rownames(RRcov$aces),names(cov))
      replicate(nrep,RRphylo(swapONE(tree,si,si2)[[1]],y,cov=cova)$rates)->Rrts
    }

    for(i in 1:dim(Rrts)[3]) Rrts[,,i][match(rownames(rts),names(Rrts[,,i]))]

    apply(Rrts,3,cbind)->multirates
    rownames(multirates)<-rownames(rts)
    apply(multirates,1,function(x) mean(abs(x)))->Rrts.distrib
    names(Rrts.distrib)<-rownames(rts)

    if(length(y)>Ntip(tree)){
      rts[match(rownames(y),rownames(rts)),]->rts
      Rrts.distrib[match(rownames(y),names(Rrts.distrib))]->Rrts.distrib
      multirates[match(rownames(y),rownames(multirates)),]->multirates
    }else{
      rts[match(names(y),rownames(rts)),]->rts
      Rrts.distrib[match(names(y),names(Rrts.distrib))]->Rrts.distrib
      multirates[match(names(y),rownames(multirates)),]->multirates
    }

    hist(log(abs(rts[match(tips(tree,node),names(rts))])))->H3
    hist(log(abs(rts[-match(tips(tree,node),names(rts))])))->H4

    hist(log(abs(Rrts.distrib[match(tips(tree,node),names(Rrts.distrib))])))->H1
    hist(log(abs(Rrts.distrib[-match(tips(tree,node),names(Rrts.distrib))])))->H2

    Rrts[,1,]->rrts
    par(mfrow=c(1,2))

    plot(H3,xlim=c(min(H4$bre)-2,max(H3$bre)+2),col=rgb(1,0,0,.5),ylim= c(0,max(H4$counts)),main="rates density distribution",xlab="rates")
    plot(H4,col=rgb(0,0,1,.5),add=TRUE)
    plot(H1,xlim=c(min(H2$bre)-2,max(H1$bre)+2),col=rgb(1,0,1,.5),ylim= c(0,max(H2$counts)),main="rates density distribution",xlab="rates")
    plot(H2,col=rgb(0,1,0,.5),add=TRUE)
    t.test(log(abs(Rrts.distrib[-match(tips(tree,node),names(Rrts.distrib))])),log(abs(Rrts.distrib[match(tips(tree,node),names(Rrts.distrib))])))->p.swap
    return(list("p.swap"=p.swap,"rates"=multirates))

  }
