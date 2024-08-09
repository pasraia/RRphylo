#' @title Test for phylogenetic clustering
#' @description The function tests the presence of phylogenetic clustering for
#'   species within a focal state.
#' @usage phyloclust(tree,state,focal,nsim=100)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric or
#'   fully dichotomous.
#' @param state the named vector of tip states.
#' @param focal the focal state to be tested for phylogenetic clustering.
#' @param nsim number of simulations to perform the phylogenetic clustering
#'   test.
#' @export
#' @return The function returns a list including the p-value ($p) for the test
#'   of phylogenetic clustering and a $declusterized object containing the
#'   declusterized versions of the original tree and state vector (i.e. tips are
#'   removed as to make p>0.05) and the vector of removed species.
#' @details To test for phylogenetic clustering, the function computes the mean
#'   cophenetic (i.e. evolutionary time) distance between all the species under
#'   the \code{focal} state. Such value is compared to a random distribution of
#'   time distances obtained by sampling \code{nsim} times as many random tips
#'   as those under the \code{focal} state. In the presence of significant
#'   phylogenetic clustering, tips under the \code{focal} state are randomly
#'   removed until the \code{p} becomes >0.05 or only 3 tips are left.
#' @author Silvia Castiglione, Pasquale Raia
#' @examples
#' data("DataFelids")
#' DataFelids$treefel->treefel
#' DataFelids$statefel->statefel
#'
#' phyloclust(tree=treefel,state=statefel,focal="saber")

phyloclust<-function(tree,state,focal,nsim=100){
  # if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
  #   data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
  #   tree$tip.label<-tips(tree,(Ntip(tree)+1))
  #   data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
  #   tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  # }

  # state <- treedata(tree, state, sort = TRUE)[[2]][,1]
  state <- treedataMatch(tree, state)[[1]][,1]
  focal->st
  cophenetic.phylo(tree)->cop
  cop[which(state==st),which(state==st)]->subcop
  mean(subcop[upper.tri(subcop)])->mds
  dim(subcop)[1]->sl

  r.mds<-array()
  for(e in 1:nsim){
    sample(tree$tip.label,sl)->test.tip
    cop[test.tip,test.tip]->r.cop
    mean(r.cop[upper.tri(r.cop)])->r.mds[e]
  }
  pval=length(which(r.mds<mds))/nsim

  remT<-c()
  length(which(state==st))->lenst
  while(pval<0.05){
    names(state)->nam
    as.numeric(as.factor(state))->rt
    data.frame(state,rt)->def
    unique(def[which(def[,1]==st),2])->stN
    data.frame(V=rle(rt)$values,L=rle(rt)$lengths)->VL
    data.frame(VL,pos=cumsum(VL[,2]))->VL
    VL[which(VL$V==stN),]->vel
    vel[which.max(vel$L),]->hit
    tree$tip.label[(hit[,3]-hit[,2]+1):hit[,3]]->hitips
    sample(hitips,ceiling(0.5*length(hitips)))->remtips
    drop.tip(tree,remtips)->tree
    state[-which(names(state)%in%remtips)]->state
    if(length(c(remtips,remT))>(lenst-3)) break else c(remtips,remT)->remT
  }
  if(length(remT)==0) remT<-NA
  return(list(p=pval,declusterized=list(tree=tree,state=state,removed.tips=remT)))
}
