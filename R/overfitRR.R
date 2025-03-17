#' @title Testing RRphylo overfit
#' @description Testing the robustness of \code{\link{RRphylo}} results to
#'   sampling effects and phylogenetic uncertainty.
#' @usage overfitRR(RR,y, phylo.list, aces=NULL,x1=NULL, aces.x1=NULL, cov=NULL,
#'   rootV=NULL, clus=0.5, s = NULL, swap.args = NULL, nsim=NULL , trend.args =
#'   NULL, shift.args = NULL, conv.args = NULL, pgls.args = NULL)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param y a named vector of phenotypes.
#' @param phylo.list a list (or \code{multiPhylo}) of alternative topologies
#'   (i.e. having the same species as the original tree arranged differently) to
#'   be tested.
#' @param aces if used to produce the \code{RR} object, the vector of those
#'   ancestral character values at nodes known in advance must be specified.
#'   Names correspond to the nodes in the tree.
#' @param x1 the additional predictor to be specified if the RR object has been
#'   created using an additional predictor (i.e. multiple version of
#'   \code{\link{RRphylo}}). \code{'x1'} vector must be as long as the number of nodes
#'   plus the number of tips of the tree, which can be obtained by running
#'   \code{\link{RRphylo}} on the predictor as well, and taking the vector of ancestral
#'   states and tip values to form the \code{x1}.
#' @param aces.x1 a named vector of ancestral character values at nodes for
#'   \code{x1}. It must be indicated if the RR object has been created using
#'   both \code{aces} and \code{x1}. Names correspond to the nodes in the tree.
#' @param cov if used to produce the \code{RR} object, the covariate must be
#'   specified. As in \code{\link{RRphylo}}, the covariate vector must be as long as
#'   the number of nodes plus the number of tips of the tree, which can be
#'   obtained by running \code{\link{RRphylo}} on the covariate as well, and taking the
#'   vector of ancestral states and tip values to form the covariate.
#' @param rootV if used to produce the \code{RR} object, the phenotypic value at
#'   the tree root must be specified.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{overfitRR} set \code{clus} = 0.
#' @param s,swap.args,nsim are deprecated. Check the function
#'   \code{\link{resampleTree}} to generate alterative phylogenies.
#' @param trend.args is deprecated. Check the function \code{\link{overfitST}}
#'   to test \code{\link{search.trend}} robustness.
#' @param shift.args is deprecated. Check the function \code{\link{overfitSS}}
#'   to test \code{\link{search.shift}}  robustness.
#' @param conv.args is deprecated. Check the function \code{\link{overfitSC}} to
#'   test \code{\link{search.conv}}  robustness.
#' @param pgls.args is deprecated. Check the function \code{\link{overfitPGLS}}
#'   to test \code{\link{PGLS_fossil}}  robustness.
#' @return The function returns a 'RRphyloList' object containing:
#' @return \strong{$RR.list} a 'RRphyloList' including the results of each
#'   \code{\link{RRphylo}} performed within \code{overfitRR}.
#' @return \strong{$root.est} the estimated root value per simulation.
#' @return \strong{$rootCI} the 95\% confidence interval around the root value.
#' @return \strong{$ace.regressions} a 'RRphyloList' including the results of
#'   linear regression between ancestral state estimates before and after the
#'   subsampling.
#' @return The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Silvia Castiglione, Carmela Serio, Giorgia Girardi, Pasquale Raia
#' @details Methods using a large number of parameters risk being overfit. This
#'   usually translates in poor fitting with data and trees other than the those
#'   originally used. With \code{\link{RRphylo}} methods this risk is usually very low.
#'   However, the user can assess how robust the results of \code{\link{RRphylo}} are
#'   by running \code{\link{resampleTree}} and \code{overfitRR}. The former is used to
#'   subsample the tree according to a \code{s} parameter (that is the
#'   proportion of tips to be removed from the tree) and to alter tree topology
#'   by means of \code{\link{swapONE}}. The list of altered topologies is fed to
#'   \code{overfitRR}, which cross-references each tree with the phenotypic data
#'   and performs \code{\link{RRphylo}} on them. Thereby, both the potential for
#'   overfit and phylogenetic uncertainty are accounted for straight away.
#'
#'   Otherwise, a list of alternative phylogenies can be supplied to
#'   \code{overfitRR}. In this case subsampling and swapping arguments are
#'   ignored, and robustness testing is performed on the alternative topologies
#'   as they are.
#' @export
#' @seealso \href{../doc/overfitRR.html}{\code{overfitRR} vignette} ;
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P. (2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M.,
#'   Profico, A., Girardi, G., & Raia, P. (2019a) Simultaneous detection of
#'   macroevolutionary patterns in phenotypic means and rate of change with and
#'   within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14:
#'   e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @examples
#' \dontrun{
#' cc<- 2/parallel::detectCores()
#' library(ape)
#'
#' ## overfitRR routine
#' # load the RRphylo example dataset including Ornithodirans tree and data
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#'
#' # extract Pterosaurs tree and data
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' # peform RRphylo on body mass
#' RRphylo(tree=treeptero,y=log(massptero),clus=cc)->RRptero
#'
#' # generate a list of subsampled and swapped phylogenies to test
#' treeptero.list<-resampleTree(RRptero$tree,s = 0.25,swap.si = 0.1,swap.si2 = 0.1,nsim=10)
#'
#' # test the robustness of RRphylo
#' ofRRptero<-overfitRR(RR = RRptero,y=log(massptero),phylo.list=treeptero.list,clus=cc)
#'
#'
#' ## overfitRR routine on multiple RRphylo
#' # load the RRphylo example dataset including Cetaceans tree and data
#' data("DataCetaceans")
#' DataCetaceans$treecet->treecet
#' DataCetaceans$masscet->masscet
#' DataCetaceans$brainmasscet->brainmasscet
#' DataCetaceans$aceMyst->aceMyst
#'
#' # cross-reference the phylogenetic tree and body and brain mass data. Remove from
#' # both the tree and vector of body sizes the species whose brain size is missing
#' drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#' masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#'
#' # peform RRphylo on the variable (body mass) to be used as additional predictor
#' RRphylo(tree=treecet.multi,y=masscet.multi,clus=cc)->RRmass.multi
#' RRmass.multi$aces[,1]->acemass.multi
#'
#' # create the predictor vector: retrieve the ancestral character estimates
#' # of body size at internal nodes from the RR object ($aces) and collate them
#' # to the vector of species' body sizes to create
#' c(acemass.multi,masscet.multi)->x1.mass
#'
#' # peform RRphylo on brain mass by using body mass as additional predictor
#' RRphylo(tree=treecet.multi,y=brainmasscet,x1=x1.mass,clus=cc)->RRmulti
#'
#' # generate a list of subsampled and swapped phylogenies to test
#' treecet.list<-resampleTree(RRmulti$tree,s = 0.25,swap.si=0.1,swap.si2=0.1,nsim=10)
#'
#' # test the robustness of multiple RRphylo
#' ofRRcet<-overfitRR(RR = RRmulti,y=brainmasscet,phylo.list=treecet.list,clus=cc,x1 =x1.mass)
#' }

overfitRR<-function(RR,y,
                    phylo.list,
                    aces=NULL,x1=NULL,
                    aces.x1=NULL,cov=NULL,
                    rootV=NULL,clus=0.5,
                    s = NULL, swap.args = NULL,nsim=NULL,
                    trend.args = NULL, shift.args = NULL,
                    conv.args = NULL, pgls.args = NULL)
{
  # require(phytools)

  remargs<-c("s","swap.args","trend.args","shift.args","conv.args","pgls.args","nsim")
  missargs<-c(missing(s),missing(swap.args),missing(trend.args),missing(shift.args),missing(conv.args),missing(pgls.args),missing(nsim))
  names(missargs)<-remargs
  if(any(!missargs)){
    warning(paste(paste(names(missargs)[which(!missargs)],collapse=", "),
                  "is/are deprecated.\n Check the 'Testing RRphylo methods overfit' vignette for the new version of overfit- functions."),immediate. = TRUE,call.=FALSE)
  }

  funcall<-match.call()
  RR$tree->tree
  y <- treedataMatch(tree, y)[[1]]
  RR$aces->y.ace
  tree$node.label<-rownames(y.ace)
  nsim<-length(phylo.list)
  phylo.size<-Ntip(phylo.list[[1]])

  pb = txtProgressBar(min = 0, max = nsim, initial = 0)

  rootlist<-acefit<-list()
  RR.list<-list()
  real.s<-array()
  for(k in 1:nsim){
    setTxtProgressBar(pb,k)

    phylo.list[[k]]->treecut
    y[match(treecut$tip.label,rownames(y)),,drop=FALSE]->ycut
    y.ace[which(rownames(y.ace)%in%treecut$node.label),,drop=FALSE]->y.acecut

    1-(Ntip(treecut)/Ntip(tree))->real.s[k]

    if(!is.null(cov)){
      treedataMatch(treecut,cov)$y->covcut
      c(RRphylo(treecut,covcut,clus=clus)$aces[,1],covcut)->covcut
      # cov[match(c(rownames(y.acecut),rownames(ycut)),names(cov))]->covcut
      # names(covcut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else covcut<-NULL

    if(!is.null(x1)) {
      as.matrix(x1)->x1
      treedataMatch(treecut,x1)$y->x1cut
      rbind(RRphylo(treecut,x1cut,clus=clus)$aces,x1cut)->x1cut
      # x1[match(c(rownames(y.acecut),rownames(ycut)),rownames(x1)),,drop=FALSE]->x1cut
      # rownames(x1cut)[1:Nnode(treecut)]<-seq((Ntip(treecut)+1),(Ntip(treecut)+Nnode(treecut)))
    }else x1cut<-NULL

    if(!is.null(aces)){
      if(is.vector(aces)) as.matrix(aces)->aces
      aces->acescut

      dropace<-c()
      for(i in 1:nrow(aces)) {
        if(length(which(tips(tree,rownames(aces)[i])%in%treecut$tip.label))>1){
          getMRCA(treecut,tips(tree,rownames(aces)[i])[which(tips(tree,rownames(aces)[i])%in%treecut$tip.label)])->newN
          if(phylo.size==Ntip(tree)){
            length(tips(treecut,newN))/length(tips(tree,rownames(aces)[i]))->sh.tips
            if(sh.tips<=1.1&sh.tips>=0.9) newN->rownames(acescut)[i] else c(dropace,i)->dropace
          }else newN->rownames(acescut)[i]
          # getMRCA(treecut,tips(tree,rownames(aces)[i])[which(tips(tree,rownames(aces)[i])%in%treecut$tip.label)])->rownames(acescut)[i]
        }else c(dropace,i)->dropace
      }
      if(length(dropace>0)) acescut[-dropace,,drop=FALSE]->acescut
      if(is.null(nrow(acescut))) acescut<-NULL
    }else acescut<-NULL

    if(!is.null(aces.x1)){
      if(is.vector(aces.x1)) as.matrix(aces.x1)->aces.x1
      aces.x1->aces.x1cut
      dropace<-c()
      for(i in 1:nrow(aces.x1)) {
        if(length(which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label))>1){
          getMRCA(treecut,tips(tree,rownames(aces.x1)[i])[which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label)])->newN1
          if(phylo.size==Ntip(tree)){
            length(tips(treecut,newN1))/length(tips(tree,rownames(aces.x1)[i]))->sh.tips
            if(sh.tips<=1.1&sh.tips>=0.9) newN1->rownames(aces.x1cut)[i] else c(dropace,i)->dropace
          }else newN1->rownames(aces.x1cut)[i]
          # getMRCA(treecut,tips(tree,rownames(aces.x1)[i])[which(tips(tree,rownames(aces.x1)[i])%in%treecut$tip.label)])->rownames(aces.x1cut)[i]
        }else c(dropace,i)->dropace
      }
      if(length(dropace>0)) aces.x1cut[-dropace,,drop=FALSE]->aces.x1cut
      if(is.null(nrow(aces.x1cut))) aces.x1cut<-NULL
    }else aces.x1cut<-NULL

    if(!is.null(rootV)) rootV->rootVcut else rootVcut<-NULL

    RRphylo(treecut,ycut,aces=acescut,x1=x1cut,aces.x1=aces.x1cut,cov=covcut,rootV = rootVcut,clus=clus)->RRcut->RR.list[[k]]
    RRcut$aces[1,,drop=FALSE]->rootlist[[k]]
    # summary(lm(y.acecut~RRcut$aces))->acefit[[k]]
    do.call(rbind,lapply(seq(1:ncol(y.acecut)),function(x) summary(lm(y.acecut[,x]~RRcut$aces[,x]))$coef[c(1,2,7,8)]))->acefit[[k]]
    if(!is.null(colnames(y))) rownames(acefit[[k]])<-colnames(y) else rownames(acefit[[k]])<-sapply(1:ncol(y),function(x) paste("y",x,sep=""))

    colnames(acefit[[k]])<-c("intercept","slope","p.intercept","p.slope")
  }
  close(pb)

  do.call(rbind,rootlist)->rootlist
  apply(rootlist,2,function(x) quantile(x,c(0.025,0.975)))->CIroot
  data.frame(root=t(y.ace)[,1],"CI 2.5"=t(CIroot)[,1],"CI 97.5"=t(CIroot)[,2])->root.conf.int
  if(!is.null(colnames(y))) rownames(root.conf.int)<-colnames(rootlist)<-colnames(y) else
    rownames(root.conf.int)<-colnames(rootlist)<-sapply(1:ncol(y),function(x) paste("y",x,sep=""))

  class(RR.list)<-class(acefit)<-"RRphyloList"

  res<-structure(list(mean.sampling = mean(real.s),
                      RR.list=RR.list,
                      root.est=rootlist,
                      rootCI=root.conf.int,
                      ace.regressions=acefit),
                 class = "RRphyloList")
  attr(res,"Call")<-funcall

  res
}
