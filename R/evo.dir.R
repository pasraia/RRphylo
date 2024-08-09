#' @title Phylogenetic vector analysis of phenotypic change
#'
#' @description This function quantifies direction, size and rate of
#'   evolutionary change of multivariate traits along node-to-tip paths and
#'   between species.
#' @usage evo.dir(RR,angle.dimension=c("rates","phenotypes"),
#'   y.type=c("original","RR"),y=NULL,pair.type=c("node","tips"),pair=NULL,
#'   node=NULL,random=c("yes","no"),nrep=100)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param angle.dimension specifies whether vectors of \code{"rates"} or
#'   \code{"phenotypes"} are used.
#' @param y.type must be indicated when \code{angle.dimension = "phenotypes"}.
#'   If \code{"original"}, it uses the phenotypes as provided by the user, if
#'   \code{"RR"} it uses \code{RR$predicted.phenotypes}.
#' @param y specifies the phenotypes to be provided if \code{y.type =
#'   "original"}.
#' @param pair.type either \code{"node"} or \code{"tips"}. Angles are computed
#'   between each possible couple of species descending from a specified node
#'   (\code{"node"}), or between a given couple of species (\code{"tips"}).
#' @param pair species pair to be specified if \code{pair.type = "tips"}. It
#'   needs to be written as in the example below.
#' @param node node number to be specified if \code{pair.type = "node"}. Notice
#'   the node number must refer to the dichotomic version of the original tree,
#'   as produced by \code{RRphylo}.
#' @param random whether to perform randomization test
#'   (\code{"yes"}/\code{"no"}).
#' @param nrep number of replications must be indicated if \code{random =
#'   "yes"}. It is set at 100 by default.
#' @importFrom ape extract.clade getMRCA
#' @importFrom stats quantile
#' @importFrom phytools getDescendants
#' @export
#' @details The way \code{evo.dir} computes vectors depends on whether
#'   phenotypes or rates are used as variables. \code{\link{RRphylo}} rates
#'   along a path are aligned along a chain of ancestor/descendant
#'   relationships. As such, each rate vector origin coincides to the tip of its
#'   ancestor, and the resultant of the path is given by vector addition. In
#'   contrast, phenotypic vectors are computed with reference to a common origin
#'   (i.e. the consensus shape in a geometric morphometrics). In this latter
#'   case, vector subtraction (rather than addition) will define the resultant
#'   of the evolutionary direction.  It is important to realize that resultants
#'   could be at any angle even if the species (the terminal vectors) have
#'   similar phenotypes, because path resultants, rather than individual
#'   phenotypes, are being contrasted. However, the function also provides the
#'   angle between individual phenotypes as 'angle.between.species'. To perform
#'   randomization test (\code{random = "yes"}), the evolutionary directions of
#'   the two species are collapsed together. Then, for each variable, the median
#'   is found, and random paths of the same size as the original paths are
#'   produced sampling at random from the 47.5th to the 52.5th percentile around
#'   the medians. This way, a random distribution of angles is obtained under
#'   the hypothesis that the two directions are actually parallel. The
#'   'angle.direction' represents the angle formed by the species phenotype and
#'   a vector of 1s (as long as the number of variables representing the
#'   phenotype). This way, each species phenotype is contrasted to the same
#'   vector. The 'angle.direction' values could be inspected to test whether
#'   individual species phenotypes evolve towards similar directions.
#' @return Under all specs, \code{evo.dir} returns a 'list' object. The length
#'   of the list is one if \code{pair.type = "tips"}.  If \code{pair.type =
#'   "node"}, the list is as long as the number of all possible species pairs
#'   descending from the node. Each element of the list contains:
#' @return \strong{angle.path.A} angle of the resultant vector of species A to
#'   MRCA
#' @return \strong{vector.size.species.A} size of the resultant vector of
#'   species A to MRCA
#' @return \strong{angle.path.B} angle of the resultant vector of species B to
#'   MRCA
#' @return \strong{vector.size.species.B} size of the resultant vector of
#'   species B to MRCA
#' @return \strong{angle.between.species.to.mrca} angle between the species
#'   paths resultant vectors to the MRCA
#' @return \strong{angle.between.species} angle between species vectors (as they
#'   are, without computing the path)
#' @return \strong{MRCA} the node identifying the most recent common ancestor of
#'   A and B
#' @return \strong{angle.direction.A} angle of the vector of species A (as it
#'   is, without computing the path) to a fixed reference vector (the same for
#'   all species)
#' @return \strong{vec.size.direction.A} size of the vector of species A
#' @return \strong{angle.direction.B} angle of the vector of species B (as it
#'   is, without computing the path) to a fixed reference vector (the same for
#'   all species)
#' @return \strong{vec.size.direction.B} size of the vector of species B
#' @return If \code{random = "yes"}, results also include p-values for the
#'   angles.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' \dontrun{
#'   data("DataApes")
#'   DataApes$PCstage->PCstage
#'   DataApes$Tstage->Tstage
#'   cc<- 2/parallel::detectCores()
#'
#'   RRphylo(tree=Tstage,y=PCstage, clus=cc)->RR
#'
#' # Case 1. Without performing randomization test
#'
#'  # Case 1.1 Computing angles between rate vectors
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="rates",pair.type="node",node=72 ,
#'     random="no")
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="rates",pair.type="tips",
#'     pair= c("Sap_1","Tro_2"),random="no")
#'
#'  # Case 1.2 computing angles between phenotypic vectors provided by the user
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="original",
#'     y=PCstage,pair.type="node",node=72,random="no")
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="original",
#'     y=PCstage,pair.type="tips",pair=c("Sap_1","Tro_2"),random="no")
#'
#'  # Case 1.3 computing angles between phenotypic vectors produced by "RRphylo"
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="RR",
#'     pair.type="node",node=72,random="no")
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="RR",
#'     pair.type="tips",pair=c("Sap_1","Tro_2"),random="no")
#'
#'
#' # Case 2. Performing randomization test
#'
#'  # Case 2.1 Computing angles between rate vectors
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="rates",pair.type="node",node=72 ,
#'     random="yes",nrep=10)
#'
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="rates",pair.type="tips",
#'     pair= c("Sap_1","Tro_2"),random="yes",nrep=10)
#'
#'  # Case 2.2 computing angles between phenotypic vectors provided by the user
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="original",
#'     y=PCstage,pair.type="node",node=72,random="yes",nrep=10)
#'
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="original",
#'     y=PCstage,pair.type="tips",pair=c("Sap_1","Tro_2"),random="yes",nrep=10)
#'
#'  # Case 2.3 computing angles between phenotypic vectors produced by "RRphylo"
#'   # for each possible couple of species descending from node 72
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="RR",
#'     pair.type="node",node=72,random="yes",nrep=10)
#'
#'   # for a given couple of species
#'     evo.dir(RR,angle.dimension="phenotypes",y.type="RR",
#'     pair.type="tips",pair=c("Sap_1","Tro_2"),random="yes",nrep=10)
#'     }

evo.dir<-function(RR,
                  angle.dimension=c("rates","phenotypes"),
                  y.type=c("original","RR"),
                  y=NULL,
                  pair.type=c("node","tips"),
                  pair=NULL,
                  node=NULL,
                  random=c("yes","no"),
                  nrep=100)
{
  #require(phytools)


  angle.btw.species<-function(tree,pair,phen,random=rdm,angle.dimension=AD)
  {
    #require(ape)
    unitV<-function(x){ sum(x^2)^.5 }
    extract.clade(tree,getMRCA(tree,pair))->Htree
    L[match(tips(tree,getMRCA(tree,pair)),rownames(L)),]->Lsub
    getMRCA(tree,pair)->subN
    getDescendants(tree,subN)->savenode
    c(subN,savenode)->savenode
    which(apply(Lsub,2,sum)==0)->cutnode
    if(length(which(names(cutnode)%in%savenode))==0) cutnode->cutnode else cutnode[-which(names(cutnode)%in%savenode)]->cutnode
    Lsub[,which(!colnames(Lsub)%in%names(cutnode))]->Lsub
    Lsub[,seq(match(getMRCA(tree,pair),colnames(Lsub)),dim(Lsub)[2])]->Lsub
    colnames(Lsub)[seq(1:(Ntip(Htree)-1))]->Htree$node.label

    phen[match(colnames(Lsub),rownames(phen)),]->subetas
    apply(subetas,1,function(x) sum(x^2)^.5)->vec.len
    as.matrix(vec.len)->vec.len

    Lsub[match(pair[1],rownames(Lsub)),]->a
    savenode[which(savenode%in%getMommy(tree,which(tree$tip.label==pair[1])))]->node.a
    names(a[c(which(names(a)%in%node.a),unname(which(a!=0)))])->a
    a[!duplicated(a)]->a
    Lsub[match(pair[2],rownames(Lsub)),]->b
    savenode[which(savenode%in%getMommy(tree,which(tree$tip.label==pair[2])))]->node.b
    names(b[c(which(names(b)%in%node.b),unname(which(b!=0)))])->b
    b[!duplicated(b)]->b

    phen[match(a,rownames(phen)),]->a.betas
    phen[match(b,rownames(phen)),]->b.betas
    as.matrix(a.betas)->a.betas
    as.matrix(b.betas)->b.betas

    a.betas[1,]->ancestor

    if(angle.dimension=="phenotypes")
    {
      if(dim(a.betas)[1]==2) apply(a.betas,2,diff)->a.resultant else apply(a.betas,2,function(x) x[1]-sum(x[2:length(x)]))->a.resultant
      if(dim(b.betas)[1]==2) apply(b.betas,2,diff)->b.resultant else apply(b.betas,2,function(x) x[1]-sum(x[2:length(x)]))->b.resultant
    } else {

      apply(a.betas,2,sum)->a.resultant
      apply(b.betas,2,sum)->b.resultant
    }

    a.betas[dim(a.betas)[1],]->a.conv
    b.betas[dim(b.betas)[1],]->b.conv
    rad2deg(acos((rep(1,dim(a.betas)[2])%*%a.conv)/(unitV(rep(1,dim(a.betas)[2]))*unitV(a.conv))))->conv.angA
    rad2deg(acos((rep(1,dim(b.betas)[2])%*%b.conv)/(unitV(rep(1,dim(b.betas)[2]))*unitV(b.conv))))->conv.angB
    unitV(a.conv)->conv.sizeA
    unitV(b.conv)->conv.sizeB

    rad2deg(acos((ancestor%*%a.resultant)/(unitV(ancestor)*unitV(a.resultant))))->angA
    rad2deg(acos((ancestor%*%b.resultant)/(unitV(ancestor)*unitV(b.resultant))))->angB
    rad2deg(acos((a.resultant%*%b.resultant)/(unitV(a.resultant)*unitV(b.resultant))))->ang.mrca
    rbind(a.betas[dim(a.betas)[1],],b.betas[dim(b.betas)[1],])->specieswise
    rownames(specieswise)<-c(pair[1],pair[2])
    rad2deg(acos((specieswise[1,]%*%specieswise[2,])/(unitV(specieswise[1,])*unitV(specieswise[2,]))))->theta.species
    unitV(a.resultant)->a.size
    unitV(b.resultant)->b.size



    if(random=="yes")
    {
      RangA<-array()
      RangB<-array()
      Rang.mrca<-array()
      theta.speciesR<-array()
      for(j in 1:(nrep-1))
      {


        phen[match(a,rownames(phen)),]->a.betas
        phen[match(b,rownames(phen)),]->b.betas
        a.betas[1,]->ancestor


        if(sum(ancestor==0)){
          ancestor<-rep(1,dim(a.betas)[2])
          names(ancestor)<-colnames(a.betas)
        }

        rbind(a.betas,b.betas)->ab
        if(dim(ab)[1]==2)
        {

          apply(ab,2,function(x) runif(1,quantile(min(x),.475),quantile(max(x),.525)))->a.betas
          apply(ab,2,function(x) runif(1,quantile(min(x),.475),quantile(max(x),.525)))->b.betas

        } else {

          apply(ab,2,function(x) runif(dim(a.betas)[1],quantile(min(x),.475),quantile(max(x),.525)))->a.betas
          apply(ab,2,function(x) runif(dim(b.betas)[1],quantile(min(x),.475),quantile(max(x),.525)))->b.betas
          rownames(a.betas)<-a
          rownames(b.betas)<-b
        }


        as.matrix(a.betas)->a.betas
        as.matrix(b.betas)->b.betas

        if(angle.dimension=="phenotypes")
        {
          if(dim(a.betas)[1]==2) apply(a.betas,2,diff)->a.resultant else apply(a.betas,2,function(x) x[1]-sum(x[2:length(x)]))->a.resultant
          if(dim(b.betas)[1]==2) apply(b.betas,2,diff)->b.resultant else apply(b.betas,2,function(x) x[1]-sum(x[2:length(x)]))->b.resultant
        } else {

          apply(a.betas,2,sum)->a.resultant
          apply(b.betas,2,sum)->b.resultant
        }

        as.numeric(ancestor)->ancestor

        rad2deg(acos((ancestor%*%a.resultant)/(unitV(ancestor)*unitV(a.resultant))))->RangA[j]
        rad2deg(acos((ancestor%*%b.resultant)/(unitV(ancestor)*unitV(b.resultant))))->RangB[j]
        rad2deg(acos((a.resultant%*%b.resultant)/(unitV(a.resultant)*unitV(b.resultant))))->Rang.mrca[j]

        rbind(a.betas[dim(a.betas)[1],],b.betas[dim(b.betas)[1],])->specieswiseR
        apply(specieswiseR,2,function(x) runif(2,min(x),max(x)))->specieswiseR
        rownames(specieswiseR)<-c(pair[1],pair[2])
        rad2deg(acos((specieswiseR[1,]%*%specieswiseR[2,])/(unitV(specieswiseR[1,])*unitV(specieswiseR[2,]))))->theta.speciesR[j]

      }
      length(which(RangA>as.numeric(angA)))/(nrep)->p.angleA
      length(which(RangB>as.numeric(angB)))/(nrep)->p.angleB
      length(which(Rang.mrca>as.numeric(ang.mrca)))/(nrep)->p.angle.mrca
      length(which(theta.speciesR>as.numeric(theta.species)))/(nrep)->p.angle.btw

      par(mfrow=c(2,2),ask=F)
      hist(RangA, xlab="random angle species A",main="path angle A")
      abline(v=angA,col="red",lwd=4)
      hist(RangB, xlab="random angle species B",main="path angle B")
      abline(v=angB,col="green",lwd=4)
      hist(Rang.mrca,xlab="random angle between paths",main="angle between species to mrca")
      abline(v=ang.mrca,col="blue",lwd=4)
      hist(theta.speciesR,xlab="random angle between species",main="angle between species")
      abline(v=theta.species,col="orange",lwd=4)




      t(data.frame("angle path A"=angA,"p angle path A"=p.angleA,"vector size species A"= a.size, "angle path B"=angB,"p angle path B"=p.angleB,"vector size species B"= b.size,"angle between species to mrca"=ang.mrca,"p.angle to mrca"=p.angle.mrca,"angle between species"=theta.species,"p.angle between species"=p.angle.btw,"MRCA"=getMRCA(tree,pair),"angle direction A"=conv.angA,"vec size direction A"=conv.sizeA,"angle direction B"=conv.angB,"vec size direction B"=conv.sizeB))->res
      paste(pair[1],pair[2],sep="/")->colnames(res)
      return(res)
    } else {
      t(data.frame("angle path A"=angA,"vector size species A"= a.size, "angle path B"=angB,"vector size species B"= b.size,"angle between species to mrca"=ang.mrca,"angle between species"=theta.species,"MRCA"=getMRCA(tree,pair),"angle direction A"=conv.angA,"vec size direction A"=conv.sizeA,"angle direction B"=conv.angB,"vec size direction B"=conv.sizeB))->res
      paste(pair[1],pair[2],sep="/")->colnames(res)
      return(res)

    }
  }

  RR$tree->tree
  RR$tip.path->L
  RR$lambda->lambda


  match.arg(angle.dimension)
  angle.dimension->AD
  if(AD=="rates")
  {
    RR$multiple.rates->phen
  } else {

    match.arg(y.type)
    if(y.type=="original"){
      # y <- treedata(RR$tree, y, sort = TRUE)[[2]]
      y <- treedataMatch(RR$tree, y)[[1]]
      y->tipP
      RR$aces->ancP
      rbind(ancP,tipP)->phen
    } else{
      RR$predicted.phenotype->tipP
      RR$aces->ancP
      rbind(ancP,tipP)->phen
    }
  }


  if (length(y) == Ntip(tree)) stop("rate vectors should be multivariate")




  deg2rad <- function(deg) {(deg * pi) / (180)}
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  match.arg(random)
  random->rdm





  match.arg(pair.type)
  if(pair.type=="node")
  {

    tips(tree,node)->leaves
    combn(leaves,2)->pairs

    results<-list()
    for(k in 1:dim(pairs)[2])
    {
      pairs[1,k]->a
      pairs[2,k]->b
      pair=c(a,b)

      angle.btw.species(tree,pair,random=rdm,phen)->results[[k]]
    }
    return(results)
  } else {

    angle.btw.species(tree,pair,phen, random=rdm)->results
  }
  return(results)
}
