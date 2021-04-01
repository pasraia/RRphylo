#'@title Extracting a user-specified subset of the evo.dir results
#'@description This function takes the result list produced by
#'  \code{\link{evo.dir}} as the input, and extracts a specific subset of it.
#'  The user may specify whether to extract the set of angles between species
#'  resultant vectors and the MRCA, the size of resultant vectors, or the set of
#'  angles between species.
#'@usage
#'  retrieve.angles(angles.res,wishlist=c("anglesMRCA","angleDir","angles.between.species"),
#'   random=c("yes","no"),focus=c("node","species","both","none"),
#'  node=NULL,species=NULL,write2csv=c("no","yes"),csvfile=NULL)
#'@param angles.res the object resulting from \code{\link{evo.dir}} function.
#'@param wishlist specifies whether to extract angles and sizes
#'  (\code{"anglesMRCA"}) of resultant vectors between individual species and
#'  the MRCA, angles and sizes (\code{"angleDir"}) of vectors between individual
#'  species and a fixed reference vector (the same for all species), or angles
#'  between species resultant vectors (\code{"angles.between.species"}).
#'@param random it needs to be \code{"yes"} if \code{'angles.res'} object
#'  contains randomization results.
#'@param focus it can be \code{"node"}, \code{"species"}, \code{"both"}, or
#'  \code{"none"}, whether the user wants the results for a focal node, or for a
#'  given species, for both, or just wants to visualize everything.
#'@param node must be indicated if \code{focus = "node"} or \code{"both"}. As
#'  for \code{evo.dir}, the node number must refer to the dichotomic version of
#'  the original tree, as produced by \code{\link{RRphylo}}.
#'@param species must be indicated if \code{focus = "species"} or \code{"both"}.
#'@param write2csv has been deprecated; please see the argument
#'   \code{csvfile} instead.
#' @param csvfile if results should be saved to a .csv file, a character
#'   indicating the name of the csv file and the path where it is to be saved.
#'   If no path is indicated the file is stored in the working directory.
#'   If left unspecified, no file will be saved.
#'@export
#'@importFrom utils write.csv
#'@details \code{retrieve.angles} allows to focalize the extraction to a
#'  particular node, species, or both. Otherwise it returns the whole dataset.
#'@return \code{retrieve.angles} outputs an object of class \code{'data.frame'}.
#'@return If \code{wishlist = "anglesMRCA"}, the data frame includes:
#'  \itemize{\item\strong{MRCA} the most recent common ancestor the angle is
#'  computed to \item\strong{species} species ID \item\strong{angle} the angle
#'  between the resultant vector of species and the MRCA
#'  \item\strong{vector.size} the size of the resultant vector computed from
#'  species to MRCA }
#'@return If \code{wishlist = "angleDir"}, the data frame includes:
#'  \itemize{\item\strong{MRCA} the most recent common ancestor the vector is
#'  computed to \item\strong{species} species ID \item\strong{angle.direction}
#'  the angle between the vector of species and a fixed reference
#'  \item\strong{vector.size} the size of the vector of species }
#'@return If \code{wishlist = "angles.between.species"}, the data frame
#'  includes: \itemize{\item\strong{MRCA} the most recent common ancestor the
#'  vector is computed from \item\strong{species} pair IDs of the species pair
#'  the "angle between species" is computed for
#'  \item\strong{angleBTWspecies2MRCA} angle between species resultant vectors
#'  to MRCA \item\strong{anglesBTWspecies} angle between species resultant
#'  vectors }
#'@author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'  Carotenuto
#' @examples
#'\dontrun{
#'     data("DataApes")
#'     DataApes$PCstage->PCstage
#'     DataApes$Tstage->Tstage
#'
#'     cc<- 2/parallel::detectCores()
#'     RRphylo(tree=Tstage,y=PCstage,clus=cc)->RR
#'
#' # Case 1. "evo.dir" without performing randomization
#'     evo.dir(RR,angle.dimension="rates",pair.type="node",
#'     node=  57,random="no")->evo.p
#'
#'  # Case 1.1 angles and sizes of resultant vectors between individual species and the MRCA:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="no",focus="node",
#'     node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="no",focus="species",
#'     species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="no",focus="both",
#'     node=68,species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="no",focus="none")
#'
#'  # Case 1.2 angles and sizes of vectors between individual species
#'  #and a fixed reference vector:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="angleDir",random="no",focus="node",
#'     node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="angleDir",random="no",focus="species",
#'     species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="angleDir",random="no",focus="both",
#'     node=68,species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="angleDir",random="no",focus="none")
#'
#'  # Case 1.3 angles between species resultant vectors:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="no",
#'     focus="node", node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="no",
#'     focus="species", species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="no",
#'     focus="both",node=68,species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="no",
#'     focus="none")
#'
#'
#' # Case 2. "evo.dir" with performing randomization
#'     evo.dir(RR,angle.dimension="rates",pair.type="node",node=57,
#'     random="yes",nrep=10)->evo.p
#'
#'  # Case 2.1 angles and sizes of resultant vectors between individual species and the MRCA:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="yes",focus="node",
#'     node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="yes", focus="species",
#'     species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="yes",focus="both",
#'     node=68,species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="anglesMRCA",random="yes",focus="none")
#'
#'  # Case 2.2 angles and sizes of vectors between individual species and a fixed reference vector:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="angleDir",random="yes",focus="node",
#'     node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="angleDir",random="yes",focus="species",
#'     species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="angleDir",random="yes",focus="both",
#'     node=68, species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="angleDir",random="yes",focus="none")
#'
#'  # Case 2.3 retrieve angles between species resultant vectors:
#'   # for a focal node
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="yes",
#'     focus="node", node=68)
#'   # for a focal species
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="yes",
#'     focus="species", species="Sap")
#'   # for both focal node and species
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="yes",
#'     focus="both",node=68,species="Sap")
#'   # without any specific requirement
#'     retrieve.angles(evo.p,wishlist="angles.between.species",random="yes",
#'     focus="none")
#'     }



retrieve.angles<-function(angles.res,
                          wishlist=c("anglesMRCA","angleDir","angles.between.species"),
                          random=c("yes","no"),
                          focus=c("node","species","both","none"),
                          node=NULL,
                          species=NULL,
                          write2csv=c("no","yes"),
                          csvfile=NULL)
{

  focus.test<-function(r1,focus=c("node","species","both","none"),node=NULL,species=NULL)
  {
    switch(focus, none={
      return(r1)
    }, node={
      r1[r1$MRCA==node,]->ab
      return(ab)
    }, species={
      r1[grep(species,r1$species),]->ac
      return(ac)
    }, both={
      r1[r1$MRCA==node,][grep(species,r1[r1$MRCA==node,]$species),]->aa
      return(aa)
    })
  }


  match.arg(wishlist)
  switch(wishlist, anglesMRCA={
    match.arg(random)
    if(random=="yes")
    {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(11,1,3)],strsplit(colnames(x),split="/")[[1]][1]),c(x[c(11,4,6)],strsplit(colnames(x),split="/")[[1]][2])))))->r1
      r1[!duplicated(r1),]->r1
      r1[order(r1[,1],r1[,4]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species","angle","vector.size")

    } else {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(7,1,2)],strsplit(colnames(x),split="/")[[1]][1]),c(x[c(7,3,4)],strsplit(colnames(x),split="/")[[1]][2])))))->r1
      r1[!duplicated(r1),]->r1
      r1[order(r1[,1],r1[,4]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species","angle","vector.size")
    }

    match.arg(focus)
    focus.test(r1,focus=focus,node=node,species=species)->r1

  }, angleDir ={
    match.arg(random)
    if(random=="yes")
    {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(11,12,13)],strsplit(colnames(x),split="/")[[1]][1]),c(x[c(11,14,15)],strsplit(colnames(x),split="/")[[1]][2])))))->r1
      r1[!duplicated(r1),]->r1
      r1[order(r1[,1],r1[,4]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species","angle.direction","vector.size")
    } else {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(7,8,9)],strsplit(colnames(x),split="/")[[1]][1]),c(x[c(7,10,11)],strsplit(colnames(x),split="/")[[1]][2])))))->r1
      r1[!duplicated(r1),]->r1
      r1[order(r1[,1],r1[,4]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species","angle.direction","vector.size")
    }

    match.arg(focus)
    focus.test(r1,focus=focus,node=node,species=species)->r1

  }, angles.between.species = {
    match.arg(random)
    if(random=="yes")
    {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(11,7,9)],colnames(x))))))->r1
      r1[order(r1[,1]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species pairs","angleBTWspecies2MRCA","anglesBTWspecies")
    } else {
      as.data.frame(do.call(rbind,lapply(angles.res,function(x) rbind(c(x[c(7,5,6)],colnames(x))))))->r1
      r1[order(r1[,1]),]->r1
      r1[,c(1,4,2,3),]->r1
      colnames(r1)<-c("MRCA","species pairs","angleBTWspecies2MRCA","anglesBTWspecies")
    }
    match.arg(focus)
    focus.test(r1,focus=focus,node=node,species=species)->r1

  })

  if(!is.null(csvfile)) write.csv(r1,file=paste(csvfile,".csv",sep=""))
  return(r1)
}
