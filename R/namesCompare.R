#' @title Checking species names for misspelling and synonyms
#' @description The function cross-references two vectors of species names
#'   checking for possible synonyms, misspelled names, and genus-species or
#'   species-subspecies correspondence.
#' @usage namesCompare(vec1,vec2,proportion=0.15,
#'   focus=c("genus","subspecies","epithet","misspelling"))
#' @param vec1,vec2 a vector of species names. Genus names only are also
#'   allowed. Generic name and specific epithet must be separated by '_'. Note
#'   that \code{vec2} is used as the reference. Incomplete or suspicious names
#'   are better placed in \code{vec1} (see example below).
#' @param proportion the maximum proportion of different characters between any
#'   \code{vec1-vec2} names pair to consider it a possible misspelling.
#' @param focus one or more of \code{"genus"}, \code{"subspecies"},
#'   \code{"epithet"}, \code{"misspelling"}. See Values for details.
#' @importFrom utils adist
#' @export
#' @return The function returns a \code{list} including (according to \code{focus}):
#' @return \strong{$genus} if \code{vec1} includes genera names which miss
#'   specific epithet, this object lists all the species in \code{vec2}
#'   belonging to each of the genera.
#' @return \strong{$subspecies} if \code{vec1} includes subspecies (i.e. two
#'   epithets after genus name), this object lists species in \code{vec2}
#'   possibly corresponding to each of the subspecies.
#' @return \strong{$epithet} lists species with matching epithets as possible
#'   synonyms.
#' @return \strong{$misspelling} lists possible misspelled names. For each
#'   proposed mismatched names pair the proportion of characters in the
#'   \code{vec1} differing from the string in \code{vec2} is returned.
#' @author Silvia Castiglione, Carmela Serio, Antonella Esposito
#' @examples
#' \dontrun{
#' names(DataFelids$statefel)->nams.fel
#' nams.fel[c(19,12,37,80,43)]<-c("Puma_yagouaroundi","Felis_manul","Catopuma",
#'                            "Pseudaelurus","Panthera_zdansky")
#' nams<-nams.fel[-81]
#'
#' namesCompare(nams,names(DataFelids$statefel))->nc1
#' namesCompare(names(DataFelids$statefel),nams)->nc2
#' }

namesCompare<-function(vec1,vec2,proportion=0.15,
                       focus=c("genus","subspecies","epithet","misspelling")){
  if(any(!focus%in%c("genus","subspecies","epithet","misspelling")))
    stop("focus must be one or more of \"genus\",\"subspecies\",\"epithet\",\"misspelling\"")

  if(!is.null(ncol(vec1))){
    if(ncol(vec1)>1) warning(paste("vec1:",ncol(vec1),"columns supplied, only the first one will be used"))
    vec1[,1]->vec1
  }

  if(!is.null(ncol(vec2))){
    if(ncol(vec2)>1) warning(paste("vec2:",ncol(vec2),"columns supplied, only the first one will be used"))
    vec2[,1]->vec2
  }

  if(all(vec1%in%vec2)) stop("vec1 and vec2 species perfectly match")

  vec1[which(!vec1%in%vec2)]->vec1
  vec2[which(!vec2%in%vec1)]->vec2

  ### GENERA CHECK ###
  if("genera"%in%focus){
    if(any(!grepl("_",vec1))){
      do.call(rbind,lapply(vec1[!grepl("_",vec1)],function(x)
        data.frame(vec1=rep(x,length(grep(x,vec2))),vec2=vec2[grep(x,vec2)])))->genera
      if(nrow(genera)<1) genera<-NULL
    }else genera<-NULL
    vec1[which(!vec1%in%genera[,1])]->vec1
  }else genera<-NULL


  ### SUBSPECIES CHECK ###
  if("subspecies"%in%focus){
    if(any(sapply(strsplit(vec1,"_"),length)>2)){
      vec1[which(sapply(strsplit(vec1,"_"),length)>2)]->suspv1
      strsplit(vec1,"_")[which(sapply(strsplit(vec1,"_"),length)>2)]->suspv1_split

      do.call(rbind,mapply(j=suspv1_split,w=suspv1,function(j,w){
        combn(j,2,function(k) paste(k,collapse="_"))[1:2]->check

        if(any(check%in%vec2)) data.frame(vec1=rep(w,length(which(vec2%in%check))),vec2=vec2[which(vec2%in%check)]) else{
          adcheck2 <- adist(check,vec2)/as.integer(nchar(check))
          adcheck2 <- data.frame(vec1=check,vec2=vec2[apply(adcheck2, 1, which.min)],proportion=apply(adcheck2, 1, min))
          if(any(adcheck2$proportion<=proportion)){
            data.frame(vec1=rep(w,length(which(adcheck2$proportion<=proportion))),
                       vec2=adcheck2[which(adcheck2$proportion<=proportion),2])
          } else NULL
        }
      },SIMPLIFY = FALSE))->subspecies
      vec1[which(!vec1%in%subspecies[,1])]->vec1
    } else subspecies<-NULL
  }else subspecies<-NULL

  ### MISSPELLING CHECK ###
  if("misspelling"%in%focus){
    ad1 <- adist(vec1,vec2)/nchar(vec1)
    ad1.dataframe <- data.frame(vec1,vec2=vec2[apply(ad1, 1, which.min)],proportion=apply(ad1, 1, min))

    ad1.dataframe[which(ad1.dataframe$proportion<=proportion),]->ad1.dataframe
    misspelling <- ad1.dataframe[order(ad1.dataframe[,3]),]
    if(nrow(misspelling)<1) misspelling<-NULL
  }else misspelling<-NULL

  ### epithet CHECK ###
  if("epithet"%in%focus){
    vec1[which(sapply(strsplit(vec1,"_"),length)>=2)]->vec1
    vec2[which(sapply(strsplit(vec2,"_"),length)>=2)]->vec2

    do.call(rbind,lapply(1:length(vec2),function(k){
      strsplit(vec2,"_")[[k]]->ep2
      if(length(ep2)>2) {
        if(any(duplicated(ep2))) data.frame(ind=k, sp=ep2[2]) else data.frame(ind=rep(k,2), sp=ep2[2:3])
      } else data.frame(ind=k, sp=ep2[2])
    }))->episp2

    do.call(rbind,lapply(1:length(strsplit(vec1,"_")),function(k){
      if(length(strsplit(vec1,"_")[[k]])>2) {
        strsplit(vec1,"_")[[k]][2:3]->ep
        if(any(duplicated(ep))) ep[1]->ep
      } else strsplit(vec1,"_")[[k]][2]->ep

      if(any(episp2[,2]%in%ep)) data.frame(vec1=rep(vec1[k],length(which(episp2[,2]%in%ep))),
                                           vec2=vec2[episp2[which(episp2[,2]%in%ep),1]]) else {
                                             adepi <- adist(ep,episp2[,2])/nchar(ep)
                                             adepi<-data.frame(ep1=ep,ep2=vec2[episp2[apply(adepi, 1, which.min),1]],proportion=apply(adepi, 1, min))
                                             if(any(adepi$proportion<=proportion)){
                                               data.frame(vec1=rep(vec1[k],length(which(adepi$proportion<=proportion))),
                                                          vec2=adepi[which(adepi$proportion<=proportion),2])
                                             } else NULL
                                           }
    }))->epithet
  }else epithet<-NULL

  res<-list(genus=genera,subspecies=subspecies,epithet=epithet,misspelling=misspelling)
  res[which(names(res)%in%focus)]
  return(res[which(names(res)%in%focus)])
}
