#' @title Run StableTraits from within R
#' @description This function runs StableTraits and StableTraitsSum
#'   (\cite{Elliot and Mooers 2014}) from within the R environment and returns
#'   its output into the workspace.
#' @usage
#'   StableTraitsR(tree,y,path,output=NULL,aces=NULL,argST=NULL,argSTS=NULL)
#' @param tree a phylogenetic tree. The tree needs not to be either ultrametric
#'   or fully dichotomous.
#' @param y a named vector of phenotypic trait.
#' @param path the folder path where the StableTraits output will be stored.
#'   Notice that the input tree and data (modified automatically if the original
#'   tree is not fully dichotomous or if \code{aces} are specified) will be
#'   stored in this folder as well.
#' @param output name of the output to be returned, if unspecified it will be
#'   named "output".
#' @param aces a named vector of ancestral character values at nodes specified
#'   in advance. Names correspond to the nodes in the tree.
#' @param argST a list of further arguments passed to StableTraits. If the
#'   argument has no value (for example "brownian") it must be specified as
#'   \code{TRUE}.
#' @param argSTS list of further arguments passed to StableTraitsSum. If the
#'   argument has no value (for example "brownian") it must be specified as
#'   \code{TRUE}.
#' @return The function returns a 'list' containing the output of StableTraits
#'   and StableTraitsSum.
#' @return \strong{$progress} a table reporting the DIC and PRSF diagnostics.
#' @return \strong{$rates_tree} a copy of the original tree with branch lengths
#'   set to the evolutionary rate imputed by the stable reconstruction.
#'   Specifically, each branch length is equal to the absolute difference in the
#'   stable reconstruction occurring on that branch divided by the square root
#'   of the input branch length.
#' @return \strong{$rates} the original branch lengths, evolutionary rates, node
#'   height and (optionally) scaled branch lengths.
#' @return \strong{$aces} the median estimates of ancestral states and stable
#'   parameters along with the 95\% credible interval.
#' @return \strong{$brownian_tree} if "brownian" is \code{TRUE} in
#'   \code{argSTS}, a copy of the original tree with branch lengths set such
#'   that the Brownian motion reconstruction of the character on this tree is
#'   approximately the same as the stable ancestral reconstruction.
#' @return \strong{$ace.prior.values} if \code{aces} is specified, the function
#'   returns a dataframe containing the corresponding node number on the
#'   \code{\link{RRphylo}} tree for each node, the original (preset) and the estimated
#'   values, and the 95\% credible interval.
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @details The StableTraits software is available at https://mickelliot.com/,
#'   along with instructions for compilation. Once it is installed, the user
#'   must set as R working directory the folder where the StableTraits software
#'   are installed. Further information about the arguments and outputs of
#'   StableTraits and StableTraitsSum can be found at https://mickelliot.com/.
#'   \code{StableTraitsR} automatically recognizes which Operating System is
#'   running on the computer (it has been tested successfully on MacOS and
#'   Windows machines).
#' @importFrom utils read.table write.table
#' @importFrom ape write.tree read.tree
#' @export
#' @references Elliot, M. G., & Mooers, A. Ø. (2014). Inferring ancestral states
#'   without assuming neutrality or gradualism using a stable model of
#'   continuous character evolution. \emph{BMC evolutionary biology}, 14: 226.
#'   doi.org/10.1186/s12862-014-0226-8
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#'
#' # Set as working directory the folder where StableTraits software are stored
#' # setwd("~/StableTraits")
#'
#' dir.create("Analyses")
#' rtree(100)->tree
#' fastBM(tree)->y
#' c(1,2,3)->acev
#' sample(Ntip(tree)+seq(1:Nnode(tree)),3)->names(acev)
#' StableTraitsR(tree,y,path="Analyses/",output="my_output",aces=acev,
#' argST=list(iterations=500000,chains=4),argSTS=list(brownian=TRUE))->STr
#' }


StableTraitsR<-function(tree,y,path,output=NULL,aces=NULL,argST=NULL,argSTS=NULL){
  # require(ape)
  # require(phytools)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }

  if (is.binary(tree)) t <- tree else t <- multi2di(tree,random=FALSE)
  y[match(t$tip.label,names(y))]->y
  t->tree.start

  if(is.null(aces)==FALSE){
    insert.at <- function(a, pos, ...){
      dots <- list(...)
      stopifnot(length(dots)==length(pos))
      result <- vector("list",2*length(pos)+1)
      result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
      result[c(FALSE,TRUE)] <- dots
      unlist(result)
    }
    aces->aceV
    L <- makeL(t)
    if(is.null(names(aceV))) stop("The vector of aces needs to be named")
    if (is.binary(tree)==FALSE){
      ac<-array()
      for(i in 1:length(aceV)){
        getMRCA(t,tips(tree,names(aceV)[i]))->ac[i]
      }
      ac->names(aceV)
    }

    aceV->P
    as.numeric(names(aceV))->N
    lapply(N, function(x) tips(t,x))->tar.tips
    names(tar.tips)<-N

    t->treeN
    y->ynew
    i=1
    while(i<=length(N)){
      getMRCA(treeN,tar.tips[[i]])->nn
      bind.tip(treeN, tip.label=paste("nod",N[i],sep=""), edge.length=0, where=nn, position=0.001)->treeN
      which(treeN$tip.label==paste("nod",N[i],sep=""))->npos
      if(npos==1) c(P[i],ynew)->ynew
      if(npos==length(ynew)+1)  c(ynew,P[i])->ynew else{
        if(npos>1 & npos<(length(ynew)+1))  insert.at(ynew,npos-1,P[i])->ynew
      }
      names(ynew)[npos]<-paste("nod",N[i],sep="")
      i=i+1
    }
    treeN->t
    ynew->y
    if(length(which(t$edge.length==0))>0) t$edge.length[which(t$edge.length==0)]<-0.0001

    paste("nod",N,sep="")->tip.rem
    nod.rem<-array()
    for(i in 1:length(N)){
      getMRCA(t,c(tip.rem[i],tar.tips[[i]]))->nod.rem[i]
    }
    # nod.rem
    # tip.rem
  }

  write.tree(t,file=paste(path,"tree",sep=""))
  write.table(y,file=paste(path,"y",sep=""),col.names = FALSE,quote=FALSE)

  if(is.null(output)) "output"->out else output->out

  if(is.null(argST)==FALSE){
    argSTnew<-list()
    for(i in 1:length(argST)){
      if(isTRUE(argST[[i]])) paste("--",names(argST)[[i]],sep="")->argSTnew[[i]] else
        paste("--",names(argST)[[i]]," ",as.integer(argST[[i]]),sep="")->argSTnew[[i]]
    }
    unlist(argSTnew)->argSTnew
    if(.Platform$OS.type=="unix"){
      c("./stabletraits", paste("--tree ", shQuote(path,"sh"),"tree",sep=""),
        paste("--data ", shQuote(path,"sh"),"y",sep=""),paste("--output ", shQuote(path,"sh"),shQuote(out,"sh"),sep=""),
        argSTnew)->tot.arg
    }else{
      c("/c","stabletraits", paste("--tree ", shQuote(path,"cmd"),"tree",sep=""),
        paste("--data ", shQuote(path,"cmd"),"y",sep=""),paste("--output ", shQuote(paste(path,out,sep=""),"cmd"),sep=""),
        argSTnew)->tot.arg
    }

  }else{
    if(.Platform$OS.type=="unix"){
      c("./stabletraits", paste("--tree ", shQuote(path,"sh"),"tree",sep=""),
        paste("--data ", shQuote(path,"sh"),"y",sep=""),paste("--output ", shQuote(path,"sh"),shQuote(out,"sh"),sep=""))->tot.arg
    }else{
      c("/c", "stabletraits", paste("--tree ", shQuote(path,"cmd"),"tree",sep=""),
        paste("--data ", shQuote(path,"cmd"),"y",sep=""),paste("--output ", shQuote(paste(path,out,sep=""),"cmd"),sep=""))->tot.arg
    }
  }


  if(is.null(argSTS)==FALSE){
    argSTSnew<-list()
    for(i in 1:length(argSTS)){
      if(isTRUE(argSTS[[i]])) paste("--",names(argSTS)[[i]],sep="")->argSTSnew[[i]] else
        paste("--",names(argSTS)[[i]]," ",as.integer(argSTS[[i]]),sep="")->argSTSnew[[i]]
    }
    unlist(argSTSnew)->argSTSnew
    if(.Platform$OS.type=="unix") c("./stabletraitssum", paste("--path ", shQuote(path,"sh"),shQuote(out,"sh"),sep=""),argSTSnew)->tot.argSTS else
      c("/c", "stabletraitssum", paste("--path ", shQuote(paste(path,out,sep=""),"cmd"),sep=""),argSTSnew)->tot.argSTS
  }else{
    if(.Platform$OS.type=="unix") c("./stabletraitssum", paste("--path ", shQuote(path,"sh"),shQuote(out,"sh"),sep=""))->tot.argSTS else
      c("/c", "stabletraitssum", paste("--path ", shQuote(paste(path,out,sep=""),"cmd"),sep=""))->tot.argSTS
  }

  if(.Platform$OS.type=="unix") {
    system2("command", args = tot.arg)
    system2("command", args = tot.argSTS)
  }else{
    system2("cmd", args = tot.arg)
    system2("cmd", args = tot.argSTS)
  }

  read.table(paste(path,out,".ancstates",sep=""),header = TRUE)->aces.list
  read.table(paste(path,out,".brlens",sep=""),header = TRUE)->rates.res
  read.tree(paste(path,out,".tree",sep=""))->STtree
  read.table(paste(path,out,".progress",sep=""),header = TRUE)->STprogress
  read.tree(paste(path,out,".rates_tree",sep=""))->rates_tree
  if("brownian"%in%names(argSTS)) read.tree(paste(path,out,".scaled_tree",sep=""))->brownian_tree else NULL->brownian_tree

  sapply(strsplit(as.character(rates.res[,1]),"->") ,"[[",1)->rates.res[,1]

  as.character(aces.list[,1])->aces.list[,1]
  as.character(rates.res[,1])->rates.res[,1]
  if(length(which(STtree$node.label=="NA"))>0){
    for(i in 1:length(which(STtree$node.label=="NA"))){
      STtree$node.label[which(STtree$node.label=="NA")[i]]<-paste("noname",i,sep="")
      aces.list[which(is.na(aces.list[,1]))[i],1]<-paste("noname",i,sep="")
      rates.res[which(rates.res[,1]=="NA")[i],1]<-paste("noname",i,sep="")
    }
  }

  data.frame(lab=STtree$node.label, num=seq((Ntip(STtree)+1), (Ntip(STtree)+Nnode(STtree))))->lab2num

  rates.res$num<-NA
  rates.res[match(lab2num[-1,1],rates.res[,1]),]$num<-lab2num[-1,2]
  rates.res[which(is.na(rates.res$num)),]$num<-rates.res[which(is.na(rates.res$num)),1]
  data.frame(num=rates.res$num,rates.res[,-which(colnames(rates.res)=="num")])->rates.res
  #subset(aces.list,Parameter%in%as.character(lab2num[,1]))->aces.res
  aces.list[match(as.character(lab2num[,1]),aces.list$Parameter),]->aces.res
  data.frame(num=lab2num[match(aces.res[,1], lab2num[,1]),2], aces.res)->aces.res


  if(is.null(aces)==FALSE){
    rates.res[-match(c(tip.rem,nod.rem),rates.res$num),]->rates.clean
    rates.clean[order(rates.clean$num),]->rates.clean
    as.character(rates.clean$num)->rates.clean$num
    rates.clean[1:(Nnode(tree.start)-1),]$num<-seq((Ntip(tree.start)+2),(Ntip(tree.start)+Nnode(tree.start)),1)
    aces.res[-match(nod.rem,aces.res$num),]->aces.clean
    aces.clean$num<-seq((Ntip(tree.start)+1),(Ntip(tree.start)+Nnode(tree.start)),1)
    data.frame(aces.clean[match(names(aces),aces.clean$num),1:2],original=aces,aces.clean[match(names(aces),aces.clean$num),3:6])->nod.sel.res
    if(is.null(brownian_tree)) res<-list(progress=STprogress,rates_tree=rates_tree,rates=rates.clean,aces=aces.clean,ace.prior.values=nod.sel.res) else
      res<-list(progress=STprogress,rates_tree=rates_tree,rates=rates.clean,aces=aces.clean,brownian_tree=brownian_tree,ace.prior.values=nod.sel.res)
  }else{
    if(is.null(brownian_tree)) res<-list(progress=STprogress,rates_tree=rates_tree,rates=rates.res,aces=aces.res) else
      res<-list(progress=STprogress,rates_tree=rates_tree,rates=rates.res,aces=aces.res,brownian_tree=brownian_tree)
  }
  return(res)
}
