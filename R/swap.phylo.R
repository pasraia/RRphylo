#' @rdname RRphylo-defunct
#' @export
swap.phylo<-function(){

  .Defunct(msg="The function swap.phylo has been removed from the package.
           Check the function overfitRR to test sampling effects and phylogenetic uncertainty on RRphylo methods")
    # maxN <- function(x, N=2){
    #   len <- length(x)
    #   if(N>len){
    #     warning('N greater than length(x).  Setting N=length(x)')
    #     N <- length(x)
    #   }
    #   sort(x,partial=len-N+1)[len-N+1]
    # }
    #
    # if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    #   data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    #   tree$tip.label<-tips(tree,(Ntip(tree)+1))
    #   data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    #   tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
    # }
    #
    # tree->tree1
    # #if(inherits(y,"data.frame")) treedata(tree,y,sort=TRUE)[[2]]->y
  # if(is.null(nrow(y))) y <- treedata(tree, y, sort = TRUE)[[2]][,1] else y <- treedata(tree, y, sort = TRUE)[[2]]
    # Rrts<-list()
    # if(is.null(cov))
    # {
    #   replicate(nrep,RRphylo(swapONE(tree,si,si2,plot.swap=FALSE)[[1]],y,clus=clus)$rates)->Rrts
    # }else{
    #   RRphylo(tree,cov)->RRcov
    #   c(RRcov$aces,cov)->cova
    #   names(cova)<-c(rownames(RRcov$aces),names(cov))
    #   replicate(nrep,RRphylo(swapONE(tree,si,si2,plot.swap=FALSE)[[1]],y,cov=cova,clus=clus)$rates)->Rrts
    # }
    #
    # for(i in 1:dim(Rrts)[3]) Rrts[,,i][match(rownames(rts),names(Rrts[,,i]))]
    #
    # apply(Rrts,3,cbind)->multirates
    # rownames(multirates)<-rownames(rts)
    # apply(multirates,1,function(x) mean(abs(x)))->Rrts.distrib
    # names(Rrts.distrib)<-rownames(rts)
    #
    # if(length(y)>Ntip(tree)){
    #   rts[match(rownames(y),rownames(rts)),]->rts
    #   Rrts.distrib[match(rownames(y),names(Rrts.distrib))]->Rrts.distrib
    #   multirates[match(rownames(y),rownames(multirates)),]->multirates
    # }else{
    #   rts[match(names(y),rownames(rts)),]->rts
    #   Rrts.distrib[match(names(y),names(Rrts.distrib))]->Rrts.distrib
    #   multirates[match(names(y),rownames(multirates)),]->multirates
    # }
    #
    # hist(log(abs(rts[match(tips(tree,node),names(rts))])))->H3
    # hist(log(abs(rts[-match(tips(tree,node),names(rts))])))->H4
    #
    # hist(log(abs(Rrts.distrib[match(tips(tree,node),names(Rrts.distrib))])))->H1
    # hist(log(abs(Rrts.distrib[-match(tips(tree,node),names(Rrts.distrib))])))->H2
    #
    # Rrts[,1,]->rrts
    # par(mfrow=c(1,2))
    #
    # plot(H3,xlim=c(min(H4$bre)-2,max(H3$bre)+2),col=rgb(1,0,0,.5),ylim= c(0,max(H4$counts)),main="rates density distribution",xlab="rates")
    # plot(H4,col=rgb(0,0,1,.5),add=TRUE)
    # plot(H1,xlim=c(min(H2$bre)-2,max(H1$bre)+2),col=rgb(1,0,1,.5),ylim= c(0,max(H2$counts)),main="rates density distribution",xlab="rates")
    # plot(H2,col=rgb(0,1,0,.5),add=TRUE)
    # t.test(log(abs(Rrts.distrib[-match(tips(tree,node),names(Rrts.distrib))])),log(abs(Rrts.distrib[match(tips(tree,node),names(Rrts.distrib))])))->p.swap
    # return(list("p.swap"=p.swap,"rates"=multirates))

  }
