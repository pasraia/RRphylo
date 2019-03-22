#' @title Plot RRphylo rates at a specified node
#' @usage plotRates(RR,node,export.tiff =TRUE,foldername)
#' @description The function \code{plotRates} plots the \code{\link{RRphylo}} rates computed for a given clade as compared to the rates computed for the rest of the tree.
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param node the node subtending the clade of interest.
#' @param export.tiff if \code{TRUE} the function save a "rate_bars.tiff" file inside the working directory. It is \code{TRUE} by default.
#' @param foldername the path of the folder where plots are to be found.
#' @export
#' @importFrom grDevices tiff
#' @importFrom graphics abline barplot legend
#' @return The function produces two barplots. On the left side, the rates (in absolute values) computed for the focal clade (blue) are plotted against the rates of the rest of the tree (green). On the right side, the absolute rates of individual branches of the focal clade are collated in increasing rate value (blue bars), and contrasted to the average rate computed over the rest of the tree branches (the vertical red line). It also returns the rate values for both nodes and species descending from the focal node.
#' @examples
#' data("DataApes")
#' DataApes$PCstage->PCstage
#' DataApes$Tstage->Tstage
#'
#'\donttest{
#' RRphylo(tree=Tstage,y=PCstage)->RR
#'
#' plotRates(RR,node=72,foldername=tempdir(),export.tiff = TRUE)
#' }




plotRates<-function(RR, node, export.tiff=TRUE,foldername){
  #require(phytools)

  RR$rates->FULLrates
  RR$tree->rr3
  getDescendants(rr3,node)->shift331

  c(FULLrates[match(shift331[shift331>Ntip(rr3)],rownames(FULLrates)),],FULLrates[match(rr3$tip.label[shift331[shift331<=Ntip(rr3)]],rownames(FULLrates)),])->shift.rates
  shift.rates[order(shift.rates,decreasing=TRUE)]->shift.rates



  if(export.tiff==TRUE)
  {
    tiff(paste(foldername,"rate_bars.tiff",sep="/"),
         width=2000,
         height=1500,
         pointsize=4,
         res=600,
         compression="lzw")
    par(mfrow=c(1,2))
    hist(log(abs(shift.rates)))->H1
    hist(log(abs(FULLrates[-match(names(shift.rates),rownames(FULLrates)),])))->H2
    log(abs(FULLrates))[(log(abs(FULLrates))!="-Inf")]->Xa
    par(mar=c(4,4,4,4))
    plot(H2,col=rgb(0,1,0,.75),xlim=c(range(Xa)[1]-diff(range(Xa))*.1,range(Xa)[2]+diff(range(Xa))*.1),ylim=c(0,max(H2$count)*1.5),xlab="log absolute rates",main="")
    plot(H1,col=rgb(0,0,1,1),add=TRUE)
    legend("topleft", c("back rates", "shift node rates"), col=c(rgb(0,1,0,.75),rgb(0,0,1,1)), box.lwd = 0,box.col = "white",bg = "white",lwd=10)
    par(las=2)
    par(mar=c(4,10,4,4))
    barplot(shift.rates,xlab="rates",horiz=TRUE,col="blue",names.arg=names(shift.rates), main="",border="red")
    abline(v=mean(FULLrates),col="red",lwd=3)

    dev.off()
    return(shift.rates)
  }else{

    par(mfrow=c(1,2))
    hist(log(abs(shift.rates)))->H1
    hist(log(abs(FULLrates[-match(names(shift.rates),rownames(FULLrates)),])))->H2
    log(abs(FULLrates))[(log(abs(FULLrates))!="-Inf")]->Xa
    par(mar=c(4,4,4,4))
    plot(H2,col=rgb(0,1,0,.75),xlim=c(range(Xa)[1]-diff(range(Xa))*.1,range(Xa)[2]+diff(range(Xa))*.1),ylim=c(0,max(H2$count)*1.3),xlab="log absolute rates",main="")
    plot(H1,col=rgb(0,0,1,1),add=TRUE)
    legend("topleft", c("back rates", "shift node rates"), col=c(rgb(0,1,0,.75),rgb(0,0,1,1)), box.lwd = 0,box.col = "white",bg = "white",lwd=10)
    par(las=2)
    par(mar=c(4,10,4,4))
    barplot(shift.rates,xlab="rates",horiz=TRUE,col="blue",names.arg=names(shift.rates), main="",border="red")
    abline(v=mean(FULLrates),col="red",lwd=3)

    return(shift.rates)
  }
}
