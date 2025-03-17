#' @export
print.RRphyloList<-function(x,...){
  if(inherits(x,"list")){
    hid <- attr(x, "hidden")
    print(x[!names(x) %in% hid])
  }else{
    if("mean.sampling"%in%attributes(x)[[1]]) cat(paste(length(x[[2]]),"overfitRR simulations",sep=" ")) else
      if("lambda"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"RRphylo outputs",sep=" ")) else
        if("single.clades"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"search.shift - clade type - outputs",sep=" "))else
          if("state.results"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"search.shift - sparse type - outputs",sep=" "))else
            if("trend.data"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"search.trend outputs",sep=" ")) else
              if("node pairs"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"search.conv - node type - outputs",sep=" "))else
                if("state.res"%in%attributes(x[[1]])[[1]]) cat("List of",paste(length(x),"search.conv - state type - outputs",sep=" "))else
                  if(class(x[[1]])[1]%in%c("gls","procD.lm")) cat("List of",paste(length(x),"PGLS_fossil outputs",sep=" ")) else
                    cat("List of",paste(length(x),"outputs",sep=" "))
  }
}
