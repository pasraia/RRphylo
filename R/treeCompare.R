#' @title Visualize the difference between phylogenetic trees
#' @description The function scans a pair of phylogenetic trees to find
#'   topological differences.
#' @usage treeCompare(tree,tree1,focal=NULL,plot=TRUE)
#' @param tree,tree1 a phylogenetic tree. The tree needs not to be ultrametric
#'   and fully dichotomous. Generic name and specific epithet must be separated
#'   by '_'.
#' @param focal a vector of focal species to search on both \code{tree} and \code{tree1}.
#' @param plot if \code{TRUE}, the function produces an interactive plotting
#'   device to check differences in species placement at the genus level.
#' @export
#' @return The function returns a data-frame indicating for each un-matching/focal
#'   species its sister species/clades on both trees.
#' @author Giorgia Girardi, Silvia Castiglione, Carmela Serio, Antonella Esposito
#' @examples
#' \dontrun{
#' DataFelids$tree->treefel
#'
#' set.seed(22)
#' drop.tip(treefel,sample(treefel$tip.label,20))->tree.red
#' swapONE(tree.red,si=0.5)[[1]]->tree.red
#'
#' sample(tree.red$tip.label,10)->focal.species
#'
#' treeCompare(treefel,tree.red)->comp
#' treeCompare(treefel,tree.red,focal=focal.species)->comp1
#'
#' }

treeCompare<-function (tree, tree1,focal=NULL, plot = TRUE){

  misspacks<-sapply(c("scales","manipulate"),requireNamespace,quietly=TRUE)
  if(any(!misspacks)){
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)],collapse=", "),
         call. = FALSE)
  }

  if (isTRUE(plot)) {
    mars <- par("mar")
    on.exit(par(mar = mars))
  }
  if (!identical(tree1$edge[tree1$edge[, 2] <= Ntip(tree1),2], seq(1, Ntip(tree1)))){
    tree1$tip.label <- tree1$tip.label[tree1$edge[tree1$edge[,2] <= Ntip(tree1), 2]]
    tree1$edge[tree1$edge[, 2] <= Ntip(tree1), 2] <- seq(1,Ntip(tree1))
  }
  if (!identical(tree$edge[tree$edge[, 2] <= Ntip(tree), 2],seq(1, Ntip(tree)))){
    tree$tip.label <- tree$tip.label[tree$edge[tree$edge[,2] <= Ntip(tree), 2]]
    tree$edge[tree$edge[, 2] <= Ntip(tree), 2] <- seq(1,Ntip(tree))
  }

  if (isFALSE(any(tree$tip.label %in% tree1$tip.label)))
    message("There are no shared species between the trees, no output is returned") else{
      if(identical(tree,tree1)) message("The  phylogenies are exactly the same, no output is returned") else {

        if(!is.null(focal)&&any(!focal%in%tree$tip.label|!focal%in%tree1$tip.label))
          stop("Focal species are missing from one of the trees")

        treeO <- tree
        tree1O <- tree1

        if(!is.null(focal)){
          tree <- drop.tip(tree, which(!tree$tip.label %in% focal))
          tree1 <- drop.tip(tree1, which(!tree1$tip.label %in%focal))
        }else{
          tree <- drop.tip(tree, which(!tree$tip.label %in% tree1$tip.label))
          tree1 <- drop.tip(tree1, which(!tree1$tip.label %in%tree$tip.label))
        }

        sispair <- data.frame(row.names = tree$tip.label)
        sispair$sist <- sapply(tree$tip.label, function(j) getSis(treeO,j, printZoom = FALSE))
        sispair$sis.list <- sapply(sispair$sist, function(ss) {
          if (!all(ss %in% tree$tip.label)) {
            c(ss[which(ss %in% treeO$tip.label)], unlist(sapply(ss[which(!ss %in%treeO$tip.label)], function(k) tips(treeO, k))))
          } else ss
        }, simplify = FALSE)
        sispair$sis.list<-lapply(sispair$sis.list,function(x) sort(x))
        sispair$poly <- sapply(sispair$sist, length)
        sispair$nsis <- sapply(sispair$sis.list, length)

        sist1 <- sapply(tree1$tip.label, function(j) getSis(tree1O,j, printZoom = FALSE))
        sispair$sist1 <- sist1[match(rownames(sispair), names(sist1))]
        sispair$sis.list1 <- sapply(sispair$sist1, function(ss) {
          if (!all(ss %in% tree1O$tip.label)) {
            c(ss[which(ss %in% tree1O$tip.label)], unlist(sapply(ss[which(!ss %in%tree1O$tip.label)], function(k) tips(tree1O, k))))
          }else ss
        }, simplify = FALSE)
        sispair$sis.list1<-lapply(sispair$sis.list1,function(x) sort(x))
        sispair$poly1 <- sapply(sispair$sist1, length)
        sispair$nsis1 <- sapply(sispair$sis.list1, length)

        if(!is.null(focal)) nomatch.sp<-sispair[, c(1, 5)] else {
          nomatch.sp <- sispair[which(sispair$poly !=sispair$poly1|sispair[,4]!=sispair[,8]),]

          sispair <- sispair[which(!rownames(sispair) %in% rownames(nomatch.sp)),]
          nomatch.sp <- rbind(nomatch.sp, sispair[which(!apply(sispair,1, function(x)
            all(x[[2]] %in% x[[6]]) & all(x[[6]] %in%x[[2]]))), ])
          nomatch.sp <- nomatch.sp[, c(1, 5)]
        }


        if (isTRUE(plot)) {
          gentree <- sapply(unique(sapply(strsplit(rownames(nomatch.sp),"_"), "[[", 1)), function(j) getGenus(treeO, j)[,3])
          gentree1 <- sapply(unique(sapply(strsplit(rownames(nomatch.sp),"_"), "[[", 1)), function(j) getGenus(tree1O, j)[,3])
          if(any(gentree<=Ntip(treeO))) gentree[which(gentree<=Ntip(treeO))]<-sapply(gentree[which(gentree<=Ntip(treeO))],function(w) getMommy(treeO,w)[1])
          if(any(gentree1<=Ntip(tree1O))) gentree1[which(gentree1<=Ntip(tree1O))]<-sapply(gentree1[which(gentree1<=Ntip(tree1O))],function(w) getMommy(tree1O,w)[1])
          gentree1 <- gentree1[order(names(gentree))]
          gentree <- gentree[order(names(gentree))]
          font.tree <- rep(3, length(tree$tip.label))
          font.tree[which(tree$tip.label %in% rownames(nomatch.sp))] <- 4
          names(font.tree) <- tree$tip.label
          plotClades <- function(cla.list, cla1.list, colo, colo1,
                                 font.cla, font.cla1) {
            par(mfrow = c(1, 2))
            plot(cla.list, tip.col = colo, no.margin = TRUE,
                 font = font.cla)
            plot(cla1.list, tip.col = colo1, no.margin = TRUE,
                 direction = "leftward", font = font.cla1)
          }
          cla.list <- cla1.list <- colo <- colo1 <- font.cla <- font.cla1 <- list()
          for (w in 1:length(gentree)) {
            cla.list[[w]] <- cla <- extract.clade(treeO, gentree[w])
            cla1.list[[w]] <- cla1 <- extract.clade(tree1O,gentree1[w])
            spcla <- match(rownames(nomatch.sp)[grep(paste(names(gentree)[w],
                                                           "_", sep = ""), rownames(nomatch.sp))], cla$tip.label)
            spcla1 <- match(rownames(nomatch.sp)[grep(paste(names(gentree1)[w],
                                                            "_", sep = ""), rownames(nomatch.sp))], cla1$tip.label)
            colo[[w]] <- rep("black", Ntip(cla))
            colo1[[w]] <- rep("black", Ntip(cla1))
            colo[[w]][spcla] <- colo1[[w]][spcla1] <- (scales::hue_pal())(length(spcla))
            font.cla[[w]] <- font.tree[match(cla$tip.label,
                                             names(font.tree))]
            font.cla1[[w]] <- font.tree[match(cla1$tip.label,
                                              names(font.tree))]
          }
          names(cla.list) <- names(cla1.list) <- names(colo) <- names(colo1) <- names(font.cla) <- names(font.cla1) <- names(gentree)
          taxon <- manipulate::picker(as.list(names(gentree)))
          manipulate::manipulate(do.call(plotClades, list(cla.list = cla.list[[which(names(cla.list)==taxon)]],
                                                          cla1.list = cla1.list[[which(names(cla1.list)==taxon)]],
                                                          colo = colo[[which(names(colo)==taxon)]],
                                                          colo1 = colo1[[which(names(colo1)==taxon)]],
                                                          font.cla = font.cla[[which(names(font.cla)==taxon)]],
                                                          font.cla1 = font.cla1[[which(names(font.cla1)==taxon)]])),taxon=taxon)
        }
        colnames(nomatch.sp) <- c("tree", "tree1")
        nomatch.sp <- nomatch.sp[order(rownames(nomatch.sp)), ]
        return(nomatch.sp)
      }
    }
}

