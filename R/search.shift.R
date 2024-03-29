#' @title Locating shifts in phenotypic evolutionary rates
#'
#' @usage search.shift(RR, status.type = c("clade", "sparse"),node = NULL, state
#'   = NULL, cov = NULL, nrep = 1000, f = NULL,foldername=NULL,filename)
#' @description The function \code{search.shift} (\cite{Castiglione et al.
#'   2018}) tests whether individual clades or isolated tips dispersed through
#'   the phylogeny evolve at different \code{\link{RRphylo}} rates as compared
#'   to the rest of the tree. Instances of rate shifts may be automatically
#'   found.
#' @param RR an object fitted by the function \code{\link{RRphylo}}.
#' @param status.type whether the \code{"clade"} or \code{"sparse"} condition
#'   must be tested.
#' @param node under the \code{"clade"} condition, the node (clade) to be tested
#'   for the rate shift. When multiple nodes are tested, they need to be written
#'   as in the example below. If \code{node} is left unspecified, the function
#'   performs under the 'auto-recognize' feature, meaning it will automatically
#'   test individual clades for deviation of their rates from the background
#'   rate of the rest of the tree (see details).
#' @param state the state of the tips specified under the \code{"sparse"}
#'   condition.
#' @param cov the covariate to be indicated if its effect on rate values must be
#'   accounted for. Contrary to \code{RRphylo}, \code{cov} needs to be as long
#'   as the number of tips of the tree.
#' @param nrep the number of simulations to be performed for the rate shift
#'   test, by default \code{nrep} is set at 1000.
#' @param f the size of the smallest clade to be tested. By default, nodes
#'   subtending to one tenth of the tree tips are tested.
#' @param foldername has been deprecated; please see the argument
#'   \code{filename} instead.
#' @param filename a character indicating the name of the pdf file and the path
#'   where it is to be saved. If no path is indicated the file is stored in the
#'   working directory
#' @importFrom graphics symbols mtext
#' @importFrom stats sd
#' @importFrom utils globalVariables
#' @export
#' @seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#' @details The function \code{search.shift} takes the object produced by
#'   \code{\link{RRphylo}}. Two different conditions of rate change can be
#'   investigated. Under the \code{"clade"} condition the vector of node or
#'   nodes subjected to the shift must be provided. Alternatively, under the
#'   \code{"sparse"} case the (named) vector of states (indicating which tips
#'   are or are not evolving under the rate shift according to the tested
#'   hypothesis) must be indicated. In the \code{"clade"} case, the function may
#'   perform an 'auto-recognize' feature. Under such specification, the function
#'   automatically tests individual clades (from clades as large as one half of
#'   the tree down to a specified clade size) for deviation of their rates from
#'   the background rate of the rest of the tree, which is identical to the
#'   \code{"clade"} case. An inclusive clade with significantly high rates is
#'   likely to include descending clades with similarly significantly high
#'   rates. Hence, with 'auto-recognize' the \code{search.shift} function is
#'   written as to scan clades individually and to select only the node
#'   subtending to the highest difference in mean absolute rates as compared to
#'   the rest of the tree. A plot of the tree highlighting nodes subtending to
#'   significantly different rates is automatically produced. Caution must be
#'   put into interpreting the 'auto-recognize' results. For instance, a clade
#'   with small rates and another with large rates could be individuated even
#'   under BM. This does not mean these clades are actual instances for rate
#'   shifts. Clades must be tested on their own without the 'auto-recognize'
#'   feature, which serves as guidance to the investigator, when no strong a
#'   priori hypothesis to be tested is advanced. The 'auto-recognize' feature is
#'   not meant to provide a test for a specific hypothesis. It serves as an
#'   optional guidance to understand whether and which clades show significantly
#'   large or small rates as compared to the rest of the tree. Individual clades
#'   are tested at once, meaning that significant instances of rate variation
#'   elsewhere on the tree are ignored. Conversely, running the \code{"clade"}
#'   condition without including the 'auto-recognize' feature, multiple clades
#'   presumed to evolve under the same shift are tested together, meaning that
#'   their rates are collectively contrasted to the rest of the tree, albeit
#'   they can additionally be compared to each other upon request. Under both
#'   the \code{"clade"} and \code{"sparse"} conditions, multiple clades could be
#'   specified at once, and optionally tested individually (for deviation of
#'   rates) against the rates of the rest of the tree and against each other.
#'   The histogram of random differences of mean rates distribution along with
#'   the position of the real difference of means is automatically generated by
#'   \code{search.shift}. Regardless of which condition is specified, the
#'   function output produces the real difference of means, and their
#'   significance value.
#' @return Under all circumstances, \code{search.shift} provides a vector of
#'   \code{$rates}. If \code{'cov'} values are provided, rates are regressed
#'   against the covariate and the residuals of such regression form the vector
#'   \strong{\code{$rates}}. Otherwise, \strong{\code{$rates}} are the same
#'   rates as with the \code{RR} argument.
#' @return Under \code{"clade"} case without specifying nodes (i.e.
#'   'auto-recognize') a list including:
#' @return \strong{$all.clades} for each detected node, the data-frame includes
#'   the average rate difference (computed as the mean rate over all branches
#'   subtended by the node minus the average rate for the rest of the tree) and
#'   the probability that it do represent a real shift. Probabilities are
#'   contrasted to simulations shuffling the rates across the tree branches for
#'   a number of replicates specified by the argument \code{nrep}. Note that the
#'   p-values refer to the number of times the real average rates are larger (or
#'   smaller) than the rates averaged over the rest of the tree, divided by the
#'   number of simulations. Hence, large rates are significantly larger than the
#'   rest of the tree (at alpha = 0.05), when the probability is > 0.975; and
#'   small rates are significantly small for p < 0.025.
#' @return \strong{$single.clades} the same as with 'all.clades' but restricted
#'   to the largest/smallest rate values along a single clade (i.e. nested
#'   clades with smaller rate shifts are excluded). Large rates are
#'   significantly larger than the rest of the tree (at alpha = 0.05), when the
#'   probability is > 0.975; and small rates are significantly small for p <
#'   0.025.
#' @return Under \code{"clade"} condition by specifying the \code{node}
#'   argument:
#' @return \strong{$all.clades.together} if more than one node is tested, this
#'   specifies the average rate difference and the significance of the rate
#'   shift, by considering all the specified nodes as evolving under a single
#'   rate. As with the 'auto-recognize' feature, large rates are significantly
#'   larger than the rest of the tree (at alpha = 0.05), when the probability is
#'   > 0.975; and small rates are significantly small for p < 0.025.
#' @return \strong{$single.clades} this gives the significance for individual
#'   clades, tested separately. As previously, large rates are significantly
#'   larger than the rest of the tree (at alpha = 0.05), when the probability is
#'   > 0.975; and small rates are significantly small for p < 0.025.
#' @return Under the \code{"sparse"} condition:
#' @return   \strong{$state.results} for each state, the data-frame includes the
#'   average rate difference (computed as the mean rate over all leaves evolving
#'   under a given state, minus the average rate for each other state or the
#'   rest of the tree) and the probability that the shift is real. Large rates
#'   are significantly larger (at alpha = 0.05), when the probability is >
#'   0.975; and small rates are significantly small for p < 0.025. States are
#'   compared pairwise.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#' Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method for
#' testing evolutionary rate variation and shifts in phenotypic evolution.
#' \emph{Methods in Ecology and Evolution}, 9:
#' 974-983.doi:10.1111/2041-210X.12954
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#'
#'
#' RRphylo(tree=treedino,y=massdino)->dinoRates
#'
#' # Case 1. Without accounting for the effect of a covariate
#'
#' # Case 1.1 "clade" condition
#' # with auto-recognize
#' search.shift(RR=dinoRates,status.type="clade",
#'              filename=paste(tempdir(), "SSauto", sep="/"))
#' # testing two hypothetical clades
#' search.shift(RR=dinoRates,status.type="clade",node=c(696,746),
#'              filename=paste(tempdir(), "SSclade", sep="/"))
#'
#' # Case 1.2 "sparse" condition
#' # testing the sparse condition.
#' search.shift(RR=dinoRates,status.type= "sparse",state=statedino,
#'              filename=paste(tempdir(), "SSsparse", sep="/"))
#'
#'
#' # Case 2. Accounting for the effect of a covariate
#'
#' # Case 2.1 "clade" condition
#' search.shift(RR=dinoRates,status.type= "clade",cov=massdino,
#'              filename=paste(tempdir(), "SSclade_cov", sep="/"))
#'
#' # Case 2.2 "sparse" condition
#' search.shift(RR=dinoRates,status.type="sparse",state=statedino,cov=massdino,
#'              filename=paste(tempdir(), "SSstate_cov", sep="/"))
#'     }



search.shift<-function(RR,
                       status.type=c("clade","sparse"),
                       node=NULL,
                       state=NULL,
                       cov=NULL,
                       nrep=1000,
                       f=NULL,
                       foldername=NULL,
                       filename)
{
  # require(phytools)
  # require(geiger)
  # require(scales)

  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!missing(foldername)){
    stop("argument foldername is deprecated; please use filename instead.",
         call. = FALSE)
  }

  tree <- RR$tree
  rates <- RR$rates
  betas<-RR$multiple.rates


  if(is.null(f)) f<-round(Ntip(tree)/10)
  if(is.null(cov)){
    rates<-rates
  }else{

    RRphylo(tree,cov)->RRcova
    abs(c(RRcova$aces,cov))->Y
    c(rownames(RRcova$aces),names(cov))->names(Y)

    if (length(betas)>length(rates)) {
      if(length(which(apply(betas,1,sum)==0))>0){
        which(apply(betas,1,sum)==0)->zeroes
        log(abs(betas))->R
        R[-zeroes,]->R
        Y[-zeroes]->Y

        residuals(lm(R~Y))->res
        which(apply(betas,1,sum)!=0)->factOut
        betas[factOut,]<-res
        betas[zeroes,]<-0

      }else {
        log(abs(betas))->R
        abs(Y)->Y
        residuals(lm(R~Y))->res
        as.matrix(res)->betas
      }
      rates<-betas
      rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
      rates <- as.matrix(rates)

    }else{

      if(length(which(betas=="0"))>0){
        which(betas=="0")->zeroes
        log(abs(betas))->R
        R[-zeroes]->R
        Y[-zeroes]->Y


        residuals(lm(R~Y))->res
        which(betas!="0")->factOut

        betas[factOut]<-res
        betas[zeroes]<-0
      } else {
        log(abs(betas))->R
        abs(Y)->Y
        residuals(lm(R~Y))->res
        as.matrix(res)->betas
      }
      betas->rates
    }

  }

  if (status.type == "clade") {
    if (is.null(node)) {
      ST <- subtrees(tree)
      len <- array()
      for (i in 1:length(ST)) {
        len[i] <- Ntip(ST[[i]])
      }
      st <- ST[which(len < (Ntip(tree)/2) & len > round(f))]
      node <- sapply(st, function(x) getMRCA(tree, x$tip.label))
      names(st) <- node
      leaf2N.diff <- array()
      p.single <- array()
      for (j in 1:length(node)) {
        Cbranch <- getDescendants(tree, node[j])
        Ctips <- tips(tree, node[j])
        Cleaf <- c(Cbranch, Ctips)
        leaf.rates <- rates[match(Cleaf, rownames(rates)),]
        leaf.rates <- na.omit(leaf.rates)
        NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
        leafR <- mean(abs(leaf.rates))
        NCR <- mean(abs(NCrates))
        leaf2N.diff[j] <- leafR - NCR
        NC <- length(rates) - length(leaf.rates)
        C <- length(leaf.rates)
        ran.diffM <- array()
        for (i in 1:nrep) {
          ran.diffM[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates),
                                                                    NC))
        }
        p.single[j] <- rank(c(leaf2N.diff[j], ran.diffM[1:(nrep -
                                                             1)]))[1]/nrep
      }
      names(leaf2N.diff) <- node
      names(p.single) <- node


      if (length(p.single[p.single >= 0.975 | p.single <=
                          0.025])==0){
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),
                                      names(leaf2N.diff))]
        p.single <- NULL
        leaf2N.diff <- NULL
      }


      if (length(p.single[p.single >= 0.975 | p.single <=
                          0.025])==1) {
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),
                                      names(leaf2N.diff))]

        p.single <- p.single[p.single >= 0.975 | p.single <=
                               0.025]
        leaf2N.diff <- leaf2N.diff[match(names(p.single),
                                         names(leaf2N.diff))]

      }

      if (length(p.single[p.single >= 0.975 | p.single <=
                          0.025]) >= 2)  {
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),
                                      names(leaf2N.diff))]

        p.single <- p.single[p.single >= 0.975 | p.single <= 0.025]

        leaf2N.diff <- leaf2N.diff[match(names(p.single), names(leaf2N.diff))]

        ups <- p.single[p.single > 0.975]
        dws <- p.single[p.single < 0.025]
        ups <- ups[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                         decreasing = FALSE)]), names(ups)))]
        dws <- dws[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                         decreasing = FALSE)]), names(dws)))]
        if (is.na(mean(dws))) {
          dws = Nnode(tree) * 2
        } else {
          s = 1
          repeat
          {
            d <- which(names(dws) %in% getDescendants(tree, names(dws)[s]))
            if (length(d) > 0) {
              leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))]->cla
              names(which.max(abs(leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))])))->IN
              dws[-match(names(cla[which(names(cla)!=IN)]),names(dws))]->dws
              s=1
            } else {
              dws <- dws
              s = s+1

            }

            if (s > length(dws))  break

          }
        }

        if (is.na(mean(ups))) {
          ups = Nnode(tree) * 2
        } else {

          z = 1
          repeat
          {
            d <- which(names(ups) %in% getDescendants(tree, names(ups)[z]))
            if (length(d) > 0) {
              leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))]->cla
              names(which.max(abs(leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))])))->IN
              ups[-match(names(cla[which(names(cla)!=IN)]),names(ups))]->ups
              z=1
            } else {
              ups <- ups
              z = z+1

            }
            if (z > length(ups))  break

          }
        }


        # p.init <- p.single
        # l2N.init <- leaf2N.diff[match(names(p.init),
        #                               names(leaf2N.diff))]
        p.single <- p.single[which(names(p.single) %in%
                                     names(c(ups, dws)))]

        leaf2N.diff <- leaf2N.diff[match(names(p.single),
                                         names(leaf2N.diff))]

        p.single[order(p.single)]->p.single
      }



      pdf(file=paste(filename,".pdf",sep=""))
      #if(length(which(p.single<=0.025|p.single>=0.975))>0){
      if(!is.null(p.single)){
        if(Ntip(tree)>100) plot(tree, show.tip.label = FALSE) else plot(tree, cex=.8)
        xy <- list()
        for (w in 1:length(p.single)) {
          xy[[w]] <- unlist(sapply(get("last_plot.phylo",
                                       envir =ape::.PlotPhyloEnv), function(x) x[as.numeric(names(p.single)[w])]))[c(21,
                                                                                                                     22)]
        }


        c(rep("red",length(which(p.single<=0.025))),rep("royalblue",length(which(p.single>=0.975))))->p.col
        symbols(lapply(xy, "[[", 1), lapply(xy, "[[", 2),
                circles = abs(leaf2N.diff[match(names(p.single),names(leaf2N.diff))])^0.5, inches = 0.25,
                add = TRUE, bg = scales::alpha(p.col, 0.5), fg = p.col)

        nodelabels(node = as.numeric(names(p.single)), adj = c(1.5,
                                                               1), text = names(p.single), frame = "none", bg = "white",
                   col = "purple")

      }else{
        if(Ntip(tree)>100) plot(tree, show.tip.label = FALSE) else plot(tree, cex=.8)
      }
      dev.off()
      data.frame("rate difference"=l2N.init[match(names(p.init),names(l2N.init))],"p-value"=p.init)->allres
      if(!is.null(p.single))
        data.frame("rate difference"=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],"p-value"=p.single)->single else single<-NULL
      res<-list(allres,single)
      names(res)<-c("all.clades","single.clades")
      if(!is.null(cov)) {
        res<-c(res,list(rates))
        names(res)[3]<-"rates"
      }
    } else {
      node = node
      Cbranch <- list()
      for (i in 1:length(node)) {
        Cbranch[[i]] <- getDescendants(tree, node[i])
      }
      Cbranch <- unlist(Cbranch)
      Cbranch <- Cbranch[-which(Cbranch < Ntip(tree))]
      Ctips <- list()
      for (i in 1:length(node)) {
        Ctips[[i]] <- tips(tree, node[i])
      }
      Ctips <- unlist(Ctips)
      Ctips <- unique(Ctips)
      Cbranch <- unique(Cbranch)
      Cleaf <- c(Cbranch, Ctips)
      leaf.rates <- rates[match(Cleaf, rownames(rates)),
                          ]
      leaf.rates <- na.omit(leaf.rates)
      NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
      leafR <- mean(abs(leaf.rates))
      NCR <- mean(abs(NCrates))
      leaf2NC.diff <- leafR - NCR
      NC <- length(rates) - length(leaf.rates)
      C <- length(leaf.rates)
      ran.diffR <- array()
      for (i in 1:nrep) {
        ran.diffR[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates),
                                                                  NC))
      }
      p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep -
                                                     1)]))[1]/nrep
      if(length(node)==1){
        pdf(file=paste(filename,".pdf",sep=""),width=8.3,height=8.3)
        par(mar = c(3, 2, 2, 1))
        hist(ran.diffR, main="",cex.lab=1.5,
             yaxt="n",xaxt="n",ylab=paste("node", node, sep = " "),xlab="",mgp=c(0.2,0,0),
             xlim = c(1.1 * min(c(leaf2NC.diff,ran.diffR)),1.1 * max(c(leaf2NC.diff,ran.diffR))))
        hist(c(leaf2NC.diff,ran.diffR),plot=FALSE)->hi
        if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
          c(hi$breaks[seq(1,length(hi$breaks),2)],
            (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
        axis(1,at=athi)
        mtext(text="random differences",1,line=1.8)
        title(main="Absolute rate difference",cex.main=2)
        abline(v = leaf2NC.diff, col = "green", lwd = 3)
        dev.off()
        # par(mar = c(1, 1, 1, 1))
        # hist(ran.diffR, xlab = "random differences",
        #      main = "selected clade", xlim = c(2.5 * range(ran.diffR)[1],
        #                                        2.5 * range(ran.diffR)[2]))
        # abline(v = leaf2NC.diff, col = "green", lwd = 3)
        names(leaf2NC.diff)<-names(p.shift)<-node
        res<-list(data.frame("rate difference"=leaf2NC.diff,"p-value"=p.shift))
        names(res)<-"single.clade"
        if(!is.null(cov)) {
          res<-c(res,list(rates))
          names(res)[2]<-"rates"
        }
      }else{
        pdf(file=paste(filename,".pdf",sep=""),width=8.3,height=11.7)
        par(mfrow = c(length(node) + 1, 1))
        par(mar = c(3, 2, 2, 1))
        hist(ran.diffR, main="",cex.lab=1.5,
             yaxt="n",ylab="All clades together",xlab="",mgp=c(0.2,0.5,0),xaxt="n",
             xlim = c(1.1 * min(c(leaf2NC.diff,ran.diffR)),1.1 * max(c(leaf2NC.diff,ran.diffR))))
        hist(c(leaf2NC.diff,ran.diffR),plot=FALSE)->hi
        if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
          c(hi$breaks[seq(1,length(hi$breaks),2)],
            (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
        axis(1,at=athi)
        mtext(text="random differences",1,line=1.8)
        title(main="Absolute rate difference",cex.main=2)
        abline(v = leaf2NC.diff, col = "green", lwd = 3)


        # par(mar = c(1, 1, 1, 1))
        # par(mfrow = c(length(node) + 1, 1))
        # hist(ran.diffR, xlab = "random differences",
        #      main = "all clades", xlim = c(2.5 * range(ran.diffR)[1],
        #                                    2.5 * range(ran.diffR)[2]))
        # abline(v = leaf2NC.diff, col = "green", lwd = 3)
        leaf2N.diff <- array()
        p.single <- array()
        ran.diff <- list()
        for (i in 1:length(node)) {
          NOD <- node[-i]
          others <- list()
          mommies <- list()
          for (j in 1:length(NOD)) {
            others[[j]] <- tips(tree, NOD[j])
            mommies[[j]] <- getDescendants(tree, NOD[j])
          }
          others <- unlist(others)
          mommies <- unlist(mommies)
          mommies <- mommies[-which(mommies < Ntip(tree) +
                                      1)]
          otmom <- c(mommies, others)
          Ctips <- tips(tree, node[i])
          Ctips <- unlist(Ctips)
          Cbranch <- getDescendants(tree, node[i])
          Cleaf <- c(Cbranch, Ctips)
          leaf.rates <- rates[match(Cleaf, rownames(rates)),
                              ]
          leaf.rates <- na.omit(leaf.rates)
          NC <- rates[-c(which(rownames(rates) %in% names(leaf.rates)),
                         which(rownames(rates) %in% otmom)), ]
          NR.r <- mean(abs(NC))
          leaf.r <- mean(abs(leaf.rates))
          leaf2N.diff[i] <- leaf.r - NR.r
          NC.l <- length(NC)
          leaf.l <- length(leaf.rates)
          tot.r <- abs(c(NC, leaf.rates))
          RAN.diff <- array()
          for (k in 1:nrep) {
            RAN.diff[k] <- mean(sample(tot.r, leaf.l)) -
              mean(sample(tot.r, NC.l))
            ran.diff[[i]] <- RAN.diff
          }
          p.single[i] <- rank(c(leaf2N.diff[i], RAN.diff[1:(nrep -
                                                              1)]))[1]/nrep
        }
        names(p.single) <- node
        names(leaf2N.diff) <- names(p.single)
        for (m in 1:length(node)) {
          par(mar = c(3, 2, 2, 1))
          hist(ran.diff[[m]], main="",cex.lab=1.5,
               yaxt="n",ylab=paste("node", node[m], sep = " "),xlab="",mgp=c(0.2,0.5,0),xaxt="n",
               xlim = c(1.1 * min(c(leaf2N.diff[m],ran.diff[[m]])),1.1 * max(c(leaf2N.diff[m],ran.diff[[m]]))))
          hist(c(leaf2N.diff[m],ran.diff[[m]]),plot=FALSE)->hi
          if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
            c(hi$breaks[seq(1,length(hi$breaks),2)],
              (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
          axis(1,at=athi)
          mtext(text="random differences",1,line=1.8)
          abline(v = leaf2N.diff[m], col = "blue", lwd = 3)

          # hist(ran.diff[[m]], xlab = "random differences",
          #      main = print(paste("Node", node[m], sep = " ")),
          #      xlim = c(2.5 * range(ran.diff[[m]])[1], 2.5 *
          #                 range(ran.diff[[m]])[2]))
          # abline(v = leaf2N.diff[m], col = "blue", lwd = 3)
        }
        dev.off()
        res<-list(data.frame("rate.difference"=leaf2NC.diff,"p.value"=p.shift),
                  data.frame("rate difference"=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],"p-value"=p.single))
        names(res)<-c("all.clades.together","single.clades")
        if(!is.null(cov)) {
          res<-c(res,list(rates))
          names(res)[3]<-"rates"
        }
      }
    }
  } else {
    state <- treedata(tree, state, sort = TRUE)[[2]][,1]
    frame <- data.frame(status = as.factor(state), rate = rates[match(names(state),
                                                           rownames(rates))])
    p.status.diff <- array()
    if (length(unique(state)) > 2) {
      status.diff <- apply(combn(tapply(abs(frame$rate),
                                        frame$status, mean), 2), 2, diff)
      sta <- tapply(abs(frame$rate), frame$status, mean)
      sta <- sta[match(unique(state), names(sta))]
      w <- array()
      for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(frame$status ==
                                                                        names(sta)[x]), 2]))
      names(w) <- names(sta)
      status.diff <- c(status.diff, w)
      names(status.diff) <- c(apply(combn(levels(frame$status),
                                          2), 2, function(x) paste(x[2], x[1], sep = "_")),
                              names(sta))
      status.diffS <- matrix(ncol = length(status.diff),
                             nrow = nrep)
      for (i in 1:nrep) {
        s.state <- frame$status
        s.ran <- sample(s.state)
        s.frame <- data.frame(s.ran, frame$rate)
        SD <- apply(combn(tapply(abs(s.frame$frame.rate),
                                 s.frame$s.ran, mean), 2), 2, diff)
        sta <- tapply(abs(frame$rate), s.ran, mean)
        sta <- sta[match(unique(state), names(sta))]
        w <- array()
        for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(s.frame$s.ran ==
                                                                          names(sta)[x]), 2]))
        status.diffS[i, ] <- c(SD, w)
      }
      colnames(status.diffS) <- names(status.diff)
      pdf(file=paste(filename,".pdf",sep=""),width=8.3,height=11.7)
      par(mfrow = c(length(unique(state)), 1))
      idx <- match(unique(state), colnames(status.diffS))
      for (i in 1:length(idx)) {
        par(mar = c(3, 2, 2, 1))
        hist(status.diffS[, idx[i]], main="",xaxt="n",cex.lab=1.5,
             yaxt="n",ylab=paste("state", colnames(status.diffS)[idx[i]], sep = " "),xlab="",mgp=c(0.2,0.5,0),
             xlim = c(1.1*min(c(status.diff[idx[i]],status.diffS[, idx[i]])),1.1*max(c(status.diff[idx[i]],status.diffS[, idx[i]]))))
        hist(c(status.diff[idx[i]],status.diffS[, idx[i]]),plot=FALSE)->hi
        if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
          c(hi$breaks[seq(1,length(hi$breaks),2)],
            (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
        axis(1,at=athi)
        mtext(text="random differences",1,line=1.8)
        if(i==1) title(main="Absolute rate difference between states",cex.main=2)
        abline(v = status.diff[idx[i]], lwd = 3, col = "green")

        #   hist(status.diffS[, idx[i]], xlab = "random differences",
        #        main = print(paste("rate difference per status",
        #                           colnames(status.diffS)[idx[i]], sep = " ")),
        #        xlim = c(min(status.diffS[, idx[i]]) - sd(status.diffS[,
        #                                                               idx[i]]), max(status.diffS[, idx[i]]) + sd(status.diffS[,
        #                                                                                                                       idx[i]])))
        #   abline(v = status.diff[idx[i]], lwd = 3, col = "green")
      }
      dev.off()
      for (i in 1:length(status.diff)) p.status.diff[i] <- rank(c(status.diff[i],
                                                                  status.diffS[1:(nrep - 1), i]))[1]/nrep
      names(p.status.diff) <- names(status.diff)
      unlist(p.status.diff)->p.status.diff
      unlist(status.diff)->status.diff

      res<-list(data.frame("rate difference"=status.diff[match(names(p.status.diff),names(status.diff))],p.status.diff))
      names(res)<-"state.results"
      if(!is.null(cov)) {
        res<-c(res,list(rates))
        names(res)[2]<-"rates"
      }
    } else {
      status.diff <- diff(tapply(abs(frame$rate), state,
                                 mean))
      status.diffS <- array()
      for (i in 1:nrep) {
        s.state <- frame$status
        s.ran <- sample(s.state)
        s.frame <- data.frame(s.ran, frame$rate)
        s.frame[, 1] <- as.factor(s.frame[, 1])
        status.diffS[i] <- diff(tapply(abs(s.frame$frame.rate),
                                       s.frame$s.ran, mean))
      }
      # hist(status.diffS, xlab = "random differences", main = "rate difference per status",
      #      xlim = c(min(status.diffS) * 2.5, max(status.diffS) *
      #                 2.5))
      # abline(v = status.diff, lwd = 3, col = "green")
      pdf(file=paste(filename,".pdf",sep=""),width=8.3,height=8.3)
      par(mar = c(3, 2, 2, 1))
      hist(status.diffS, main="",xaxt="n",
           yaxt="n",ylab="",xlab="random differences",
           cex.lab=1,mgp=c(1.5,0.8,0),
           xlim = c(1.1*min(c(status.diff,status.diffS)), 1.1*max(c(status.diff,status.diffS))))
      hist(c(status.diff,status.diffS),plot=FALSE)->hi
      if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
        c(hi$breaks[seq(1,length(hi$breaks),2)],
          (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
      axis(1,at=athi,mgp=c(0.2,0.5,0))
      title(main="Absolute rate difference between states",cex.main=2)
      abline(v = status.diff, lwd = 3, col = "green")
      dev.off()
      p.status.diff <- rank(c(status.diff, status.diffS[1:(nrep -
                                                             1)]))[1]/nrep
      state.results<-data.frame("rate difference"=status.diff[match(names(p.status.diff),names(status.diff))],"p.value"=p.status.diff)
      rownames(state.results)<-paste(names(p.status.diff),unique(state)[which(unique(state)!=names(p.status.diff))],sep="-")
      res<-list(state.results)
      names(res)<-c("state.results")
      if(!is.null(cov)) {
        res<-c(res,list(rates))
        names(res)[2]<-"rates"
      }
    }
  }
  return(res)
}
