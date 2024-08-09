#' @title Locating shifts in phenotypic evolutionary rates
#' @usage search.shift(RR, status.type = c("clade", "sparse"),node = NULL, state
#'   = NULL, cov = NULL, nrep = 1000, f = NULL,filename=NULL)
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
#' @param cov the covariate vector to be indicated if its effect on rate values must be
#'   accounted for. Contrary to \code{RRphylo}, \code{cov} needs to be as long
#'   as the number of tips of the tree.
#' @param nrep the number of simulations to be performed for the rate shift
#'   test, by default \code{nrep} is set at 1000.
#' @param f the size of the smallest clade to be tested. By default, nodes
#'   subtending to one tenth of the tree tips are tested.
#' @param filename is deprecated. \code{search.shift} does not return plots
#'   anymore, check the function \code{\link{plotShift}} instead.
#' @importFrom graphics symbols mtext
#' @importFrom stats sd
#' @importFrom utils globalVariables
#' @importFrom grDevices	pdf	dev.off
#' @export
#' @seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#' @seealso \href{../doc/Plotting-tools.html}{\code{plotShift} vignette}
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
#'   the rest of the tree. Caution must be put into interpreting the
#'   'auto-recognize' results. For instance, a clade with small rates and
#'   another with large rates could be individuated even under BM. This does not
#'   mean these clades are actual instances for rate shifts. Clades must be
#'   tested on their own without the 'auto-recognize' feature, which serves as
#'   guidance to the investigator, when no strong a priori hypothesis to be
#'   tested is advanced. The 'auto-recognize' feature is not meant to provide a
#'   test for a specific hypothesis. It serves as an optional guidance to
#'   understand whether and which clades show significantly large or small rates
#'   as compared to the rest of the tree. Individual clades are tested at once,
#'   meaning that significant instances of rate variation elsewhere on the tree
#'   are ignored. Conversely, running the \code{"clade"} condition without
#'   including the 'auto-recognize' feature, multiple clades presumed to evolve
#'   under the same shift are tested together, meaning that their rates are
#'   collectively contrasted to the rest of the tree, albeit they can
#'   additionally be compared to each other upon request. Under both the
#'   \code{"clade"} and \code{"sparse"} conditions, multiple clades could be
#'   specified at once, and optionally tested individually (for deviation of
#'   rates) against the rates of the rest of the tree and against each other.
#'   Regardless of which condition is specified, the function output produces
#'   the real difference of means, and their significance value.
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
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#'
#' # Case 1. Without accounting for the effect of a covariate
#'
#' # Case 1.1 "clade" condition
#' # with auto-recognize
#' search.shift(RR=dinoRates,status.type="clade")
#' # testing two hypothetical clades
#' search.shift(RR=dinoRates,status.type="clade",node=c(696,746))
#'
#' # Case 1.2 "sparse" condition
#' # testing the sparse condition.
#' search.shift(RR=dinoRates,status.type= "sparse",state=statedino)
#'
#'
#' # Case 2. Accounting for the effect of a covariate
#'
#' # Case 2.1 "clade" condition
#' search.shift(RR=dinoRates,status.type= "clade",cov=massdino)
#'
#' # Case 2.2 "sparse" condition
#' search.shift(RR=dinoRates,status.type="sparse",state=statedino,cov=massdino)
#'     }



search.shift<-function(RR,
                       status.type=c("clade","sparse"),
                       node=NULL,
                       state=NULL,
                       cov=NULL,
                       nrep=1000,
                       f=NULL,
                       filename=NULL)
{
  # require(phytools)
  # require(scales)

  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!missing(filename))
    warning("The argument filename is deprecated. Check the function plotShift to plot results",immediate. = TRUE)

  tree <- RR$tree
  rates <- RR$rates[,,drop=FALSE]
  betas<-RR$multiple.rates[,,drop=FALSE]


  if(is.null(f)) f<-round(Ntip(tree)/10)

  if(!is.null(cov)){
    RRphylo(tree,cov,clus=0)->RRcova
    abs(c(RRcova$aces,cov))->Y
    c(rownames(RRcova$aces),names(cov))->names(Y)

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
      residuals(lm(R~Y))->res
      as.matrix(res)->betas
    }

    if(ncol(betas)>1) rates <- as.matrix(apply(betas, 1, function(x) sqrt(sum(x^2)))) else rates<-betas
    # rates <- apply(betas, 1, function(x) sqrt(sum(x^2)))
    # rates <- as.matrix(rates)
  }

  if (status.type == "clade") {
    if (is.null(node)) {
      st <- subtrees(tree)
      len<-sapply(st,Ntip)
      st <- st[which(len < (Ntip(tree)/2) & len > round(f))]
      node <- sapply(st, function(x) getMRCA(tree, x$tip.label))
      leaf2N.diff <- p.single <- array()
      for (j in 1:length(node)) {
        Cleaf <- c(getDescendants(tree, node[j]), tips(tree, node[j]))
        leaf.rates <-  na.omit(rates[match(Cleaf, rownames(rates)),])
        NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
        leaf2N.diff[j] <- mean(abs(leaf.rates))- mean(abs(NCrates))
        C <- length(leaf.rates)
        NC <- length(rates) - C
        ran.diffM <- replicate(nrep,mean(sample(abs(rates), C)) - mean(sample(abs(rates),NC)))
        p.single[j] <- rank(c(leaf2N.diff[j], ran.diffM[-nrep]))[1]/nrep
      }
      names(leaf2N.diff) <- names(p.single) <- node


      if (length(p.single[p.single>=0.975|p.single<=0.025])==0){
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),names(leaf2N.diff))]
        p.single <- leaf2N.diff <- NULL
      }

      if (length(p.single[p.single>=0.975|p.single<=0.025])==1){
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),names(leaf2N.diff))]

        p.single <- p.single[p.single>=0.975|p.single<=0.025]
        leaf2N.diff <- leaf2N.diff[match(names(p.single),names(leaf2N.diff))]

      }

      if (length(p.single[p.single>=0.975|p.single<=0.025])>= 2){
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init),names(leaf2N.diff))]

        p.single <- p.single[p.single>=0.975|p.single<=0.025]

        leaf2N.diff <- leaf2N.diff[match(names(p.single), names(leaf2N.diff))]

        ups <- p.single[p.single >= 0.975]
        dws <- p.single[p.single <= 0.025]
        ups <- ups[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                         decreasing = FALSE)]), names(ups)))]
        dws <- dws[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                         decreasing = FALSE)]), names(dws)))]
        if (is.na(mean(dws))) {
          dws = Nnode(tree) * 2
        } else {
          s = 1
          repeat{
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
          repeat{
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
        p.single <- p.single[which(names(p.single)%in%names(c(ups, dws)))]

        leaf2N.diff <- leaf2N.diff[match(names(p.single),names(leaf2N.diff))]

        p.single[order(p.single)]->p.single
      }

      data.frame(rate.difference=l2N.init[match(names(p.init),names(l2N.init))],p.value=p.init)->allres
      if(!is.null(p.single))
        data.frame(rate.difference=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],
                   p.value=p.single)->single else single<-NULL
      res<-list(allres,single)
      names(res)<-c("all.clades","single.clades")
    } else {
      Cbranch<-unlist(lapply(node,function(k) getDescendants(tree,k)))
      Cbranch <- unique(Cbranch[-which(Cbranch < Ntip(tree))])
      Ctips<-unique(unlist(lapply(node,function(k) tips(tree,k))))
      Cleaf <- c(Cbranch, Ctips)
      leaf.rates <- na.omit(rates[match(Cleaf, rownames(rates)),])
      NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
      leaf2NC.diff <- mean(abs(leaf.rates))-mean(abs(NCrates))
      C <- length(leaf.rates)
      NC <- length(rates) - C
      ran.diffR <- replicate(nrep,mean(sample(abs(rates), C)) -
                               mean(sample(abs(rates),NC)))
      p.shift <- rank(c(leaf2NC.diff, ran.diffR[-nrep]))[1]/nrep
      #names(leaf2NC.diff)<-names(p.shift)<-node
      res<-list(data.frame(rate.difference=leaf2NC.diff,p.value=p.shift))

      if(length(node)>1){
        leaf2N.diff <- p.single <- array()
        for (i in 1:length(node)) {
          NOD <- node[-i]
          others <- unlist(lapply(NOD,function(k) tips(tree, k)))
          mommies <- unlist(lapply(NOD,function(k) getDescendants(tree, k)))
          mommies <- mommies[-which(mommies< Ntip(tree)+1)]
          otmom <- c(mommies, others)
          Cleaf <- c(getDescendants(tree, node[i]),unlist(tips(tree, node[i])))
          leaf.rates <- na.omit(rates[match(Cleaf, rownames(rates)),])
          NC <- rates[-c(which(rownames(rates) %in% names(leaf.rates)),
                         which(rownames(rates) %in% otmom)), ]
          leaf2N.diff[i] <- mean(abs(leaf.rates))-mean(abs(NC))
          NC.l <- length(NC)
          leaf.l <- length(leaf.rates)
          tot.r <- abs(c(NC, leaf.rates))
          RAN.diff<-replicate(nrep,mean(sample(tot.r, leaf.l))-mean(sample(tot.r, NC.l)))
          p.single[i] <- rank(c(leaf2N.diff[i], RAN.diff[-nrep]))[1]/nrep
        }
        names(p.single) <- names(leaf2N.diff) <-node

        res<-c(res,list(data.frame(rate.difference=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],
                                 p.value=p.single)))

      }

      if(length(res)>1) names(res)<-c("all.clades.together","single.clades") else{
        rownames(res[[1]])<-node
        names(res)<-"single.clades"
      }
    }
  } else {
    state<-as.matrix(state)
    state <- treedataMatch(tree, state)[[1]][,1]
    frame <- data.frame(status = as.factor(state),
                        rate = rates[match(names(state),rownames(rates))])
    p.status.diff <- array()
    if (length(unique(state)) > 2) {
      status.diff <- apply(combn(tapply(abs(frame$rate),
                                        frame$status, mean), 2), 2, diff)
      sta <- tapply(abs(frame$rate), frame$status, mean)
      sta <- sta[match(unique(state), names(sta))]

      w <- sapply(1:length(sta),function(x) sta[x]-
                    mean(abs(frame[-which(frame$status ==names(sta)[x]), 2])))
      status.diff <- c(status.diff, w)
      names(status.diff) <- c(apply(combn(levels(frame$status),
                                          2), 2, function(x) paste(x[2], x[1], sep = "_")),
                              names(sta))
      status.diffS <- matrix(ncol = length(status.diff),
                             nrow = nrep)
      for (i in 1:nrep) {
        s.ran <- sample(frame$status)
        s.frame <- data.frame(s.ran, frame$rate)
        SD <- apply(combn(tapply(abs(s.frame$frame.rate),
                                 s.frame$s.ran, mean), 2), 2, diff)
        sta <- tapply(abs(frame$rate), s.ran, mean)
        sta <- sta[match(unique(state), names(sta))]

        w <- sapply(1:length(sta),function(x) sta[x]-
                      mean(abs(frame[-which(s.frame$s.ran ==names(sta)[x]), 2])))
        status.diffS[i, ] <- c(SD, w)
      }
      colnames(status.diffS) <- names(status.diff)

      p.status.diff<-sapply(1:length(status.diff),function(i)
        rank(c(status.diff[i],status.diffS[-nrep, i]))[1]/nrep)
      names(p.status.diff) <- names(status.diff)
      unlist(p.status.diff)->p.status.diff
      unlist(status.diff)->status.diff

      res<-list(data.frame(rate.difference=status.diff[match(names(p.status.diff),names(status.diff))],p.status.diff))
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
        s.frame <- data.frame(sample(s.state), frame$rate)
        s.frame[, 1] <- as.factor(s.frame[, 1])
        status.diffS[i] <- diff(tapply(abs(s.frame$frame.rate),
                                       s.frame[,1], mean))
      }

      p.status.diff <- rank(c(status.diff, status.diffS[-nrep]))[1]/nrep
      state.results<-data.frame(rate.difference=status.diff[match(names(p.status.diff),names(status.diff))],p.value=p.status.diff)
      rownames(state.results)<-paste(names(p.status.diff),unique(state)[which(unique(state)!=names(p.status.diff))],sep="-")
      res<-list(state.results)
      names(res)<-c("state.results")
    }
  }

  if(!is.null(cov)) {
    res<-c(res,list(rates))
    names(res)[length(res)]<-"rates"
  }

  return(res)
}
