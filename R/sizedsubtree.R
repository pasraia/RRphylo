#' @title Find a node subtending to a clade of desired size
#'
#' @description The function \code{sizedsubtree} scans a phylogentic tree to
#'   randomly find nodes subtending to a subtree of desired minimum size, up to
#'   one half of the tree size (number of tips).
#' @usage sizedsubtree(tree,Size=NULL,time.limit=10)
#' @param tree a phylogenetic tree.
#' @param Size the desired size of the tree subtending to the extracted node. By
#'   default, the minimum tree size is set at one tenth of the tree size (i.e.
#'   number of tips).
#' @param time.limit specifies a limit to the searching time, a warning message
#'   is thrown if the limit is reached.
#' @export
#' @importFrom stats na.omit
#' @details The argument \code{time.limit} sets the searching time. The
#'   algorithm stops if that limit is reached, avoiding recursive search when no
#'   solution is in fact possible.
#' @return A node subtending to a subtree of desired minimum size.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#'
#' sizedsubtree(tree=treedino,Size=40)

sizedsubtree <-function(tree,Size=NULL,time.limit=10)
{
  #require(R.utils)
  #require(ape)
  #require(geiger)

  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package \"R.utils\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  szs<-function(tree,size=Size){
    if(is.null(size)) size<-Ntip(tree)/10
    nod <- array()
    repeat {
      i = 1
      node <- sample(seq(Ntip(tree)+2, dim(tree$edge)[1]-1), 1)
      a <- length(tips(tree, node))
      i = i + 1
      if (a >= size & a <= (Ntip(tree)-1)/2)
        nod[i] <- node
      else nod[i] <- 1
      if (nod[i] > 1) {
        break
      }
      if (is.na(sum(nod))) nod <- na.omit(nod) else nod <- nod
    }
    return(nod[2])
  }
  R.utils::withTimeout({szs(tree,Size)},timeout=time.limit,onTimeout="silent")->a
  if(is.null(a)) print("searching time exceeded the time limit. It is possible there is no clade large enough to satisfy the condition") else return(a)
}

