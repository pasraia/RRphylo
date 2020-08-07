#' @title Finding distance between nodes and tips
#' @description The function computes the distance between pairs of nodes, pairs
#'   of tips, or between nodes and tips. The distance is meant as both patristic
#'   distance and the number of nodes intervening between the pair.
#' @usage distNodes(tree,node=NULL,clus=0.5)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param node either a single node/tip or a pair of nodes/tips.
#' @param clus the proportion of clusters to be used in parallel computing. To
#'   run the single-threaded version of \code{distNodes} set \code{clus} = 0.
#' @export
#' @return If \code{node} is specified, the function returns a data frame with
#'   distances between the focal node/tip and the other nodes/tips on the tree
#'   (or for the selected pair only). Otherwise, the function returns a matrix
#'   containing the number of nodes intervening between each pair of nodes and
#'   tips.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' \donttest{
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' cc<- 2/parallel::detectCores()
#' distNodes(tree=Tstage,clus=cc)
#' distNodes(tree=Tstage,node=64,clus=cc)
#' distNodes(tree=Tstage,node="Tro_2",clus=cc)
#' distNodes(tree=Tstage,node=c(64,48),clus=cc)
#' distNodes(tree=Tstage,node=c(64,"Tro_2"),clus=cc)
#' }

distNodes<-function(tree,node=NULL,clus=0.5){
  #require(ape)
  #require(phytools)
  #require(geiger)
  #require(doParallel)
  #require(parallel)

  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }

  makeL(tree)->L
  makeL1(tree)->L1
  colnames(L)->nam
  nam[(Nnode(tree)+1):(Nnode(tree)+Ntip(tree))]<-which(tree$tip.label%in%colnames(L)[(Nnode(tree)+1):(Nnode(tree)+Ntip(tree))])
  as.numeric(nam)->nam
  if(is.null(node)){
    matrix(ncol=ncol(L),nrow=ncol(L))->mat
    i=1
    res<-list()
    if(round((detectCores() * clus), 0)==0) cl<-makeCluster(1) else cl <- makeCluster(round((detectCores() * clus), 0))
    registerDoParallel(cl)
    res <- foreach(i = 1:length(nam),
                   .packages = c("RRphylo","ape", "geiger", "phytools", "doParallel")) %dopar%
                   {
                     gc()
                     nam[i]->n1
                     mata<-array()
                     for(j in i:length(nam)){
                       nam[j]->n2
                       if(n1==n2) ll<-0 else{
                         if(n1>Ntip(tree)&n2>Ntip(tree)){
                           if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
                             c(n1,n2)->nn
                             nn[which.min(nn)]->mrca
                             nn[which.max(nn)]->n
                             getDescendants(tree,mrca)->des
                             des[which(des>Ntip(tree))]->ndes
                             length(which(names(which(L1[which(rownames(L1)==n),]>0))%in%ndes))->ll
                           }else{
                             getMRCA(tree,c(n1,n2))->mrca
                             getDescendants(tree,mrca)->des
                             des[which(des>Ntip(tree))]->ndes
                             length(which(names(which(L1[which(rownames(L1)==n1),]>0))%in%ndes))->l1
                             length(which(names(which(L1[which(rownames(L1)==n2),]>0))%in%ndes))->l2
                             l1+l2->ll
                           }
                         }


                         if(n1<=Ntip(tree)&n2<=Ntip(tree)){
                           getMRCA(tree,c(tree$tip.label[n1],tree$tip.label[n2]))->mrca
                           getDescendants(tree,mrca)->des
                           des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
                           length(which(names(which(L[which(rownames(L)==tree$tip.label[n1]),]>0))%in%des))->l1
                           length(which(names(which(L[which(rownames(L)==tree$tip.label[n2]),]>0))%in%des))->l2
                           l1+l2-1->ll
                         }

                         if(n1>Ntip(tree)&n2<=Ntip(tree)|n1<=Ntip(tree)&n2>Ntip(tree)){
                           if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
                             c(n1,n2)->nn
                             nn[which.max(nn)]->mrca
                             nn[which.min(nn)]->n
                             getDescendants(tree,mrca)->des
                             des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
                             length(which(names(which(L[which(rownames(L)==tree$tip.label[n]),]>0))%in%des))->ll
                           }else{
                             getMRCA(tree,c(n1,n2))->mrca
                             getDescendants(tree,mrca)->des
                             des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
                             length(which(names(which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0))%in%des))->l1
                             length(which(names(which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0))%in%des))->l2
                             l1+l2->ll
                           }
                         }
                       }

                       ll->mata[j]
                     }
                     mata->res[[i]]
                   }
    stopCluster(cl)

    do.call(rbind,res)->mat
    t(mat)[lower.tri(t(mat))]->mat[lower.tri(mat)]
    colnames(mat)<-rownames(mat)<-colnames(L)
  }else{
    if(length(node)==1){
      node->n1
      if(n1%in%tree$tip.label)  which(tree$tip.label==n1)->n1 else as.numeric(n1)->n1
      if(node<=Ntip(tree)) tree$tip.label[n1]->node

      mat<-matrix(ncol=2,nrow=length(nam))
      for(j in 1:length(nam)){
        nam[j]->n2
        if(n1==n2){
          ll<-0
          lt<-0
        }else{
          if(n1>Ntip(tree)&n2>Ntip(tree)){
            if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
              c(n1,n2)->nn
              nn[which.min(nn)]->mrca
              nn[which.max(nn)]->n
              getDescendants(tree,mrca)->des
              des[which(des>Ntip(tree))]->ndes
              length(which(names(which(L1[which(rownames(L1)==n),]>0))%in%ndes))->ll
              sum(L1[which(rownames(L1)==n),][which(L1[which(rownames(L1)==n),]>0)][which(names(which(L1[which(rownames(L1)==n),]>0))%in%ndes)])->lt
            }else{
              getMRCA(tree,c(n1,n2))->mrca
              getDescendants(tree,mrca)->des
              des[which(des>Ntip(tree))]->ndes
              length(which(names(which(L1[which(rownames(L1)==n1),]>0))%in%ndes))->l1
              length(which(names(which(L1[which(rownames(L1)==n2),]>0))%in%ndes))->l2
              sum(L1[which(rownames(L1)==n1),][which(L1[which(rownames(L1)==n1),]>0)][which(names(which(L1[which(rownames(L1)==n1),]>0))%in%ndes)])->lt1
              sum(L1[which(rownames(L1)==n2),][which(L1[which(rownames(L1)==n2),]>0)][which(names(which(L1[which(rownames(L1)==n2),]>0))%in%ndes)])->lt2
              l1+l2->ll
              lt1+lt2->lt
            }
          }


          if(n1<=Ntip(tree)&n2<=Ntip(tree)){
            getMRCA(tree,c(tree$tip.label[n1],tree$tip.label[n2]))->mrca
            getDescendants(tree,mrca)->des
            des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
            length(which(names(which(L[which(rownames(L)==tree$tip.label[n1]),]>0))%in%des))->l1
            length(which(names(which(L[which(rownames(L)==tree$tip.label[n2]),]>0))%in%des))->l2
            sum(L[which(rownames(L)==tree$tip.label[n1]),][which(L[which(rownames(L)==tree$tip.label[n1]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n1]),]>0))%in%des)])->lt1
            sum(L[which(rownames(L)==tree$tip.label[n2]),][which(L[which(rownames(L)==tree$tip.label[n2]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n2]),]>0))%in%des)])->lt2
            l1+l2-1->ll
            lt1+lt2->lt
          }

          if(n1>Ntip(tree)&n2<=Ntip(tree)|n1<=Ntip(tree)&n2>Ntip(tree)){
            if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
              c(n1,n2)->nn
              nn[which.max(nn)]->mrca
              nn[which.min(nn)]->n
              getDescendants(tree,mrca)->des
              des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
              length(which(names(which(L[which(rownames(L)==tree$tip.label[n]),]>0))%in%des))->ll
              sum(L[which(rownames(L)==tree$tip.label[n]),][which(L[which(rownames(L)==tree$tip.label[n]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n]),]>0))%in%des)])->lt
            }else{
              getMRCA(tree,c(n1,n2))->mrca
              getDescendants(tree,mrca)->des
              des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
              length(which(names(which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0))%in%des))->l1
              length(which(names(which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0))%in%des))->l2
              sum(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),][which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0))%in%des)])->lt1
              sum(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),][which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0)][which(names(which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0))%in%des)])->lt2
              l1+l2->ll
              lt1+lt2->lt
            }
          }
        }
        ll->mat[j,1]
        lt->mat[j,2]
      }
      rownames(mat)<-colnames(L)
      colnames(mat)<-c("node","time")
    }else{
      node[1]->n1
      node[2]->n2
      if(n1%in%tree$tip.label)  which(tree$tip.label==n1)->n1 else as.numeric(n1)->n1
      if(n2%in%tree$tip.label)  which(tree$tip.label==n2)->n2 else as.numeric(n2)->n2
      if(node[1]<=Ntip(tree)) tree$tip.label[n1]->node[1]
      if(node[2]<=Ntip(tree)) tree$tip.label[n2]->node[2]


      if(n1==n2) {
        ll<-0
        lt<-0
      }else{
        if(n1>Ntip(tree)&n2>Ntip(tree)){
          if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
            c(n1,n2)->nn
            nn[which.min(nn)]->mrca
            nn[which.max(nn)]->n
            getDescendants(tree,mrca)->des
            des[which(des>Ntip(tree))]->ndes
            length(which(names(which(L1[which(rownames(L1)==n),]>0))%in%ndes))->ll
            sum(L1[which(rownames(L1)==n),][which(L1[which(rownames(L1)==n),]>0)][which(names(which(L1[which(rownames(L1)==n),]>0))%in%ndes)])->lt
          }else{
            getMRCA(tree,c(n1,n2))->mrca
            getDescendants(tree,mrca)->des
            des[which(des>Ntip(tree))]->ndes
            length(which(names(which(L1[which(rownames(L1)==n1),]>0))%in%ndes))->l1
            length(which(names(which(L1[which(rownames(L1)==n2),]>0))%in%ndes))->l2
            sum(L1[which(rownames(L1)==n1),][which(L1[which(rownames(L1)==n1),]>0)][which(names(which(L1[which(rownames(L1)==n1),]>0))%in%ndes)])->lt1
            sum(L1[which(rownames(L1)==n2),][which(L1[which(rownames(L1)==n2),]>0)][which(names(which(L1[which(rownames(L1)==n2),]>0))%in%ndes)])->lt2
            l1+l2->ll
            lt1+lt2->lt
          }
        }


        if(n1<=Ntip(tree)&n2<=Ntip(tree)){
          getMRCA(tree,c(tree$tip.label[n1],tree$tip.label[n2]))->mrca
          getDescendants(tree,mrca)->des
          des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
          length(which(names(which(L[which(rownames(L)==tree$tip.label[n1]),]>0))%in%des))->l1
          length(which(names(which(L[which(rownames(L)==tree$tip.label[n2]),]>0))%in%des))->l2
          sum(L[which(rownames(L)==tree$tip.label[n1]),][which(L[which(rownames(L)==tree$tip.label[n1]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n1]),]>0))%in%des)])->lt1
          sum(L[which(rownames(L)==tree$tip.label[n2]),][which(L[which(rownames(L)==tree$tip.label[n2]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n2]),]>0))%in%des)])->lt2
          l1+l2-1->ll
          lt1+lt2->lt
        }

        if(n1>Ntip(tree)&n2<=Ntip(tree)|n1<=Ntip(tree)&n2>Ntip(tree)){
          if(n2%in%c(getMommy(tree,n1),getDescendants(tree,n1))){
            c(n1,n2)->nn
            nn[which.max(nn)]->mrca
            nn[which.min(nn)]->n
            getDescendants(tree,mrca)->des
            des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
            length(which(names(which(L[which(rownames(L)==tree$tip.label[n]),]>0))%in%des))->ll
            sum(L[which(rownames(L)==tree$tip.label[n]),][which(L[which(rownames(L)==tree$tip.label[n]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[n]),]>0))%in%des)])->lt
          }else{
            getMRCA(tree,c(n1,n2))->mrca
            getDescendants(tree,mrca)->des
            des[which(des<=Ntip(tree))]<-tree$tip.label[des[which(des<=Ntip(tree))]]
            length(which(names(which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0))%in%des))->l1
            length(which(names(which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0))%in%des))->l2
            sum(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),][which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0)][which(names(which(L[which(rownames(L)==tree$tip.label[c(n1,n2)[which.min(c(n1,n2))]]),]>0))%in%des)])->lt1
            sum(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),][which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0)][which(names(which(L1[which(rownames(L1)==c(n1,n2)[which.max(c(n1,n2))]),]>0))%in%des)])->lt2
            l1+l2->ll
            lt1+lt2->lt
          }
        }
      }
      data.frame(node=ll,time=lt)->mat
    }
  }
  return(mat)
}


