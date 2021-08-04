#' @title Fast addition of tips and clades on an existing tree
#' @description The function attaches new tips and/or clades derived from a
#'   source phylogeny to a pre-existing backbone tree.
#' @usage tree.merger(backbone,data,source.tree=NULL,age.offset=NULL,tip.ages =
#'   NULL, node.ages = NULL,min.branch=NULL,plot=TRUE,filename=NULL)
#' @param backbone the backbone tree to attach tips/clades on.
#' @param data a dataset including as columns:\enumerate{\item bind = the
#'   tips/clades to be attached; \item reference = the reference tip or clade
#'   where 'bind' must be attached; \item poly = logical specifying if 'bind'
#'   and 'reference' should form a polytomous clade.} See details for further
#'   explanations.
#' @param source.tree the tree where 'bind' clades are to be extracted from. If
#'   no clade has to be attached, it can be left unspecified.
#' @param age.offset if the most recent age (i.e. the maximum distance from the
#'   tree root) differs between the source and the backbone trees, the
#'   “age.offset” is the difference between them in this exact order (source
#'   minus backbone). It is positive when the backbone tree attains younger age
#'   than the source tree, and vice-versa.
#' @param tip.ages as in \code{\link{scaleTree}}, a named vector including the
#'   ages (i.e. the time distance from the youngest tip within the tree) of the
#'   tips. If unspecified, the function assumes all the tips on the backbone
#'   tree are correctly placed and places all the new tips at the maximum
#'   distance from the tree root (i.e. the present if the tips are extant).
#' @param node.ages as in \code{scaleTree}, a named vector including the ages
#'   (i.e. the time distance from the youngest tip within the tree) of the
#'   nodes. The nodes must be defined by collating the names of the two
#'   phylogenetically furthest tips it subtends to, separated by a "-" (see
#'   examples). If no calibration date for nodes is supplied, the function may
#'   shift the node position back in time as to place new tips/clades and to fit
#'   tip ages.
#' @param min.branch as in \code{scaleTree}, the minimum branch length that will
#'   be imposed for shifted nodes.
#' @param plot if \code{TRUE}, the function produces an interactive plotting
#'   device to check the placing of each \code{bind}.
#' @param filename if \code{plot=TRUE} and provided a \code{filename} (with or
#'   without the path), the function stores a pdf file showing the plot of the
#'   entire phylogeny.
#' @importFrom ape bind.tree
#' @return Merged phylogenetic tree.
#' @author Silvia Castiglione, Carmela Serio, Pasquale Raia
#' @details The function attaches tips and/or clades from the \code{source} tree
#'   to the \code{backbone} tree according to the \code{data} object. Within the
#'   latter, a clade, either to be bound or to be the reference, must be
#'   indicated by collating the names of the two phylogenetically furthest tips
#'   belonging to it, separated by a "-". Duplicated 'bind' produce error.
#'   Tips/clades set to be attached to the same 'reference' are considered to
#'   represent a polytomy. Tips set as 'bind' which are already on the backbone
#'   tree are removed from the latter and placed according to the 'reference'.
#'   See examples and
#'   \href{../doc/Tree-Manipulation.html#tree.merger.html}{vignette} for
#'   clarifications.
#' @export
#' @seealso
#' \href{../doc/Tree-Manipulation.html#tree.merger.html}{\code{tree.merger}
#' vignette}; \href{../doc/Tree-Manipulation.html#scaleTree}{\code{scaleTree}
#' vignette};
#' @references aaa
#' @examples
#'  \dontrun{
#'  require(ape)
#'  require(geiger)
#'  DataCetaceans$treecet->tree
#'  data.frame(bind=c("Balaena_mysticetus-Caperea_marginata",
#'                    "Aetiocetus_weltoni",
#'                    "Saghacetus_osiris",
#'                    "Zygorhiza_kochii",
#'                    "Ambulocetus_natans",
#'                    "Kentriodon_pernix",
#'                    "Kentriodon_schneideri",
#'                    "Kentriodon_obscurus",
#'                    "Tursiops_truncatus-Delphinus_delphis",
#'                    "Kogia_sima",
#'                    "Grampus_griseus"),
#'             reference=c("Fucaia_buelli-Aetiocetus_weltoni",
#'                         "Aetiocetus_cotylalveus",
#'                         "Fucaia_buelli-Tursiops_truncatus",
#'                         "Saghacetus_osiris-Fucaia_buelli",
#'                         "Dalanistes_ahmedi-Fucaia_buelli",
#'                         "Kentriodon_schneideri",
#'                         "Phocoena_phocoena-Delphinus_delphis",
#'                         "Kentriodon_schneideri",
#'                         "Stenella_attenuata-Stenella_longirostris",
#'                         "Kogia_breviceps",
#'                         "Globicephala_melas-Pseudorca_crassidens"),
#'             poly=c(FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE,
#'                    FALSE))->dato
#'
#'  c(Aetiocetus_weltoni=28.0,
#'    Saghacetus_osiris=33.9,
#'    Zygorhiza_kochii=34.0,
#'    Ambulocetus_natans=40.4,
#'    Kentriodon_pernix=15.9,
#'    Kentriodon_schneideri=11.61,
#'    Kentriodon_obscurus=13.65)->tip.ages
#'  c("Ambulocetus_natans-Fucaia_buelli"=52.6,
#'    "Balaena_mysticetus-Caperea_marginata"=21.5)->node.ages
#'
#'  # remove some tips from the original tree and create a source tree
#'  drop.tip(tree,c(names(tip.ages),
#'                  tips(tree,131)[-which(tips(tree,131)%in%
#'                                c("Caperea_marginata","Eubalaena_australis"))],
#'                  tips(tree,195)[-which(tips(tree,195)=="Tursiops_aduncus")]))->backtree
#'  drop.tip(tree,which(!tree$tip.label%in%c(names(tip.ages),
#'                                           tips(tree,131),
#'                                           tips(tree,195))))->sourcetree
#'
#'  plot(backtree,cex=.6)
#'  plot(sourcetree,cex=.6)
#'
#'  tree.merger(backbone=backtree,data=dato,source.tree=sourcetree,
#'              tip.ages=tip.ages,node.ages = node.ages, plot=TRUE)->treeM
#'    }

tree.merger<-function(backbone,data,source.tree=NULL,age.offset=NULL,
                      tip.ages = NULL, node.ages = NULL,min.branch=NULL,plot=TRUE,filename=NULL){
  # require(ape)
  # require(phytools)
  # require(geiger)
  # require(manipulate)

  if (!requireNamespace("manipulate", quietly = TRUE)) {
    stop("Package \"manipulate\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  backbone->tree
  data->dat
  source.tree->tree2
  max(diag(vcv(tree)))->H
  H-diag(vcv(tree))->ages

  if(!is.null(age.offset)&&age.offset<0){
    ages+abs(age.offset)->ages
    H+abs(age.offset)->Hset
  }else H->Hset

  if(is.null(min.branch)) min(tree$edge.length)->min.branch

  ### Check data ###
  if(!all(colnames(dat)%in%c("bind","reference","poly"))) {
    if(any(is.na(as.logical(dat[,3])))) stop("Check columns order: it should be 'bind', 'reference', 'poly'")
    warning("Colnames not matching: columns assumed to be ordered as 'bind','reference','poly'",immediate. = TRUE)
    colnames(dat)<-c("bind","reference","poly")
  }

  ### Check for duplicated bind ###
  if(any(duplicated(dat$bind))) stop(paste(paste(dat$bind[duplicated(dat$bind)],collapse = ", "),"names duplicated in supplied tips"))

  if(!is.logical(dat$poly)) as.logical(dat$poly)->dat$poly
  data.frame(dat,bind.type=sapply(strsplit(dat[,1],"-"),length))->dat

  if(all(dat$bind.type==1)&(!is.null(tree2))) tree2<-NULL

  if(any(dat$bind%in%tree2$tip.label)) drop.tip(tree2,dat[which(dat$bind%in%tree2$tip.label),1])->tree2

  # if(all(dat$bind.type==1))
  #   data.frame(bind=dat$bind,cbind(unlist(strsplit(dat[,1],"-")),unlist(strsplit(dat[,1],"-"))),dat[,2:4])->dat else
  #     data.frame(bind=dat$bind,do.call(rbind,strsplit(dat[,1],"-")),dat[,2:4])->dat
  # colnames(dat)[2:3]<-paste("bind",1:2,sep="")

  dat$bind.tips<-NA
  if(all(dat$bind.type==1)) dat$bind.tips<-dat$bind else{
    dat[which(dat$bind.type==2),]$bind.tips<-lapply(dat$bind[which(dat$bind.type==2)],function(x) tips(tree2,getMRCA(tree2,strsplit(x,"-")[[1]])))
    if(any(is.na(dat$bind.tips))) dat[which(is.na(dat$bind.tips)),]$bind.tips<-dat[which(is.na(dat$bind.tips)),]$bind
  }

  ### bind missing from tree2 ###
  # if(!all(c(dat[which(dat$bind.type==2),2],dat[which(dat$bind.type==2),3])%in%tree2$tip.label)){
  #   stop(paste(paste(c(dat[which(dat$bind.type==2),2],
  #                      dat[which(dat$bind.type==2),3])[which(!c(dat[which(dat$bind.type==2),2],
  #                                                               dat[which(dat$bind.type==2),3])%in%tree2$tip.label)],
  #                    collapse=", "),
  #              "not in source.tree"))
  # }
  #

  if(!all(unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)){
    stop(paste(paste(unlist(dat[which(dat$bind.type==2),]$bind.tips)[which(
      unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)],
      collapse=", "),"not in source.tree"))
  }

  ### reference missing ###
  # if(!all(unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,dat$bind1,dat$bind2,tree2$tip.label))){
  #   stop(paste(paste(unlist(strsplit(dat$reference,"-"))[which(!unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,dat$bind1,dat$bind2))],
  #                    collapse=","),"missing from the backbone, the source and the tips to be attached"))
  # }

  if(!all(unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))){
    stop(paste(paste(unlist(strsplit(dat$reference,"-"))[which(!unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))],
                     collapse=","),"missing from the backbone, the source and the tips to be attached"))
  }

  ### bind already on the backbone ###
  if(any(dat$bind%in%tree$tip.label)){
    warning(paste(paste(dat[which(dat$bind%in%tree$tip.label),1],collapse=", "),"removed from the backbone tree"),immediate. = TRUE)
    drop.tip(tree,dat[which(dat$bind%in%tree$tip.label),1])->tree
  }

  ### duplicated reference ###
  table(dat$reference)->tab.ref
  if(any(tab.ref>1)){
    for(j in 1:length(which(tab.ref>1))){
      dat[which(dat$reference==names(which(tab.ref>1)[j])),]->ref.mult
      if(any(isTRUE(ref.mult$poly))) ref.mult[-which(isTRUE(ref.mult$poly)),]->ref.mult
      paste(strsplit(ref.mult$reference[1],"-")[[1]][1],strsplit(ref.mult$bind[1],"-")[[1]][1],sep="-")->ref.mult$reference[-1]
      ref.mult[-1,]$poly<-TRUE
      dat[match(ref.mult$bind,dat$bind),]<-ref.mult
    }
  }

  ### ordering ###
  strsplit(dat$reference,"-")->refs
  dat$ref.tree1<-NA
  dat$ref.tree2<-NA

  lapply(refs,function(x){
    if(length(x[which(!x%in%tree$tip.label)])<1) NA else x[which(!x%in%tree$tip.label)]
  })->ref.tree

  dat[which(is.na(ref.tree)),][order(dat[which(is.na(ref.tree)),]$bind.type),]$ref.tree1<-seq(1,length(which(is.na(ref.tree))))

  if(any(which(!is.na(ref.tree)))){
    dat[which(!is.na(ref.tree)),]->dat.new
    # dat.new[,7:8]<-do.call(rbind,ref.tree[-which(is.na(ref.tree))])
    #
    # if(any(dat.new[,7]%in%c(dat.new[,2],dat.new[,3])|dat.new[,8]%in%c(dat.new[,2],dat.new[,3]))){
    #   while(nrow(dat.new>0)){
    #     which(!(dat.new[,7]%in%c(dat.new[,2],dat.new[,3])|dat.new[,8]%in%c(dat.new[,2],dat.new[,3])))->outs
    #     dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-max(dat$ref.tree1,na.rm=TRUE)+1:length(outs)
    #     dat.new[-outs,]->dat.new
    #   }
    # }else{
    #   which(!(dat.new[,7]%in%c(dat.new[,2],dat.new[,3])|dat.new[,8]%in%c(dat.new[,2],dat.new[,3])))->outs
    #   dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-
    #     max(dat$ref.tree1,na.rm=TRUE)+1:length(which(!dat.new$ref.tree1%in%dat.new[,1]))
    # }


    dat.new[,6:7]<-do.call(rbind,ref.tree[-which(is.na(ref.tree))])

    if(any(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips))){
      while(nrow(dat.new)>0){
        which(!(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips)))->outs
        dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-max(dat$ref.tree1,na.rm=TRUE)+1:length(outs)
        dat.new[-outs,]->dat.new
      }
    }else{
      which(!(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips)))->outs
      dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-
        max(dat$ref.tree1,na.rm=TRUE)+1:length(which(!dat.new$ref.tree1%in%dat.new[,1]))
    }
  }

  dat[order(dat$ref.tree1),]->dat

  if(!is.null(tree2)){
    # if(any(dat$bind%in%tree2$tip.label)) drop.tip(tree2,dat[which(dat$bind%in%tree2$tip.label),1])->tree2
    dat$MRCAbind<-NA
    # apply(dat[which(dat$bind.type==2),2:3],1,function(x) getMRCA(tree2,x))->dat$MRCAbind[which(dat$bind.type==2)]
    sapply(dat[which(dat$bind.type==2),]$bind.tips,function(x) getMRCA(tree2,x))->dat$MRCAbind[which(dat$bind.type==2)]

    if(any(dat$bind.type==2)){
      unlist(sapply(dat[which(dat$bind.type==2),]$MRCAbind,function(x) tree$tip.label[which(tree$tip.label%in%tips(tree2,x))]))->remt
      if(length(remt)>0){
        drop.tip(tree,remt)->tree
        ages[-match(remt,names(ages))]->ages
        warning(paste(paste(remt,collapse=", "),"already on the source tree: removed from the backbone tree"),immediate. = TRUE)
      }
    }
  }

  ### binding ###
  for(k in 1:nrow(dat)){
    if(length(strsplit(dat$reference,"-")[[k]])>1)
      getMRCA(tree,strsplit(dat$reference,"-")[[k]])->where.ref else
        which(tree$tip.label==strsplit(dat$reference,"-")[[k]])->where.ref

    if(where.ref!=(Ntip(tree)+1)) tree$edge.length[which(tree$edge[,2]==where.ref)]->br.len else{
      if(is.null(tree$root.edge)) tree$root.edge<-mean(tree$edge.length)
      tree$root.edge->br.len
    }

    if(dat$bind.type[k]==1){
      if(isTRUE(dat$poly[k])) 0->pos.ref else br.len/2->pos.ref
      bind.tip(tree,dat$bind[k],where.ref,position = pos.ref,edge.length =br.len/2)->tree
    }else{
      extract.clade(tree2,dat$MRCAbind[k])->cla

      if(isTRUE(dat$poly[k])){
        0->pos.ref
        if(where.ref==(Ntip(tree)+1)&&(max(diag(vcv(cla)))+max(diag(vcv(cla)))/10)>H)
          rescale(cla,"depth",(H+max(diag(vcv(cla)))/10))->cla else
            if(where.ref!=(Ntip(tree)+1)&(max(diag(vcv(cla)))+max(diag(vcv(cla)))/10)>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))
              rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]+max(diag(vcv(cla)))/10))->cla
        cla$root.edge<-max(diag(vcv(cla)))/10
      }else {
        if(where.ref==(Ntip(tree)+1)) br.len/2->pos.ref else {
          max(diag(vcv(cla)))+br.len/2->Hcla
          if((H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2])<H/1000){
            rescale(cla,"depth",br.len/4)->cla
            pos.ref<-br.len-br.len/4
          }else{
            if(Hcla>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),1])){
              rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))->cla
              pos.ref<-br.len/2
            }else if(Hcla>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))
              (max(diag(vcv(cla)))+br.len/2)-(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2])->pos.ref else br.len/2->pos.ref
          }
        }
        cla$root.edge<-br.len/2
      }
      bind.tree(tree,cla,where=where.ref,position = pos.ref)->tree
    }
  }

  ### tips calibration ages ###
  if(is.null(tip.ages)){
    rep(0,length(tree$tip.label))->tip.ages
    names(tip.ages)<-tree$tip.label
    tip.ages[match(names(ages),names(tip.ages))]<-ages
  }else{
    if(!all(names(ages)%in%names(tip.ages))) c(ages[which(!names(ages)%in%names(tip.ages))],tip.ages)->tip.ages
    if(!all(tree$tip.label%in%names(tip.ages))){
      rep(0,length(which(!tree$tip.label%in%names(tip.ages))))->tip.add
      names(tip.add)<-tree$tip.label[which(!tree$tip.label%in%names(tip.ages))]
      c(tip.ages,tip.add)->tip.ages
    }
  }

  ### nodes calibration ages ###
  if(!is.null(node.ages)){
    sapply(names(node.ages),function(x){
      getMRCA(tree,unlist(strsplit(x,"-")))
    })->names(node.ages)
  }else{
    node.ages<-c()
  }

  ### age original root ###
  if(!getMRCA(tree,names(ages))%in%names(node.ages)){
    node.ages<-c(node.ages,Hset)
    names(node.ages)[length(node.ages)]<-getMRCA(tree,names(ages))
  }

  ### time distances inside attached clades ###
  if(any(dat$bind.type==2)){
    unlist(lapply(dat$MRCAbind[which(dat$bind.type==2)],function(x){
      c(x,getDescendants(tree2,x)[which(getDescendants(tree2,x)>Ntip(tree2))])->des
      dist.nodes(tree2)->dn
      max(diag(vcv(tree2)))-dn[which(rownames(dn)==(Ntip(tree2)+1)),match(des,rownames(dn))]->dndes
      names(dndes)<-des
      dndes
    }))->ages.fix

    if(!is.null(age.offset)&&age.offset>0) ages.fix+age.offset->ages.fix
    sapply(names(ages.fix),function(x) getMRCA(tree,tips(tree2,as.numeric(x))))->names(ages.fix)

    if(any(!names(ages.fix)%in%names(node.ages)))
      c(ages.fix[which(!names(ages.fix)%in%names(node.ages))],node.ages)->node.ages
  }

  if(max(diag(vcv(tree)))>H&&(!(Ntip(tree)+1)%in%names(node.ages)))
    warning(paste("Root age not indicated: the tree root arbitrarily set at",round(max(diag(vcv(tree))),2)),immediate.=TRUE)

  scaleTree(tree,node.ages=node.ages,tip.ages =tip.ages,min.branch=min.branch)->tree.final->tree.plot

  if(isTRUE(plot)){
    if(any(dat$bind.type==2)) lapply(dat$MRCAbind[which(dat$bind.type==2)],function(x)
      c(getMRCA(tree.plot,tips(tree2,x)),getDescendants(tree.plot,getMRCA(tree.plot,tips(tree2,x)))))->cla.plot else cla.plot<-c()

      c(which(tree.plot$tip.label%in%dat$bind[which(dat$bind.type==1)]),unlist(cla.plot))->all.plot
      colo<-rep(scales::hue_pal()(2)[2],nrow(tree.plot$edge))
      colo[match(all.plot,tree.plot$edge[,2])]<-scales::hue_pal()(2)[1]
      names(colo)<-tree.plot$edge[,2]
      colo[which(as.numeric(names(colo))<=Ntip(tree.plot))]->colo.tips
      tree.plot$tip.label[which(!tree.plot$tip.label%in%unique(c(unlist(strsplit(dat$bind,"-")),unlist(strsplit(dat$reference,"-")))))]<-" "
      names(colo.tips)<-tree.plot$tip.label[as.numeric(names(colo.tips))]

      if(!is.null(filename)){
        pdf(file=paste(filename,".pdf",sep=""))
        if(Ntip(tree.plot)<100)
          plot(tree.plot,edge.color=colo,cex=.6,tip.color=colo.tips[match(tree.plot$tip.label,names(colo.tips))]) else
            plot(tree.plot,edge.color=colo,type="fan",cex=.6,tip.color=colo[which(tree.plot$edge[,2]<=Ntip(tree.plot))])
        dev.off()
      }

      lapply(1:nrow(dat),function(x){
        # getMRCA(tree.final,c(as.matrix(dat[x,2:3]),unlist(strsplit(dat[x,4],"-"))))->MRCAplot
        getMRCA(tree.final,c(unlist(dat$bind.tips[x]),unlist(strsplit(dat$reference[x],"-"))))->MRCAplot
        extract.clade(tree.final,getMRCA(tree.final,c(unlist(dat$bind.tips[x]),unlist(strsplit(dat$reference[x],"-")))))->cla
        colo[which(names(colo)%in%getDescendants(tree.final,MRCAplot))]->colo.cla
        as.numeric(names(colo.cla))->nam.colo

        which(cla$tip.label%in%tree.final$tip.label[nam.colo[which(nam.colo<=Ntip(tree.final))]])->nam.colo[which(nam.colo<=Ntip(tree.final))]
        if(any(nam.colo>Ntip(tree.final))) sapply(nam.colo[which(nam.colo>Ntip(tree.final))],function(r) getMRCA(cla,tips(tree.final,r)))->nam.colo[which(nam.colo>Ntip(tree.final))]
        names(colo.cla)<-nam.colo

        colo.cla<-colo.cla[match(cla$edge[,2],names(colo.cla))]
        list(cla,colo.cla)
      })->clades.plot

      names(clades.plot)<-dat[,1]


      taxon<-manipulate::picker(as.list(dat[,1]))
      manipulate::manipulate(plot(clades.plot[[which(names(clades.plot)==taxon)]][[1]],edge.color=clades.plot[[which(names(clades.plot)==taxon)]][[2]],cex=.6),
                             taxon=taxon)

  }

  return(tree.final)
}

