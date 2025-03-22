# search.conv
unitV <- function(x)  sum(x^2)^0.5
deg2rad <- function(deg) (deg * pi)/(180)
rad2deg <- function(rad)  (rad * 180)/(pi)
theta<-function(a,b) rad2deg(acos((a%*%b)/(unitV(a) *unitV(b))))
angle.vecs<-function(vec1,vec2){
  ((vec1%*%vec2)/(RRphylo:::unitV(vec1) *RRphylo:::unitV(vec2)))->ppA
  if((ppA-1)>0) ppA<-1 else if((1+ppA)<0) ppA<-(-1)
  RRphylo:::rad2deg(acos(ppA))
}

veclen <- function(v) sqrt(sum(v^2))
xprod <- function(v, w) c( v[2]*w[3] - v[3]*w[2],
                           v[3]*w[1] - v[1]*w[3],
                           v[1]*w[2] - v[2]*w[1] )

arc3d <- function(from, to, center, radius, n, circle = 50, base = 0, plot = TRUE, ...) {
  fixarg <- function(arg) {
    if (is.matrix(arg))
      arg[, 1:3, drop = FALSE]
    else
      matrix(arg, 1, 3)
  }
  normalize <- function(v)
    v / veclen(v)
  getrow <- function(arg, i) {
    arg[1 + (i - 1) %% nrow(arg),]
  }
  from <- fixarg(from)
  to <- fixarg(to)
  center <- fixarg(center)

  m <- max(nrow(from), nrow(to), nrow(center), length(base))
  base <- rep_len(base, m)

  result <- matrix(NA_real_, nrow = 1, ncol = 3)

  for (j in seq_len(m)) {
    from1 <- getrow(from, j)
    to1 <- getrow(to, j)
    center1 <- getrow(center, j)
    base1 <- base[j]
    logr1 <- log(veclen(from1 - center1))
    logr2 <- log(veclen(to1 - center1))
    A <- normalize(from1 - center1)
    B <- normalize(to1 - center1)
    steps <- if (base1 <= 0) 4*abs(base1) + 1 else 4*base1 - 1
    for (k in seq_len(steps)) {
      if (k %% 2) {
        A1 <- A * (-1)^(k %/% 2)
        B1 <- B * (-1)^(k %/% 2 + (base1 > 0))
      } else {
        A1 <- B * (-1)^(k %/% 2 + (base1 <= 0))
        B1 <- A * (-1)^(k %/% 2)
      }
      theta <- acos(sum(A1*B1))
      if (isTRUE(all.equal(theta, pi)))
        warning("Arc ", j, " points are opposite each other!  Arc is not well defined.")
      if (missing(n))
        n1 <- ceiling(circle*theta/(2*pi))
      else
        n1 <- n

      if (missing(radius)) {
        pretheta <- (k %/% 2)*pi - (k %% 2 == 0)*theta
        if (k == 1)
          totaltheta <- (steps %/% 2)*pi - (steps %% 2 == 0)*theta + theta
        p1 <- pretheta/totaltheta
        p2 <- (pretheta + theta)/totaltheta
        radius1 <- exp(seq(from = (1 - p1)*logr1 + p1*logr2,
                           to   = (1 - p2)*logr1 + p2*logr2,
                           length.out = n1 + 1))
      } else
        radius1 <- rep_len(radius, n1)
      arc <- matrix(NA_real_, nrow = n1 + 1, ncol = 3)
      p <- seq(0, 1, length.out = n1 + 1)
      arc[1,] <- center1 + radius1[1]*A1
      arc[n1 + 1,] <- center1 + radius1[n1 + 1]*B1
      AB <- veclen(A1 - B1)
      for (i in seq_len(n1)[-1]) {
        ptheta <- p[i]*theta
        phi <- pi/2 + (0.5 - p[i])*theta
        q <- (sin(ptheta) / sin(phi))/AB
        D <- (1-q)*A1 + q*B1
        arc[i,] <- center1 + radius1[i] * normalize(D)
      }
      if (k == 1)
        result <- rbind(result, arc)
      else
        result <- rbind(result[-nrow(result), ,drop = FALSE], arc)
    }
    result <- rbind(result, result[1,])
  }
  if (plot)
    lines3d(result[c(-1, -nrow(result)), , drop = FALSE], ...)
  else
    result[c(-1, -nrow(result)), , drop = FALSE]
}

conv3dplot<-function(){
  sapply(1:2,function(w) apply(y[tips(tree1,nod.par[w]),],2,mean)->>y.vec[w,])
  matrix(c(0,0,0),1,3)->rV

  colt<-c("firebrick1","deepskyblue1")
  coltx<-c("red4","blue3")
  rep("gray80",length(y))->coltips
  sapply(1:2,function(w) coltips[match(tips(tree1,nod.par[w]),rownames(y))]<<-colt[w])
  plot3d(y,col=coltips,size=1,bbox=FALSE,axes=FALSE,box=FALSE,type="s")
  box3d()
  title3d(xlab="y[,1]",ylab="y[,2]",zlab="y[,3]")
  plot3d(rV,col="green",add=TRUE,size=2,bbox=FALSE,axes=FALSE,box=FALSE,type="s")
  plot3d(ant,col="gray16",size=1,bbox=FALSE,axes=TRUE,box=FALSE,add=TRUE,type="s")
  lapply(1:2,function(w){
    text3d(ant[match(nod.par[w],rownames(ant)),,drop=FALSE],texts="mrca1",add=TRUE,adj=1,cex=1,col=coltx[w])
    segments3d(rbind(rV,ant[match(nod.par[w],rownames(ant)),,drop=FALSE]),col=coltx[w],lwd=4)
    segments3d(rbind(rV,y.vec[w,]),col=colt[w],lwd=4)
  })


  ant[match(nod.par[1],rownames(ant)),]->a
  ant[match(nod.par[2],rownames(ant)),]->b
  as.numeric(round(angle.vecs(a,b),1))->theta.ace
  veclen(a)->la
  veclen(b)->lb
  if(la>lb) (la/lb)*b->b else (lb/la)*a->a
  if(la<0.5*lb)
    segments3d(rbind(ant[match(nod.par[1],rownames(ant)),,drop=FALSE],(lb/la)*ant[match(nod.par[1],rownames(ant)),,drop=FALSE]))
  if(la>0.5*lb)
    segments3d(rbind(ant[match(nod.par[2],rownames(ant)),,drop=FALSE],(la/lb)*ant[match(nod.par[2],rownames(ant)),,drop=FALSE]))


  veclen(a)->la
  arc3d(a,b,rV,add=TRUE,col="red",radius=0.5*la)
  arc3d(b,a,rV,add=TRUE,col="red",radius=0.5*la)

  text3d(rV,texts=paste("theta.ace",theta.ace),adj=c(-0.5,0.5),col="red")

  y.vec[1,]->at1
  y.vec[2,]->at2
  as.numeric(round(angle.vecs(at1,at2),1))->theta.real
  veclen(at1)->lat1
  veclen(at2)->lat2
  if(lat1>lat2) (lat1/lat2)*at2->at2 else (lat2/lat1)*at1->at1
  if(lat1<0.5*lat2) segments3d(rbind(y.vec[1,,drop=FALSE],(lat2/lat1)*y.vec[1,,drop=FALSE]))
  if(lat2<0.5*lat1) segments3d(rbind(y.vec[2,,drop=FALSE],(lat1/lat2)*y.vec[2,,drop=FALSE]))
  veclen(at1)->lat1
  arc3d(at1,at2,rV,add=TRUE,col="blue",radius=0.5*lat1)
  arc3d(at2,at1,rV,add=TRUE,col="blue",radius=0.5*lat1)
  text3d(rV,texts=paste("theta.tips",theta.real),adj=c(-0.5,lat1),col="blue")

  range(y[,3])->ran
  ifelse(ran[1]<0,ran[1]*0.9,ran[1]*1.1)->ran[1]
  ifelse(ran[2]<0,ran[2]*1.1,ran[2]*0.9)->ran[2]


  points3d(matrix(c(rep(range(y[,1])[2]*1.5,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=FALSE),
           col=c("gray80","gray16","green","firebrick1","deepskyblue1","blue3","red4"),size=8)
  text3d(matrix(c(rep(range(y[,1])[2]*1.6,5),rep(range(y[,3])[1]*1.2,5),seq(ran[1],ran[2],length.out=5)),5,3,byrow=FALSE),
         texts=c("phenotypes at tips","phenotypes at nodes","origin","clade 2","clade 1"),adj=c(0,0.5),cex=1.2)
}

# tree.merger
tm<-function(backbone=NULL,data,source.tree=NULL,
             tip.ages = NULL, node.ages = NULL,age.offset=NULL,title){

  data->dat
  if(is.null(backbone)){
    if(dat[1,]$reference%in%dat$bind)
      stop("Please indicate as first row of data the first pair of tips to build the tree on.")
    startsp<-unlist(dat[1,1:2])
    dat<-dat[-1,]
    if(any(dat$reference==startsp[2])){
      dat$poly[which(dat$reference==startsp[2])]<-TRUE
      dat$reference[which(dat$reference==startsp[2])]<-paste(startsp,collapse="-")
    }

    tree<-list()
    tree$edge<-matrix(c(3,1,3,2),ncol=2,byrow=TRUE)
    tree$tip.label<-startsp
    tree$Nnode<-1
    tree$edge.length<-rep(1,2)
    class(tree)<-"phylo"
  } else backbone->tree

  source.tree->tree2
  max(diag(vcv(tree)))->H
  H-diag(vcv(tree))->ages

  if(!is.null(age.offset)&&age.offset<0){
    ages+abs(age.offset)->ages
    H+abs(age.offset)->Hset
  }else H->Hset

  # if(is.null(min.branch)) min(tree$edge.length)->min.branch

  #### DATA CHECK ####
  if(!all(colnames(dat)%in%c("bind","reference","poly"))) {
    if(any(is.na(as.logical(dat[,3])))) stop("Check columns order: it should be 'bind', 'reference', 'poly'")
    warning("Colnames not matching: columns assumed to be ordered as 'bind','reference','poly\n'",immediate. = TRUE)
    colnames(dat)<-c("bind","reference","poly")
  }
  if(!is.logical(dat$poly)) as.logical(dat$poly)->dat$poly

  ### Check for wrong poly assignment ###
  if(any(apply(dat,1,function(k) (!grepl("-",k[2]))&!grepl("Genus",k[2])&!grepl("Clade",k[2])&as.logical(k[3])))){
    warning("Found poly=TRUE at binding two tips, automatically changed to FALSE\n",immediate.=TRUE)
    lapply(1:nrow(dat),function(k){
      if((!grepl("-",dat[k,2]))&!grepl("Genus",k[k,2])&!grepl("Clade",k[k,2])&dat[k,3]){
        dat[k,3]<<-FALSE
        message(paste("  Please check",dat[k,1],"&",dat[k,2],"in your dataset",sep=" "))
      }
    })
  }

  ### Check for genera and clades as bind ###
  if(any(grepl("Genus",dat$bind))){
    if(is.null(tree2)) stop("Please provide the source.tree")
    sapply(1:nrow(dat),function(k){
      if(grepl("Genus",dat[k,]$bind)){
        trimws(gsub("Genus ","",dat[k,]$bind))->genref
        getGenus(tree2,genref)->getgenref
        if(getgenref[2]==1) dat[k,]$bind<<-grep(genref,tree2$tip.label,value=TRUE) else
          dat[k,]$bind<<-paste(tips(tree2,getgenref[,3])[c(1,length(tips(tree2,getgenref[,3])))],collapse="-")
      }
    })
  }

  if(any(grepl("Clade",dat$bind))){
    sapply(1:nrow(dat),function(k){
      if(grepl("Clade",dat[k,]$bind)){
        if(!gsub("Clade ","",dat[k,]$bind)%in%tree2$node.label) stop("Required node.label not indicated on the source.tree")
        tips(tree2,Ntip(tree2)+which(tree2$node.label==gsub("Clade ","",dat[k,]$bind)))->claref
        dat[k,]$bind<<-paste(claref[c(1,length(claref))],collapse="-")
      }
    })
  }

  data.frame(dat,bind.type=sapply(strsplit(dat[,1],"-"),length))->dat
  if(any(dat$bind.type==2)&is.null(tree2)) stop("Please provide the source.tree")

  bind.all<-unlist(sapply(1:nrow(dat),function(x) {
    ifelse(dat[x,4]==2,des<-tips(tree2,getMRCA(tree2,strsplit(dat[x,1],"-")[[1]])),des<-dat[x,1])
    des
  }))
  if(any(duplicated(bind.all))) stop(paste(paste(bind.all[duplicated(bind.all)],collapse = ", "),"names duplicated in supplied tips"))

  if(all(dat$bind.type==1)&(!is.null(tree2))) tree2<-NULL
  # if(any(dat$bind%in%tree2$tip.label)) drop.tip(tree2,dat[which(dat$bind%in%tree2$tip.label),1])->tree2

  # dat$bind.tips<-NA
  # if(all(dat$bind.type==1)) dat$bind.tips<-dat$bind else{
  #   dat[which(dat$bind.type==2),]$bind.tips<-lapply(dat$bind[which(dat$bind.type==2)],function(x) tips(tree2,getMRCA(tree2,strsplit(x,"-")[[1]])))
  #   if(any(is.na(dat$bind.tips))) dat[which(is.na(dat$bind.tips)),]$bind.tips<-dat[which(is.na(dat$bind.tips)),]$bind
  # }

  dat$bind.tips<-dat$bind
  if(any(dat$bind.type==2))
    dat[which(dat$bind.type==2),]$bind.tips<-lapply(dat$bind[which(dat$bind.type==2)],function(x) tips(tree2,getMRCA(tree2,strsplit(x,"-")[[1]])))


  ### bind missing from tree2 ###
  if(!all(unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)){
    stop(paste(paste(unlist(dat[which(dat$bind.type==2),]$bind.tips)[which(
      unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)],
      collapse=", "),"not in source.tree"))
  }

  ### reference missing ###
  # if(!all(unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))){
  #   stop(paste(paste(unlist(strsplit(dat$reference,"-"))[which(!unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))],
  #                    collapse=","),"missing from the backbone, the source and the tips to be attached"))
  # }

  ### bind already on the backbone ###
  if(any(dat$bind%in%tree$tip.label)){
    warning(paste(paste(dat[which(dat$bind%in%tree$tip.label),1],collapse=", "),"removed from the backbone tree\n"),immediate. = TRUE)
    drop.tip(tree,dat[which(dat$bind%in%tree$tip.label),1])->tree
  }

  ### Check for genera and clades as references ###
  if(any(grepl("Genus",dat$reference))){
    dat$bindgen<-sapply(dat$bind.tips,function(j) unique(sapply(j,function(x) strsplit(x,"_")[[1]][1])))
    sapply(1:nrow(dat),function(k){
      if(grepl("Genus",dat[k,]$reference)){
        trimws(gsub("Genus ","",dat[k,]$reference))->genref
        try(getGenus(tree,genref),silent = TRUE)->getgenref
        if(!inherits(getgenref,"try-error")){
          if(getgenref[2]==1) gentree<-grep(genref,tree$tip.label,value=TRUE) else
            gentree<-paste(tips(tree,getgenref[,3])[c(1,length(tips(tree,getgenref[,3])))],collapse="-")
          genbind<-dat$bindgen[[k]]
          dat[k,]$reference<<-gentree

          if(genref%in%genbind){
            if(length(genbind)>1&any(sapply(dat$bindgen,function(w) all(w%in%genref))))
              dat[k,]$reference<<-paste(c(gentree,sapply(dat[which(sapply(dat$bindgen,function(w) all(w%in%genref))),]$bind.tips,"[[",1)),collapse="-")
          }else{
            if(any(sapply(dat$bindgen,function(w) any(w%in%genref))))
              dat[k,]$reference<<-paste(c(gentree,
                                          sapply(dat[which(sapply(dat$bindgen,function(w) any(w%in%genref))),]$bind.tips,"[[",1)),collapse="-")
          }
        }else stop(paste("Genus",genref,"missing from the backbone tree"))
      }
    })
    dat$bindgen<-NULL
  }


  if(any(grepl("Clade",dat$reference))){
    sapply(1:nrow(dat),function(k){
      if(grepl("Clade",dat[k,]$reference)){
        if(!gsub("Clade ","",dat[k,]$reference)%in%tree$node.label) stop("Required node.label not indicated on the tree")
        tips(tree,Ntip(tree)+which(tree$node.label==gsub("Clade ","",dat[k,]$reference)))->claref
        dat[k,]$reference<<-paste(claref[c(1,length(claref))],collapse="-")

      }
    })
  }


  ### duplicated reference ###
  table(dat$reference)->tab.ref
  if(any(tab.ref>1)){
    for(j in 1:length(which(tab.ref>1))){
      dat[which(dat$reference==names(which(tab.ref>1)[j])),]->ref.mult
      # if(any(isTRUE(ref.mult$poly))) ref.mult[-which(isTRUE(ref.mult$poly)),]->ref.mult
      if(any(ref.mult$poly)) ref.mult[which(!ref.mult$poly),]->ref.mult
      paste(strsplit(ref.mult$reference[1],"-")[[1]][1],strsplit(ref.mult$bind[1],"-")[[1]][1],sep="-")->ref.mult$reference[-1]
      # ref.mult[-1,]$poly<-TRUE
      ref.mult$poly[-1]<-TRUE
      dat[match(ref.mult$bind,dat$bind),]<-ref.mult
    }
  }

  if(!is.null(tree2)){
    dat$MRCAbind<-NA
    sapply(dat[which(dat$bind.type==2),]$bind.tips,function(x) getMRCA(tree2,x))->dat$MRCAbind[which(dat$bind.type==2)]

    if(any(dat$bind.type==2)){
      unlist(sapply(dat[which(dat$bind.type==2),]$MRCAbind,function(x) tree$tip.label[which(tree$tip.label%in%tips(tree2,x))]))->remt
      if(length(remt)>0){
        drop.tip(tree,remt)->tree
        ages[-match(remt,names(ages))]->ages
        warning(paste(paste(remt,collapse=", "),"already on the source tree: removed from the backbone tree\n"),immediate. = TRUE)
      }
    }
  }

  ### ordering ###
  strsplit(dat$reference,"-")->refs
  dat$ref.tree1<-NA

  lapply(refs,function(x){
    if(all(x%in%tree$tip.label)) NA else x[which(!x%in%tree$tip.label)]
  })->ref.tree

  dat[which(is.na(ref.tree)),][order(dat[which(is.na(ref.tree)),]$bind.type),]$ref.tree1<-seq(1,length(which(is.na(ref.tree))))

  if(any(!is.na(ref.tree))){
    dat[which(!is.na(ref.tree)),]->dat.new
    dat.new$ref.tree1<-ref.tree[-which(is.na(ref.tree))]

    if(any(unlist(dat.new$ref.tree1)%in%unlist(dat.new$bind.tips))){
      while(nrow(dat.new)>0){
        which(!sapply(dat.new$ref.tree1,function(w) any(w%in%unlist(dat.new$bind.tips))))->outs
        # which(!(dat.new[,7]%in%unlist(dat.new$bind.tips)|dat.new[,8]%in%unlist(dat.new$bind.tips)))
        if(length(outs)<1) stop("Recursive species attachment: check rows ",paste(rownames(dat.new),collapse = ", ")," in data")
        dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-max(dat$ref.tree1,na.rm=TRUE)+1:length(outs)
        dat.new[-outs,]->dat.new
      }
    }else{
      # which(!(dat.new[,7]%in%unlist(dat.new$bind.tips)|dat.new[,8]%in%unlist(dat.new$bind.tips)))->outs
      which(!sapply(dat.new$ref.tree1,function(w) any(w%in%unlist(dat.new$bind.tips))))->outs
      dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-
        max(dat$ref.tree1,na.rm=TRUE)+1:length(which(!dat.new$ref.tree1%in%dat.new[,1]))
    }
  }

  dat[order(dat$ref.tree1),]->dat

  ### binding ###
  pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
  for(k in 1:nrow(dat)){
    setTxtProgressBar(pb,k)
    if(length(strsplit(dat$reference,"-")[[k]])>1)
      getMRCA(tree,strsplit(dat$reference,"-")[[k]])->where.ref else
        which(tree$tip.label==strsplit(dat$reference,"-")[[k]])->where.ref

    if(where.ref!=(Ntip(tree)+1))
      tree$edge.length[which(tree$edge[,2]==where.ref)]->br.len else{
        if(is.null(tree$root.edge)||tree$root.edge==0) tree$root.edge<-mean(tree$edge.length)
        tree$root.edge->br.len
      }

    if(dat$bind.type[k]==1){
      if(dat$poly[k]){
        0->pos.ref
        max(tree$edge.length[which(tree$edge[,1]%in%where.ref)])->br.len
      } else br.len/2->pos.ref
      bind.tip(tree,dat$bind[k],where.ref,position = pos.ref,edge.length =br.len/2)->tree
    }else{
      extract.clade(tree2,dat$MRCAbind[k])->cla

      if(dat$poly[k]){
        0->pos.ref
        if(where.ref==(Ntip(tree)+1)&&(max(diag(vcv(cla)))+max(diag(vcv(cla)))/10)>H)
          # rescale(cla,"depth",(H+max(diag(vcv(cla)))/10))->cla else
          rescaleRR(cla,height=(H+max(diag(vcv(cla)))/10))->cla else
            if(where.ref!=(Ntip(tree)+1)&(max(diag(vcv(cla)))+max(diag(vcv(cla)))/10)>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))
              # rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]+max(diag(vcv(cla)))/10))->cla
              rescaleRR(cla,height=(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]+max(diag(vcv(cla)))/10))->cla
        cla$root.edge<-max(diag(vcv(cla)))/10
      }else{
        if(where.ref==(Ntip(tree)+1)) br.len/2+(max(diag(vcv(cla)))-max(diag(vcv(tree))))->pos.ref else {
          max(diag(vcv(cla)))+br.len/2->Hcla
          if((H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2])<H/1000){
            # rescale(cla,"depth",br.len/4)->cla
            rescaleRR(cla,height=br.len/4)->cla
            pos.ref<-br.len-br.len/4
          }else{
            if(Hcla>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),1])){
              # rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))->cla
              rescaleRR(cla,height=(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))->cla
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
  close(pb)
  message("Binding done!")

  if(is.null(backbone)){
    tree$edge.length<-rep(1,length(tree$edge.length))
    if(!is.null(tip.ages)) tree<-rescaleRR(tree,height=tip.ages[which.max(tip.ages)]+tip.ages[which.max(tip.ages)]/10)
  }

  trycal<-try({
    ### tips calibration ages ###
    if(is.null(tip.ages)){
      rep(0,length(tree$tip.label))->tip.ages
      names(tip.ages)<-tree$tip.label
      tip.ages[match(names(ages),names(tip.ages))]<-ages
    }else{
      if(any(!names(tip.ages)%in%tree$tip.label)){
        warning(paste(paste(names(tip.ages)[which(!names(tip.ages)%in%tree$tip.label)],collapse=", "),
                      "not on the final tree: removed from the vector of tip.ages\n"),immediate. = TRUE)
        tip.ages<-tip.ages[which(names(tip.ages)%in%tree$tip.label)]

      }
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
    }else node.ages<-c()

    ### age original root ###
    if(!is.null(backbone)&!getMRCA(tree,names(ages))%in%names(node.ages))
      node.ages<-setNames(c(node.ages,Hset),c(names(node.ages),getMRCA(tree,names(ages))))
    # {
    #   node.ages<-c(node.ages,Hset)
    #   names(node.ages)[length(node.ages)]<-getMRCA(tree,names(ages))
    #   }

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
      # sapply(names(ages.fix),function(x) getMRCA(tree,tips(tree2,as.numeric(x))))->names(ages.fix)
      sapply(1:length(ages.fix),function(x){
        getMRCA(tree,tips(tree2,as.numeric(names(ages.fix)[x])))->nodex
        tip.ages[match(tips(tree,nodex),names(tip.ages))]->tipx.ages
        if(any(tipx.ages>ages.fix[x])) NA else nodex
      })->names(ages.fix)
      ages.fix<-ages.fix[which(!is.na(names(ages.fix)))]

      if(any(!names(ages.fix)%in%names(node.ages)))
        c(ages.fix[which(!names(ages.fix)%in%names(node.ages))],node.ages)->node.ages
    }

    if(max(diag(vcv(tree)))>H&&(!(Ntip(tree)+1)%in%names(node.ages)))
      warning(paste("Root age not indicated: the tree root arbitrarily set at\n",round(max(diag(vcv(tree))),2)),immediate.=TRUE)

    scaleTree(tree,node.ages=node.ages,tip.ages =tip.ages)->tree.final->tree.plot
    message("Age Calibration done!")
  },silent=TRUE)
  if(inherits(trycal,"try-error")){
    warning("Age Calibration failed with the error: \n",trycal[1],
            "\n Returning the uncalibrated version of the tree",call.=FALSE)
    tree->tree.final->tree.plot
  }


  if(any(dat$bind.type==2)) lapply(dat$MRCAbind[which(dat$bind.type==2)],function(x)
    c(getMRCA(tree.plot,tips(tree2,x)),getDescendants(tree.plot,getMRCA(tree.plot,tips(tree2,x)))))->cla.plot else cla.plot<-c()

  c(which(tree.plot$tip.label%in%dat$bind[which(dat$bind.type==1)]),unlist(cla.plot))->all.plot
  # colo<-rep(scales::hue_pal()(2)[2],nrow(tree.plot$edge))
  # colo[which(tree.plot$edge[,2]%in%all.plot)]<-scales::hue_pal()(2)[1]
  colo<-rep("gray60",nrow(tree.plot$edge))
  colo[which(tree.plot$edge[,2]%in%all.plot)]<-"red"
  names(colo)<-tree.plot$edge[,2]

  par(mar=c(2,0,1,0))
  plot(tree.plot,edge.color=colo,edge.width=1.5,tip.color=colo[which(tree.plot$edge[,2]<=Ntip(tree.plot))])
  title(title)
  axisPhylo()


  return(tree.final)
}

# tree.merger, search.shift
plotTab<-function(tab,col.names=TRUE,row.names=FALSE,
                  text.cex=1,text.font=1,text.col="black",
                  text.hadj="center",text.vadj="center",
                  text.box=NULL,text.grid=NULL,text.highlight=NULL,
                  colN.cex=1,colN.font=1,colN.col="black",
                  colN.hadj="center",colN.vadj="center",
                  colN.box=NULL,colN.grid=NULL,colN.highlight=NULL,
                  rowN.cex=1,rowN.font=1,rowN.col="black",
                  rowN.hadj="center",rowN.vadj="center",
                  rowN.box=NULL,rowN.grid=NULL,rowN.highlight=NULL,
                  main=NULL,main.cex=2,main.font=2,main.col="black",
                  main.hadj="center",main.vadj="center",
                  main.box=NULL,main.highlight=NULL){
  if(any(is.na(tab))) tab[which(is.na(tab))]<-"NA"

  y0=0
  if(col.names) y0=y0+1
  if(!is.null(main)) y0=y0+1
  ifelse(row.names,1,0)->x0
  if(!is.matrix(tab)) tab<-as.matrix(tab)
  tab.start<-tab

  # text  -------------------------------------------------------------------
  # text adj  ---------------------------------------------------------------
  if(all(text.hadj=="center")){
    xc=0.5
    xad=0.5
  }else if(all(text.hadj=="left")){
    xc=0.05
    xad=0
  }else if(all(text.hadj=="right")){
    xc=0.95
    xad=1
  }else{
    if(!is.null(text.hadj["xcoord"])) xc<-text.hadj["xcoord"] else xc=0.5
    if(!is.null(text.hadj["xadj"])) xad<-text.hadj["xadj"] else xad=0.5
  }

  if(all(text.vadj=="center")){
    yc=0.5
    yad=0.5
  }else if(all(text.vadj=="top")){
    yc=0.05
    yad=1
  }else if(all(text.vadj=="bottom")){
    yc=0.95
    yad=0
  }else{
    if(!is.null(text.vadj["ycoord"])) yc<-text.vadj["ycoord"] else yc=0.5
    if(!is.null(text.vadj["yadj"])) yad<-text.vadj["yadj"] else yad=0.5
  }

  xcoords<-matrix(rep(seq(xc,ncol(tab)),each=nrow(tab)),ncol=ncol(tab))
  ycoords<-matrix(rep(seq(yc,nrow(tab)),ncol(tab)),ncol=ncol(tab))
  xadj<-matrix(rep(xad,each=nrow(tab)*ncol(tab)),ncol=ncol(tab))
  yadj<-matrix(rep(yad,each=nrow(tab)*ncol(tab)),ncol=ncol(tab))

  cex<-matrix(rep(text.cex,length.out=length(tab)),ncol = ncol(tab))
  font<-matrix(rep(text.font,length.out=length(tab)),ncol = ncol(tab))
  col<-matrix(rep(text.col,length.out=length(tab)),ncol = ncol(tab))

  # text box  ---------------------------------------------------------------
  if(!is.list(text.box)) text.box<-as.list(text.box)
  if(is.null(text.box$box)) text.bbox<-"o" else text.bbox<-text.box$box

  if(text.bbox!="n"){
    text.bbox%in%c("o","]","c", "7")->btop
    text.bbox%in%c("o","c", "u", "l")->bleft
    text.bbox%in%c("o","]", "u", "7")->bright
    text.bbox%in%c("o","]","c", "u", "l")->bbottom
    bframe<-c(bbottom,bleft,btop,bright)

    text.boxlty<-rep(1,length.out=4)
    text.boxlwd<-rep(1,length.out=4)
    text.boxcol<-rep("black",length.out=4)

    if(!is.null(text.box$lty)) text.boxlty[which(bframe)]<-rep(text.box$lty,length.out=sum(bframe))[1:sum(bframe)]
    if(!is.null(text.box$lwd)) text.boxlwd[which(bframe)]<-rep(text.box$lwd,length.out=sum(bframe))[1:sum(bframe)]
    if(!is.null(text.box$col)) text.boxcol[which(bframe)]<-rep(text.box$col,length.out=sum(bframe))[1:sum(bframe)]

    if(btop){
      box.top<-expression(segments(x0=x0,y0=y0,x1=ncol(tab),y1=y0,lty=text.boxlty[3],lwd=text.boxlwd[3],col=text.boxcol[3]))
    } else box.top<-NULL
    if(bleft){
      box.left<-expression(segments(x0=x0,y0=y0,x1=x0,y1=nrow(tab),lty=text.boxlty[2],lwd=text.boxlwd[2],col=text.boxcol[2]))
    } else box.left<-NULL
    if(bright){
      box.right<-expression(segments(x0=ncol(tab),y0=y0,x1=ncol(tab),y1=nrow(tab),lty=text.boxlty[4],lwd=text.boxlwd[4],col=text.boxcol[4]))
    } else box.right<-NULL
    if(bbottom){
      box.bottom<-expression(segments(x0=x0,y0=nrow(tab),x1=ncol(tab),y1=nrow(tab),lty=text.boxlty[3],lwd=text.boxlwd[3],col=text.boxcol[3]))
    } else box.bottom<-NULL
  }

  # text grid  --------------------------------------------------------------
  if(is.null(text.grid))
    text.grid<-data.frame(rown=c(1:nrow(tab),rep(0,ncol(tab))),
                          coln=c(rep(0,nrow(tab)),1:ncol(tab)),
                          col="black",lwd=1,lty=1)
  if(all(text.grid!="n")){
    if(is.null(ncol(text.grid))) text.grid<-as.data.frame(matrix(text.grid,nrow=1,dimnames=list(NULL,names(text.grid))))
    if(any(!colnames(text.grid)%in%c("rown","coln","col","lwd","lty")))
      stop("text.grid columns should include 'rown','coln','col','lwd','lty'")
    if("col"%in%colnames(text.grid)) colnames(text.grid)[which(colnames(text.grid)=="col")]<-"color"

    if(all(!c("rown","coln")%in%colnames(text.grid))){
      text.grid<-data.frame(rown=c(1:nrow(tab),rep(0,ncol(tab))),
                            coln=c(rep(0,nrow(tab)),1:ncol(tab)),text.grid)
    }
    text.grid<-as.data.frame(text.grid)
    if(is.null(text.grid$color)) text.grid<-data.frame(text.grid,color="black")
    if(is.null(text.grid$lwd)) text.grid<-data.frame(text.grid,lwd=1)
    if(is.null(text.grid$lty)) text.grid<-data.frame(text.grid,lty=1)
    text.grid[,which(colnames(text.grid)!="color")]<-apply(text.grid[,which(colnames(text.grid)!="color"),drop=FALSE],2,as.numeric)

    if(any(text.grid$rown==0)){
      gridcol<-text.grid[which(text.grid$rown==0),,drop=FALSE]
      text.grid<-text.grid[which(text.grid$rown!=0),,drop=FALSE]
      coltext.grid<-expression(segments(x0=x0+gridcol$coln[j],y0=y0,x1=x0+gridcol$coln[j],y1=nrow(tab),
                                        col=gridcol$color[j],lwd=gridcol$lwd[j],lty=gridcol$lty[j]))

    }else coltext.grid<-NULL

    if(any(text.grid$coln==0)){
      gridrow<-text.grid[which(text.grid$coln==0),,drop=FALSE]
      text.grid<-text.grid[which(text.grid$coln!=0),,drop=FALSE]
      rowtext.grid<-expression(segments(x0=x0,y0=y0+gridrow$rown[j],x1=ncol(tab),y1=y0+gridrow$rown[j],
                                        col=gridrow$color[j],lwd=gridrow$lwd[j],lty=gridrow$lty[j]))
    }else rowtext.grid<-NULL

    if(nrow(text.grid)>0){
      celltext.grid<-expression(polygon(x=c(x0+text.grid$coln[j]-1,x0+text.grid$coln[j],x0+text.grid$coln[j],x0+text.grid$coln[j]-1),
                                        y=c(y0+text.grid$rown[j]-1,y0+text.grid$rown[j]-1,y0+text.grid$rown[j],y0+text.grid$rown[j]),
                                        border=text.grid$color[j],lwd=text.grid$lwd[j],lty=text.grid$lty[j]))
    }else celltext.grid<-NULL

  }

  # text high  --------------------------------------------------------------
  if(is.null(text.highlight)) text.hig<-NULL else text.hig<-text.highlight
  if(!is.null(text.hig)){
    if(all(!c("rown","coln")%in%colnames(text.hig)))
      stop("text.highlight: please specify at least one rown and coln")
    if(any(!colnames(text.hig)%in%c("rown","coln","col","alpha")))
      stop("text.highlight: columns should include 'rown','coln','col','alpha'")
    if("col"%in%colnames(text.hig)) colnames(text.hig)[which(colnames(text.hig)=="col")]<-"color"

    text.hig<-as.data.frame(text.hig)
    if(is.null(text.hig$color)) text.hig<-data.frame(text.hig,color="gray80")
    if(is.null(text.hig$alpha)) text.hig<-data.frame(text.hig,alpha=0.8)
    text.hig[,which(colnames(text.hig)!="color")]<-apply(text.hig[,which(colnames(text.hig)!="color"),drop=FALSE],2,as.numeric)

    if(any(text.hig$rown==0)){
      higcol<-text.hig[which(text.hig$rown==0),,drop=FALSE]
      text.hig<-text.hig[which(text.hig$rown!=0),,drop=FALSE]
      coltext.hig<-expression(polygon(x=c(x0+higcol$coln[j]-1,x0+higcol$coln[j],
                                          x0+higcol$coln[j],x0+higcol$coln[j]-1),
                                      y=c(y0,y0,nrow(tab),nrow(tab)),
                                      border=NA,col=scales::alpha(higcol$color[j],higcol$alpha[j])))

    }else coltext.hig<-NULL

    if(any(text.hig$coln==0)){
      higrow<-text.hig[which(text.hig$coln==0),,drop=FALSE]
      text.hig<-text.hig[which(text.hig$coln!=0),,drop=FALSE]
      rowtext.hig<-expression(polygon(x=c(x0,ncol(tab),ncol(tab),x0),
                                      y=c(y0+higrow$rown[j]-1,y0+higrow$rown[j]-1,y0+higrow$rown[j],y0+higrow$rown[j]),
                                      border=NA,col=scales::alpha(higrow$color[j],higrow$alpha[j])))
    }else rowtext.hig<-NULL

    if(nrow(text.hig)>0){
      celltext.hig<-expression(polygon(x=c(x0+text.hig$coln[j]-1,x0+text.hig$coln[j],x0+text.hig$coln[j],x0+text.hig$coln[j]-1),
                                       y=c(y0+text.hig$rown[j]-1,y0+text.hig$rown[j]-1,y0+text.hig$rown[j],y0+text.hig$rown[j]),
                                       border=NA,col=scales::alpha(text.hig$color[j],text.hig$alpha[j])))
    }else celltext.hig<-NULL

  }

  # colN --------------------------------------------------------------------
  if(col.names){
    tab<-as.matrix(rbind(colnames(tab),tab))
    # colN adj --------------------------------------------------------------
    if(all(colN.hadj=="center")){
      xc.colN=0.5
      xad.colN=0.5
    }else if(all(colN.hadj=="left")){
      xc.colN=0.05
      xad.colN=0
    }else if(all(colN.hadj=="right")){
      xc.colN=0.95
      xad.colN=1
    }else{
      if(!is.null(colN.hadj["xcoord"])) xc.colN<-colN.hadj["xcoord"] else xc.colN=0.5
      if(!is.null(colN.hadj["xadj"])) xad.colN<-colN.hadj["xadj"] else xad.colN=0.5
    }

    if(all(colN.vadj=="center")){
      yc.colN=0.5
      yad.colN=0.5
    }else if(all(colN.vadj=="top")){
      yc.colN=0.05
      yad.colN=1
    }else if(all(colN.vadj=="bottom")){
      yc.colN=0.95
      yad.colN=0
    }else{
      if(!is.null(colN.vadj["ycoord"])) yc.colN<-colN.vadj["ycoord"] else yc.colN=0.5
      if(!is.null(colN.vadj["yadj"])) yad.colN<-colN.vadj["yadj"] else yad.colN=0.5
    }

    xcoords<-rbind(seq(xc.colN,ncol(xcoords)),xcoords)
    ycoords<-rbind(rep(yc.colN,ncol(ycoords)),ycoords+1)
    xadj<-rbind(rep(xad.colN,ncol(xadj)),xadj)
    yadj<-rbind(rep(yad.colN,ncol(yadj)),yadj)

    # colN box --------------------------------------------------------------
    if(!is.list(colN.box)) colN.box<-as.list(colN.box)
    if(is.null(colN.box$box)) colN.bbox<-"o" else colN.bbox<-colN.box$box

    if(colN.bbox!="n"){
      colN.bbox%in%c("o","]","c", "7")->colN.btop
      colN.bbox%in%c("o","c", "u", "l")->colN.bleft
      colN.bbox%in%c("o","]", "u", "7")->colN.bright
      colN.bbox%in%c("o","]","c", "u", "l")->colN.bbottom
      colN.bframe<-c(colN.bbottom,colN.bleft,colN.btop,colN.bright)

      colN.boxlty<-rep(1,length.out=4)
      colN.boxlwd<-rep(1,length.out=4)
      colN.boxcol<-rep("black",length.out=4)

      if(!is.null(colN.box$lty)) colN.boxlty[which(colN.bframe)]<-rep(colN.box$lty,length.out=sum(colN.bframe))[1:sum(colN.bframe)]
      if(!is.null(colN.box$lwd)) colN.boxlwd[which(colN.bframe)]<-rep(colN.box$lwd,length.out=sum(colN.bframe))[1:sum(colN.bframe)]
      if(!is.null(colN.box$col)) colN.boxcol[which(colN.bframe)]<-rep(colN.box$col,length.out=sum(colN.bframe))[1:sum(colN.bframe)]

      if(colN.btop){
        colNbox.top<-expression(segments(x0=x0,y0=y0-1,x1=ncol(tab),y1=y0-1,lty=colN.boxlty[3],lwd=colN.boxlwd[3],col=colN.boxcol[3]))
      } else colNbox.top<-NULL
      if(colN.bleft){
        colNbox.left<-expression(segments(x0=x0,y0=y0-1,x1=x0,y1=y0,lty=colN.boxlty[2],lwd=colN.boxlwd[2],col=colN.boxcol[2]))
      } else colNbox.left<-NULL
      if(colN.bright){
        colNbox.right<-expression(segments(x0=ncol(tab),y0=y0-1,x1=ncol(tab),y1=y0,lty=colN.boxlty[4],lwd=colN.boxlwd[4],col=colN.boxcol[4]))
      } else colNbox.right<-NULL
      if(colN.bbottom){
        colNbox.bottom<-expression(segments(x0=x0,y0=y0,x1=ncol(tab),y1=y0,lty=colN.boxlty[1],lwd=colN.boxlwd[1],col=colN.boxcol[1]))
      } else colNbox.bottom<-NULL
    }


    # colN grid --------------------------------------------------------------
    if(is.null(colN.grid))
      colN.grid<-data.frame(coln=1:ncol(tab),col="black",lwd=1,lty=1)

    if(all(colN.grid!="n")){
      if(is.null(ncol(colN.grid)))
        colN.grid<-as.data.frame(matrix(colN.grid,byrow=TRUE,nrow=1,dimnames=list(NULL,names(colN.grid))))
      if(any(!colnames(colN.grid)%in%c("coln","col","lwd","lty")))
        stop("colN.grid columns should include 'coln','col','lwd','lty'")
      if("col"%in%colnames(colN.grid)) colnames(colN.grid)[which(colnames(colN.grid)=="col")]<-"color"
      colN.grid<-as.data.frame(colN.grid)
      if(is.null(colN.grid$coln)) colN.grid<-data.frame(coln=1:ncol(tab),colN.grid)
      if(is.null(colN.grid$color)) colN.grid<-data.frame(colN.grid,color="black")
      if(is.null(colN.grid$lwd)) colN.grid<-data.frame(colN.grid,lwd=1)
      if(is.null(colN.grid$lty)) colN.grid<-data.frame(colN.grid,lty=1)
      colN.grid[,which(colnames(colN.grid)!="color")]<-apply(colN.grid[,which(colnames(colN.grid)!="color"),drop=FALSE],2,as.numeric)

      col.grid<-expression(segments(x0=x0+colN.grid$coln[j],x1=x0+colN.grid$coln[j],y0=y0-1,y1=y0,
                                    col=colN.grid$color[j],lwd=colN.grid$lwd[j],lty=colN.grid$lty[j]))

    }

    # colN high -------------------------------------------------------------
    if(is.null(colN.highlight))
      colN.hig<-data.frame(coln=1:ncol(tab),col="gray80",alpha=1) else
        colN.hig<-colN.highlight

    if(all(colN.hig!="n")){
      if(is.null(ncol(colN.hig)))
        colN.hig<-as.data.frame(matrix(colN.hig,byrow=TRUE,nrow=1,dimnames=list(NULL,names(colN.hig))))
      if(any(!colnames(colN.hig)%in%c("coln","col","alpha")))
        stop("colN.highlight: columns should include 'coln','col','alpha'")
      if("col"%in%colnames(colN.hig)) colnames(colN.hig)[which(colnames(colN.hig)=="col")]<-"color"
      colN.hig<-as.data.frame(colN.hig)
      if(is.null(colN.hig$coln)) colN.hig<-data.frame(coln=1:ncol(tab),colN.hig)
      if(is.null(colN.hig$color)) colN.hig<-data.frame(colN.hig,color="gray80")
      if(is.null(colN.hig$alpha)) colN.hig<-data.frame(colN.hig,alpha=0.8)
      colN.hig[,which(colnames(colN.hig)!="color")]<-apply(colN.hig[,which(colnames(colN.hig)!="color"),drop=FALSE],2,as.numeric)

      col.hig<-expression(polygon(x=x0+c(colN.hig$coln[j]-1,colN.hig$coln[j],colN.hig$coln[j],colN.hig$coln[j]-1),
                                  y=c(y0-1,y0-1,y0,y0),
                                  border=NA,col=scales::alpha(colN.hig$color[j],colN.hig$alpha[j])))

    }

    cex<-rbind(colN.cex,cex)
    font<-rbind(colN.font,font)
    col<-rbind(colN.col,col)

  }

  # rowN --------------------------------------------------------------------
  if(row.names){
    # rowN adj --------------------------------------------------------------
    tab<-data.frame(row=rownames(tab),tab)
    if(all(rowN.hadj=="center")){
      xc.rowN=0.5
      xad.rowN=0.5
    }else if(all(rowN.hadj=="left")){
      xc.rowN=0.05
      xad.rowN=0
    }else if(all(rowN.hadj=="right")){
      xc.rowN=0.95
      xad.rowN=1
    }else{
      if(!is.null(rowN.hadj["xcoord"])) xc.rowN<-rowN.hadj["xcoord"] else xc.rowN=0.5
      if(!is.null(rowN.hadj["xadj"])) xad.rowN<-rowN.hadj["xadj"] else xad.rowN=0.5
    }

    if(all(rowN.vadj=="center")){
      yc.rowN=0.5
      yad.rowN=0.5
    }else if(all(rowN.vadj=="top")){
      yc.rowN=0.05
      yad.rowN=1
    }else if(all(rowN.vadj=="bottom")){
      yc.rowN=0.95
      yad.rowN=0
    }else{
      if(!is.null(rowN.vadj["ycoord"])) yc.rowN<-rowN.vadj["ycoord"] else yc.rowN=0.5
      if(!is.null(rowN.vadj["yadj"])) yad.rowN<-rowN.vadj["yadj"] else yad.rowN=0.5
    }

    xcoords<-cbind(rep(xc.rowN,nrow(xcoords)),xcoords+1)
    ycoords<-cbind(seq(yc.rowN,nrow(ycoords)),ycoords)
    xadj<-cbind(rep(xad.rowN,nrow(xadj)),xadj)
    yadj<-cbind(rep(yad.rowN,nrow(xadj)),yadj)

    # rowN box --------------------------------------------------------------
    if(!is.list(rowN.box)) rowN.box<-as.list(rowN.box)
    if(is.null(rowN.box$box)) rowN.bbox<-"o" else rowN.bbox<-rowN.box$box
    if(rowN.bbox!="n"){
      rowN.bbox%in%c("o","]","c", "7")->rowN.btop
      rowN.bbox%in%c("o","c", "u", "l")->rowN.bleft
      rowN.bbox%in%c("o","]", "u", "7")->rowN.bright
      rowN.bbox%in%c("o","]","c", "u", "l")->rowN.bbottom
      rowN.bframe<-c(rowN.bbottom,rowN.bleft,rowN.btop,rowN.bright)

      rowN.boxlty<-rep(1,length.out=4)
      rowN.boxlwd<-rep(1,length.out=4)
      rowN.boxcol<-rep("black",length.out=4)

      if(!is.null(rowN.box$lty)) rowN.boxlty[which(rowN.bframe)]<-rep(rowN.box$lty,length.out=sum(rowN.bframe))[1:sum(rowN.bframe)]
      if(!is.null(rowN.box$lwd)) rowN.boxlwd[which(rowN.bframe)]<-rep(rowN.box$lwd,length.out=sum(rowN.bframe))[1:sum(rowN.bframe)]
      if(!is.null(rowN.box$col)) rowN.boxcol[which(rowN.bframe)]<-rep(rowN.box$col,length.out=sum(rowN.bframe))[1:sum(rowN.bframe)]

      if(rowN.btop){
        rowNbox.top<-expression(segments(x0=0,y0=y0,x1=x0,y1=y0,lty=rowN.boxlty[3],lwd=rowN.boxlwd[3],col=rowN.boxcol[3]))
      } else rowNbox.top<-NULL
      if(rowN.bleft){
        rowNbox.left<-expression(segments(x0=0,y0=y0,x1=0,y1=nrow(tab),lty=rowN.boxlty[2],lwd=rowN.boxlwd[2],col=rowN.boxcol[2]))
      } else rowNbox.left<-NULL
      if(rowN.bright){
        rowNbox.right<-expression(segments(x0=x0,y0=y0,x1=x0,y1=nrow(tab),lty=rowN.boxlty[4],lwd=rowN.boxlwd[4],col=rowN.boxcol[4]))
      } else rowNbox.right<-NULL
      if(rowN.bbottom){
        rowNbox.bottom<-expression(segments(x0=0,y0=nrow(tab),x1=x0,y1=nrow(tab),lty=rowN.boxlty[1],lwd=rowN.boxlwd[1],col=rowN.boxcol[1]))
      } else rowNbox.bottom<-NULL
    }

    # rowN grid  --------------------------------------------------------------
    if(is.null(rowN.grid))
      rowN.grid<-data.frame(rown=1:nrow(tab),col="black",lwd=1,lty=1)

    if(all(rowN.grid!="n")){
      if(is.null(ncol(rowN.grid)))
        rowN.grid<-as.data.frame(matrix(rowN.grid,byrow=TRUE,nrow=1,dimnames=list(NULL,names(rowN.grid))))
      if(any(!colnames(rowN.grid)%in%c("rown","col","lwd","lty")))
        stop("rowN.grid columns should include 'rown','col','lwd','lty'")
      if("col"%in%colnames(rowN.grid)) colnames(rowN.grid)[which(colnames(rowN.grid)=="col")]<-"color"
      rowN.grid<-as.data.frame(rowN.grid)
      if(is.null(rowN.grid$rown)) rowN.grid<-data.frame(rown=1:nrow(tab),rowN.grid)
      if(is.null(rowN.grid$color)) rowN.grid<-data.frame(rowN.grid,color="black")
      if(is.null(rowN.grid$lwd)) rowN.grid<-data.frame(rowN.grid,lwd=1)
      if(is.null(rowN.grid$lty)) rowN.grid<-data.frame(rowN.grid,lty=1)
      rowN.grid[,which(colnames(rowN.grid)!="color")]<-apply(rowN.grid[,which(colnames(rowN.grid)!="color"),drop=FALSE],2,as.numeric)

      row.grid<-expression(segments(x0=x0-1,x1=x0,y0=y0+rowN.grid$rown[j],y1=y0+rowN.grid$rown[j],
                                    col=rowN.grid$color[j],lwd=rowN.grid$lwd[j],lty=rowN.grid$lty[j]))

    }

    # rowN high -------------------------------------------------------------
    if(is.null(rowN.highlight)) rowN.hig<-data.frame(rown=1:nrow(tab.start),col="gray80",alpha=1) else
      rowN.hig<-rowN.highlight
    if(all(rowN.hig!="n")){
      if(is.null(ncol(rowN.hig)))
        rowN.hig<-as.data.frame(matrix(rowN.hig,byrow=TRUE,nrow=1,dimnames=list(NULL,names(rowN.hig))))
      if(any(!colnames(rowN.hig)%in%c("rown","col","alpha")))
        stop("rowN.highlight: columns should include 'rown','col','alpha'")
      if("col"%in%colnames(rowN.hig)) colnames(rowN.hig)[which(colnames(rowN.hig)=="col")]<-"color"
      rowN.hig<-as.data.frame(rowN.hig)
      if(is.null(rowN.hig$rown)) rowN.hig<-data.frame(rown=1:nrow(tab.start),rowN.hig)
      if(is.null(rowN.hig$color)) rowN.hig<-data.frame(rowN.hig,color="gray80")
      if(is.null(rowN.hig$alpha)) rowN.hig<-data.frame(rowN.hig,alpha=0.8)
      rowN.hig[,which(colnames(rowN.hig)!="color")]<-apply(rowN.hig[,which(colnames(rowN.hig)!="color"),drop=FALSE],2,as.numeric)

      row.hig<-expression(polygon(x=c(0,1,1,0),
                                  y=y0+c(rowN.hig$rown[j]-1,rowN.hig$rown[j]-1,rowN.hig$rown[j],rowN.hig$rown[j]),
                                  border=NA,col=scales::alpha(rowN.hig$color[j],rowN.hig$alpha[j])))
    }

    if(col.names){
      rncex<-rep(rowN.cex,length.out=nrow(cex))
      rowN.cex<-c(rncex[length(rncex)],rncex[-length(rncex)])
      rnfont<-rep(rowN.font,length.out=nrow(font))
      rowN.font<-c(rnfont[length(rnfont)],rnfont[-length(rnfont)])
      rncol<-rep(rowN.col,length.out=nrow(col))
      rowN.col<-c(rncol[length(rncol)],rncol[-length(rncol)])
    }
    cex<-cbind(rowN.cex,cex)
    font<-cbind(rowN.font,font)
    col<-cbind(rowN.col,col)
  }
  if(row.names&col.names) tab[1,1]<-""

  # main  -------------------------------------------------------------------
  if(!is.null(main)){
    # main adj --------------------------------------------------------------
    if(all(main.hadj=="center")){
      xc.main=ncol(tab)/2
      xad.main=0.5
    }else if(all(main.hadj=="left")){
      xc.main=0.05
      xad.main=0
    }else if(all(main.hadj=="right")){
      xc.main=ncol(tab)-0.05
      xad.main=1
    }else{
      if(!is.null(main.hadj["xcoord"])) xc.main<-main.hadj["xcoord"] else xc.main=0.5
      if(!is.null(main.hadj["xadj"])) xad.main<-main.hadj["xadj"] else xad.main=0.5
    }

    if(all(main.vadj=="center")){
      yc.main=0.5
      yad.main=0.5
    }else if(all(main.vadj=="top")){
      yc.main=0.05
      yad.main=1
    }else if(all(main.vadj=="bottom")){
      yc.main=0.95
      yad.main=0
    }else{
      if(!is.null(main.vadj["ycoord"])) yc.main<-main.vadj["ycoord"] else yc.main=0.5
      if(!is.null(main.vadj["yadj"])) yad.main<-main.vadj["yadj"] else yad.main=0.5
    }

    main.x<-xc.main
    main.y<-yc.main
    main.xadj<-xad.main
    main.yadj<-yad.main
    tab<-as.matrix(rbind("",tab))
    ycoords<-ycoords+1

    # main box --------------------------------------------------------------
    if(!is.null(main.box)){
      if(!is.list(main.box)) main.box<-as.list(main.box)
      main.bbox<-main.box$box
      main.bbox%in%c("o","]","c", "7")->main.btop
      main.bbox%in%c("o","c", "u", "l")->main.bleft
      main.bbox%in%c("o","]", "u", "7")->main.bright
      main.bbox%in%c("o","]","c", "u", "l")->main.bbottom
      main.bframe<-c(main.bbottom,main.bleft,main.btop,main.bright)

      main.boxlty<-rep(1,length.out=4)
      main.boxlwd<-rep(1,length.out=4)
      main.boxcol<-rep("black",length.out=4)

      if(!is.null(main.box$lty)) main.boxlty[which(main.bframe)]<-rep(main.box$lty,length.out=sum(main.bframe))[1:sum(main.bframe)]
      if(!is.null(main.box$lwd)) main.boxlwd[which(main.bframe)]<-rep(main.box$lwd,length.out=sum(main.bframe))[1:sum(main.bframe)]
      if(!is.null(main.box$col)) main.boxcol[which(main.bframe)]<-rep(main.box$col,length.out=sum(main.bframe))[1:sum(main.bframe)]

      if(main.btop){
        mainbox.top<-expression(segments(x0=x0,y0=y0-2,x1=ncol(tab),y1=y0-2,lty=main.boxlty[3],lwd=main.boxlwd[3],col=main.boxcol[3]))
      } else mainbox.top<-NULL
      if(main.bleft){
        mainbox.left<-expression(segments(x0=x0,y0=y0-2,x1=x0,y1=y0-1,lty=main.boxlty[2],lwd=main.boxlwd[2],col=main.boxcol[2]))
      } else mainbox.left<-NULL
      if(main.bright){
        mainbox.right<-expression(segments(x0=ncol(tab),y0=y0-2,x1=ncol(tab),y1=y0-1,lty=main.boxlty[4],lwd=main.boxlwd[4],col=main.boxcol[4]))
      } else mainbox.right<-NULL
      if(main.bbottom){
        mainbox.bottom<-expression(segments(x0=x0,y0=y0-1,x1=ncol(tab),y1=y0-1,lty=main.boxlty[1],lwd=main.boxlwd[1],col=main.boxcol[1]))
      } else mainbox.bottom<-NULL
    }

    # main high -------------------------------------------------------------
    if(!is.null(main.highlight)){
      if(is.null(ncol(main.highlight)))
        main.highlight<-as.data.frame(matrix(main.highlight,byrow=TRUE,nrow=1,dimnames=list(NULL,names(main.highlight))))
      if(any(!colnames(main.highlight)%in%c("col","alpha")))
        stop("main.highlight: columns should include 'col','alpha'")

      if(is.null(main.highlight$col)) main.highlight<-c(main.highlight,col="gray80")
      if(is.null(main.highlight$alpha)) main.highlight<-c(main.highlight,alpha=0.8)
      main.hig<-expression(polygon(x=c(x0,ncol(tab),ncol(tab),x0),y=c(0,0,y0-1,y0-1),border=NA,
                                   col=scales::alpha(main.highlight$col,as.numeric(main.highlight$alpha))))
    }

  }

  # plot  -------------------------------------------------------------------
  as.matrix(tab)->tab
  plot(NA,xlim=c(-0.1,ncol(tab)+0.1),ylim=c(nrow(tab)+0.1,-0.1),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

  if(col.names&&all(colN.hig!="n")) lapply(1:nrow(colN.hig),function(j) eval(col.hig))
  if(row.names&&all(rowN.hig!="n")) lapply(1:nrow(rowN.hig),function(j) eval(row.hig))
  if(col.names&&all(colN.grid!="n")) lapply(1:nrow(colN.grid),function(j) eval(col.grid))
  if(row.names&&all(rowN.grid!="n")) lapply(1:nrow(rowN.grid),function(j) eval(row.grid))

  if(!is.null(text.hig)&&!is.null(coltext.hig)) lapply(1:nrow(higcol),function(j) eval(coltext.hig))
  if(!is.null(text.hig)&&!is.null(rowtext.hig)) lapply(1:nrow(higrow),function(j) eval(rowtext.hig))
  if(!is.null(text.hig)&&!is.null(celltext.hig)) lapply(1:nrow(text.hig),function(j) eval(celltext.hig))

  if(all(text.grid!="n")&&!is.null(coltext.grid)) lapply(1:nrow(gridcol),function(j) eval(coltext.grid))
  if(all(text.grid!="n")&&!is.null(rowtext.grid)) lapply(1:nrow(gridrow),function(j) eval(rowtext.grid))
  if(all(text.grid!="n")&&!is.null(celltext.grid)) lapply(1:nrow(text.grid),function(j) eval(celltext.grid))

  if(text.bbox!="n"){
    eval(box.bottom)
    eval(box.top)
    eval(box.left)
    eval(box.right)
  }

  if(col.names&&colN.bbox!="n"){
    eval(colNbox.bottom)
    eval(colNbox.top)
    eval(colNbox.left)
    eval(colNbox.right)
  }

  if(row.names&&rowN.bbox!="n"){
    eval(rowNbox.bottom)
    eval(rowNbox.top)
    eval(rowNbox.left)
    eval(rowNbox.right)
  }

  if(!is.null(main)){
    if(!is.null(main.highlight)) eval(main.hig)
    if(!is.null(main.box)){
      eval(mainbox.bottom)
      eval(mainbox.top)
      eval(mainbox.left)
      eval(mainbox.right)
    }
    text(x=main.x,y=main.y,labels = main,font=main.font,
         cex=main.cex,col=main.col,adj=c(main.xadj,main.yadj))
    tab<-tab[-1,]
  }
  invisible(mapply(x=xcoords,y=ycoords,lab = tab,adj1=xadj,adj2=yadj,col=col,font=font,cex=cex,
         function(x,y,lab,adj1,adj2,col,font,cex){
           text(x=x,y=y,labels = lab,adj=c(adj1,adj2),col=col,font=font,cex=cex)
         }))

}


