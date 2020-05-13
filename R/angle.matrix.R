#' @title Ontogenetic shape vectors analysis
#'
#' @description This function computes and compares ontogenetic vectors among species in a tree.
#' @usage angle.matrix(RR,node,Y=NULL,select.axes=c("no","yes"),
#' type=c("phenotypes","rates"),cova=NULL,clus=0.5)
#' @param RR an object produced by \code{\link{RRphylo}}.
#' @param node the number identifying the most recent common ancestor to all the species the user wants ontogenetic vectors be computed.
#' @param Y multivariate trait values at tips.
#' @param type specifies weather to perform the analysis on phenotypic (\code{"phenotypes"}) or rate (\code{"rates"}) vectors.
#' @param select.axes if \code{"yes"}, \code{Y} variables are individually regressed against developmental stages and only significant variables are retained to compute ontogenetic vectors. All variables are retained otherwise.
#' @param cova the covariate to be indicated if its effect on rate values must be accounted for. Contrary to \code{RRphylo}, \code{cova} needs to be as long as the number of tips in the tree. As the covariate only affects rates computation, there is no covariate to provide when \code{type = "phenotypes"}.
#' @param clus the proportion of clusters to be used in parallel computing.
#' @importFrom stats confint logLik var xtabs
#' @importFrom smatr sma
#' @importFrom rlist list.append
#' @export
#' @details The \code{angle.matrix} function takes as objects a phylogenetic tree (retrieved directly from an \code{\link{RRphylo}} object), including the different ontogenetic stages of each species as polytomies. Names at tips must be written as species ID and stage number separated by the underscore. The \code{RRphylo} object \code{angle.matrix} is fed with is just used to extract the dichotomized version of the phylogeny. This is necessary because node numbers change randomly at dichotomizing non-binary trees. However, when performing \code{angle.matrix} with the covariate the \code{RRphylo } object must be produced without accounting for the covariate. Furthermore, as the covariate only affects the rates computation, it makes no sense to use it when computing vectors for phenotypic variables. Once angles and vectors are computed, \code{angle.matrix} performs two tests by means of standard major axis (SMA) regression. For each species pair, the "biogenetic test" verifies whether the angle between species grows during development, meaning that the two species becomes less similar to each other during growth. The "paedomorphosis test" tells whether there is heterochronic shape change in the data. Under paedomorphosis, the adult stages of one (paedomorphic) species will resemble the juvenile stages of the other (peramorphic) species. The test regresses the angles formed by the shapes at different ontogenetic stages of a species to the shape at the youngest stage of the other in the pair, against age. Then, it tests whether the two regression lines (one per species) have different slopes, and whether they have different signs. If the regression lines point to different directions, it means that one of the two species in the pair resembles, with age, the juveniles of the other, indicating paedomorphosis. Ontogenetic vectors of individual species are further computed, in reference to the MRCA of the pair, and to the first stage of each species (i.e. intraspecifically). Importantly, the size of the ontogenetic vectors of rates tell whether the two species differ in terms of developmental rate, which is crucial to understand which process is behind paedomorphosis, where it applies.While performing the analysis, the function prints messages on-screen informing about tests results. If \code{select.axes = "yes"}, informs the user about which phenotypic variables are used. Secondly, it specifies whether ontogenetic vectors to MRCA, and intraspecific ontogenetic vectors significantly differ in angle or size between species pairs. Then, for each species pair, it indicates if the biogenetic law and paedomorphosis apply.
#' @return A list containing 4 objects:
#' @return \enumerate{\item \strong{$regression.matrix} a 'list' including 'angles between species' and 'angles between species to MRCA' matrices for all possible combinations of species pairs from the two sides descending from the MRCA. For each matrix, corresponding biogenetic and paedomorphosis tests are reported.
#' \item \strong{$angles.2.MRCA.and.vector.size} a 'data.frame' including angles between the resultant vector of species and the MRCA and the size of the resultant vector computed from species to MRCA, per stage per species.
#' \item \strong{$ontogenetic.vectors2MRCA} a 'data.frame' including angle, size, and corresponding x and y components, of ontogenetic vectors computed between each species and the MRCA. For both angle and size, the p-value for the difference between species pairs is reported.
#' \item \strong{$ontogenetic.vectors.to.1st.stage} a 'list' containing:
#' \itemize{\item$matrices: for all possible combinations of species pairs from the two sides descending form the MRCA, the upper triangle of the matrix contains the angles between different ontogenetic stages for the first species. The same applies to the lower triangle, but for the second species.
#' \item$vectors: for all possible combinations of species pairs from the two sides descending form the MRCA, angles and sizes of ontogenetic vectors computed to the first stage of each species. For both, the p-value for the difference between the species pair is reported.
#' }
#' }
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#'   \dontrun{
#'   data("DataApes")
#'   DataApes$PCstage->PCstage
#'   DataApes$Tstage->Tstage
#'   DataApes$CentroidSize->CS
#'
#'   cc<- 2/parallel::detectCores()
#'   RRphylo(tree=Tstage,y=PCstage,clus=cc)->RR
#' # Case 1. without accounting for the effect of a covariate
#'
#'  # Case 1.1 selecting shape variables that show significant relationship with age
#'   # on phenotypic vectors
#'     angle.matrix(RR,node=72,Y=PCstage,select.axes="yes",type="phenotypes",clus=cc)
#'   # on rates vectors
#'     angle.matrix(RR,node=72,Y=PCstage,select.axes="yes",type="rates",clus=cc)
#'
#'  # Case 1.2 using all shape variables
#'   # on phenotypic vectors
#'     angle.matrix(RR,node=72,Y=PCstage,select.axes="no",type="phenotypes",clus=cc)
#'   # on rates vectors
#'     angle.matrix(RR,node=72,Y=PCstage,select.axes="no",type="rates",clus=cc)
#'
#'
#' # Case 2. accounting for the effect of a covariate (on rates vectors only)
#'
#'  # Case 2.1 selecting shape variables that show significant relationship with age
#'    angle.matrix(RR,node=72,Y=PCstage,select.axes="yes",type="rates", cova=CS,clus=cc)
#'
#'
#'  # Case 2.2 using all shape variables
#'    angle.matrix(RR,node=72,Y=PCstage,select.axes="no",type="rates",cova=CS,clus=cc)
#'   }


angle.matrix<-function(RR,node,Y=NULL,select.axes=c("no","yes"),type=c("phenotypes","rates"),cova=NULL,clus=0.5)
{
  #require(smatr)
  #require(rlist)
  #require(geiger)

  unitV<-function(x){ sum(x^2)^.5 }

  ### Choose whether to select only significant axes or not ###
  match.arg(select.axes)
  node->n
  RR$tree->tree
  if(select.axes=="yes")
  {
    data.frame(species=unlist(lapply(strsplit(rownames(Y),"_"),"[[",1)),stage=unlist(lapply(strsplit(rownames(Y),"_"),"[[",2)))->df
    data.frame(df,Y)->fulldata
    as.numeric(fulldata[,2])->fulldata[,2]
    fulldata[which(rownames(fulldata)%in%tips(multi2di(tree),node)),]->fd
    which(apply(sapply(summary(lm(as.matrix(fd[,-c(1:2)])~fd[,1]:fd[,2])), function(x) x$coef[,4])[-1,],2,function(x) length(which(x<0.05)))==length(unique(fd$species)))->sel
    if(length(sel)>1){
      Y[,sel]->y.sel
      print(paste("I am using", paste(colnames(y.sel),collapse=" & "),"as variables"))
    }else{
      Y->y.sel
      print(paste("I am using the whole dataset"))
    }


    if(is.null(cova)){

      RRphylo(tree,y.sel,clus=clus)->rr.sel

    }else{

      RRphylo(tree,cova)->RRcova
      c(RRcova$aces,cova)->covari
      names(covari)<-c(rownames(RRcova$aces),names(cova))
      RRphylo(tree,y.sel,cov=covari,clus=clus)->rr.sel
    }
    match.arg(type)
    if(type=="rates") {
      rr.sel$multiple.rates->y.onto
      evo.dir(rr.sel,angle.dimension = "rates",pair.type="node",node=n,random="no")->evo.p
    }else{
      y.sel->y.onto
      evo.dir(rr.sel,y.type = "original",angle.dimension = "phenotypes",y=y.sel,pair.type="node",node=n,random="no")->evo.p
    }
  } else {


    if(is.null(cova)){

      RRphylo(tree,Y,clus=clus)->RR

    }else{

      RRphylo(tree,cova)->RRcova
      c(RRcova$aces,cova)->covari
      names(covari)<-c(rownames(RRcova$aces),names(cova))
      RRphylo(tree,Y,cov=covari,clus=clus)->RR
    }

    if(type=="rates") {
      RR$multiple.rates->y.onto
      evo.dir(RR,angle.dimension = "rates",pair.type="node",node=n,random="no")->evo.p
    }else{
      Y->y.onto
      evo.dir(RR,y.type = "original",angle.dimension = "phenotypes",y=Y,pair.type="node",node=n,random="no")->evo.p
    }
  }
  ### End select axes ###




  deg2rad <- function(deg) {(deg * pi) / (180)}
  rad2deg <- function(rad) {(rad * 180) / (pi)}

  node->n
  retrieve.angles(evo.p,wishlist="angles.between.species",focus="node",write2csv="no",node=n,random="no")->ret.ang
  retrieve.angles(evo.p,wishlist="anglesMRCA",focus="node",write2csv="no",node=72,random="no")[,c(1,2,4)]->vecS
  retrieve.angles(evo.p,wishlist="anglesMRCA",focus="node",write2csv="no",node=n,random="no")->angMRCA
  ret.ang[,c(2,4,3)]->btwsp
  data.frame(do.call(rbind,strsplit(as.character(btwsp[,1]), "/")),btwsp[,c(2,3)])->ansp
  colnames(ansp)<-c("sp1","sp2","anglesBTWspecies","angleBTWspecies2MRCA")
  ansp[order(ansp[,1],ansp[,2]),]->ansp
  matrix(ncol=length(unique(ansp[,1])),nrow=length(unique(ansp[,2])))->anglesBTW.mat
  matrix(ncol=length(unique(ansp[,1])),nrow=length(unique(ansp[,2])))->anglesBTW2MRCA.mat
  colnames(anglesBTW.mat)<-unique(ansp[,1])
  rownames(anglesBTW.mat)<-unique(ansp[,2])
  colnames(anglesBTW2MRCA.mat)<-unique(ansp[,1])
  rownames(anglesBTW2MRCA.mat)<-unique(ansp[,2])

  for (i in 1:length(colnames(anglesBTW.mat))){
    as.numeric(as.character(ansp[which(ansp[,1]%in%colnames(anglesBTW.mat)[i]),3]))->anglesBTW.mat[,i]
    as.numeric(as.character(ansp[which(ansp[,1]%in%colnames(anglesBTW2MRCA.mat)[i]),4]))->anglesBTW2MRCA.mat[,i]
  }

  vecS[order(vecS[,2]),]->vecS
  angMRCA[order(angMRCA[,2]),]->angMRCA
  data.frame(angMRCA,"vector.size"=vecS[,3])->anglesMRCA


  data.frame(species=anglesMRCA$species,angle=as.numeric(as.character(anglesMRCA$angle)),size=as.numeric(as.character(anglesMRCA$vector.size)))->VS
  cos(deg2rad(VS$angle))*VS$size->xx
  sin(deg2rad(VS$angle))*VS$size->yy
  data.frame(VS,XX=xx,YY=yy)->VS
  unique(unlist(lapply(strsplit(as.character(VS$species),split="_"),"[[",1)))->group
  if(length(group)==1) stop("Check group names or attempt to compare single group")


  ### Ontogenetic vectors to MRCA ####

  gg<-list()
  for (i in 1:length(group)){
    VS[grep(group[i],VS$species),]->gg[[i]]
  }
  names(gg)<-group

  lapply(gg,function(x) apply(x[,4:5],2,function(z) z[1]-sum(z[2:length(z)])))->x.and.y
  lapply(x.and.y,function(x) rad2deg(atan(x[2]/x[1])))->vector.angle
  lapply(x.and.y,function(x) (x[1]^2+x[2]^2)^.5)->vector.length
  rbind(ontogenetic.angle=as.data.frame(vector.angle),ontogenetic.size=as.data.frame(vector.length),as.data.frame(x.and.y))->onto.vector
  abs(apply(onto.vector,1,diff))[1:2]->onto.real
  combn(colnames(onto.vector),2)->couples

  onto.real<-matrix(ncol=2,nrow=dim(couples)[2])
  for (i in 1 :dim(couples)[2]){
    apply(onto.vector[1:2,which(colnames(onto.vector)%in%couples[,i])],1,diff)->onto.real[i,]
  }
  rownames(onto.real)<-apply(couples,2,function(x) paste(x,collapse = "/"))
  colnames(onto.real)<-c("ontogenetic.angle.diff","ontogenetic.size.diff")

  ##### Randomization #####

  onto.random<-list()
  for(w in 1:1000)
  {
    data.frame(species=VS$species,VS[sample(as.numeric(rownames(VS))),2:5])->VS.r
    unique(unlist(lapply(strsplit(as.character(VS.r$species),split="_"),"[[",1)))->group
    gg<-list()
    for (i in 1:length(group)){
      VS.r[grep(group[i],VS.r$species),]->gg[[i]]
    }
    names(gg)<-group

    lapply(gg,function(x) apply(x[,4:5],2,function(z) z[1]-sum(z[2:length(z)])))->x.and.y
    lapply(x.and.y,function(x) rad2deg(atan(x[2]/x[1])))->vector.angle
    lapply(x.and.y,function(x) (x[1]^2+x[2]^2)^.5)->vector.length
    rbind(ontogenetic.angle=as.data.frame(vector.angle),ontogenetic.size=as.data.frame(vector.length),as.data.frame(x.and.y))->onto.vectorR


    onto.r<-matrix(ncol=2,nrow=dim(couples)[2])
    for (i in 1 :dim(couples)[2]){
      apply(onto.vectorR[1:2,which(colnames(onto.vectorR)%in%couples[,i])],1,diff)->onto.r[i,]
    }
    rownames(onto.r)<-apply(couples,2,function(x) paste(x,collapse = "/"))
    colnames(onto.r)<-c("ontogenetic.angle.diff","ontogenetic.size.diff")
    onto.r->onto.random[[w]]
  }

  p.ontogenetic.vector.angle<-array()
  p.ontogenetic.vector.size<-array()
  for(m in 1:dim(couples)[2]){
    length(which(unlist(lapply(onto.random,function(x) x[m,1]))>onto.real[m,1]))/1000->p.ontogenetic.vector.angle[m]
    length(which(unlist(lapply(onto.random,function(x) x[m,2]))>onto.real[m,2]))/1000->p.ontogenetic.vector.size[m]
  }
  names(p.ontogenetic.vector.angle)<-rownames(onto.real)
  names(p.ontogenetic.vector.size)<-rownames(onto.real)
  if(length(which(p.ontogenetic.vector.angle>0.950))>0) p.ontogenetic.vector.angle[which(p.ontogenetic.vector.angle>0.950)]<-1-p.ontogenetic.vector.angle[which(p.ontogenetic.vector.angle>0.950)]
  if(length(which(p.ontogenetic.vector.size>0.950))>0) p.ontogenetic.vector.size[which(p.ontogenetic.vector.size>0.950)]<-1-p.ontogenetic.vector.size[which(p.ontogenetic.vector.size>0.950)]


  ##### end Randomization #####

  cbind(onto.vector,"p"=rbind(p.ontogenetic.vector.angle,p.ontogenetic.vector.size,rep("",length(p.ontogenetic.vector.angle)),rep("",length(p.ontogenetic.vector.angle))))->onto.vector
  if(length(which(p.ontogenetic.vector.angle<0.05))>0) print(paste("ontogenetic vectors to MRCA differ in angle for",paste(names(p.ontogenetic.vector.angle[which(p.ontogenetic.vector.angle<0.05)]),collapse=" and ")))
  if(length(which(p.ontogenetic.vector.size<0.05))>0) print(paste("ontogenetic vectors to MRCA differ in size for",paste(names(p.ontogenetic.vector.size[which(p.ontogenetic.vector.size<0.05)]),collapse=" and ")))

  ### end Ontogenetic vectors to MRCA ###

  list("angles between species"=anglesBTW.mat,"angles between species 2 MRCA"=anglesBTW2MRCA.mat)->matt

  mats<-list()
  for(f in 1:length(matt)){
    matt[[f]]->aaa
    unique(unlist(lapply(strsplit(colnames(aaa),"_"),"[[",1)))->sp.col
    unique(unlist(lapply(strsplit(rownames(aaa),"_"),"[[",1)))->sp.row

    colmat<-list()
    for (i in 1:length(sp.col)){
      aaa[,grep(sp.col[i],colnames(aaa))]->colmat[[i]]
    }


    couple.mat<-list()
    for (i in 1:length(colmat)){
      colmat[[i]]->bbb
      rowmat<-list()
      for (k in 1:length(sp.row)){
        bbb[grep(sp.row[k],rownames(bbb)),]->rowmat[[k]]

      }
      for(w in 1:length(rowmat))
      {
        list.append(couple.mat,rowmat[[w]])->couple.mat
      }

    }
    couple.mat->mats[[f]]
  }

  names(matt)->names(mats)

  ### Ontogenetic vectors from 1st stage###

  mats$`angles between species`->matte
  t(sapply(matte,function(x) unique(sapply(strsplit(c(rownames(x),colnames(x)),"_"),"[[",1))))->pair.group


  ontogenesis<-list()
  ontogenesis.mats<-list()
  for(g in 1:dim(pair.group)[1]){
    pair.group[g,]->group

    group.mats<-list()
    for(i in 1:length(group)){
      retrieve.angles(evo.p,wishlist="angles.between.species",focus="species",species=group[i],random="no",write2csv = "no")->g1
      g1[grep(group[i],sapply(strsplit(as.character(g1[,2]),split="/"),"[[",1))[which(grep(group[i],sapply(strsplit(as.character(g1[,2]),split="/"),"[[",1))%in%grep(group[i],sapply(strsplit(as.character(g1[,2]),split="/"),"[[",2)))],]->g1
      data.frame(do.call(rbind,strsplit(as.character(g1[,2]), "/")),g1[,c(3,4)])->g1
      unique(c(as.character(g1[,1]),as.character(g1[,2])))[order(unique(c(as.character(g1[,1]),as.character(g1[,2]))))]->group1
      g1[,c(2,1,3,4)]->g1.inv
      colnames(g1.inv)<-colnames(g1)
      rbind(g1,g1.inv)->g11
      xtabs(as.numeric(as.character(anglesBTWspecies))~as.character(X1)+as.character(X2),data=g11)->group.mat
      names(dimnames(group.mat))<-NULL
      as.matrix(group.mat)->group.mat
      group.mat->group.mats[[i]]
    }

    group.matt<-list()
    which(sapply(strsplit(colnames(group.mats[[1]]),"_"),"[[",2)%in%sapply(strsplit(colnames(group.mats[[2]]),"_"),"[[",2))->in1
    which(sapply(strsplit(colnames(group.mats[[2]]),"_"),"[[",2)%in%sapply(strsplit(colnames(group.mats[[1]]),"_"),"[[",2))->in2
    group.mats[[1]][in1,in1]->s1
    group.mats[[2]][in2,in2]->s2
    s1->group.matt[[1]]
    s2->group.matt[[2]]
    s1[which(lower.tri(s1))]<-s2[which(lower.tri(s2))]
    rownames(s2)->rownames(s1)
    diag(s1) <- rep("", length(diag(s1)))
    as.data.frame(rbind(s1,rownames(s1)))->s3
    data.frame(s3,c(colnames(s3),""))->s3
    colnames(s3)[length(colnames(s1))+1]<-""

    s3->ontogenesis.mats[[g]]

    PCgroup<-list()
    ontog<-list()
    for (i in 1:length(group.matt)){
      group.matt[[i]]->group.mat
      y.onto[match(rownames(group.mat),rownames(y.onto)),]->groupY
      groupY->PCgroup[[i]]
      as.matrix(groupY)->groupY
      apply(groupY,1,unitV)->group.len
      apply(groupY,2,function(x) x[1]-sum(x[2:length(x)]))->Y.resultant
      rad2deg(acos((groupY[1,]%*%Y.resultant)/(unitV(groupY[1,])*unitV(Y.resultant))))->onto.angle
      unitV(Y.resultant)->onto.vector.size
      cbind(onto.angle,onto.vector.size)->onto.tot
      colnames(onto.tot)<-c("angle","size")
      rownames(onto.tot)<-group[i]
      onto.tot->ontog[[i]]
    }
    do.call(rbind,ontog)->ontog


    #### Randomization ####
    onto1.ran<-matrix(ncol=2,nrow=100)
    for(h in 1:100){
      seq(1,dim(PCgroup[[1]])[1],1)->ss
      sample(c(rep(1,length(ss)),rep(2,length(ss))),length(ss))->sam
      while(sum(sam)==length(ss)|sum(sam)==2*length(ss)) sample(c(rep(1,length(ss)),rep(2,length(ss))),length(ss))->sam
      sam2<-array()
      for(k in 1:length(sam)){ if (sam[k]==2) sam2[k]<-1 else sam2[k]<-2}
      c(sam,sam2)->sam3
      c(ss,ss)->ss
      matss<-matrix(ncol=dim(PCgroup[[1]])[2],nrow = length(ss))
      for(i in 1:length(sam3)){
        PCgroup[[sam3[i]]][ss[i],]->tem
        as.matrix(tem)->tem
        tem->matss[i,]
      }
      unlist(lapply(PCgroup,rownames))->rownames(matss)
      unlist(lapply(PCgroup,colnames)[[1]])->colnames(matss)

      l<-list()
      matss[1:(length(ss)/2),]->l[[1]]
      matss[((length(ss)/2)+1):length(ss),]->l[[2]]

      matrix(nrow=2,ncol=2)->ontog.ran
      for (i in 1:2){
        l[[i]]->groupY
        as.matrix(groupY)->groupY
        apply(groupY,1,unitV)->group.len
        apply(groupY,2,function(x) x[1]-sum(x[2:length(x)]))->Y.resultant
        rad2deg(acos((groupY[1,]%*%Y.resultant)/(unitV(groupY[1,])*unitV(Y.resultant))))->onto.angle
        unitV(Y.resultant)->onto.vector.size
        cbind(onto.angle,onto.vector.size)->onto.tot
        colnames(onto.tot)<-c("angle","size")
        rownames(onto.tot)<-group[i]
        onto.tot->ontog.ran[i,]
      }
      apply(ontog.ran,2,diff)->onto1.ran[h,]
    }
    colnames(onto1.ran)<-c("angleDiff","sizeDiff")
    apply(ontog,2,diff)->ontodiffs
    length(which(onto1.ran[,1]>ontodiffs[1]))/100->p.onto1.angle
    length(which(onto1.ran[,2]>ontodiffs[2]))/100->p.onto1.size
    #### end Randomization ####
    if(p.onto1.angle<0.05|p.onto1.angle>0.95) print(paste("ontogenetic vectors to 1st stage differ in angle for",paste(rownames(ontog),collapse=" and ")))
    if(p.onto1.size<0.05|p.onto1.size>0.95) print(paste("ontogenetic vectors to 1st stage differ in size for",paste(rownames(ontog),collapse=" and ")))


    rbind(ontog,p=c(p.onto1.angle,p.onto1.size))->ontogenesis[[g]]

  }
  apply(pair.group,1, function(x) paste(x[1],x[2],sep="/"))->names(ontogenesis.mats)
  apply(pair.group,1, function(x) paste(x[1],x[2],sep="/"))->names(ontogenesis)
  ontogenesis.list<-list()
  ontogenesis.mats->ontogenesis.list[[1]]
  ontogenesis->ontogenesis.list[[2]]
  names(ontogenesis.list)<-c("matrices","vectors")

  ### End Ontogenetic vectors from 1st stage ###

  SMAT.res<-list()
  for (h in 1:length(mats))
  {
    mats[[h]]->matta
    smat.res<-list()
    nam<-array()
    for (i in 1:length(matta)){
      matta[[i]]->mat
      which(is.na(match(unlist(lapply(strsplit(colnames(mat),"_"),"[[",2)),unlist(lapply(strsplit(rownames(mat),"_"),"[[",2)))))->colout
      which(is.na(match(unlist(lapply(strsplit(rownames(mat),"_"),"[[",2)),unlist(lapply(strsplit(colnames(mat),"_"),"[[",2)))))->rowout

      if (length(colout)>0) mat[,-colout]->mat1 else mat->mat1
      if (length(rowout)>0) mat1[-rowout,]->mat1

      if(dim(mat1)[1]<3) {
        list(matrix=mat,paedomorphosis.test=NULL,biogenetic.test=NULL)->smat.res[[i]]
        c(colnames(mat),rownames(mat))->totnam
        unlist(lapply(strsplit(totnam,split="_"),"[[",1))->group
        paste(unique(group),collapse="/")->nam[i]
        warning("too few ontogenetic stages for meaningful test for ",paste(nam[i]))
      } else {

        diag(mat1)->biogen.vector

        as.data.frame(c(mat1[1,],mat1[,1]))->mat2

        unlist(lapply(strsplit(rownames(mat2),split="_"),"[[",1))->group
        unlist(lapply(strsplit(rownames(mat2),split="_"),"[[",2))->age
        data.frame(group=group,age=age,angle=mat2[,1])->mat3
        as.numeric(as.character(mat3$age))->mat3$age

        if(class(try(sma(angle~age*group,data=mat3)->res.slope,silent=TRUE))=="try-error") {

          summary(lm(angle~group/age-1,data=mat3))->a
          t(a$coef[c(3,4),c(1,4)])->a1
          confint(lm(angle~group/age-1,data=mat3))->a2
          t(a2[3:4,])->a2
          summary(lm(angle~age,data=subset(mat3,mat3$group==levels(mat3$group)[1])))$r.squared->r1
          summary(lm(angle~age,data=subset(mat3,mat3$group==levels(mat3$group)[2])))$r.squared->r2
          c(r1,r2)->rr
          names(rr)<-colnames(a1)
          rbind(a1[1,],a2,rr,a1[2,])->A
          rownames(A)<-c("slope","lower.CI.lim","upper.CI.lim","R2","p")
          unique(group)->colnames(A)

          data.frame(LogL=logLik(lm(angle~group/age-1,data=mat3))[1],p=summary(lm(angle~group*age,data=mat3))$coef[4,4])->B
          list("Results of comparing lines among groups by lm"=B,"Coefficients by group in variable 'group' by lm"=A)->sma.res

        } else {
          res.slope=res.slope
          res.slope$commoncoef[c(1,2,7)]->a
          rbind(a[[3]],rbind(t(res.slope$r2),t(res.slope$pval)))->b
          rownames(b)[4:5]<-c("R2","p")
          list(as.data.frame(a[1:2]),b)->sma.res
          names(sma.res)<-c("Results of comparing lines among groups","Coefficients by group in variable 'group'")
        }
        data.frame(age=as.numeric(unique(age)),biogen.vector)->bg.data
        sma(biogen.vector~age,data=bg.data)->bg.res
        list(as.data.frame(bg.res$coef),data.frame(p=unlist(bg.res$p),R2=unlist(bg.res$r2)))->bg.test
        names(bg.test)<-c("biogenetic law test coefficients","biogenetic law test significance")
        list(matrix=mat,paedomorphosis.test=sma.res,biogenetic.test=bg.test)->smat.res[[i]]
        paste(unique(group),collapse="/")->nam[i]
      }
    }

    names(smat.res)<-nam
    if(h==1){
      for (i in 1:length(smat.res)){
        if(length(smat.res[[i]]$paedomorphosis.test)>0){
          if(as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][1,1])*as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][1,2])<0&smat.res[[i]]$paedomorphosis.test[[1]]$p<.05 |     length(which(as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][5,])>.05))==1 &smat.res[[i]]$paedomorphosis.test[[1]]$p<.05){
            print(paste("It looks there is paedomorphosis between",nam[i]))
          }
          if(smat.res[[i]]$biogenetic.test$`biogenetic law test coefficients`$coef.SMA.[2]>0&smat.res[[i]]$biogenetic.test$`biogenetic law test significance`$p<.05){
            print(paste("Biogenetic law is confirmed for",nam[i]))
          }
        }
      }
    }else{
      for (i in 1:length(smat.res)){
        if(length(smat.res[[i]]$paedomorphosis.test)>0){
          if(as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][1,1])*as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][1,2])<0&smat.res[[i]]$paedomorphosis.test[[1]]$p<.05 |     length(which(as.data.frame(smat.res[[i]]$paedomorphosis.test[[2]][5,])>.05))==1 &smat.res[[i]]$paedomorphosis.test[[1]]$p<.05){
            print(paste("It looks there is paedomorphosis between",nam[i],"2 MRCA"))
          }
          if(smat.res[[i]]$biogenetic.test$`biogenetic law test coefficients`$coef.SMA.[2]>0&smat.res[[i]]$biogenetic.test$`biogenetic law test significance`$p<.05){
            print(paste("Biogenetic law is confirmed for",nam[i],"2 MRCA"))
          }
        }
      }
    }
    smat.res->SMAT.res[[h]]
  }
  names(SMAT.res)<-names(matt)

  return(list(regression.matrix=SMAT.res,angles.2.MRCA.and.vector.size=anglesMRCA,ontogenetic.vectors2MRCA=onto.vector,ontogenetic.vectors.to.1st.stage=ontogenesis.list))

}
