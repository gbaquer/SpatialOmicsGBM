#'@author Gerard Baquer, \email{gbaquer@@bwh.harvard}, \email{gerard.baquer@@alumni.urv.cat},\email{baquer.gomez@@gmail.com}
#'@keywords Mass Spectrometry Imaging, Ion Mobility, Tandem Mass Spectrometry, Bioinformatics, Cheminformatics, Image Registration, Data Fusion, Quantification, Metabolomics, Lipidomics, Proteomics.


# REQUIRED LIBRARIES
#'@import "XML"
#'@import "sp"
#'@import "Rcpp"
#'@import "ggplot2"
#'@import "ggrepel"
#'@import "jsonlite"
#'@import "raster"
#'@import "magrittr"
#'@import "tictoc"
#'@import "grDevices"
#'@import "mzR"
#'@import "DBI"
#'@import "dplyr"

pkg.env <- new.env()


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
rearrange <- function(pks,v){
  i<-order(v)
  tmp<-utils.subsetPks(pks,i)
  return(tmp)
}


row.match <-  function(a, b, nomatch=NA){
  ca<-paste(a[,1],a[,2])
  cb<-paste(b[,1],b[,2])
  match(ca,cb,nomatch=nomatch)
}
na.rm<-function(x){
  x[!is.na(x)]
}

# Load imzml & regions (rsd)
forceUUID <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ))
{
  xmlRes <- rMSI:::CimzMLParse(path.expand(imzML_File))
  xmlRes$continuous_mode=F;
  bincon <- file(description = path.expand(ibd_File), open = "rb")
  binUUID <- paste(sprintf("%.2X", readBin(bincon, integer(), 16, size = 1, signed = F)), collapse = "")
  close(bincon)
  xmlRes$UUID=binUUID;
  rMSI:::CimzMLStore(path.expand(imzML_File),xmlRes)
}

#' Load Raw Centroid data
#'
#' @param files (.imzml, the same folder should contain the .ibd and .rsd counterparts)
#' @param ppm (binning tolerance. Dedfaults to 10 ppm)
#' @param a (binning tolerance. Dedfaults to 10 ppm)
#'
#' @returns pks (output data)
#'
#' @export
load<-function(files,ppm=10,a=0.5,mzmin=0,mzmax=Inf,imin=0){
  #Load raw data
  tictoc::tic()
  tmp<-lapply(files,parse)
  data<-list()
  data$pos<-do.call(rbind,lapply(tmp,function(x)x$pos))
  data$peakList<-do.call(c,lapply(tmp,function(x)x$peakList))
  pks<-list()
  pks$pos<-data$pos
  tictoc::toc()
  filelengths<-sapply(tmp,function(x)nrow(x$pos))
  rm(tmp)
  gc()
  print(pryr::mem_used())

  tictoc::tic()
  m<-unlist(lapply(data$peakList,function(x)x$mass))
  i<-unlist(lapply(data$peakList,function(x)x$intensity))
  id<-rep(1:nrow(pks$pos),(sapply(data$peakList,function(x)length(x$mass))))
  is<-(m>=mzmin)&(m<=mzmax)&(i>=imin)
  m<-m[is]
  i<-i[is]
  id<-id[is]
  rm(is)
  gc()
  tictoc::toc()
  print(pryr::mem_used())


  tictoc::tic()
  f<-100000
  m<-as.integer(m*f)
  tictoc::toc()
  print(pryr::mem_used())

  tictoc::tic()
  rm(data)
  gc()
  d<-density(m,bw=0.0001,weights = i/sum(i),n = 2^24)
  b<-pracma::findpeaks(d$y,minpeakheight = median(d$y)+a*sd(d$y),nups=1,minpeakdistance = 1)
  pks$mz<-sort(d$x[b[,2]])/f
  pks$intensity<-matrix(0,nrow(pks$pos),length(pks$mz))
  tictoc::toc()
  rm(b,d)
  gc()
  print(pryr::mem_used())

  tictoc::tic()
  o<-sort(m,index.return=T,method="radix")
  tictoc::toc()
  print(pryr::mem_used())

  tictoc::tic()
  hi<-sapply(as.integer((pks$mz+pks$mz*ppm*1e-6)*f),function(x)Rfast::binary_search(o$x,x,T))
  hi[hi>length(m)]<-length(m)
  li<-sapply(as.integer((pks$mz-pks$mz*ppm*1e-6)*f),function(x)Rfast::binary_search(o$x,x,T))
  tictoc::toc()
  print(pryr::mem_used())

  tictoc::tic()
  ind1<-o$ix[unlist(lapply(1:length(pks$mz),function(x)li[x]:hi[x]))]
  ind2<-rep(1:length(pks$mz),hi-li+1)
  pks$intensity[id[ind1]+(ind2-1)*nrow(pks$pos)]<-i[ind1]
  tictoc::toc()
  rm(m,i,id,o,hi,li,ind1,ind2)
  gc()
  print(pryr::mem_used())

  pks$dfmz=data.frame(name=rep("",length(pks$mz)))
  pks$dfmz$mean<-apply(pks$intensity,2,mean)

  pks$df=data.frame(name=rep("",nrow(pks$pos)))
  pks$df$file<-rep(files,filelengths)
  pks$df$batch<-as.numeric(as.factor(pks$df$file))
  pks$df$tic<-apply(pks$intensity,1,sum)

  tryCatch({
    for(file in unique(pks$df$file)){
      r<-loadRegions(gsub(".imzML",".srd",file))
      is<-pks$df$file==file
      pks$df$name[is]<-mergePksRegions(utils.subsetPks(pks,is),r)
    }
  },
  error=function(e){})

  tryCatch({
    pks<-utils.arrange(pks,pks$df$batch)
  },
  error=function(e){})


  class(pks)<-"peakMatrix"
  rm(list=ls()[! ls() %in% c("pks")])
  gc()
  print(pryr::mem_used())

  return(pks)
}

saveimzml<-function(pks,filename){
  tmp<-list()
  tmp$pos<-pks$pos
  tmp$pixel_size_um<-100
  tmp$peakList<-list()
  for(i in 1:nrow(pks$intensity)){
    tmp$peakList[[i]]<-list()
    tmp$peakList[[i]]$mass<-pks$mz
    tmp$peakList[[i]]$intensity<-pks$intensity[i,]
  }
  rMSIproc::export_imzMLpeakList(tmp$peakList,tmp$pos,tmp$pixel_size_um,filename)
}


loadRegions<-function(fileName){
  regions<-jsonlite::fromJSON(txt=readChar(fileName, file.info(fileName)$size))
}


mergePksRegions<-function(pks,r,subset=rep(T,length(r$Regions$Name))){
  names<-r$Regions$Name
  run<-as.numeric(as.factor(sapply(r$Regions$Sources,function(x)x$Path[1])))
  children<-sapply(r$Regions$Sources,function(x)length(x$Path))==1
  pks$df=data.frame(name=rep("",nrow(pks$pos)))
  for (k in unique(run)){
    for (i in unique(run)){
      a<-do.call(rbind,lapply(which(run==i&children),function(j)r$Regions$Sources[[j]]$Spots[[1]]))
      a$region<-rep(r$Regions$Name[run==i&children],sapply(which(run==i&children),function(j)nrow(r$Regions$Sources[[j]]$Spots[[1]])))
      a$c<-as.numeric(as.factor(a$region))
      a.pos<-data.matrix(a[,c("X","Y")])
      b2<-t(t(a.pos)-getReferencePoint(a.pos)+getReferencePoint(pks$pos[pks$df$name=="",]))
      b2<-b2[a$region%in%names[subset],]
      b<-rMSIworkflows:::row.match(b2,pks$pos)
      if(sum(is.na(b))==0){
        pks$df$name[b]<-a$region[a$region%in%names[subset]]
        break;

      }
    }
  }
  return(pks$df$name)
}

## Tests
fc<-function(pks,fixed=i~ GBM,mag="intensity"){
  a<-pks$df[[fixed[[3]]]]
  #sapply(1:length(pks$mz),function(i)mean(pks[[mag]][(a==unique(a)[1]),i])/mean(pks[[mag]][(a==unique(a)[2]),i]))
  sapply(1:length(pks$mz),function(i)mean(sapply(unique(pks$df$id[(a==unique(a)[1])]),function(j)mean(pks[[mag]][pks$df$id==j,i])))/mean(sapply(unique(pks$df$id[(a==unique(a)[2])]),function(j)mean(pks[[mag]][pks$df$id==j,i]))))
}
#1. t.test
test.ttest<-function(pks,fixed=i~ GBM,mag="intensity",dropmodels=F){
  start<-Sys.time()

  a<-pks$df[[fixed[[3]]]]

  data<-list()
  data$mz<-pks$mz
  data$model<-lapply(1:length(pks$mz),function(i)t.test(pks[[mag]][(a==unique(a)[1]),i],pks[[mag]][(a==unique(a)[2]),i],alternative="two.sided"))
  data$pval<-sapply(data$model,function(x)x$p.value)
  data$fc<- fc(pks,fixed,mag)
  data$name<-"ttest"
  data$time<-Sys.time()-start
  if(dropmodels)
    data$model<-NULL

  return(data)
}

test.wilcox<-function(pks,fixed=i~ GBM,mag="intensity",dropmodels=F){
  start<-Sys.time()

  a<-pks$df[[fixed[[3]]]]

  data<-list()
  data$mz<-pks$mz
  data$model<-lapply(1:length(pks$mz),function(i)wilcox.test(pks[[mag]][(a==unique(a)[1]),i],pks[[mag]][(a==unique(a)[2]),i],alternative="two.sided"))
  data$pval<-sapply(data$model,function(x)x$p.value)
  data$fc<-fc(pks,fixed,mag)
  data$name<-"Wilcox rank sum test"
  data$time<-Sys.time()-start

  if(dropmodels)
    data$model<-NULL

  return(data)
}
test.multiple<-function(pks,test,fixed=i~ GBM,random=~1|id,mag="intensity",other=NULL,dedicated=F){
  # no_cores <- parallel::detectCores()
  # cl<-parallel::makeCluster(no_cores-dedicated*1)
  # doParallel::registerDoParallel(cl)

  data<-list()
  data$mz<-pks$mz
  data$fc<-fc(pks,fixed,mag)
  data$model<-foreach::foreach(i=1:length(pks$mz))%dopar%
    {
      test(pks,i,fixed,random,mag,other)
    }
  # parallel::stopCluster(cl)
  gc()
  return(data)
}
test.single.nlme<-function(pks,i,fixed=i~ GBM,random=~1|id,mag="intensity",other=NULL){
  tryCatch(
    expr={
      d<-as.data.frame(pks$df)
      d$i<-pks$intensity[,i]

      if(is.null(other$pixels))
        return(nlme::lme(fixed=fixed,data=d,random=random,method="ML"))
      else{
        d<-cbind(pks$pos,d)
        n<-(d$x%%2==1)&(d$y%%2==1)
        d<-d[n,]
        return(nlme::lme(fixed=fixed,data=d,random=random,method="ML",correlation = nlme::Initialize(nlme::corLin(other$pixels,~x+y|id),d)))
      }

    },
    error=function(e){
      return(NULL)
    }
  )

}
test.nlme<-function(pks,fixed=i~ GBM,random=~1|id,mag="intensity",dropmodels=F,other=NULL){
  start<-Sys.time()
  data<-test.multiple(pks,test.single.nlme,fixed,random,mag,other=other)
  data$pval<-sapply(data$model,function(x)if(is.null(x)) NA else coef(summary(x))[2,"p-value"])
  data$AIC<-sapply(data$model,function(x)if(is.null(x)) NA else AIC(x))
  data$name<-if(is.null(other$pixels)) "LME (nlme)" else paste("LME (nlme)",other$pixels,"pixels")
  data$time<-Sys.time()-start

  if(dropmodels)
    data$model<-NULL

  return(data)
}

test.single.lmerTest<-function(pks,i,fixed=i~ GBM,random=~1|id,mag="intensity",other=NULL){
  tryCatch(
    expr={
      d<-as.data.frame(pks$df)
      d$i<-pks$intensity[,i]
      return(lmerTest::lmer(i ~ GBM+(1|id), data=d))
    },
    error=function(e){
      return(NULL)
    }
  )

}
test.lmerTest<-function(pks,fixed=i~ GBM,random=~1|id,mag="intensity",dropmodels=F,other=NULL){
  start<-Sys.time()
  data<-test.multiple(pks,test.single.lmerTest,fixed,random,mag,other=other)
  data$pval<-sapply(data$model,function(x)if(is.null(x)) NA else coef(summary(X))[2,5])
  data$name<-if(is.null(other$pixels)) "LME (lmerTest)" else paste("LME (lmerTest)",other$pixels,"pixels")
  data$time<-Sys.time()-start

  if(dropmodels)
    data$model<-NULL

  return(data)
}

test.Cardinal<-function(pks,fixed=i~ GBM,random=~1|id,mag="intensity",dropmodels=F,other=NULL){
  start<-Sys.time()
  card<-pks2Cardinal(pks)

  data<-list()
  data$mz<-pks$mz

  if(other$mode=="means"|is.null(other$mode)){
    data$model<-Cardinal::meansTest(card,fixed=~ GBM)
    data$pval<-summary(data$model)[,3]
  }
  else if(other$mode=="segmentation"){
    dgmm2 <- Cardinal::spatialDGMM(card, r=1, k=5, groups=Cardinal::run(card))
    data$model <- Cardinal::segmentationTest(dgmm2, ~ GBM, classControl="Ymax")
    data$pval<-summary(stest)[,'PValue']
  }
  else{
    data$model<-NULL
    data$pval<-NULL
  }

  data$fc<-fc(pks,fixed,mag)
  data$name<-paste("LME (Cardinal ",other$mode,")",sep="")
  data$time<-Sys.time()-start

  if(dropmodels)
    data$model<-NULL

  return(data)
}

plotNames<-function(pks,v=pks$df$name){
  df<-cbind(pks$df,pks$pos)
  df$v<-v
  label=unique(v)
  x.t=sapply(label,function(n)mean(pks$pos[v==n,1]))
  y.t=sapply(label,function(n)mean(pks$pos[v==n,2]))
  ggplot(df,aes(x=x,y=-y,fill=as.factor(as.numeric(as.factor(v)))))+geom_tile()+coord_fixed(ratio=1)+theme_bw()+
    annotate(geom="text",x=x.t,y=-y.t,label=label,size=2)+theme_void()+
    theme(legend.position="none")
}


## UTILS MODULE ##
utils.pad<-function(a,p=c(10,10),v=NA){
  b<-array(v,dim=dim(a)+2*p)
  b[p[1]+(1:dim(a)[1]),p[2]+(1:dim(a)[2])]<-a
  return(b)
}
utils.unpad<-function(a,p=c(10,10)){
  if(length(dim(a))==3)
    return(a[p[1]+(1:(dim(a)[1]-2*p[1])),p[2]+(1:(dim(a)[2]-2*p[2])),1])
  else
    return(a[p[1]+(1:(dim(a)[1]-2*p[1])),p[2]+(1:(dim(a)[2]-2*p[2]))])
}
utils.mergePks<-function(pks_list){
  pks<-pks_list[[1]]
  pks$pos<-unique(do.call(rbind,lapply(pks_list,function(x)x$pos)))
  pks$mz<-unique(unlist(lapply(pks_list,function(x)x$mz)))
  pks$intensity<-matrix(0,nrow(pks$pos),length(pks$mz))
  for(x in pks_list){
    i<-row.match(x$pos,pks$pos)
    j<-match(x$mz,pks$mz)
    pks$intensity[i,j]<-x$intensity
  }
  return(pks)

}
utils.subsetPks<-function(pks,i=rep(T,nrow(pks$pos)),j=rep(T,length(pks$mz))){
  # ISSUE: It messes up the pks$dfmz data frame.
  pks$mz<-pks$mz[j]
  pks$intensity<-pks$intensity[,j]
  for(x in c("pos","intensity","df")){
    pks[[x]]<-pks[[x]][i,]
  }
  for(x in c("dfmz")){
    pks[[x]]<-pks[[x]][j,]
  }
  if(!is.null(pks$cells)){
    is<-pks$cells$pixel%in%which(i)
    pks$cells<-pks$cells[is,]
  }
  if(!is.null(pks$pathway)){
    pks$pathway$enrichment<-pks$pathway$enrichment[i,]
    pks$pathway$p<-pks$pathway$p[i,]
  }
  return(pks)
}
#' @export
utils.geti<-function(pks,mz){
  which.min(abs(pks$mz-mz))
}
#' @export
utils.rotate<-function(p,a){
  newp<-p
  newp[,1]<-p[,1]*cos(a)-p[,2]*sin(a)
  newp[,2]<-p[,1]*sin(a)+p[,2]*cos(a)
  return(newp)
}
#' @export
utils.flip<-function(p,h=T,v=T){
  p[,1]<-((!h)*p[,1])+h*(1+abs(p[,1]-max(p[,1])))
  p[,2]<-((!v)*p[,2])+v*(1+abs(p[,2]-max(p[,2])))
  return(p)
}
#' @export
utils.morans<-function(pks){
  sapply(1:length(pks$mz),function(i)tryCatch(expr = {with(pks,moranfast::moranfast(intensity[,i],pos[,1],pos[,2])$observed)},error=function(e){return(NA)}))
}
#' @export
utils.coolColors<-function(){
  cols<-c("#8DD3C4","#FDB462","#BEBADA","#FB8072","#B3DE69","#FAA9D2","#80B1D3","#FFFF69")
  return(rep(cols,10))
}
#' @export
utils.getNormalization<-function(v){
  m<-mean(v)
  v[v==0]<-min(v[v!=0])
  return(m/v)
}

## EXPLORATORY MODULE ##
#' @export
exp.UMAP<-function(pks,name="",is=1:nrow(pks$pos),js=1:length(pks$mz),n=2){
  tictoc::tic()
  tmp<-apply(pks$intensity[is,js],2,rMSIworkflows:::range01)
  tmp[is.na(tmp)]<-0
  umap<-uwot::umap(tmp,n_components = n)
  pks$df[is,paste(paste("umap",1:n,sep=""),name,sep=".")]<-umap
  tictoc::toc()
  return(pks)
}
#' @export
exp.DBSCAN<-function(pks,name="",eps=0.2,minPts=5){
  tictoc::tic()
  df<-pks$df
  df$umap1<-pks$df[,paste("umap1.",name,sep="")]
  df$umap2<-pks$df[,paste("umap2.",name,sep="")]
  is<-!is.na(df$umap1)
  df<-df[is,]
  pks$df[is,paste(c("c"),name,sep=".")]<-as.factor(dbscan::dbscan(cbind(df$umap1,df$umap2),eps=eps,minPts = minPts)$cluster)
  tictoc::toc()
  return(pks)
}
exp.plotUMAP<-function(pks,name="",n=NULL,v=pks$df[,paste("c.",name,sep="")],t="Cluster ID",a=0.3){
  palette <- colorRampPalette(RColorBrewer::brewer.pal(9,"RdYlBu"))
  df<-pks$df
  umaps<-names(df)[grepl(sprintf("(umap[0-9]+\\.%s$)",name),names(df),perl = T)]
  umaps<-umaps[order(umaps)]
  umaps.short<-sapply(strsplit(umaps,".",fixed=T),function(x)x[1])
  for(i in 1:length(umaps))
    df[[umaps.short[i]]]<-df[[umaps[i]]]
  df$c<-v

  if(is.null(n)){
    df<-df[sample(1:nrow(df)),]
    return(ggplot(df,aes(x=umap1,y=umap2,color=c))+geom_point(alpha=a,stroke=0,size=2)+theme_classic()+#scale_color_distiller(name="Cluster ID",palette="RdYlBu")
      scale_color_manual("Cluster ID",values=palette(length(unique(df$c))))+
      xlab("UMAP 1")+ylab("UMAP 2")+
      guides(colour = guide_legend(title=t,override.aes = list(alpha = 1)))+
      theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "transparent"),legend.background =element_rect(fill = "transparent"), text = element_text(color="white"),axis.text = element_text(color="white"),axis.line = element_line(color="white"),axis.ticks = element_line(color="white"))+
      coord_fixed())
  }
  else{
    df<-cbind(df,pks$pos)
    if(n==-1){
      #INCLUDE
    }
    else{
      df$umap<-df[[paste("umap",n,sep="")]]
      p<-ggplot(df,aes(x=x,y=-y,fill=umap))+geom_raster()+coord_fixed()+ scale_fill_distiller(name=paste("UMAP",n),palette="RdYlBu")+theme_void()+coord_fixed()+
        theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "transparent"),legend.background =element_rect(fill = "transparent"), text = element_text(color="white"))
    }
    return(p)

  }
}
exp.plotDBSCAN<-function(pks,name="",annotate=F){
  palette <- colorRampPalette(RColorBrewer::brewer.pal(9,"RdYlBu"))
  df<-pks$df
  df$c<-pks$df[,paste("c.",name,sep="")]
  df<-cbind(df,pks$pos)
  df2<-data.frame(c=unique(df$c))
  df2$x<-sapply(df2$c,function(i)median(df$x[df$c==i]))
  df2$y<-sapply(df2$c,function(i)median(df$y[df$c==i]))
  p<-ggplot(df,aes(x=x,y=-y,fill=c))+geom_raster()+coord_fixed()+ scale_fill_manual("Cluster ID",values=palette(length(unique(df$c))))+theme_void()+coord_fixed()+
    theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "transparent"),legend.background =element_rect(fill = "transparent"), text = element_text(color="white"))
  if(annotate)
    p<-p+annotate("text",x=df2$x,y=-df2$y,label=df2$c)
  return(p)
}
exp.loadAnnotations<-function(pks,file,ppm=3,type="metaspace"){
  if(type=="metaspace"){
    ann<-read.csv2(file,comment.char = "#",sep=",")
    ppms<-outer(pks$mz,as.numeric(ann$mz),function(x,y)1e6*abs(x-y)/max(x,y))
    hits<-which(ppms<ppm,arr.ind = T)

    pks$dfmz<-data.frame(molecule=rep(NA,length(pks$mz)),id=rep(NA,length(pks$mz)))
    pks$dfmz$molecule[hits[,1]]<-sapply(strsplit(ann$moleculeNames[hits[,2]],", ",fixed=T),function(x)x[1])
    pks$dfmz$id[hits[,1]]<-sapply(strsplit(ann$moleculeIds[hits[,2]],", ",fixed=T),function(x)x[1])
    for(n in c("ion","adduct","formula","isomers","isobars","msm","fdr","rhoSpectral","rhoSpatial","moleculeNames","moleculeIds")){
      if(n %in% names(ann)){
        pks$dfmz[[n]]<-NA
        pks$dfmz[[n]][hits[,1]]<-ann[[n]][hits[,2]]
      }
    }
    return(pks)
  }
  else{
    if(type=="metaboscape"){
      ann<-read.csv2(file,comment.char = "#",sep=",")
      ppms<-outer(pks$mz,as.numeric(ann$m.z.meas.),function(x,y)1e6*abs(x-y)/max(x,y))
      hits<-which(ppms<ppm,arr.ind = T)

      pks$dfmz<-data.frame(molecule=rep(NA,length(pks$mz)),id=rep(NA,length(pks$mz)))
      pks$dfmz$molecule[hits[,1]]<-sapply(strsplit(ann$Name[hits[,2]],",",fixed=T),function(x)x[1])
      pks$dfmz$id[hits[,1]]<-sapply(strsplit(ann$Name[hits[,2]],",",fixed=T),function(x)x[2])
      pks$dfmz$formula<-NA
      pks$dfmz$formula[hits[,1]]<-ann$Molecular.Formula[hits[,2]]
      pks$dfmz$adduct<-NA
      pks$dfmz$adduct[hits[,1]]<-ann$Ions[hits[,2]]

      return(pks)
    }
    else{
      stop("Unknown type")
    }
  }

}

plot.ion<-function(pks,v,name="",q=0.99){
  max=quantile(v,q,na.rm = T)
  v[v>max]<-max
  df<-as.data.frame(pks$pos)
  df$v<-v
  return(ggplot(df,aes(x=x,y=-y,fill=v))+geom_raster()+coord_fixed()+scale_fill_gradientn(name=name,colours = viridisLite::viridis(100),na.value="transparent")+theme_void()+coord_fixed()+
           theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "transparent"),legend.background =element_rect(fill = "transparent"), text = element_text(color="white"))+ labs(fill='') )

}
plot.contour<-function(pks,v,z= pks$df$z,hi=quantile(v,0.99,na.rm=T)){
  if(is.null(z))
    z=rep(1,nrow(pks$pos))
  v[v>hi]<-hi
  a<-reg.get2Dvolume(pks,v,NA)
  b<-reg.get2Dvolume(pks,z,0)
  df<-as.data.frame(expand.grid(x=a$x,y=b$y))
  complete<-as.data.frame(expand.grid(x=(min(a$x)-5):(max(a$x)+5),y=(min(a$y)-5):(max(a$y)+5)))
  is<-row.match(complete,df)
  complete<-complete[is.na(is),]
  complete$v<-NA
  complete$z<-0
  df$v<-c(a$a)
  df$z<-c(b$a)
  df<-rbind(df,complete)
  return(ggplot(df,aes(x=x,y=-y,fill=v,z=z))+geom_raster()+coord_fixed()+scale_fill_gradientn(colours = viridisLite::viridis(100),limits=c(0,hi),na.value="transparent")+geom_contour(col="white",size=0.1)+theme_void()+coord_fixed()+
           theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "transparent"),legend.background =element_rect(fill = "transparent"), text = element_text(color="white"))+ labs(fill=''))

}

#' Get 2D volume
#' Return a 2D array for a given vector. The 3D volume is smoothed out using a moving average window of size n to improve 2D visualization.
#'
#' @param pks peak matrix
#' @param v vector of intensities to plot in the volume (pks$intensity[,*])
#' @param f number of pixels between sections (distance between sections / step size)
#' @param n kernel size for the moving averaging window (smooth out features for better 3D visualization)
#'
#' @returns list containing x,y,z coordinates and the 3D array (a)
#'
#' @export
reg.get2Dvolume<-function(pks,v,d=0){
  #Arrange array
  x<-sort(unique(pks$pos[,1]))
  y<-sort(unique(pks$pos[,2]))

  i<-match(pks$pos[,1],x)
  j<-match(pks$pos[,2],y)

  a<-rep(d,length=length(x)*length(y))
  a[i+j*length(x)]<-v
  a<-array(a,dim=c(length(x), length(y)))
  na<-is.na(a)
  a[na]<-d

  return(list(x=x,y=y,a=a))
}

utils.blur<-function(a,n=3){
  r <- raster::raster(as.matrix(a))
  return(raster::as.matrix(raster::focal(r, matrix(1, n, n), mean, pad = T, padValue = 0)))
}
utils.gaussianblur<-function(a,s=2,n=3){
  r <- raster::raster(as.matrix(a))
  return(raster::as.matrix(raster::focal(r, w = utils.gaussian.kernel(sigma=s, n=n), fun = mean, na.rm=TRUE, pad=FALSE)))
}
utils.blurvec<-function(pks,v,n=3,st=NA){
  img<-reg.get2Dvolume(pks,v,st)
  a<-utils.blur(img$a,n)
  b<-reg.set2Dvolume(img$x,img$y,a)
  js<-row.match(b$pos,pks$pos)
  kp<-!is.na(js)
  v[js[kp]]<-b$v[kp]
  return(v)
}
utils.gaussianblurvec<-function(pks,v,n=3,st=NA){
  img<-reg.get2Dvolume(pks,v,st)
  a<-utils.gaussianblur(img$a,n)
  b<-reg.set2Dvolume(img$x,img$y,a)
  js<-row.match(b$pos,pks$pos)
  kp<-!is.na(js)
  v[js[kp]]<-b$v[kp]
  return(v)
}
utils.gaussian.kernel <- function(sigma=2, n=5)
{
  m <- matrix(ncol=n, nrow=n)
  mcol <- rep(1:n, n)
  mrow <- rep(1:n, each=n)
  x <- mcol - ceiling(n/2)
  y <- mrow - ceiling(n/2)
  m[cbind(mrow, mcol)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
  m / sum(m)
}

#' @export
reg.exportimages<-function(pks,v,folder,n=NULL,individual=T,padValue = 0){
  dir.create(folder,recursive=T)
  if(individual){
    for(n in unique(pks$df$name)){
      is<-pks$df$name==n
      pks.small<-utils.subsetPks(pks,is)
      img<-reg.get2Dvolume(pks.small,v[is],padValue)$a
      #img<-img-min(img)
      img<-utils.pad(0.9*img/max(img),v = 0.9*padValue/max(img))
      png::writePNG(t(img),paste(folder,n,".png",sep=""))
    }
  }
  else{
    img<-reg.get2Dvolume(pks,v,padValue)$a
    #img<-img-min(img)
    img<-utils.pad(0.9*img/max(img),v = 0.9*padValue/max(img))
    png::writePNG(t(img),paste(folder,n,".png",sep=""))
  }
}
#' @export
reg.importimages<-function(pks,name,folder){
  pks$df[[name]]<-rep(NA,nrow(pks$pos))
  f<-list.files(folder,full.names = T)
  ns<-sapply(strsplit(basename(f),".",fixed=T),function(x)x[1])
  for(i in 1:length(ns)){
    is<-pks$df$name==ns[i]
    #Read
    img<-png::readPNG(f[i])
    img[img>0.95]<-0
    img<-utils.unpad(img)
    #Sort out indexing
    x<-sort(unique(pks$pos[is,1]))
    y<-sort(unique(pks$pos[is,2]))
    vol<-reg.set2Dvolume(x,y,t(img))
    js<-row.match(vol$pos,pks$pos)
    kp<-!is.na(js)
    pks$df[[name]][js[kp]]<-vol$v[kp]
  }
  return(pks)
}

#' Get 2D volume
#' Return a 2D array for a given vector. The 3D volume is smoothed out using a moving average window of size n to improve 2D visualization.
#'
#' @param pks peak matrix
#' @param v vector of intensities to plot in the volume (pks$intensity[,*])
#' @param f number of pixels between sections (distance between sections / step size)
#' @param n kernel size for the moving averaging window (smooth out features for better 3D visualization)
#'
#' @returns list containing x,y,z coordinates and the 3D array (a)
#'
#' @export
reg.set2Dvolume<-function(x,y,a){
  return(list(pos=expand.grid(x=x,y=y),v=c(a)))
}

utils.rotateflip<-function(pks,rot=rep(0,nrow(pks$pos)),hflip=rep(F,nrow(pks$pos)),vflip=rep(F,nrow(pks$pos))){
  for(n in unique(pks$df$name)){
    is<-which(pks$df$name == n)
    pks$pos[is,]=utils.flip(pks$pos[is,],hflip[is][1],vflip[is][1])
    pks$pos[is,]=utils.rotate(pks$pos[is,],rot[is][1])
  }
  return(pks)
}
utils.rotate<-function(p,a){
  newp<-p
  newp[,1]<-p[,1]*cos(a)-p[,2]*sin(a)
  newp[,2]<-p[,1]*sin(a)+p[,2]*cos(a)
  return(newp)
}
#' @export
utils.flip<-function(p,h=T,v=T){
  p[,1]<-((!h)*p[,1])+h*(min(p[,1])+abs(p[,1]-max(p[,1])))
  p[,2]<-((!v)*p[,2])+v*(min(p[,2])+abs(p[,2]-max(p[,2])))
  return(p)
}

pathway.initializeEnrichment<-function(){
  s.file="kegg/metpa/hsa.qs"
  l.file=tempfile()
  lib.url <- paste0("https://www.metaboanalyst.ca/resources/libs/", s.file);
  download.file(lib.url, destfile=l.file, method="curl")
  pkg.env$lib <- qs::qread(l.file)
}
pathway.completeEnrichment<-function(pks,q=0.5,is=1:nrow(pks$pos),type="set"){
  if(is.null(pkg.env$lib))
    pathway.initializeEnrichment()
  js<-which(!is.na(pks$dfmz$moleculeIds))
  ids<-pks$dfmz$moleculeIds[js]
  qtop<-apply(pks$intensity[is,js],2,quantile,q)
  enrichment<-lapply(1:nrow(pks$pos),function(i) pathway.enrich(pks$intensity[i,js],qtop,ids,type))
  pks$pathway$enrichment<-do.call(rbind,lapply(enrichment,function(x)x[,"Enrichment"]))
  pks$pathway$p<-do.call(rbind,lapply(enrichment,function(x)x[,"Raw p"]))
  return(pks)
}
pathway.enrich<-function(v,vref,ids,type="set"){
  #Prepare data
  ora.vec<-unique(unlist(strsplit(ids[v>vref],", ")))
  q.size <- length(ora.vec);
  current.mset <- pkg.env$lib$mset.list;
  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);
  set.size <- length(current.mset);
  uniq.count <- pkg.env$lib$uniq.count;
  imp.list <- pkg.env$lib$rbc.list;
  gd.sets <- names(current.mset);
  imp.list <- imp.list[gd.sets];
  hits<-lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num<-unlist(lapply(hits, function(x) length(x)), use.names = FALSE);

  #Compute results
  res.mat<-matrix(NA, nrow=set.size, ncol=7);
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR p","Enrichment");
  res.mat[,1]<-set.num;
  res.mat[,2]<-q.size*(set.num/uniq.count);
  res.mat[,3]<-hit.num;
  res.mat[,4]<-phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F) #GetFisherPvalue(hit.num, q.size, set.num, uniq.count);
  res.mat[,5] <-p.adjust(res.mat[,4], "holm");
  res.mat[,6] <-p.adjust(res.mat[,4], "fdr");
  if(type=="set")
    res.mat[,7] <-res.mat[,3]/res.mat[,2]
  else
    res.mat[,7] <- mapply(function(x, y){sum(x[y])}, imp.list, hits);
  row.names(res.mat)<-names(pkg.env$lib$path.ids)
  return(res.mat)
}

pathway.getCoverage<-function(pks,p){
  js<-which(!is.na(pks$dfmz$id))
  ids<-pks$dfmz$id[js]
  setids<-pkg.env$lib$mset.list[[match(p,names(pkg.env$lib$path.ids))]]
  return(list(matched=setids[setids%in%ids],notmatched=setids[!setids%in%ids]))
}

pathway.plotSpatialEnrichment<-function(pks,ref=as.numeirc(as.factor(pks$df$name)),n=ncol(pks$pathway$enrichment)){
  o<-apply(pks$pathway$enrichment,2,cor,ref,use="complete.obs")
  for(i in 1:n){
    x=sort(o,decreasing = T)[i]
    print(rMSIworkflows:::plot.arrange(pks,pks$pathway$enrichment[,names(x)])+ggtitle(paste0(names(x)," r=",round(x,2))))
  }
}

reg.importCells<-function(pks,path,name="all"){
  df<-read.csv2(path,comment.char = "#",sep=",")
  df$x<-as.numeric(df$x)
  df$y<-as.numeric(df$y)
  df$numcells<-ceiling(df$SizeInPixels/median(df$SizeInPixels))
  df$pixel<-rMSIworkflows:::row.match(cbind(round(df$x),round(df$y)),pks$pos)
  df<-df[!is.na(df$pixel),]
  s=sapply(df$pixel,function(x)sum(df$numcells[which(x==df$pixel)]))
  pks$df[[paste("cells.",name,sep="")]]=rep(0,nrow(pks$pos))
  pks$df[[paste("cells.",name,sep="")]][df$pixel]<-s
  pks$cells[[name]]<-df
  return(pks)
}

utils.updatecells<-function(pks){
  cids<-rMSIworkflows:::row.match(cbind(round(pks$merfish$xtrans),round(pks$merfish$ytrans)),pks$pos)
  c<-lapply(1:nrow(pks$pos),function(i)table(pks$merfish$leiden_class[which(cids==i)])[unique(pks$merfish$leiden_class)])
  df2<-do.call(rbind,c)
  colnames(df2)<-unique(pks$merfish$leiden_class)
  df2[is.na(df2)]<-0
  pks$cells<-df2
  return(pks)
}




utils.arrange<-function(pks,x=rep(1,nrow(pks$pos)),y=rep(1,nrow(pks$pos)),margin=5,id=pks$df$name,n=NULL,h=T,m=1){
  utils.grid(utils.group(pks,paste(x,y),n=n,h=h,margin=margin,m=m),x=x,y=y,margin=5*margin)
}
utils.group<-function(pks,v,id=pks$df$name,n=NULL,h=T,margin=5,m=NULL){
  groups<-lapply(unique(v),function(i)as.factor(unique(id[v==i])))
  dict<-unlist(lapply(groups,as.numeric))
  names(dict)<-unlist(lapply(groups,as.character))
  i<-dict[id]-1
  if(is.null(n)){
    n<-ceiling(sqrt(max(sapply(groups,length))))
  }
  if(is.null(m))
    m<-ceiling((max(sapply(groups,length)))/n)
  if(h){
    x<-1+i%%m
    y<-1+floor(i/m)
  }
  else {
    y<-1+i%%m
    x<-1+floor(i/m)
  }
  utils.grid(pks,x,y,margin,id)
}

utils.grid<-function(pks,x=rep(1,nrow(pks$pos)),y=rep(1,nrow(pks$pos)),margin=5,id=paste(x,y)){
  is<-(!is.na(x))&(!is.na(y))
  x<-x[is]
  y<-y[is]
  id<-id[is]
  pks<-utils.subsetPks(pks,is)

  xlabs<-factor(sort(unique(x)))
  ylabs<-factor(sort(unique(y)))

  grid<-data.frame(id=unique(id))
  centroids<-sapply(grid$id,function(n)apply(pks$pos[id==n,],2,function(v)round(min(v)+(max(v)-min(v))/2)))
  grid$x<-centroids[1,]
  grid$y<-centroids[2,]
  size<-sapply(grid$id,function(n)apply(pks$pos[id==n,],2,max)-apply(pks$pos[id==n,],2,min))
  grid$width<-size[1,]
  grid$height<-size[2,]
  grid$cell.x<-sapply(grid$id,function(n)match(x[match(n,id)],xlabs))*(max(grid$width)+margin)
  grid$cell.y<-sapply(grid$id,function(n)match(y[match(n,id)],ylabs))*(max(grid$height)+margin)
  grid$offset.x<-grid$cell.x-grid$x
  grid$offset.y<-grid$cell.y-grid$y

  pks$pos[,1]<-pks$pos[,1]+grid$offset.x[match(id,grid$id)]
  pks$pos[,2]<-pks$pos[,2]+grid$offset.y[match(id,grid$id)]

  pks$plot<-list()
  pks$plot$grid<-grid
  pks$plot$xlabs<-xlabs
  pks$plot$ylabs<-ylabs
  pks$plot$margin<-margin

  return(pks)
}

sc.import<-function(pks,path,exclude=c()){
  pks$cells<-NULL
  for(n in unique(pks$df$name)){
    f<-paste0(path,n,"_reg.csv")
    if(file.exists(f)){
      cells<-read.table(f,sep=",",header=T)
      cells$x<-(cells$x)+min(pks$pos[pks$df$name==n,1])
      cells$y<-(cells$y)+min(pks$pos[pks$df$name==n,2])
      pks$cells<-rbind(pks$cells,cells)
    }
  }
  return(utils.updatecellsdf(pks,exclude))
}
utils.updatecellsdf<-function(pks,exclude=c()){
  pks$cells$pixel<-rMSIworkflows:::row.match(cbind(round(pks$cells$x),round(pks$cells$y)),pks$pos)
  pks<-utils.summarizecolumn(pks,pks$cells$pixel,rep(T,length(pks$cells$pixel)),"c.density")
  ns<-colnames(pks$cells)
  ns<-ns[!ns%in%exclude]
  for(n in ns)
    pks<-utils.summarizecolumn(pks,pks$cells$pixel,pks$cells[[n]],paste0("c.",n))

  return(pks)
}
utils.summarizecolumn<-function(pks,p,v,n,f=mean){
  if(is.logical(v)){
    pks$df[[n]]=rep(0,nrow(pks$pos))
    w<-table(pks$cells$pixel[v])
    pks$df[[n]][as.numeric(names(w))]<-w
  }

  if(is.numeric(v)){
    pks$df[[n]]=rep(0,nrow(pks$pos))
    df<-data.frame(p=p,v=v)
    df<-df %>% group_by(p) %>%  summarise(f = f(v))
    df<-df[!is.na(df$p),]
    pks$df[[n]][df$p]<-df$f
  }

  if(is.character(v)||is.factor(v)){
    v=as.factor(v)
    for(m in unique(v)){
      pks$df[[paste0(n,".",m)]]=rep(0,nrow(pks$pos))
      w<-table(pks$cells$pixel[v==m])
      pks$df[[paste0(n,".",m)]][as.numeric(names(w))]=w
    }
  }
  return(pks)
}
