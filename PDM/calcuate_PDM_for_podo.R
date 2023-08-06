library(dplyr);library(tidyverse);library(MASS);library(Seurat);library(rrcov)

##functions for calculating KL
multivar_KLscore<-function(expmat1,expmat2){
  Mu1<-rowMeans(expmat1)
  Mu2<-rowMeans(expmat2)
  Sigma1<-cov(t(expmat1))
  Sigma2<-cov(t(expmat2))
  ln_ratio<-log((det(Sigma2)+10e-7)/(det(Sigma1)+10e-7))
  n<-length(Mu1)
  tr<-sum(diag(ginv(Sigma2)%*%Sigma1))
  d_z<-(t(Mu2-Mu1))%*%(ginv(Sigma2))%*%(Mu2-Mu1)
  kl=0.5 * (ln_ratio - n + tr + d_z)
  score=kl/log(ncol(expmat1)*ncol(expmat2))
  return(score)
}

pca_process<-function(expmat){
  x<-t(expmat)
  new<-prcomp(x)$x[,1:2]
  new<-t(new)
  return(new)
}

# load files
s0<-read.table('./podo_state0.txt',
               header = T,row.names = 1) # Ctrl
s1<-read.table('./podo_3state1.txt',
                          header = T,row.names = 1) # cell state 1
s2<-read.table('./podo_3state2.txt',
                          header = T,row.names = 1) # cell state 2
s3<-read.table('./podo_3state3.txt',
               header = T,row.names = 1) # cell state 3

#construct regulon
load(file='./MARAregulon.RData') 
tfs<-regulon$tf[which(regulon$target=='BMP2')] %>% as.character()
tfs<-tfs[tfs%in%rownames(s1)]
regulon.list<-regulon.list[tfs]

#
subregklres<-sapply(regulon.list,function(x){
  n=sum(rownames(s1)%in%x)
  if(n>=5){
    d0<-as.matrix(s0[which(rownames(s0)%in%x),]) %>% pca_process()
    d1<-as.matrix(s1[which(rownames(s1)%in%x),]) %>% pca_process()
    d2<-as.matrix(s2[which(rownames(s2)%in%x),]) %>% pca_process()
    d3<-as.matrix(s3[which(rownames(s3)%in%x),]) %>% pca_process()
    
    kls01<-tryCatch(1/2*(multivar_KLscore(d0,d1)+multivar_KLscore(d1,d0))/T2.test(t(d1),t(d0))$statistic[1],error=function(e){return(NA)})
    kls02<-tryCatch(1/2*(multivar_KLscore(d0,d2)+multivar_KLscore(d2,d0))/T2.test(t(d2),t(d0))$statistic[1],error=function(e){return(NA)})
    kls03<-tryCatch(1/2*(multivar_KLscore(d0,d3)+multivar_KLscore(d3,d0))/T2.test(t(d3),t(d0))$statistic[1],error=function(e){return(NA)})
  
  }else{
    kls01=kls02=kls03=NA
  }
  print(paste0(x,' is calculated'))
  return(c(kls01,kls02,kls03))
})
# 
# names(regklres)<-names(regulon.list)
# 
ex=unique(which(is.na(subregklres),arr.ind = T)[,2])
if(length(ex)>0){
  subregklres<-subregklres[,-ex]
}

row.names(subregklres)<-c('kls01','kls02','kls03','kls12','kls23')
