####enrichment analysis in intercellular crosstalk

net.dif <- subset(net,pval<0.05) #significant signals predicted between PMN and Ctrl

result1 <- data.frame()
a <- nrow(subset(net.dif,datasets == 'MN'))
b <- nrow(subset(net.dif,datasets == 'ctrl'))
for (i in unique(net.dif$interaction_name)) {
  c <- nrow(subset(net.dif,interaction_name == i & datasets == 'MN'))
  d <- nrow(subset(net.dif,interaction_name == i & datasets == 'ctrl'))
  test<-matrix(c(c,d,a-c,b-d),nrow=2,ncol=2,byrow=TRUE)
  qq<-fisher.test(test) # p-val

  fe=ifelse(c!=0 & a!=0,ifelse(d!=0,(c/a)/(d/b),1000),0.001) # fold enrichment 
  fe<-log2(fe)
  result1 <- rbind(result1,c(qq$p.value,fe))
}

colnames(result1) <- c('enrich_pval','logfc_enrich')
result1$enrich_fdr<-p.adjust(result1$enrich_pval,method = 'fdr')
rownames(result1) <- unique(net.dif$interaction_name)


result1$sig <- "Stable"
result1$sig[result1$enrich_pval < 0.05 & result1$logfc_enrich > 0.3]<- "MN"
result1$sig[result1$enrich_pval < 0.05 & result1$logfc_enrich < -0.3] <- "ctrl"
result1$sig<-factor(result1$sig,levels=c('ctrl','Stable','MN'))
result1$interaction_name<-rownames(result1)
table(result1$sig)

EnUp<-intersect(net.up$interaction_name,result1$interaction_name[result1$sig=='MN']) #signals upregulated and enriched in PMN
EnDown<-intersect(net.down$interaction_name,result1$interaction_name[result1$sig=='ctrl']) #signals downregulated and un-enriched in PMN
