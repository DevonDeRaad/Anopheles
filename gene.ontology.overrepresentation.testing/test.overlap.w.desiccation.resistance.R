
lades<-read.csv("~/Downloads/2la.dess.csv")
lades[,1]<-as.character(lades[,1])
lades<-strsplit(lades[,1], split=" ")
gene.2la<-unlist(lapply(lades, '[[', 3))
gene.2la

rbdes<-read.csv("~/Downloads/2rb.dess.csv")
rbdes[,1]<-as.character(rbdes[,1])
rbdes<-strsplit(rbdes[,1], split=" ")
gene.2rb<-unlist(lapply(rbdes, '[[', 3))
gene.2rb

df<-data.frame(gene=c(unique(gene.2la),unique(gene.2rb)), inversion=c(rep("la", times=length(unique(gene.2la))),
                                                                      rep("rb", times=length(unique(gene.2rb)))))
#write.csv(df, "~/Downloads/mapping.dessication.resistance.agam.csv")

## read in lfmm hum outlier genes
lfmm.hum<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/lfmm.hum.gene.snp.match.txt")[-1,]
lfmm.hum<-as.vector(unique(lfmm.hum[,3]))
## read in baye hum outlier genes
baye.hum<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/bayescan.hum.gene.snp.match.txt")[-1,]
baye.hum<-as.vector(unique(baye.hum[,3]))
## read in r2vim hum outlier genes
r2vim.hum<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/r2vim.hum.gene.snp.match.txt")[-1,]
r2vim.hum<-as.vector(unique(r2vim.hum[,3]))

## read in lfmm prec outlier genes
lfmm.prec<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/lfmm.prec.gene.snp.match.txt")[-1,]
lfmm.prec<-as.vector(unique(lfmm.prec[,3]))
## read in baye prec outlier genes
baye.prec<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/bayescan.prec.gene.snp.match.txt")[-1,]
baye.prec<-as.vector(unique(baye.prec[,3]))
## read in r2vim prec outlier genes
r2vim.prec<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/r2vim.prec.gene.snp.match.txt")[-1,]
r2vim.prec<-as.vector(unique(r2vim.prec[,3]))

#all dessication genes
dess<-unique(c(lfmm.hum,lfmm.prec,baye.hum,baye.prec,r2vim.hum,r2vim.prec))
sum(df$gene %in% dess)
#calculate proportion of genes predicted by each relevant dataset that were also found in Ayala 2018 candidate gene set
sum(df$gene %in% lfmm.hum)/length(lfmm.hum)
sum(df$gene %in% lfmm.prec)/length(lfmm.prec)
sum(df$gene %in% baye.hum)/length(baye.hum)
sum(df$gene %in% baye.prec)/length(baye.prec)
sum(df$gene %in% r2vim.hum)/length(r2vim.hum)
sum(df$gene %in% r2vim.prec)/length(r2vim.prec)

#calculate significance of each value
phyper(q=sum(df$gene %in% lfmm.hum),m=194,n=12576-194,k=length(lfmm.hum),lower.tail=FALSE)
phyper(q=sum(df$gene %in% lfmm.prec),m=194,n=12576-194,k=length(lfmm.prec),lower.tail=FALSE)
phyper(q=sum(df$gene %in% baye.hum),m=194,n=12576-194,k=length(baye.hum),lower.tail=FALSE)
phyper(q=sum(df$gene %in% baye.prec),m=194,n=12576-194,k=length(baye.prec),lower.tail=FALSE)
phyper(q=sum(df$gene %in% r2vim.hum),m=194,n=12576-194,k=length(r2vim.hum),lower.tail=FALSE)
phyper(q=sum(df$gene %in% r2vim.prec),m=194,n=12576-194,k=length(r2vim.prec),lower.tail=FALSE)

32/194
phyper(152,900,10000-900,1500,lower.tail=FALSE)
phyper(q=32,m=194,n=12576-194,k=813,lower.tail=FALSE)
20/813
194/12576


#check whether top genes are in each dataset
"AGAP006026" %in% lfmm.hum
"AGAP006026" %in% lfmm.prec
"AGAP006026" %in% baye.hum
"AGAP006026" %in% baye.prec
"AGAP006026" %in% r2vim.hum
"AGAP006026" %in% r2vim.prec

"AGAP002578" %in% lfmm.hum
"AGAP002578" %in% lfmm.prec
"AGAP002578" %in% baye.hum
"AGAP002578" %in% baye.prec
"AGAP002578" %in% r2vim.hum
"AGAP002578" %in% r2vim.prec


## read in lfmm temp outlier genes
lfmm.temp<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/lfmm.temp.gene.snp.match.txt")[-1,]
lfmm.temp<-as.vector(unique(lfmm.temp[,3]))
## read in baye temp outlier genes
baye.temp<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/bayescan.temp.gene.snp.match.txt")[-1,]
baye.temp<-as.vector(unique(baye.temp[,3]))
## read in r2vim temp outlier genes
r2vim.temp<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/r2vim.temp.gene.snp.match.txt")[-1,]
r2vim.temp<-as.vector(unique(r2vim.temp[,3]))

#read in differentially expressed candidate genes in response to heat treatment
heat<-read.csv("~/Downloads/heat.tolerance.outliers.csv")
heat<-as.data.frame(heat[-319,])
heat<-as.data.frame(heat[-132,])

heat[,1]<-gsub("\\+ |\\+|# ","",heat[,1])

969/12576
#calculate proportion of genes predicted by each relevant dataset that were also found in heat exposure set
sum(heat[,1] %in% lfmm.temp)/length(lfmm.temp)
sum(heat[,1] %in% baye.temp)/length(baye.temp)
sum(heat[,1] %in% r2vim.temp)/length(r2vim.temp)

#calculate significance of each value
phyper(q=sum(df$gene %in% lfmm.hum),m=194,n=12576-194,k=length(lfmm.hum),lower.tail=FALSE)
phyper(q=sum(df$gene %in% lfmm.prec),m=194,n=12576-194,k=length(lfmm.prec),lower.tail=FALSE)
phyper(q=sum(df$gene %in% baye.hum),m=194,n=12576-194,k=length(baye.hum),lower.tail=FALSE)



