library(ape)
library(clusterProfiler)
library(ggplot2)

#bring in all significance values
all.sigs<-read.csv("~/Downloads/all.sigs.csv")
#bring in phase1 allele frequencies
phase1.all.freqs<-read.csv(file = "~/Desktop/anoph.phase2/phase1.all.freqs.csv")
#make column for id
phase1.all.freqs$id<-paste(phase1.all.freqs$chrom,phase1.all.freqs$POS)
#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Downloads/phase1_all_vairables_last.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","p_annual","t_mean","h_mean","autocorrelated","random")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

## read in genomic regions and attributes file
genomic_location <- read.gff("~/Downloads/Anopheles_gambiae.AgamP4.47.chr.gff3.gz")

#keep only regions identified as "gene"
gen_regions <- genomic_location[genomic_location$type == "gene", ] # subset to only gene regions

#pull all genes queried
all.gene.id<-data.frame(do.call('rbind', strsplit(as.character(gen_regions$attributes),';',fixed=TRUE)))[,1]
all.genes<-data.frame(do.call('rbind', strsplit(as.character(all.gene.id),':',fixed=TRUE)))[,2]
#write.table(all.genes, "~/Downloads/all.genes.anopheles.gambiae.txt", row.names = F, col.names = F, quote = F)


###
##
#bring in vetted outliers for all GEA tests
out<-read.csv("~/Desktop/anoph.phase2/vetted.outliers.csv")
#retain only vetted outliers
out<-out[out$vet== "pass",]

#split by test
X <- split(out, out$test)
str(X)

#open empty list to hold output gene lists
gene.dfs <- vector(mode = "list", length = 0)
#define chroms
chroms <- c("2R", "2L", "3R", "3L", "X") # chromosomes of interest
#define the tests conducted
tests<-c("lfmm.hum","lfmm.prec","lfmm.temp","lfmm.rand","lfmm.corr",
         "baye.hum","baye.prec","baye.temp","baye.rand","baye.corr",
         "vim.hum","vim.prec","vim.temp","vim.rand","vim.corr",
         "cors.hum","cors.prec","cors.temp","cors.rand","cors.corr")
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in tests){
  set<-X[[i]] #use the ith GEA candidate SNP dataset from the list of all 15 stored in 'X'
  gene.list<-data.frame() #initialize empty gene.list
  for (j in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[j], ] #subet genome-wide gene info into only the relevant chromosome
  pos <- as.numeric(as.character(set[set$chrom == chroms[j],2])) #subet list of GEA outlier SNPs to only the relevant chromosome
    for (w in 1:length(pos)) {
      #retain all genes within an outlier SNP within a 10kb distance of the start or end of the gene
      if (nrow(genloc[pos[w] >= genloc$start-1000 & pos[w] <= genloc$end+1000,]) !=0){
      gene.list<-rbind(gene.list, cbind(genloc[pos[w] >= genloc$start-1000 & pos[w] <= genloc$end+1000,], pos[w]))
      }
    }
  cat("dataset", i, "in progress\n\n") #print the chromosome
  }
  gene.dfs[[i]]<-gene.list #store the set of associated genes as a distinct df in this list
}

#for each of the 20 datasets
for (i in tests){
  #pull out gene ID as its own column
  gene.dfs[[i]]$gene.id<-data.frame(do.call('rbind', strsplit(as.character(gene.dfs[[i]]$attributes),';',fixed=TRUE)))[,1]
  gene.dfs[[i]]$gene<-data.frame(do.call('rbind', strsplit(as.character(gene.dfs[[i]]$gene.id),':',fixed=TRUE)))[,2]
  #make SNP ID column
  gene.dfs[[i]]$SNPid<-paste(gene.dfs[[i]]$seqid,gene.dfs[[i]]$`pos[w]`)
}

#calc number of unique genes in each outlier SNP dataset
j=1
for (i in tests){
  print(tests[j])
  print(length(unique(gene.dfs[[i]]$gene)))
  j=j+1
  }

#write out each candidate gene dataset to publish on the github page
write.csv(cbind(gene.dfs[["baye.corr"]][,c(1,10,12)]), "~/Downloads/baye.corr.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["baye.hum"]][,c(1,10,12)]), "~/Downloads/baye.hum.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["baye.prec"]][,c(1,10,12)]), "~/Downloads/baye.prec.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["baye.rand"]][,c(1,10,12)]), "~/Downloads/baye.rand.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["baye.temp"]][,c(1,10,12)]), "~/Downloads/baye.temp.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["lfmm.corr"]][,c(1,10,12)]), "~/Downloads/lfmm.corr.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["lfmm.hum"]][,c(1,10,12)]), "~/Downloads/lfmm.hum.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["lfmm.prec"]][,c(1,10,12)]), "~/Downloads/lfmm.prec.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["lfmm.rand"]][,c(1,10,12)]), "~/Downloads/lfmm.rand.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["lfmm.temp"]][,c(1,10,12)]), "~/Downloads/lfmm.temp.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["vim.corr"]][,c(1,10,12)]), "~/Downloads/r2VIM.corr.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["vim.hum"]][,c(1,10,12)]), "~/Downloads/r2VIM.hum.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["vim.prec"]][,c(1,10,12)]), "~/Downloads/r2VIM.prec.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["vim.rand"]][,c(1,10,12)]), "~/Downloads/r2VIM.rand.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["vim.temp"]][,c(1,10,12)]), "~/Downloads/r2VIM.temp.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["cors.corr"]][,c(1,10,12)]), "~/Downloads/cors.corr.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["cors.hum"]][,c(1,10,12)]), "~/Downloads/cors.hum.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["cors.prec"]][,c(1,10,12)]), "~/Downloads/cors.prec.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["cors.rand"]][,c(1,10,12)]), "~/Downloads/cors.rand.gene.list.csv", row.names = F, col.names = T, quote = F)
write.csv(cbind(gene.dfs[["cors.temp"]][,c(1,10,12)]), "~/Downloads/cors.temp.gene.list.csv", row.names = F, col.names = T, quote = F)

#read in 2La dessication-tolerance associated genes from Ayala et al. 2019
lades<-read.csv("~/Downloads/2la.dess.csv")
lades[,1]<-as.character(lades[,1])
lades<-strsplit(lades[,1], split=" ")
gene.2la<-unlist(lapply(lades, '[[', 3))
gene.2la

#read in 2Rb dessication-tolerance associated genes from Ayala et al. 2019
rbdes<-read.csv("~/Downloads/2rb.dess.csv")
rbdes[,1]<-as.character(rbdes[,1])
rbdes<-strsplit(rbdes[,1], split=" ")
gene.2rb<-unlist(lapply(rbdes, '[[', 3))
gene.2rb

#combine in a single dataframe
df<-data.frame(gene=c(unique(gene.2la),unique(gene.2rb)), inversion=c(rep("la", times=length(unique(gene.2la))),
                                                                      rep("rb", times=length(unique(gene.2rb)))))

#calc total number of unique dessication associated genes identified in Ayala et al. 2019
length(unique(df$gene)) #194

#calculate number of genes predicted by each dataset that were also found in Ayala candidate gene set
for (i in tests){
  print(i)
  print(sum(df$gene %in% unique(gene.dfs[[i]]$gene)))}

#calculate significance
for (i in tests){
  print(i)
  print(phyper(q=sum(df$gene %in% unique(gene.dfs[[i]]$gene)),m=194,n=12576-194,
               k=length(unique(gene.dfs[[i]]$gene)),lower.tail=FALSE))
}

#Is top candidate gene from 2La found in each dataset
for (i in tests){
  print(i)
  print("AGAP006026" %in% gene.dfs[[i]]$gene)
}

#find significant SNPs in these genes
gene.dfs[["cors.prec"]][gene.dfs[["cors.prec"]]$gene == "AGAP006026",]

#plot AGAP006026
frq1<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2L 25188870",3:16]))
df<-data.frame(frq=frq1,
               prec=sort.pops.phase1$p_annual,
               SNP=c(rep("2L 25188870",times=14)))
prec.6026<-ggplot(df, aes(x=prec, y=frq, color=SNP)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  labs(x="Precipitation (mm/year)", y="Allele Frequency")+
  ggtitle("AGAP006026")

#find significant SNPs in these genes
gene.dfs[["lfmm.hum"]][gene.dfs[["lfmm.hum"]]$gene == "AGAP006026",]

#plot AGAP006026
frq1<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2L 25194338",3:16]))
df<-data.frame(frq=frq1,
               hum=sort.pops.phase1$h_mean,
               SNP=c(rep("2L 25194338",times=14)))
hum.6026<-ggplot(df, aes(x=hum, y=frq, color=SNP)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  labs(x="Humidity (10^5 Kg H20 / Kg air)", y="Allele Frequency")



for (i in tests){
  print(i)
  print("AGAP002578" %in% gene.dfs[[i]]$gene)
}

gene.dfs[["cors.prec"]][gene.dfs[["cors.prec"]]$gene == "AGAP002578",]

frq1<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2L 25188870",3:16]))
frq2<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23163309",3:16]))
frq3<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23183183",3:16]))
frq4<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23192809",3:16]))
frq5<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23203396",3:16]))

df<-data.frame(frq=c(frq1,frq2,frq3,frq4,frq5),
               prec=rep(sort.pops.phase1$p_annual, times=5),
               SNP=c(rep("2L 25188870",times=14),rep("2R 23163309",times=14),rep("2R 23183183",times=14),rep("2R 23192809",times=14),
                     rep("2R 23203396",times=14)))
prec.2578<-ggplot(df, aes(x=prec, y=frq, color=SNP)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  labs(x="Precipitation (mm/year)", y="Allele Frequency")
p<-arrangeGrob(prec.2578, hum.6026, nrow = 1)
ggsave("~/Downloads/allele.freq.fig.pdf", p, width=8, height=2.75)







