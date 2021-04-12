

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

#
rs.hum<-c()
for (i in 1:nrow(phase1.all.freqs)){
  a<-cor.test(sort.pops.phase1$h_mean, as.numeric(as.vector(phase1.all.freqs[i,c(3:16)])))
  rs.hum[i]<-a$p.value
}

rs.prec<-c()
for (i in 1:nrow(phase1.all.freqs)){
  a<-cor.test(sort.pops.phase1$p_annual, as.numeric(as.vector(phase1.all.freqs[i,c(3:16)])))
rs.prec[i]<-a$p.value
}

rs.temp<-c()
for (i in 1:nrow(phase1.all.freqs)){
  a<-cor.test(sort.pops.phase1$t_mean, as.numeric(as.vector(phase1.all.freqs[i,c(3:16)])))
rs.temp[i]<-a$p.value
}

rs.corr<-c()
for (i in 1:nrow(phase1.all.freqs)){
  a<-cor.test(sort.pops.phase1$autocorrelated, as.numeric(as.vector(phase1.all.freqs[i,c(3:16)])))
rs.corr[i]<-a$p.value
}

rs.rand<-c()
for (i in 1:nrow(phase1.all.freqs)){
  a<-cor.test(sort.pops.phase1$random, as.numeric(as.vector(phase1.all.freqs[i,c(3:16)])))
rs.rand[i]<-a$p.value
}


dff<-data.frame(chrom=phase1.all.freqs$chrom,
               pos=phase1.all.freqs$POS,
               hum=rs.hum,
               prec=rs.prec,
               temp=rs.temp,
               corr=rs.corr,
               rand=rs.rand)

table(dff$hum < .0001)
table(dff$prec < .0001)
table(dff$temp < .0001)
table(dff$rand < .0001)
table(dff$corr < .0001)


write.csv(dff, "~/Desktop/anoph.phase2/repeated.corrs.csv", row.names = F)







## read in genomic regions and attributes file
genomic_location <- read.gff("~/Downloads/Anopheles_gambiae.AgamP4.47.chr.gff3.gz")
#keep only regions identified as "gene"
gen_regions <- genomic_location[genomic_location$type == "gene", ] # subset to only gene regions

set<-dff[dff$corr < .0001,c(1,2)]
set<-dff[order(dff$temp),c(1,2)][1:100,]
#open empty list to hold output gene lists
gene.dfs <- vector(mode = "list", length = 15)
#define chroms
chroms <- c("2R", "2L", "3R", "3L", "X") # chromosomes of interest
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
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
  }

#pull out gene ID as its own column
gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(gene.list$attributes),';',fixed=TRUE)))[,1]
gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
gene.list$SNPid<-paste(gene.list$seqid,gene.list$`pos[w]`)

unique(gene.list$gene)
length(unique(gene.list$gene))
sum(df$gene %in% unique(gene.list$gene))
phyper(q=sum(df$gene %in% unique(gene.list$gene)),m=194,n=12576-194,
             k=length(unique(gene.list$gene)),lower.tail=FALSE)


par(mfrow=c(3,4))
for (i in paste(set$chrom,set$pos)){
  plot(sort.pops.phase1$h_mean, as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == i,c(3:16)])))
}

cor(sort.pops.phase1$h_mean, as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2L 10224126",c(3:16)])))

