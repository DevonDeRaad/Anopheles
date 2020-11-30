#Do GO term enrichment with all SNPs identified as outliers regardless of vetting

library(ape)
library(clusterProfiler)

## read in genomic regions and attributes file
genomic_location <- read.gff("~/Downloads/Anopheles_gambiae.AgamP4.47.chr.gff3.gz")
#keep only regions identified as "gene"
gen_regions <- genomic_location[genomic_location$type == "gene", ] # subset to only gene regions

#pull all genes queried
all.gene.id<-data.frame(do.call('rbind', strsplit(as.character(gen_regions$attributes),';',fixed=TRUE)))[,1]
all.genes<-data.frame(do.call('rbind', strsplit(as.character(all.gene.id),':',fixed=TRUE)))[,2]

#read in each GO term to gene map (col1=geneID,col2=GOterm)
mart.export<-read.table("~/Downloads/mart_export (3).txt", sep = "\t", header = T)

#read in all sig values
all.sigs<-read.table(file="~/Desktop/anoph.phase2/all.sigs.adjusted.txt", header=T)
head(all.sigs)

###LFMM HUM
lfmm.hum.outlier.table<-all.sigs[,c(1,2)][all.sigs$p_value_hannual < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for lfmm
lfmm.hum.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(lfmm.hum.outlier.table[lfmm.hum.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      lfmm.hum.gene.list<-rbind(lfmm.hum.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
lfmm.hum.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(lfmm.hum.gene.list$attributes),';',fixed=TRUE)))[,1]
lfmm.hum.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(lfmm.hum.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
lfmm.hum.gene.list$SNPid<-paste(lfmm.hum.gene.list$seqid,lfmm.hum.gene.list$`pos[w]`)
#subset
lfmm.hum.gene.list<-lfmm.hum.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(lfmm.hum.outlier.table)
#calc number of unique genes
length(unique(lfmm.hum.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for lfmm prec
lfmm.hum.enriched<-enricher(gene = lfmm.hum.gene.list$gene,
                             universe = all.genes,
                             pAdjustMethod = "fdr",
                             TERM2GENE= mart.export[,c(2,1)]
)
result<-lfmm.hum.enriched@result

lfmm.hum.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(lfmm.hum.enriched.sig)


###LFMM prec
lfmm.prec.outlier.table<-all.sigs[,c(1,2)][all.sigs$p_value_pannual < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for lfmm
lfmm.prec.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(lfmm.prec.outlier.table[lfmm.prec.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      lfmm.prec.gene.list<-rbind(lfmm.prec.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
lfmm.prec.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(lfmm.prec.gene.list$attributes),';',fixed=TRUE)))[,1]
lfmm.prec.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(lfmm.prec.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
lfmm.prec.gene.list$SNPid<-paste(lfmm.prec.gene.list$seqid,lfmm.prec.gene.list$`pos[w]`)
#subset
lfmm.prec.gene.list<-lfmm.prec.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(lfmm.prec.outlier.table)
#calc number of unique genes
length(unique(lfmm.prec.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for lfmm prec
lfmm.prec.enriched<-enricher(gene = lfmm.prec.gene.list$gene,
                            universe = all.genes,
                            pAdjustMethod = "fdr",
                            TERM2GENE= mart.export[,c(2,1)]
)
result<-lfmm.prec.enriched@result

lfmm.prec.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(lfmm.prec.enriched.sig)


###LFMM temp
lfmm.temp.outlier.table<-all.sigs[,c(1,2)][all.sigs$p_value_tmean < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for lfmm
lfmm.temp.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(lfmm.temp.outlier.table[lfmm.temp.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      lfmm.temp.gene.list<-rbind(lfmm.temp.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
lfmm.temp.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(lfmm.temp.gene.list$attributes),';',fixed=TRUE)))[,1]
lfmm.temp.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(lfmm.temp.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
lfmm.temp.gene.list$SNPid<-paste(lfmm.temp.gene.list$seqid,lfmm.temp.gene.list$`pos[w]`)
#subset
lfmm.temp.gene.list<-lfmm.temp.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(lfmm.temp.outlier.table)
#calc number of unique genes
length(unique(lfmm.temp.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for lfmm prec
lfmm.temp.enriched<-enricher(gene = lfmm.temp.gene.list$gene,
                            universe = all.genes,
                            pAdjustMethod = "fdr",
                            TERM2GENE= mart.export[,c(2,1)]
)
result<-lfmm.temp.enriched@result

lfmm.temp.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(lfmm.temp.enriched.sig)


###baye HUM
baye.hum.outlier.table<-all.sigs[,c(1,2)][all.sigs$hum.q < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for baye
baye.hum.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(baye.hum.outlier.table[baye.hum.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      baye.hum.gene.list<-rbind(baye.hum.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
baye.hum.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(baye.hum.gene.list$attributes),';',fixed=TRUE)))[,1]
baye.hum.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(baye.hum.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
baye.hum.gene.list$SNPid<-paste(baye.hum.gene.list$seqid,baye.hum.gene.list$`pos[w]`)
#subset
baye.hum.gene.list<-baye.hum.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(baye.hum.outlier.table)
#calc number of unique genes
length(unique(baye.hum.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for baye prec
baye.hum.enriched<-enricher(gene = baye.hum.gene.list$gene,
                            universe = all.genes,
                            pAdjustMethod = "fdr",
                            TERM2GENE= mart.export[,c(2,1)]
)
result<-baye.hum.enriched@result

baye.hum.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(baye.hum.enriched.sig)


###baye prec
baye.prec.outlier.table<-all.sigs[,c(1,2)][all.sigs$prec.q < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for baye
baye.prec.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(baye.prec.outlier.table[baye.prec.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      baye.prec.gene.list<-rbind(baye.prec.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
baye.prec.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(baye.prec.gene.list$attributes),';',fixed=TRUE)))[,1]
baye.prec.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(baye.prec.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
baye.prec.gene.list$SNPid<-paste(baye.prec.gene.list$seqid,baye.prec.gene.list$`pos[w]`)
#subset
baye.prec.gene.list<-baye.prec.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(baye.prec.outlier.table)
#calc number of unique genes
length(unique(baye.prec.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for baye prec
baye.prec.enriched<-enricher(gene = baye.prec.gene.list$gene,
                             universe = all.genes,
                             pAdjustMethod = "fdr",
                             TERM2GENE= mart.export[,c(2,1)]
)
result<-baye.prec.enriched@result

baye.prec.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(baye.prec.enriched.sig)


###baye temp
baye.temp.outlier.table<-all.sigs[,c(1,2)][all.sigs$temp.q < .01,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for baye
baye.temp.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(baye.temp.outlier.table[baye.temp.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      baye.temp.gene.list<-rbind(baye.temp.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
baye.temp.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(baye.temp.gene.list$attributes),';',fixed=TRUE)))[,1]
baye.temp.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(baye.temp.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
baye.temp.gene.list$SNPid<-paste(baye.temp.gene.list$seqid,baye.temp.gene.list$`pos[w]`)
#subset
baye.temp.gene.list<-baye.temp.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(baye.temp.outlier.table)
#calc number of unique genes
length(unique(baye.temp.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for baye prec
baye.temp.enriched<-enricher(gene = baye.temp.gene.list$gene,
                             universe = all.genes,
                             pAdjustMethod = "fdr",
                             TERM2GENE= mart.export[,c(2,1)]
)
result<-baye.temp.enriched@result

baye.temp.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(baye.temp.enriched.sig)


###r2vim HUM
r2vim.hum.outlier.table<-all.sigs[,c(1,2)][all.sigs$vim.hum > 1,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for r2vim
r2vim.hum.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(r2vim.hum.outlier.table[r2vim.hum.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      r2vim.hum.gene.list<-rbind(r2vim.hum.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
r2vim.hum.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(r2vim.hum.gene.list$attributes),';',fixed=TRUE)))[,1]
r2vim.hum.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(r2vim.hum.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
r2vim.hum.gene.list$SNPid<-paste(r2vim.hum.gene.list$seqid,r2vim.hum.gene.list$`pos[w]`)
#subset
r2vim.hum.gene.list<-r2vim.hum.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(r2vim.hum.outlier.table)
#calc number of unique genes
length(unique(r2vim.hum.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for r2vim prec
r2vim.hum.enriched<-enricher(gene = r2vim.hum.gene.list$gene,
                            universe = all.genes,
                            pAdjustMethod = "fdr",
                            TERM2GENE= mart.export[,c(2,1)]
)
result<-r2vim.hum.enriched@result

r2vim.hum.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(r2vim.hum.enriched.sig)


###r2vim prec
r2vim.prec.outlier.table<-all.sigs[,c(1,2)][all.sigs$vim.prec > 1,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for r2vim
r2vim.prec.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(r2vim.prec.outlier.table[r2vim.prec.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      r2vim.prec.gene.list<-rbind(r2vim.prec.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
r2vim.prec.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(r2vim.prec.gene.list$attributes),';',fixed=TRUE)))[,1]
r2vim.prec.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(r2vim.prec.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
r2vim.prec.gene.list$SNPid<-paste(r2vim.prec.gene.list$seqid,r2vim.prec.gene.list$`pos[w]`)
#subset
r2vim.prec.gene.list<-r2vim.prec.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(r2vim.prec.outlier.table)
#calc number of unique genes
length(unique(r2vim.prec.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for r2vim prec
r2vim.prec.enriched<-enricher(gene = r2vim.prec.gene.list$gene,
                             universe = all.genes,
                             pAdjustMethod = "fdr",
                             TERM2GENE= mart.export[,c(2,1)]
)
result<-r2vim.prec.enriched@result

r2vim.prec.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(r2vim.prec.enriched.sig)


###r2vim temp
r2vim.temp.outlier.table<-all.sigs[,c(1,2)][all.sigs$vim.temp > 1,]
#define chroms
chroms <- c("2R", "2L", "3R", "3L","X") # chromosomes of interest
#initialize gene.list for r2vim
r2vim.temp.gene.list<-data.frame()
#fill gene.list by matching location of genes of interest in gene feature regions by chromosomes
for (i in 1:length(chroms)){
  genloc <- gen_regions[gen_regions$seqid == chroms[i], ] # genome info per chromosome
  pos <- as.numeric(as.character(r2vim.temp.outlier.table[r2vim.temp.outlier.table$chrom == chroms[i],2])) # outlierFST position per chromosome
  for (w in 1:length(pos)) {
    if (nrow(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,]) !=0){
      r2vim.temp.gene.list<-rbind(r2vim.temp.gene.list, cbind(genloc[pos[w] >= genloc$start & pos[w] <= genloc$end,], pos[w]))
    }
  }
  cat("chromosome", i, "of", length(chroms), "finished\n\n")
}
#pull out gene ID as its own column
r2vim.temp.gene.list$gene.id<-data.frame(do.call('rbind', strsplit(as.character(r2vim.temp.gene.list$attributes),';',fixed=TRUE)))[,1]
r2vim.temp.gene.list$gene<-data.frame(do.call('rbind', strsplit(as.character(r2vim.temp.gene.list$gene.id),':',fixed=TRUE)))[,2]
#make SNP ID column
r2vim.temp.gene.list$SNPid<-paste(r2vim.temp.gene.list$seqid,r2vim.temp.gene.list$`pos[w]`)
#subset
r2vim.temp.gene.list<-r2vim.temp.gene.list[,c("SNPid","gene")]
#number of candidate snps
nrow(r2vim.temp.outlier.table)
#calc number of unique genes
length(unique(r2vim.temp.gene.list$gene))
#overrepresentation testing
#calculate overrepresentation for r2vim prec
r2vim.temp.enriched<-enricher(gene = r2vim.temp.gene.list$gene,
                             universe = all.genes,
                             pAdjustMethod = "fdr",
                             TERM2GENE= mart.export[,c(2,1)]
)
result<-r2vim.temp.enriched@result

r2vim.temp.enriched.sig<-result[result$p.adjust <= .05,]
#calc number of enriched GO terms
nrow(r2vim.temp.enriched.sig)





