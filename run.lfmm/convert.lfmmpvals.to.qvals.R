#May 17 2020
#plot results from phase1 and verify with phase2 allele freq patterns
#devtools::install_github('exaexa/scattermore')
library(scattermore)
library(gridExtra)
library(VennDiagram)
library(qvalue)
library(outliers)

#read in all phase1 significance values:
all.sigs<-read.csv(file="~/Desktop/anoph.3.march.2020/all.sigs.csv")
table(all.sigs$chrom)

#convert 0 bayescenv q-vals to 1e-4, to avoid infinite values when plotting with -log10()
all.sigs$prec.q[all.sigs$prec.q == 0]<-.0001
all.sigs$temp.q[all.sigs$temp.q == 0]<-.0001
all.sigs$hum.q[all.sigs$hum.q == 0]<-.0001

#convert lfmm prec p values to genome-wide FDR adjusted q values
for (i in levels(all.sigs$chrom)){
  qvals<-qvalue(all.sigs$prec.p[all.sigs$chrom== i])
  all.sigs$prec.p[all.sigs$chrom== i]<-qvals[["qvalues"]]
}
#convert lfmm hum p values to genome-wide FDR adjusted q values
for (i in levels(all.sigs$chrom)){
  qvals<-qvalue(all.sigs$hum.p[all.sigs$chrom== i])
  all.sigs$hum.p[all.sigs$chrom== i]<-qvals[["qvalues"]]
}
#convert lfmm temp p values to genome-wide FDR adjusted q values
for (i in levels(all.sigs$chrom)){
  qvals<-qvalue(all.sigs$temp.p[all.sigs$chrom== i])
  all.sigs$temp.p[all.sigs$chrom== i]<-qvals[["qvalues"]]
}
rm(qvals) #space saver

#give order for plotting manhattan plots
all.sigs$order<-c(all.sigs$pos[all.sigs$chrom== "2R"],
                  all.sigs$pos[all.sigs$chrom== "2L"]+61545105,
                  all.sigs$pos[all.sigs$chrom== "3R"]+61545105+49364325,
                  all.sigs$pos[all.sigs$chrom== "3L"]+61545105+49364325+53200684,
                  all.sigs$pos[all.sigs$chrom== "X"]+61545105+49364325+53200684+41963435)

write.csv(allsigs, file="fdr.adj.all.sigs.csv")
