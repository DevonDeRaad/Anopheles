#15 May 2020
#Script to turn genotype matrices into allele frequency csv files for making plots

calc.allele.freqs <- function(x, pos, pop.id){
  #calculate the number of haplotypes (2 * # of non-missing genotypes)
  #haplotypes<-(colSums(x >=0))*2
  #calculate non-reference haplotypes, 1 means het, 2 means hom alt, so add the # of 1 genotypes to 2* # of 2 genotypes
  #count.alt.haplotypes<-(colSums(x == 1) + (colSums(x ==2)*2))
  ##calculate allele frequency
  ##allele.freq<- count.alt.haplotypes / haplotypes
  ##allele.freq<-((colSums(x == 1) + (colSums(x ==2)*2)))/((colSums(x >=0))*2)
  
  #pop.id must be a vector with length=nrow(chrom.matrix) with a pop identifier for each row (sample) in the matrix
  #pos must be a vector with length=ncol(chrom.matrix) with a snp ID for each column (snp) in the matrix
  
  #initialize allele.freq df with pos
  allele.freq<-data.frame(pos)
  
  #initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(levels(as.factor(pop.id))), style = 3)
  j<-1 #loop tracker
  
  #loop calculates alternate allele frequency for each pop, stores vector in a pop-specific column of allele freq df
  for (i in levels(as.factor(pop.id))){
    #calc alt allele freq and store in pop.freq vector
    pop.freq<-(colSums(x[pop.id == i,] ==1) +(colSums(x[pop.id == i,] ==2)*2)) / (colSums(x[pop.id == i,] >=0)*2)
    #assign pop freq vector to column in df
    allele.freq<-cbind(allele.freq, pop.freq)
    #update progress bar
    setTxtProgressBar(pb, j)
    j<-j+1 #update loop tracker
  }
  
  close(pb)
  
  #fix rownames
  colnames(allele.freq) <- c("POS", levels(as.factor(pop.id)))
  #export
  return(allele.freq)
  
}

#bring in dataframe with locality and environmental data for each individual
old.data<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
phase1.alldata<-read.csv(file = "~/Downloads/phase1_all_vairables_last.csv")
#reorder new df based on old df
phase1.alldata<-phase1.alldata[match(old.data$ox_code,phase1.alldata$ox_code),]
rownames(phase1.alldata)<-1:nrow(phase1.alldata)

#identify pop.id (vector equal in length to number of rows in chrom matrix w/ pop identifier for each sample)
pops<-phase1.alldata$latitude

#identify pos from matrix output for X
pos.X<-read.table(file = "~/Downloads/maf05.ld100.matrix.X.012.pos")
pos.X<-pos.X$V2
#identify pos from 2L
pos.2L<-read.table(file = "~/Downloads/maf05.ld100.matrix.2L.012.pos")
pos.2L<-pos.2L$V2
#identify pos from 2R
pos.2R<-read.table(file = "~/Downloads/maf05.ld100.matrix.2R.012.pos")
pos.2R<-pos.2R$V2
#identify pos from 3L
pos.3L<-read.table(file = "~/Downloads/maf05.ld100.matrix.3L.012.pos")
pos.3L<-pos.3L$V2
#identify pos from 3R
pos.3R<-read.table(file = "~/Downloads/maf05.ld100.matrix.3R.012.pos")
pos.3R<-pos.3R$V2

##
##

#bring in X chrom SNP matrix
phase1.X <- data.matrix(read.csv("~/Downloads/zip.maf05.ld100.matrix.X.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 2L chrom SNP matrix
phase1.2L <- data.matrix(read.csv("~/Downloads/zip.maf05.ld100.matrix.2L.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 2R chrom SNP matrix
phase1.2R <- data.matrix(read.csv("~/Downloads/zip.maf05.ld100.matrix.2R.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 3L chrom SNP matrix
phase1.3L <- data.matrix(read.csv("~/Downloads/zip.maf05.ld100.matrix.3L.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 3R chrom SNP matrix
phase1.3R <- data.matrix(read.csv("~/Downloads/zip.maf05.ld100.matrix.3R.012", header = FALSE, row.names = 1, sep = "\t"))

#calc allele freqs and write them out as a csv for X chrom phase1
phase1.X.freqs<-calc.allele.freqs(phase1.X, pos=pos.X, pop.id=pops)
#add chrom info
phase1.X.freqs<-data.frame(chrom="X",phase1.X.freqs)
write.csv(phase1.X.freqs, file = "~/Desktop/anoph.phase2/phase1.X.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase1.2L.freqs<-calc.allele.freqs(phase1.2L, pos=pos.2L, pop.id = pops)
phase1.2L.freqs<-data.frame(chrom="2L",phase1.2L.freqs)
write.csv(phase1.2L.freqs, file = "~/Desktop/anoph.phase2/phase1.2L.freqs.csv")

#calc allele freqs and write them out as a csv for 2R chrom phase1
phase1.2R.freqs<-calc.allele.freqs(phase1.2R, pos=pos.2R, pop.id = pops)
phase1.2R.freqs<-data.frame(chrom="2R",phase1.2R.freqs)
write.csv(phase1.2R.freqs, file = "~/Desktop/anoph.phase2/phase1.2R.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase1.3L.freqs<-calc.allele.freqs(phase1.3L, pos=pos.3L, pop.id = pops)
phase1.3L.freqs<-data.frame(chrom="3L",phase1.3L.freqs)
write.csv(phase1.3L.freqs, file = "~/Desktop/anoph.phase2/phase1.3L.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase1.3R.freqs<-calc.allele.freqs(phase1.3R, pos=pos.3R, pop.id = pops)
phase1.3R.freqs<-data.frame(chrom="3R",phase1.3R.freqs)
write.csv(phase1.3R.freqs, file = "~/Desktop/anoph.phase2/phase1.3R.freqs.csv")

#write full phase1 allele freq file
write.csv(as.data.frame(rbind(phase1.X.freqs, phase1.2L.freqs, phase1.2R.freqs, phase1.3L.freqs, phase1.3R.freqs)),
          file = "~/Desktop/anoph.phase2/phase1.all.freqs.csv", row.names = F)

head(as.data.frame(rbind(phase1.X.freqs, phase1.2L.freqs, phase1.2R.freqs, phase1.3L.freqs, phase1.3R.freqs)))
dim(as.data.frame(rbind(phase1.X.freqs, phase1.2L.freqs, phase1.2R.freqs, phase1.3L.freqs, phase1.3R.freqs)))

###
###Phase2
###

#bring in dataframe with locality and environmental data for each individual
phase2.alldata<-read.csv(file="~/Downloads/phase2.unique.sampling.locs.csv")
#identify pop.id (vector equal in length to number of rows in chrom matrix w/ pop identifier for each sample)
pops<-phase2.alldata$latitude

#identify pos from matrix output for X
pos.X<-read.table(file = "~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.X.geno.matrix.012.pos")
pos.X<-pos.X$V2
#identify pos from 2L
pos.2L<-read.table(file = "~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2L.geno.matrix.012.pos")
pos.2L<-pos.2L$V2
#identify pos from 2R
pos.2R<-read.table(file = "~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2R.geno.matrix.012.pos")
pos.2R<-pos.2R$V2
#identify pos from 3L
pos.3L<-read.table(file = "~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3L.geno.matrix.012.pos")
pos.3L<-pos.3L$V2
#identify pos from 3R
pos.3R<-read.table(file = "~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3R.geno.matrix.012.pos")
pos.3R<-pos.3R$V2

##
##

#bring in X chrom SNP matrix
phase2.X <- data.matrix(read.csv("~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.X.geno.matrix.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 2L chrom SNP matrix
phase2.2L <- data.matrix(read.csv("~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2L.geno.matrix.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 2R chrom SNP matrix
phase2.2R <- data.matrix(read.csv("~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2R.geno.matrix.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 3L chrom SNP matrix
phase2.3L <- data.matrix(read.csv("~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3L.geno.matrix.012", header = FALSE, row.names = 1, sep = "\t"))
#bring in 3R chrom SNP matrix
phase2.3R <- data.matrix(read.csv("~/Downloads/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3R.geno.matrix.012", header = FALSE, row.names = 1, sep = "\t"))

#calc allele freqs and write them out as a csv for X chrom phase1
phase2.X.freqs<-calc.allele.freqs(phase2.X, pos=pos.X, pop.id = pops)
phase2.X.freqs<-data.frame(chrom="X",phase2.X.freqs)
write.csv(phase2.X.freqs, file = "~/Desktop/anoph.phase2/phase2.X.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase2.2L.freqs<-calc.allele.freqs(phase2.2L, pos=pos.2L, pop.id = pops)
phase2.2L.freqs<-data.frame(chrom="2L",phase2.2L.freqs)
write.csv(phase2.2L.freqs, file = "~/Desktop/anoph.phase2/phase2.2L.freqs.csv")

#calc allele freqs and write them out as a csv for 2R chrom phase1
phase2.2R.freqs<-calc.allele.freqs(phase2.2R, pos=pos.2R, pop.id = pops)
phase2.2R.freqs<-data.frame(chrom="2R",phase2.2R.freqs)
write.csv(phase2.2R.freqs, file = "~/Desktop/anoph.phase2/phase2.2R.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase2.3L.freqs<-calc.allele.freqs(phase2.3L, pos=pos.3L, pop.id = pops)
phase2.3L.freqs<-data.frame(chrom="3L",phase2.3L.freqs)
write.csv(phase2.3L.freqs, file = "~/Desktop/anoph.phase2/phase2.3L.freqs.csv")

#calc allele freqs and write them out as a csv for 2L chrom phase1
phase2.3R.freqs<-calc.allele.freqs(phase2.3R, pos=pos.3R, pop.id = pops)
phase2.3R.freqs<-data.frame(chrom="3R",phase2.3R.freqs)
write.csv(phase2.3R.freqs, file = "~/Desktop/anoph.phase2/phase2.3R.freqs.csv")

#write out all chroms to single outfile
#write full phase1 allele freq file
write.csv(as.data.frame(rbind(phase2.X.freqs, phase2.2L.freqs, phase2.2R.freqs, phase2.3L.freqs, phase2.3R.freqs)),
          file = "~/Desktop/anoph.phase2/phase2.all.freqs.csv", row.names = F)
