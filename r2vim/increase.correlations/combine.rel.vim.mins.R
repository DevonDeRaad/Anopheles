
#script to manually take the min of each SNP across runs that were arrayed due to 
#multithreading issues, and create a single csv with genome wide importance values for each env variable 
setwd("~/Downloads/r2vim.corrs/")

#read in prec
prec.2L.vim<-read.csv("2L.prec.rel.vim.min.csv")
hist(prec.2L.vim$min.rel.vim)
prec.2R.vim<-read.csv("2R.hum.rel.vim.min.csv")
hist(prec.2R.vim$min.rel.vim)

#Start with prec 3L
prec.3L.1<-read.csv("3L.prec.1.rel.vim.min.csv")
prec.3L.2<-read.csv("3L.prec.2.rel.vim.min.csv")
prec.3L.3<-read.csv("3L.prec.3.rel.vim.min.csv")
prec.3L.4<-read.csv("3L.prec.4.rel.vim.min.csv")
prec.3L.5<-read.csv("3L.prec.5.rel.vim.min.csv")
prec.3L.6<-read.csv("3L.prec.6.rel.vim.min.csv")
prec.3L.7<-read.csv("3L.prec.7.rel.vim.min.csv")
prec.3L.8<-read.csv("3L.prec.8.rel.vim.min.csv")
prec.3L.9<-read.csv("3L.prec.9.rel.vim.min.csv")
prec.3L.10<-read.csv("3L.prec.10.rel.vim.min.csv")

prec.3L<-data.frame(chrom=prec.3L.1$chrom,
                    snp=prec.3L.1$snp,
                    run1=prec.3L.1$run.1,run2=prec.3L.2$run.2,
                    run3=prec.3L.3$run.3,run4=prec.3L.4$run.4,
                    run5=prec.3L.5$run.5,run6=prec.3L.6$run.6,
                    run7=prec.3L.7$run.7,run8=prec.3L.8$run.8,
                    run9=prec.3L.9$run.9,run10=prec.3L.10$run.10)

#check convergence
cor(prec.3L[,3:12])
hist(cor(prec.3L[,3:12]))

#calc relvim
for (i in 3:12){
  prec.3L[,i]<- prec.3L[,i] / abs(min(prec.3L[,i]))
}

#calc min
prec.3L.vim<-data.frame(chrom=prec.3L$chrom,
                        snp=prec.3L$snp,
                        min.rel.vim=apply(prec.3L[,3:12], 1, min))
hist(prec.3L.vim$min.rel.vim)

#prec 3R
prec.3R.1<-read.csv("3R.prec.1.rel.vim.min.csv")
prec.3R.2<-read.csv("3R.prec.2.rel.vim.min.csv")
prec.3R.3<-read.csv("3R.prec.3.rel.vim.min.csv")
prec.3R.4<-read.csv("3R.prec.4.rel.vim.min.csv")
prec.3R.5<-read.csv("3R.prec.5.rel.vim.min.csv")
prec.3R.6<-read.csv("3R.prec.6.rel.vim.min.csv")
prec.3R.7<-read.csv("3R.prec.7.rel.vim.min.csv")
prec.3R.8<-read.csv("3R.prec.8.rel.vim.min.csv")
prec.3R.9<-read.csv("3R.prec.9.rel.vim.min.csv")
prec.3R.10<-read.csv("3R.prec.10.rel.vim.min.csv")

prec.3R<-data.frame(chrom=prec.3R.1$chrom,
                    snp=prec.3R.1$snp,
                    run1=prec.3R.1$run.1,run2=prec.3R.2$run.2,
                    run3=prec.3R.3$run.3,run4=prec.3R.4$run.4,
                    run5=prec.3R.5$run.5,run6=prec.3R.6$run.6,
                    run7=prec.3R.7$run.7,run8=prec.3R.8$run.8,
                    run9=prec.3R.9$run.9,run10=prec.3R.10$run.10)

cor(prec.3R[,3:12])
hist(cor(prec.3R[,3:12]))

#calc relvim
for (i in 3:12){
  prec.3R[,i]<- prec.3R[,i] / abs(min(prec.3R[,i]))
}

#calc min
prec.3R.vim<-data.frame(chrom=prec.3R$chrom,
                        snp=prec.3R$snp,
                        min.rel.vim=apply(prec.3R[,3:12], 1, min))
hist(prec.3R.vim$min.rel.vim)

#
prec.X.vim<-read.csv("X.hum.rel.vim.min.csv")
hist(prec.X.vim$min.rel.vim)

#combine all of prec
prec<-as.data.frame(rbind(prec.2L.vim,prec.2R.vim,
                          prec.3L.vim,prec.3R.vim,
                          prec.X.vim))

#move on to temp
temp.2L.vim<-read.csv("2L.temp.rel.vim.min.csv")
hist(temp.2L.vim$min.rel.vim)

#temp 2R
temp.2R.1<-read.csv("2R.temp.6.rel.vim.min.csv")
temp.2R.2<-read.csv("2R.temp.7.rel.vim.min.csv")
temp.2R.3<-read.csv("2R.temp.8.rel.vim.min.csv")
temp.2R.4<-read.csv("2R.temp.9.rel.vim.min.csv")
temp.2R.5<-read.csv("2R.temp.10.rel.vim.min.csv")
temp.2R.6<-read.csv("2R.temp.11.rel.vim.min.csv")
temp.2R.7<-read.csv("2R.temp.12.rel.vim.min.csv")
temp.2R.8<-read.csv("2R.temp.13.rel.vim.min.csv")
temp.2R.9<-read.csv("2R.temp.14.rel.vim.min.csv")
temp.2R.10<-read.csv("2R.temp.15.rel.vim.min.csv")

temp.2R<-data.frame(chrom=temp.2R.1$chrom,
                    snp=temp.2R.1$snp,
                    run1=temp.2R.1$run.1,run2=temp.2R.2$run.2,
                    run3=temp.2R.3$run.3,run4=temp.2R.4$run.4,
                    run5=temp.2R.5$run.5,run6=temp.2R.6$run.6,
                    run7=temp.2R.7$run.7,run8=temp.2R.8$run.8,
                    run9=temp.2R.9$run.9,run10=temp.2R.10$run.10)
cor(temp.2R[,3:12])
hist(cor(temp.2R[,3:12]))

#calc relvim
for (i in 3:12){
  temp.2R[,i]<- temp.2R[,i] / abs(min(temp.2R[,i]))
}

#calc min
temp.2R.vim<-data.frame(chrom=temp.2R$chrom,
                        snp=temp.2R$snp,
                        min.rel.vim=apply(temp.2R[,3:12], 1, min))
hist(temp.2R.vim$min.rel.vim)

#temp 3L
temp.3L.1<-read.csv("3L.temp.1.rel.vim.min.csv")
temp.3L.2<-read.csv("3L.temp.2.rel.vim.min.csv")
temp.3L.3<-read.csv("3L.temp.3.rel.vim.min.csv")
temp.3L.4<-read.csv("3L.temp.4.rel.vim.min.csv")
temp.3L.5<-read.csv("3L.temp.5.rel.vim.min.csv")
temp.3L.6<-read.csv("3L.temp.6.rel.vim.min.csv")
temp.3L.7<-read.csv("3L.temp.7.rel.vim.min.csv")
temp.3L.8<-read.csv("3L.temp.8.rel.vim.min.csv")
temp.3L.9<-read.csv("3L.temp.9.rel.vim.min.csv")
temp.3L.10<-read.csv("3L.temp.10.rel.vim.min.csv")

temp.3L<-data.frame(chrom=temp.3L.1$chrom,
                    snp=temp.3L.1$snp,
                    run1=temp.3L.1$vim.run.1,run2=temp.3L.2$vim.run.2,
                    run3=temp.3L.3$vim.run.3,run4=temp.3L.4$vim.run.4,
                    run5=temp.3L.5$vim.run.5,run6=temp.3L.6$vim.run.6,
                    run7=temp.3L.7$vim.run.7,run8=temp.3L.8$vim.run.8,
                    run9=temp.3L.9$vim.run.9,run10=temp.3L.10$vim.run.10)
cor(temp.3L[,3:12])
hist(temp.3L[,3:12]))

#calc relvim
for (i in 3:12){
  temp.3L[,i]<- temp.3L[,i] / abs(min(temp.3L[,i]))
}

#calc min
temp.3L.vim<-data.frame(chrom=temp.3L$chrom,
                        snp=temp.3L$snp,
                        min.rel.vim=apply(temp.3L[,3:12], 1, min))
hist(temp.3L.vim$min.rel.vim)

#move on to temp 3R
temp.3R.1<-read.csv("3R.temp.1.rel.vim.min.csv")
temp.3R.2<-read.csv("3R.temp.2.rel.vim.min.csv")
temp.3R.3<-read.csv("3R.temp.3.rel.vim.min.csv")
temp.3R.4<-read.csv("3R.temp.4.rel.vim.min.csv")
temp.3R.5<-read.csv("3R.temp.5.rel.vim.min.csv")

temp.3R<-data.frame(chrom=temp.3R.1$chrom,
                    snp=temp.3R.1$snp,
                    run1=temp.3R.1$vim.run.1,run2=temp.3R.1$vim.run.2,
                    run3=temp.3R.2$vim.run.3,run4=temp.3R.2$vim.run.4,
                    run5=temp.3R.3$vim.run.5,run6=temp.3R.3$vim.run.6,
                    run7=temp.3R.4$vim.run.7,run8=temp.3R.4$vim.run.8,
                    run9=temp.3R.5$vim.run.9,run10=temp.3R.5$vim.run.10)
cor(temp.3R[,3:12])
hist(cor(temp.3R[,3:12]))

#calc relvim
for (i in 3:12){
 temp.3R[,i]<- temp.3R[,i] / abs(min(temp.3R[,i]))
}

#calc min
temp.3R.vim<-data.frame(chrom=temp.3R$chrom,
                        snp=temp.3R$snp,
                        min.rel.vim=apply(temp.3R[,3:12], 1, min))
hist(temp.3R.vim$min.rel.vim)

#X
temp.X.vim<-read.csv("X.temp.rel.vim.min.csv")
hist(temp.X.vim$min.rel.vim)

#combine all of temp
temp<-as.data.frame(rbind(temp.2L.vim,temp.2R.vim,
                          temp.3L.vim,temp.3R.vim,
                          temp.X.vim))


#humidity
hum.2L.vim<-read.csv("2L.hum.rel.vim.min.csv")
hist(hum.2L.vim$min.rel.vim)

#2R
hum.2R.vim<-read.csv("2R.hum.rel.vim.min.csv")
hist(hum.2R.vim$min.rel.vim)

#3L
hum.3L.vim<-read.csv("3L.hum.rel.vim.min.csv")
hist(hum.2L.vim$min.rel.vim)

#3R
hum.3R.vim<-read.csv("3R.hum.rel.vim.min.csv")
hist(hum.3R.vim$min.rel.vim)

#X
hum.X.vim<-read.csv("X.hum.rel.vim.min.csv")
hist(hum.X.vim$min.rel.vim)

#combine all of hum
hum<-as.data.frame(rbind(hum.2L.vim,hum.2R.vim,
                          hum.3L.vim,hum.3R.vim,
                          hum.X.vim))

#combine all three variables
all.vars.rel.vim.mins<-as.data.frame(cbind(temp,prec,hum))
#confirm that the SNPs are in the same order for each variable
table(all.vars.rel.vim.mins[,2] == all.vars.rel.vim.mins[,5] & all.vars.rel.vim.mins[,2] == all.vars.rel.vim.mins[,8])
#subset relevant columns, rename appropriately, and write out
all.vars.rel.vim.mins<-all.vars.rel.vim.mins[,c(1,2,3,6,9)]
colnames(all.vars.rel.vim.mins)<-c("chrom","snp","vim.temp","vim.prec","vim.hum")
head(all.vars.rel.vim.mins)
write.csv(all.vars.rel.vim.mins, file="~/Desktop/anoph.phase2/r2vim/all.vars.rel.vim.mins.csv",
          quote=F, row.names = F)


