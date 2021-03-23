
#script to manually take the min of each SNP across runs that were arrayed due to 
#multithreading issues, and create a single csv with genome wide importance values for each env variable 
setwd("~/Desktop/anoph.phase2/r2vim/")

#read in humidity

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


#read in prec
#prec
prec.2L.vim<-read.csv("2L.prec.rel.vim.min.csv")
hist(prec.2L.vim$min.rel.vim)

#2R
prec.2R.vim<-read.csv("2R.prec.rel.vim.min.csv")
hist(prec.2R.vim$min.rel.vim)

#3L
prec.3L.vim<-read.csv("3L.prec.rel.vim.min.csv")
hist(prec.2L.vim$min.rel.vim)

#3R
prec.3R.vim<-read.csv("3R.prec.rel.vim.min.csv")
hist(prec.3R.vim$min.rel.vim)

#X
prec.X.vim<-read.csv("X.prec.rel.vim.min.csv")
hist(prec.X.vim$min.rel.vim)

#combine all of prec
prec<-as.data.frame(rbind(prec.2L.vim,prec.2R.vim,
                         prec.3L.vim,prec.3R.vim,
                         prec.X.vim))

#read in temp
#temp
temp.2L.vim<-read.csv("2L.temp.rel.vim.min.csv")
hist(temp.2L.vim$min.rel.vim)

#2R
temp.2R.vim<-read.csv("2R.temp.rel.vim.min.csv")
hist(temp.2R.vim$min.rel.vim)

#3L
temp.3L.vim<-read.csv("3L.temp.rel.vim.min.csv")
hist(temp.2L.vim$min.rel.vim)

#3R
temp.3R.vim<-read.csv("3R.temp.rel.vim.min.csv")
hist(temp.3R.vim$min.rel.vim)

#X
temp.X.vim<-read.csv("X.temp.rel.vim.min.csv")
hist(temp.X.vim$min.rel.vim)

#combine all of temp
temp<-as.data.frame(rbind(temp.2L.vim,temp.2R.vim,
                          temp.3L.vim,temp.3R.vim,
                          temp.X.vim))

#read in rand
#rand
rand.2L.vim<-read.csv("2L.rand.rel.vim.min.csv")
hist(rand.2L.vim$min.rel.vim)

#2R
rand.2R.vim<-read.csv("2R.rand.rel.vim.min.csv")
hist(rand.2R.vim$min.rel.vim)

#3L
rand.3L.vim<-read.csv("3L.rand.rel.vim.min.csv")
hist(rand.2L.vim$min.rel.vim)

#3R
rand.3R.vim<-read.csv("3R.rand.rel.vim.min.csv")
hist(rand.3R.vim$min.rel.vim)

#X
rand.X.vim<-read.csv("X.rand.rel.vim.min.csv")
hist(rand.X.vim$min.rel.vim)

#combine all of rand
rand<-as.data.frame(rbind(rand.2L.vim,rand.2R.vim,
                          rand.3L.vim,rand.3R.vim,
                          rand.X.vim))

#read in corr
#corr
corr.2L.vim<-read.csv("2L.corr.rel.vim.min.csv")
hist(corr.2L.vim$min.rel.vim)

#2R
corr.2R.vim<-read.csv("2R.corr.rel.vim.min.csv")
hist(corr.2R.vim$min.rel.vim)

#3L
corr.3L.vim<-read.csv("3L.corr.rel.vim.min.csv")
hist(corr.2L.vim$min.rel.vim)

#3R
corr.3R.vim<-read.csv("3R.corr.rel.vim.min.csv")
hist(corr.3R.vim$min.rel.vim)

#X
corr.X.vim<-read.csv("X.corr.rel.vim.min.csv")
hist(corr.X.vim$min.rel.vim)

#combine all of corr
corr<-as.data.frame(rbind(corr.2L.vim,corr.2R.vim,
                          corr.3L.vim,corr.3R.vim,
                          corr.X.vim))


#combine all variables
all.vars.rel.vim.mins<-as.data.frame(cbind(temp,prec,hum,rand,corr))
#confirm that the SNPs are in the same order for each variable
table(all.vars.rel.vim.mins[,2] == all.vars.rel.vim.mins[,5] & all.vars.rel.vim.mins[,2] == all.vars.rel.vim.mins[,8])
#subset relevant columns, rename appropriately, and write out
all.vars.rel.vim.mins<-all.vars.rel.vim.mins[,c(1,2,3,6,9,12,15)]
colnames(all.vars.rel.vim.mins)<-c("chrom","snp","vim.temp","vim.prec","vim.hum","vim.rand","vim.corr")
head(all.vars.rel.vim.mins)
write.csv(all.vars.rel.vim.mins, file="~/Desktop/anoph.phase2/r2vim/all.vars.rel.vim.mins.csv",
          quote=F, row.names = F)


