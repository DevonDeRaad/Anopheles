#plot results from phase1 and verify with phase2 allele freq patterns
library(scattermore)
library(gridExtra)
library(VennDiagram)
library(qvalue)
library(outliers)
library(e1071)

#read in all sig values
all.sigs<-read.csv("~/Downloads/all.sigs.csv")
head(all.sigs)

#bring in phase1 allele frequencies
phase1.all.freqs<-read.csv(file = "~/Desktop/anoph.phase2/phase1.all.freqs.csv")
#make column for id
phase1.all.freqs$id<-paste(phase1.all.freqs$chrom,phase1.all.freqs$POS)
head(phase1.all.freqs)

#bring in phase2 allele frequencies
phase2.all.freqs<-read.csv(file = "~/Desktop/anoph.phase2/phase2.all.freqs.csv")
#make id column
phase2.all.freqs$id<-paste(phase2.all.freqs$chrom,phase2.all.freqs$POS)

#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file="~/Downloads/phase1_all_vairables_last.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","p_annual","t_mean","h_mean","autocorrelated","random")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

#bring in dataframe with locality and environmental data for each individual phase2
phase2.alldata<-read.csv(file = "~/Downloads/phase2_all_vairables_last.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase2.alldata[,c("latitude","longitude","p_annual","t_mean","h_mean","autocorrelated","random")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase2<-metapops[order(metapops$latitude),]
sort.pops.phase2<-sort.pops.phase2[c(3,8,10:28),]

#pull all significant outlier SNP IDs for LFMM humidity
out<- all.sigs[all.sigs$lfmm.hum < .01,c(1,2)]
#add the allele frequency at all 14 sites for each outlier SNP to the dataframe
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
#initialize empty dataframe
z<-data.frame()
#make matching dataframe with the value of the environmental variable used for testing for each SNP
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$h_mean)}
#fix column names to be consistent
colnames(z)<-c(sort.pops.phase1$latitude)
#combine into a single data.frame
lfmm.hum.outliers<-cbind(x,z,test=rep("lfmm.hum", times=nrow(x)))
#each row contains the identifier for an outlier SNP, the allele frequency of that SNP across phase1,
#the environmental values used for testing that SNP, and the name of the test

#lfmm precipitation
out<- all.sigs[all.sigs$lfmm.prec < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$p_annual)}
colnames(z)<-c(sort.pops.phase1$latitude)
lfmm.prec.outliers<-cbind(x,z,test=rep("lfmm.prec", times=nrow(x)))

#lfmm temperature
out<- all.sigs[all.sigs$lfmm.temp < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$t_mean)}
colnames(z)<-c(sort.pops.phase1$latitude)
lfmm.temp.outliers<-cbind(x,z,test=rep("lfmm.temp", times=nrow(x)))

#lfmm spatially autocorrelated variable
out<- all.sigs[all.sigs$lfmm.corr < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$autocorrelated)}
colnames(z)<-c(sort.pops.phase1$latitude)
lfmm.corr.outliers<-cbind(x,z,test=rep("lfmm.corr", times=nrow(x)))

#lfmm random variable
out<- all.sigs[all.sigs$lfmm.rand < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$random)}
colnames(z)<-c(sort.pops.phase1$latitude)
lfmm.rand.outliers<-cbind(x,z,test=rep("lfmm.rand", times=nrow(x)))

#bayescenv prec 
out<- all.sigs[all.sigs$baye.prec < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$p_annual)}
colnames(z)<-c(sort.pops.phase1$latitude)
baye.prec.outliers<-cbind(x,z,test=rep("baye.prec", times=nrow(x)))

#bayescenv hum
out<- all.sigs[all.sigs$baye.hum < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$h_mean)}
colnames(z)<-c(sort.pops.phase1$latitude)
baye.hum.outliers<-cbind(x,z,test=rep("baye.hum", times=nrow(x)))

#bayescenv temp
out<- all.sigs[all.sigs$baye.temp < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$t_mean)}
colnames(z)<-c(sort.pops.phase1$latitude)
baye.temp.outliers<-cbind(x,z,test=rep("baye.temp", times=nrow(x)))

#bayescenv corr
out<- all.sigs[all.sigs$baye.corr < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$autocorrelated)}
colnames(z)<-c(sort.pops.phase1$latitude)
baye.corr.outliers<-cbind(x,z,test=rep("baye.corr", times=nrow(x)))

#bayescenv rand 
out<- all.sigs[all.sigs$baye.rand < .01,c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$random)}
colnames(z)<-c(sort.pops.phase1$latitude)
baye.rand.outliers<-cbind(x,z,test=rep("baye.rand", times=nrow(x)))

#vim prec
out<- all.sigs[all.sigs$vim.prec > 1, c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$p_annual)}
colnames(z)<-c(sort.pops.phase1$latitude)
vim.prec.outliers<-cbind(x,z,test=rep("vim.prec", times=nrow(x)))

#vim hum
out<- all.sigs[all.sigs$vim.hum > 1, c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$h_mean)}
colnames(z)<-c(sort.pops.phase1$latitude)
vim.hum.outliers<-cbind(x,z,test=rep("vim.hum", times=nrow(x)))

#vim temp
out<- all.sigs[all.sigs$vim.temp > 1, c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$t_mean)}
colnames(z)<-c(sort.pops.phase1$latitude)
vim.temp.outliers<-cbind(x,z,test=rep("vim.temp", times=nrow(x)))

#vim corr
out<- all.sigs[all.sigs$vim.corr > 1, c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$autocorrelated)}
colnames(z)<-c(sort.pops.phase1$latitude)
vim.corr.outliers<-cbind(x,z,test=rep("vim.corr", times=nrow(x)))

#vim rand
out<- all.sigs[all.sigs$vim.rand > 1, c(1,2)]
x <- phase1.all.freqs[phase1.all.freqs$id %in% paste(out[,1],out[,2]),c(1:16)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase1$random)}
colnames(z)<-c(sort.pops.phase1$latitude)
vim.rand.outliers<-cbind(x,z,test=rep("vim.rand", times=nrow(x)))

#combine all the dataframes
full<-rbind(lfmm.hum.outliers,lfmm.prec.outliers,lfmm.temp.outliers,lfmm.rand.outliers,lfmm.corr.outliers,
            baye.hum.outliers,baye.prec.outliers,baye.temp.outliers,baye.rand.outliers,baye.corr.outliers,
            vim.hum.outliers,vim.prec.outliers,vim.temp.outliers,vim.rand.outliers,vim.corr.outliers)


#build correlation dataframe for phase1 data
store<-c() #init vector
for (i in 1:nrow(full)){
  store[i]<-cor(as.numeric(as.vector(full[i,17:30])), as.numeric(as.vector(full[i,3:16]))) #perform linear reg
}

#
full$cor<-store

#####step2 repeat this process with phase2 data
####
###
##
#

#pull all significant outlier SNP IDs for LFMM humidity
out<- all.sigs[all.sigs$lfmm.hum < .01,c(1,2)]
#add the allele frequency at all 14 sites for each outlier SNP to the dataframe
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
#initialize empty dataframe
z<-data.frame()
#make matching dataframe with the value of the environmental variable used for testing for each SNP
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$h_mean)}
#fix column names to be consistent
colnames(z)<-c(sort.pops.phase2$latitude)
#combine into a single data.frame
lfmm.hum.outliers<-cbind(x,z,test=rep("lfmm.hum", times=nrow(x)))
#each row contains the identifier for an outlier SNP, the allele frequency of that SNP across phase1,
#the environmental values used for testing that SNP, and the name of the test

#lfmm precipitation
out<- all.sigs[all.sigs$lfmm.prec < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$p_annual)}
colnames(z)<-c(sort.pops.phase2$latitude)
lfmm.prec.outliers<-cbind(x,z,test=rep("lfmm.prec", times=nrow(x)))

#lfmm temperature
out<- all.sigs[all.sigs$lfmm.temp < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$t_mean)}
colnames(z)<-c(sort.pops.phase2$latitude)
lfmm.temp.outliers<-cbind(x,z,test=rep("lfmm.temp", times=nrow(x)))

#lfmm spatially autocorrelated variable
out<- all.sigs[all.sigs$lfmm.corr < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$autocorrelated)}
colnames(z)<-c(sort.pops.phase2$latitude)
lfmm.corr.outliers<-cbind(x,z,test=rep("lfmm.corr", times=nrow(x)))

#lfmm random variable
out<- all.sigs[all.sigs$lfmm.rand < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$random)}
colnames(z)<-c(sort.pops.phase2$latitude)
lfmm.rand.outliers<-cbind(x,z,test=rep("lfmm.rand", times=nrow(x)))

#bayescenv prec 
out<- all.sigs[all.sigs$baye.prec < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$p_annual)}
colnames(z)<-c(sort.pops.phase2$latitude)
baye.prec.outliers<-cbind(x,z,test=rep("baye.prec", times=nrow(x)))

#bayescenv hum
out<- all.sigs[all.sigs$baye.hum < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$h_mean)}
colnames(z)<-c(sort.pops.phase2$latitude)
baye.hum.outliers<-cbind(x,z,test=rep("baye.hum", times=nrow(x)))

#bayescenv temp
out<- all.sigs[all.sigs$baye.temp < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$t_mean)}
colnames(z)<-c(sort.pops.phase2$latitude)
baye.temp.outliers<-cbind(x,z,test=rep("baye.temp", times=nrow(x)))

#bayescenv corr
out<- all.sigs[all.sigs$baye.corr < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$autocorrelated)}
colnames(z)<-c(sort.pops.phase2$latitude)
baye.corr.outliers<-cbind(x,z,test=rep("baye.corr", times=nrow(x)))

#bayescenv rand 
out<- all.sigs[all.sigs$baye.rand < .01,c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$random)}
colnames(z)<-c(sort.pops.phase2$latitude)
baye.rand.outliers<-cbind(x,z,test=rep("baye.rand", times=nrow(x)))

#vim prec
out<- all.sigs[all.sigs$vim.prec > 1, c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$p_annual)}
colnames(z)<-c(sort.pops.phase2$latitude)
vim.prec.outliers<-cbind(x,z,test=rep("vim.prec", times=nrow(x)))

#vim hum
out<- all.sigs[all.sigs$vim.hum > 1, c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$h_mean)}
colnames(z)<-c(sort.pops.phase2$latitude)
vim.hum.outliers<-cbind(x,z,test=rep("vim.hum", times=nrow(x)))

#vim temp
out<- all.sigs[all.sigs$vim.temp > 1, c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$t_mean)}
colnames(z)<-c(sort.pops.phase2$latitude)
vim.temp.outliers<-cbind(x,z,test=rep("vim.temp", times=nrow(x)))

#vim corr
out<- all.sigs[all.sigs$vim.corr > 1, c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$autocorrelated)}
colnames(z)<-c(sort.pops.phase2$latitude)
vim.corr.outliers<-cbind(x,z,test=rep("vim.corr", times=nrow(x)))

#vim rand
out<- all.sigs[all.sigs$vim.rand > 1, c(1,2)]
x <- phase2.all.freqs[phase2.all.freqs$id %in% paste(out[,1],out[,2]),c(1:23)]
z<-data.frame()
for (i in 1:nrow(x)){z<-rbind(z,sort.pops.phase2$random)}
colnames(z)<-c(sort.pops.phase2$latitude)
vim.rand.outliers<-cbind(x,z,test=rep("vim.rand", times=nrow(x)))

#combine all the dataframes
full.2<-rbind(lfmm.hum.outliers,lfmm.prec.outliers,lfmm.temp.outliers,lfmm.rand.outliers,lfmm.corr.outliers,
              baye.hum.outliers,baye.prec.outliers,baye.temp.outliers,baye.rand.outliers,baye.corr.outliers,
              vim.hum.outliers,vim.prec.outliers,vim.temp.outliers,vim.rand.outliers,vim.corr.outliers)

#build correlation dataframe for phase1 data
store<-c() #init vector
for (i in 1:nrow(full.2)){
  store[i]<-cor(as.numeric(as.vector(full.2[i,24:44])), as.numeric(as.vector(full.2[i,3:23]))) #perform linear reg
}

#add as column in full.2
full.2$cor<-store

#subset corframe 1 SNPs found in cor.frame 2
hom.full <- full[paste(full[,1],full[,2],full[,31]) %in% paste(full.2[,1],full.2[,2],full.2[,45]),]

#check to make sure all cells are homologous between these data frames
table(paste(hom.full$chrom,hom.full$pos,hom.full$test) == paste(full.2$chrom,full.2$pos,full.2$test))

#absolute value difference in correlation coefficient
hist(abs(hom.full$cor-full.2$cor))

#
table(hom.full$cor > 0 & full.2$cor < 0 | hom.full$cor < 0 & full.2$cor > 0)

vetted<-data.frame(chrom=hom.full$chrom,
                   pos=hom.full$POS,
                   test=hom.full$test,
                   pass.direction.test=hom.full$cor < 0 & full.2$cor < 0 | hom.full$cor > 0 & full.2$cor > 0,
                   pass.strength.test=abs(hom.full$cor-full.2$cor) < .5)


#calculate overall false positive rate
1-sum(vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/nrow(vetted) #total tested SNPs

#lfmm calculate number of successfully vetted SNPs and the false positive rate for each test
#lfmm hum
sum(vetted$test == "lfmm.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "lfmm.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "lfmm.hum") #false positive rate

#lfmm prec
sum(vetted$test == "lfmm.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "lfmm.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "lfmm.prec") #false positive rate

#lfmm temp
sum(vetted$test == "lfmm.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "lfmm.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "lfmm.temp") #false positive rate

#lfmm rand
sum(vetted$test == "lfmm.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "lfmm.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "lfmm.rand") #false positive rate

#lfmm corr
sum(vetted$test == "lfmm.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "lfmm.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "lfmm.corr") #false positive rate

#baye calculate number of successfully vetted SNPs and the false positive rate for each test
#baye hum
sum(vetted$test == "baye.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "baye.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "baye.hum") #false positive rate

#baye prec
sum(vetted$test == "baye.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "baye.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "baye.prec") #false positive rate

#baye temp
sum(vetted$test == "baye.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "baye.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "baye.temp") #false positive rate

#baye rand
sum(vetted$test == "baye.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "baye.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "baye.rand") #false positive rate

#baye corr
sum(vetted$test == "baye.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "baye.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "baye.corr") #false positive rate


#vim calculate number of successfully vetted SNPs and the false positive rate for each test
#vim hum
sum(vetted$test == "vim.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "vim.hum" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "vim.hum") #false positive rate

#vim prec
sum(vetted$test == "vim.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "vim.prec" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "vim.prec") #false positive rate

#vim temp
sum(vetted$test == "vim.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "vim.temp" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "vim.temp") #false positive rate

#vim rand
sum(vetted$test == "vim.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "vim.rand" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "vim.rand") #false positive rate

#vim corr
sum(vetted$test == "vim.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE) #nSNPs
1-sum(vetted$test == "vim.corr" & vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE)/sum(vetted$test == "vim.corr") #false positive rate

#add column saying whether passed both
vetted$vet<-vetted$pass.direction.test == TRUE & vetted$pass.strength.test == TRUE
vetted$vet[vetted$vet == TRUE]<-"pass"
vetted$vet[vetted$vet == FALSE]<-"fail"

#plot 100 random points to visually confirm that the classification worked
set.seed(62399)
plot.samples<-vetted[sample(nrow(vetted),100),]
par(mfrow=c(3,3))
j<-1 #initialize loop tracker
pb <- txtProgressBar(min = 0, max = length(sig.pos), style = 3) #initialize progress bar
for (i in paste(plot.samples$chrom,plot.samples$pos,plot.samples$test)){
  plot(as.numeric(as.vector(full[paste(full$chrom,full$POS,full$test) == i,c(17:30)])),
       as.numeric(as.vector(full[paste(full$chrom,full$POS,full$test) == i,c(3:16)])), main = i, ylim =c(0,1), xlab=plot.samples$vet[paste(plot.samples$chrom,plot.samples$pos,plot.samples$test) == i],ylab="allele frequency")
  abline(lm(as.numeric(as.vector(full[paste(full$chrom,full$POS,full$test) == i,c(3:16)]))~
              as.numeric(as.vector(full[paste(full$chrom,full$POS,full$test) == i,c(17:30)]))))
  points(as.numeric(as.vector(full.2[paste(full.2$chrom,full.2$POS,full.2$test) == i,c(24:44)])),
         as.numeric(as.vector(full.2[paste(full.2$chrom,full.2$POS,full.2$test) == i,c(3:23)])), col="red")
  abline(lm(as.numeric(as.vector(full.2[paste(full.2$chrom,full.2$POS,full.2$test) == i,c(3:23)]))~
              as.numeric(as.vector(full.2[paste(full.2$chrom,full.2$POS,full.2$test) == i,c(24:44)]))), col="red")
  setTxtProgressBar(pb, j)
  j<-j+1
}

write.csv(vetted, file = "~/Desktop/anoph.phase2/vetted.outliers.csv", quote = F, row.names = F)


