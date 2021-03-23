#Devon DeRaad
#sticking all of our p-values together for analysis

#read in bayescenv p vals for 2R
baye.2R.hum<-read.table("~/Desktop/anoph.phase2/bayescenv.results/1.hum_fst.txt")[,2]
baye.2R.precip<-read.table("~/Desktop/anoph.phase2/bayescenv.results/1.prec_fst.txt")[,2]
baye.2R.temp<-read.table("~/Desktop/anoph.phase2/bayescenv.results/1.temp_fst.txt")[,2]
baye.2R.rand<-read.table("~/Desktop/anoph.phase2/bayescenv.results/1.rand_fst.txt")[,2]
baye.2R.corr<-read.table("~/Desktop/anoph.phase2/bayescenv.results/1.corr_fst.txt")[,2]

#read in bayescenv p vals for 2L
baye.2L.hum<-read.table("~/Desktop/anoph.phase2/bayescenv.results/2.hum_fst.txt")[,2]
baye.2L.precip<-read.table("~/Desktop/anoph.phase2/bayescenv.results/2.prec_fst.txt")[,2]
baye.2L.temp<-read.table("~/Desktop/anoph.phase2/bayescenv.results/2.temp_fst.txt")[,2]
baye.2L.rand<-read.table("~/Desktop/anoph.phase2/bayescenv.results/2.rand_fst.txt")[,2]
baye.2L.corr<-read.table("~/Desktop/anoph.phase2/bayescenv.results/2.corr_fst.txt")[,2]

#read in bayescenv p vals for 3R
baye.3R.hum<-read.table("~/Desktop/anoph.phase2/bayescenv.results/3.hum_fst.txt")[,2]
baye.3R.precip<-read.table("~/Desktop/anoph.phase2/bayescenv.results/3.prec_fst.txt")[,2]
baye.3R.temp<-read.table("~/Desktop/anoph.phase2/bayescenv.results/3.temp_fst.txt")[,2]
baye.3R.rand<-read.table("~/Desktop/anoph.phase2/bayescenv.results/3.rand_fst.txt")[,2]
baye.3R.corr<-read.table("~/Desktop/anoph.phase2/bayescenv.results/3.corr_fst.txt")[,2]

#read in bayescenv p vals for 3L
baye.3L.hum<-read.table("~/Desktop/anoph.phase2/bayescenv.results/4.hum_fst.txt")[,2]
baye.3L.precip<-read.table("~/Desktop/anoph.phase2/bayescenv.results/4.prec_fst.txt")[,2]
baye.3L.temp<-read.table("~/Desktop/anoph.phase2/bayescenv.results/4.temp_fst.txt")[,2]
baye.3L.rand<-read.table("~/Desktop/anoph.phase2/bayescenv.results/4.rand_fst.txt")[,2]
baye.3L.corr<-read.table("~/Desktop/anoph.phase2/bayescenv.results/4.corr_fst.txt")[,2]

#read in bayescenv p vals for X
baye.X.hum<-read.table("~/Desktop/anoph.phase2/bayescenv.results/5.hum_fst.txt")[,2]
baye.X.precip<-read.table("~/Desktop/anoph.phase2/bayescenv.results/5.prec_fst.txt")[,2]
baye.X.temp<-read.table("~/Desktop/anoph.phase2/bayescenv.results/5.temp_fst.txt")[,2]
baye.X.rand<-read.table("~/Desktop/anoph.phase2/bayescenv.results/5.rand_fst.txt")[,2]
baye.X.corr<-read.table("~/Desktop/anoph.phase2/bayescenv.results/5.corr_fst.txt")[,2]


#make a single vector with qvals for temp
baye.temp<-c(baye.2R.temp,baye.2L.temp,baye.3R.temp,baye.3L.temp,baye.X.temp)

#make a single df with qvals for precip
baye.precip<-c(baye.2R.precip,baye.2L.precip,baye.3R.precip,baye.3L.precip,baye.X.precip)

#make a single df with qvals for temp
baye.hum<-c(baye.2R.hum,baye.2L.hum,baye.3R.hum,baye.3L.hum,baye.X.hum)

#make a single df with qvals for rand
baye.rand<-c(baye.2R.rand,baye.2L.rand,baye.3R.rand,baye.3L.rand,baye.X.rand)

#make a single df with qvals for corr
baye.corr<-c(baye.2R.corr,baye.2L.corr,baye.3R.corr,baye.3L.corr,baye.X.corr)

#read in position info
pos.2R<-read.table("~/Downloads/maf05.ld100.matrix.2R.012.pos")
pos.2L<-read.table("~/Downloads/maf05.ld100.matrix.2L.012.pos")
pos.3R<-read.table("~/Downloads/maf05.ld100.matrix.3R.012.pos")
pos.3L<-read.table("~/Downloads/maf05.ld100.matrix.3L.012.pos")
pos.X<-read.table("~/Downloads/maf05.ld100.matrix.X.012.pos")

chrom<-c(rep("2R", times=nrow(pos.2R)),rep("2L", times=nrow(pos.2L)),
         rep("3R", times=nrow(pos.3R)),rep("3L", times=nrow(pos.3L)),
         rep("X", times=nrow(pos.X)))
pos<-c(pos.2R[,2],pos.2L[,2],pos.3R[,2],pos.3L[,2],pos.X[,2])

#make df
baye.df<-data.frame(chrom=chrom, pos=pos,
                    temp.q=baye.temp, hum.q=baye.hum, prec.q=baye.precip, rand.q=baye.rand, corr.q=baye.corr)

dim(baye.df)
head(baye.df)
write.csv(baye.df, file="~/Desktop/anoph.phase2/bayescenv.results/baye.qvalues.output.csv", quote=FALSE, row.names = FALSE)

