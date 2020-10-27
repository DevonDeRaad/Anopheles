#Devon DeRaad
#sticking all of our p-values together for analysis

#read in bayescenv p vals for 2R
baye.2R.hum<-read.table("~/Desktop/bayescenv.results/1.hum_fst.txt")[,2]
baye.2R.precip<-read.table("~/Desktop/bayescenv.results/1.precip_fst.txt")[,2]
baye.2R.temp<-read.table("~/Desktop/bayescenv.results/1.temp_fst.txt")[,2]

#read in bayescenv p vals for 2L
baye.2L.hum<-read.table("~/Desktop/bayescenv.results/2.hum_fst.txt")[,2]
baye.2L.precip<-read.table("~/Desktop/bayescenv.results/2.precip_fst.txt")[,2]
baye.2L.temp<-read.table("~/Desktop/bayescenv.results/2.temp_fst.txt")[,2]

#read in bayescenv p vals for 3R
baye.3R.hum<-read.table("~/Desktop/bayescenv.results/3.hum_fst.txt")[,2]
baye.3R.precip<-read.table("~/Desktop/bayescenv.results/3.precip_fst.txt")[,2]
baye.3R.temp<-read.table("~/Desktop/bayescenv.results/3.temp_fst.txt")[,2]

#read in bayescenv p vals for 3L
baye.3L.hum<-read.table("~/Desktop/bayescenv.results/4.hum_fst.txt")[,2]
baye.3L.precip<-read.table("~/Desktop/bayescenv.results/4.precip_fst.txt")[,2]
baye.3L.temp<-read.table("~/Desktop/bayescenv.results/4.temp_fst.txt")[,2]

#read in bayescenv p vals for X
baye.X.hum<-read.table("~/Desktop/bayescenv.results/5.hum_fst.txt")[,2]
baye.X.precip<-read.table("~/Desktop/bayescenv.results/5.precip_fst.txt")[,2]
baye.X.temp<-read.table("~/Desktop/bayescenv.results/5.temp_fst.txt")[,2]


#make a single vector with qvals for temp
baye.temp<-c(baye.2R.temp,baye.2L.temp,baye.3R.temp,baye.3L.temp,baye.X.temp)

#make a single df with qvals for precip
baye.precip<-c(baye.2R.precip,baye.2L.precip,baye.3R.precip,baye.3L.precip,baye.X.precip)

#make a single df with qvals for temp
baye.hum<-c(baye.2R.hum,baye.2L.hum,baye.3R.hum,baye.3L.hum,baye.X.hum)

#read in position info
pos.2R<-read.table("~/Desktop/anoph.3.march.2020/ag1000g.phase1.ar3.pass.biallelic.maf05.2R.geno.matrix.012.pos")
pos.2L<-read.table("~/Desktop/anoph.3.march.2020/ag1000g.phase1.ar3.pass.biallelic.maf05.2L.geno.matrix.012.pos")
pos.3R<-read.table("~/Desktop/anoph.3.march.2020/ag1000g.phase1.ar3.pass.biallelic.maf05.3R.geno.matrix.012.pos")
pos.3L<-read.table("~/Desktop/anoph.3.march.2020/ag1000g.phase1.ar3.pass.biallelic.maf05.3L.geno.matrix.012.pos")
pos.X<-read.table("~/Desktop/anoph.3.march.2020/ag1000g.phase1.ar3.pass.biallelic.maf05.X.geno.matrix.012.pos")

chrom<-c(rep("2R", times=nrow(pos.2R)),rep("2L", times=nrow(pos.2L)),
         rep("3R", times=nrow(pos.3R)),rep("3L", times=nrow(pos.3L)),
         rep("X", times=nrow(pos.X)))
pos<-c(pos.2R[,2],pos.2L[,2],pos.3R[,2],pos.3L[,2],pos.X[,2])

#make df
baye.df<-data.frame(chrom=chrom, pos=pos,
                    temp.q=baye.temp, hum.q=baye.hum, prec.q=baye.precip)

write.csv(baye.df, file="~/Desktop/bayescenv.results/baye.qvalues.output.csv", quote=FALSE, row.names = FALSE)

