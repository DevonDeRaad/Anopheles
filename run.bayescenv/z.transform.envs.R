#Z transform each env variable
env.sampling<-read.csv("~/Downloads/phase1_all_vairables_last.csv")
unique(env.sampling$latitude)

#transform hum
z.hum.pop.list<-scale(unique(env.sampling$h_mean))

#transform precip
z.precip.pop.list<-scale(unique(env.sampling$p_annual))

#transform temp
z.temp.pop.list<-scale(unique(env.sampling$t_mean))

#transform corr
z.corr.pop.list<-scale(unique(env.sampling$autocorrelated))

#transform rand
z.rand.pop.list<-scale(unique(env.sampling$random))

#write them out
write.table(t(z.temp.pop.list), file ="~/Desktop/anoph.phase2/phase1.ztemp.txt", col.names=F, row.names = F)
write.table(t(z.precip.pop.list), file ="~/Desktop/anoph.phase2/phase1.zprecip.txt", col.names=F, row.names = F)
write.table(t(z.hum.pop.list), file ="~/Desktop/anoph.phase2/phase1.zhum.txt", col.names=F, row.names = F)
write.table(t(z.corr.pop.list), file ="~/Desktop/anoph.phase2/phase1.zcorr.txt", col.names=F, row.names = F)
write.table(t(z.rand.pop.list), file ="~/Desktop/anoph.phase2/phase1.zrand.txt", col.names=F, row.names = F)