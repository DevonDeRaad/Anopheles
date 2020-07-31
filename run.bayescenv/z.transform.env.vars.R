
#Z transform each env variable
env.sampling<-read.csv("~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")

#transform hum
z.hum.pop.list<-scale(unique(env.sampling$hannual))

#transform precip
z.precip.pop.list<-scale(unique(env.sampling$pannual))

#transform temp
z.temp.pop.list<-scale(unique(env.sampling$tmean))

#write them out
write.table(t(z.temp.pop.list), file ="~/Desktop/anoph.phase2/phase1.ztemp.txt", col.names=F, row.names = F)
write.table(t(z.precip.pop.list), file ="~/Desktop/anoph.phase2/phase1.zprecip.txt", col.names=F, row.names = F)
write.table(t(z.hum.pop.list), file ="~/Desktop/anoph.phase2/phase1.zhum.txt", col.names=F, row.names = F)
