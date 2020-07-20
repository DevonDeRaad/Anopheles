
#Z transform each env variable
env.sampling<-read.csv("~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")

#transform hum
z.hum.pop.list<-unique(scale(env.sampling$hannual))

#transform precip
z.precip.pop.list<-unique(scale(env.sampling$pannual))

#transform temp
z.temp.pop.list<-unique(scale(env.sampling$tmean))

#write them out
write.table(t(z.temp.pop.list), file ="~/Downloads/z.temp.txt", col.names=F, row.names = F)
write.table(t(z.precip.pop.list), file ="~/Downloads/z.precip.txt", col.names=F, row.names = F)
write.table(t(z.hum.pop.list), file ="~/Downloads/z.hum.txt", col.names=F, row.names = F)
