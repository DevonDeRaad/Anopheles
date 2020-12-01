#make figure 6 searching for patterns in AGAP006026 and AGAP002578
library(ggplot2)
library(gridExtra)

#bring in phase1 allele frequencies
phase1.all.freqs<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.all.freqs.csv")
#make column for id
phase1.all.freqs$id<-paste(phase1.all.freqs$chrom,phase1.all.freqs$POS)

#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","pannual","tmean","hannual")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

#bring in outlier datsets identifying AGAP002578 as a vetted outlier
lfmm.prec<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/lfmm.prec.gene.snp.match.txt")[-1,]
lfmm.prec<-lfmm.prec[lfmm.prec$V3 == "AGAP002578",]
lfmm.hum<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/lfmm.hum.gene.snp.match.txt")[-1,]
lfmm.hum<-lfmm.hum[lfmm.hum$V3 == "AGAP002578",]
baye.hum<-read.table("/Users/devder/Desktop/anoph.phase2/overrep/bayescan.hum.gene.snp.match.txt")[-1,]
baye.hum<-baye.hum[baye.hum$V3 == "AGAP002578",]

#plot allele frequency versus prec for the only SNP associated with precipitation in AGAP002578
frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2L 25194357",3:16]))
frq1<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23203682",3:16]))
d.f<-data.frame(frq=c(frq,frq1), prec=rep(sort.pops.phase1$pannual, times=2),
               SNP=c(rep("2L 25194357", times=14), rep("2R 23203682", times=14)))
rb.prec<-ggplot(d.f, aes(x=prec, y=frq, color=SNP)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  #scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Precipitation",
       y="Allele Frequency")


#plot allele frequency versus hum for the SNPs associated with humidity in AGAP002578
#baye hum: 2R 23168254, 2R 23168258
#lfmm hum: 2R 23168254, 2R 23168258, 2R 23170542, 2R 23170622, 2R 23170902, 2R 23171091, 2R 23181995
frq1<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23168254",3:16]))
frq2<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23168258",3:16]))
frq3<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23170542",3:16]))
frq4<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23170622",3:16]))
frq5<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23170902",3:16]))
frq6<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23171091",3:16]))
frq7<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 23181995",3:16]))
df<-data.frame(frq=c(frq1,frq2,frq3,frq4,frq5,frq6,frq7),
               hum=rep(sort.pops.phase1$hannual, times=7),
               SNP=c(rep("2R 23168254",times=14),rep("2R 23168258",times=14),rep("2R 23170542",times=14),
                     rep("2R 23170622",times=14),rep("2R 23170902",times=14),rep("2R 23171091",times=14),
                     rep("2R 23181995",times=14)))
rb.hum<-ggplot(df, aes(x=hum, y=frq, color=SNP)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Humidity",
       y="Allele Frequency")

g<-arrangeGrob(rb.prec,rb.hum, ncol=2)
ggsave("~/Desktop/anoph.phase2/Fig.6.pdf", g, width=8, height = 2.5, units = "in")





