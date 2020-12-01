
#script to plot Figure S3
library(ggplot2)
library(gridExtra)
#bring in the hand classified SNP set for environment vs outlier associated with lfmm humidity
lfmm.hum.vet1<-read.csv("~/Downloads/ml.outputs/lfmm.hum.cor.training.csv")

#bring in hand classifications for matching phase 2, lfmm hum
lfmm.hum.vet2<-read.csv("~/Downloads/lfmm.hum.validate.env.csv")

#bring in phase1 allele frequencies
phase1.all.freqs<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.all.freqs.csv")
#make column for id
phase1.all.freqs$id<-paste(phase1.all.freqs$chrom,phase1.all.freqs$POS)

#bring in phase2 allele frequencies
phase2.all.freqs<-read.csv(file = "~/Desktop/anoph.phase2/phase2.all.freqs.csv")
#make id column
phase2.all.freqs$id<-paste(phase2.all.freqs$chrom,phase2.all.freqs$POS)

#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","pannual","tmean","hannual")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase1 <-metapops[order(metapops$latitude),]

#bring in dataframe with locality and environmental data for each individual phase2
phase2.alldata<-read.csv(file = "~/Downloads/phase2_localities_all_vairables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase2.alldata[,c("latitude","longitude","pannual","tmean","hannual")]
#subset only the unique lat/longs
metapops<-unique(meta)
#sort by latitude to match the order of the allele frequency files
sort.pops.phase2 <-metapops[order(metapops$latitude),]

#Vetting step 1
#plot allele frequency versus humidity for a SNP that was tossed due to lack of general signal
frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "2R 49154088",3:16]))
df<-data.frame(frq=frq, hum=sort.pops.phase1$hannual)
fail<-ggplot(df, aes(x=hum, y=frq)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Humidity",
       y="Allele Frequency")+
  ggtitle("SNP 2R 49154088")

#plot allele frequency versus humidity for a SNP that was retained as environmentally associated
frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "3R 28605030",3:16]))
df<-data.frame(frq=frq, hum=sort.pops.phase1$hannual)
kept<-ggplot(df, aes(x=hum, y=frq)) + 
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Humidity",
       y="Allele Frequency")+
  ggtitle("SNP 3R 28605030")

g<-arrangeGrob(fail, kept, ncol=2)
ggsave("~/Desktop/anoph.phase2/Fig.S4.pdf", g, width=5, height = 2.5, units = "in")


#vetting step 2
#plot allele frequency versus humidity for a SNP that was tossed due to lack of concordance
frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "3L 5904736",3:16]))
frq2<-as.numeric(as.vector(phase2.all.freqs[phase2.all.freqs$id == "3L 5904736",3:23]))
df<-data.frame(frq=c(frq,frq2), phase=c(rep("Phase 1",times=14), rep("Phase 2",times=21)),hum=c(sort.pops.phase1$hannual,sort.pops.phase2$hannual))
fail<-ggplot(df, aes(x=hum, y=frq, colour=phase, fill=phase)) +
  scale_color_manual(values=c("black","red"))+
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Humidity",
       y="Allele Frequency")+
  ggtitle("SNP 3L 5904736")+
  theme(legend.title = element_blank())

#plot allele frequency versus humidity for a SNP that was retained
frq<-as.numeric(as.vector(phase1.all.freqs[phase1.all.freqs$id == "X 15477814",3:16]))
frq2<-as.numeric(as.vector(phase2.all.freqs[phase2.all.freqs$id == "X 15477814",3:23]))
df<-data.frame(frq=c(frq,frq2), phase=c(rep("Phase 1",times=14), rep("Phase 2",times=21)),hum=c(sort.pops.phase1$hannual,sort.pops.phase2$hannual))
kept<-ggplot(df, aes(x=hum, y=frq, colour=phase, fill=phase)) +
  scale_color_manual(values=c("black","red"))+
  geom_point(alpha=.4, size=2.5)+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1200, 1800, by = 200))+
  labs(x="Humidity",
       y="Allele Frequency")+
  ggtitle("SNP X 15477814")+
  theme(legend.title = element_blank())

g<-arrangeGrob(fail, kept, ncol=2)
ggsave("~/Desktop/anoph.phase2/Fig.S3.pdf", g, width=8, height = 2.5, units = "in")

