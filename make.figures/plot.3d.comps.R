library(scattermore)
library(gridExtra)
library(VennDiagram)
library(qvalue)
library(outliers)
library(e1071)
library(dplyr)
library(ggplot2)
library(rgl)

#read in all sig values
all.sigs<-read.table(file="~/Desktop/anoph.phase2/all.sigs.adjusted.txt", header=T)
head(all.sigs)

#bring in vetted outlier SNP list by variable and method
#hum
#lfmm
lfmm.hum <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.hum.env.csv")
#baye
baye.hum <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.hum.env.csv")
#r2vim
r2vim.hum <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.hum.env.csv")
colnames(r2vim.hum)<-c("","chrom","pos")

#prec
#lfmm
lfmm.prec <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.prec.env.csv")
#baye
baye.prec <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.prec.env.csv")
#r2vim
r2vim.prec <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.prec.env.csv")
colnames(r2vim.prec)<-c("","chrom","pos")

#temp
#lfmm
lfmm.temp <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.temp.env.csv")
#baye
baye.temp <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.temp.env.csv")
#r2vim
r2vim.temp <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.temp.env.csv")
colnames(r2vim.temp)<-c("","chrom","pos")

#build venn hum
venn.diagram(x = list(paste(lfmm.hum$chrom,lfmm.hum$pos),
                      paste(baye.hum$chrom,baye.hum$pos),
                      paste(r2vim.hum$chrom,r2vim.hum$pos)),
             category.names = c("LFMM", "BayeScEnv", "r2VIM"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.3,
             cat.just=list(c(-.3,-.35),c(.75,-.35),c(.05,.29)),
             filename = "~/Downloads/hum.png")

#build venn prec
venn.diagram(x = list(paste(lfmm.prec$chrom,lfmm.prec$pos),
                      paste(baye.prec$chrom,baye.prec$pos),
                      paste(r2vim.prec$chrom,r2vim.prec$pos)),
             category.names = c("LFMM", "BayeScEnv", "r2VIM"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.3,
             cat.just=list(c(-.3,-.35),c(.75,-.35),c(.05,.29)),
             filename = "~/Downloads/prec.png")

#build venn temp
venn.diagram(x = list(paste(lfmm.temp$chrom,lfmm.temp$pos),
                      paste(baye.temp$chrom,baye.temp$pos),
                      paste(r2vim.temp$chrom,r2vim.temp$pos)),
             category.names = c("LFMM", "BayeScEnv", "r2VIM"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.3,
             cat.just=list(c(-.3,-.35),c(.75,-.35),c(.05,.29)),
             filename = "~/Downloads/temp.png")


#build 3D plots
#humidity
all.sigs$bayescenv.vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.hum$chrom,lfmm.hum$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.hum$chrom,baye.hum$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.hum$chrom,r2vim.hum$pos)]<-"red"
table(all.sigs$bayescenv.vett)

plot3d(x=-log10(all.sigs$p_value_hannual),
       y=-log10(all.sigs$hum.q),
       z=all.sigs$vim.hum,
       col = all.sigs$bayescenv.vett,
       xlab="", ylab="", zlab="")
#rgl.snapshot('~/Desktop/anoph.phase2/hum.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/hum.3D.nolabs.pdf', fmt = 'pdf')

#precipitation
all.sigs$bayescenv.vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.prec$chrom,lfmm.prec$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.prec$chrom,baye.prec$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.prec$chrom,r2vim.prec$pos)]<-"red"
table(all.sigs$bayescenv.vett)

plot3d(x=-log10(all.sigs$p_value_pannual),
       y=-log10(all.sigs$prec.q),
       z=all.sigs$vim.prec,
       col = all.sigs$bayescenv.vett,
       xlab="", ylab="", zlab="", lit=FALSE)
#rgl.snapshot('~/Desktop/anoph.phase2/prec.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/prec.3D.nolabs.pdf', fmt = 'pdf')

#temp
all.sigs$bayescenv.vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.temp$chrom,lfmm.temp$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.temp$chrom,baye.temp$pos)]<-"red"
all.sigs$bayescenv.vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.temp$chrom,r2vim.temp$pos)]<-"red"
table(all.sigs$bayescenv.vett)

plot3d(x=-log10(all.sigs$p_value_tmean),
       y=-log10(all.sigs$temp.q),
       z=all.sigs$vim.temp,
       col = all.sigs$bayescenv.vett,
       xlab="", ylab="", zlab="")
#rgl.snapshot('~/Desktop/anoph.phase2/temp.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/temp.3D.nolabs.pdf', fmt = 'pdf')

#make histograms of each variable
#bring in dataframe with locality and environmental data for each individual phase1
phase1.alldata<-read.csv(file = "~/Desktop/anoph.3.march.2020/phase1.allvariables.csv")
#subset to only the variables we need (lat,long,hum,temp,precip)
meta<-phase1.alldata[,c("latitude","longitude","pannual","tmean","hannual")]
#subset only the unique lat/longs
metapops<-unique(meta)

#plot
ggplot(metapops, aes(x=pannual)) + 
  geom_histogram(color="black", fill="white", binwidth=500)+
  labs(x="precipitation (mm/year)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=9))+
ggsave("~/Desktop/anoph.phase2/prec.hist.pdf", width=2.5, height=2.5, units="in")

ggplot(metapops, aes(x=hannual)) + 
  geom_histogram(color="black", fill="white", binwidth=100)+
  labs(x="humidity (100000*kg of water/kg of air)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=9))+
ggsave("~/Desktop/anoph.phase2/hum.hist.pdf", width=2.5, height=2.5, units="in")

ggplot(metapops, aes(x=tmean)) + 
  geom_histogram(color="black", fill="white", binwidth=10)+
  labs(x="mean annual temperature (celsius*100)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=9))+
ggsave("~/Desktop/anoph.phase2/temp.hist.pdf", width=2.5, height=2.5, units="in")



#Now plot the same venn diagrams and 3D plots but by method rather than env variable, 
#to search for signature of genes acting in multiple dimensions

#build venn lfmm
venn.diagram(x = list(paste(lfmm.hum$chrom,lfmm.hum$pos),
                      paste(lfmm.prec$chrom,lfmm.prec$pos),
                      paste(lfmm.temp$chrom,lfmm.temp$pos)),
             category.names = c("humidity", "precipitation", "temperature"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.2,
             cat.just=list(c(-.1,-.35),c(.77,-.35),c(.4,.1)),
             filename = "~/Downloads/lfmm.png")

#build venn bayescenv
venn.diagram(x = list(paste(baye.hum$chrom,baye.hum$pos),
                      paste(baye.prec$chrom,baye.prec$pos),
                      paste(baye.temp$chrom,baye.temp$pos)),
             category.names = c("humidity", "precipitation", "temperature"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.2,
             cat.just=list(c(-.1,-.35),c(.77,-.35),c(.4,.1)),
             filename = "~/Downloads/baye.png")

#build venn r2VIM
venn.diagram(x = list(paste(r2vim.hum$chrom,r2vim.hum$pos),
                      paste(r2vim.prec$chrom,r2vim.prec$pos),
                      paste(r2vim.temp$chrom,r2vim.temp$pos)),
             category.names = c("humidity", "precipitation", "temperature"),
             main = "",
             euler.d=FALSE,
             scaled=FALSE,
             width=2.5, height=2.5, units = "in",
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             cex = 1.3,
             cat.cex=1.2,
             cat.just=list(c(-.1,-.35),c(.77,-.35),c(.4,.1)),
             filename = "~/Downloads/r2vim.png")

#build 3D plots
#lfmm
all.sigs$vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.hum$chrom,lfmm.hum$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.prec$chrom,lfmm.prec$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(lfmm.temp$chrom,lfmm.temp$pos)]<-"red"
table(all.sigs$vett)

plot3d(x=-log10(all.sigs$p_value_hannual),
       y=-log10(all.sigs$p_value_pannual),
       z=-log10(all.sigs$p_value_tmean),
       col = all.sigs$vett,
       xlab="", ylab="", zlab="")
#rgl.snapshot('~/Desktop/anoph.phase2/lfmm.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/lfmm.3D.nolabs.pdf', fmt = 'pdf')


#baye
all.sigs$vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.hum$chrom,baye.hum$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.prec$chrom,baye.prec$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(baye.temp$chrom,baye.temp$pos)]<-"red"
table(all.sigs$vett)

plot3d(x=-log10(all.sigs$hum.q),
       y=-log10(all.sigs$prec.q),
       z=-log10(all.sigs$temp.q),
       col = all.sigs$vett,
       xlab="", ylab="", zlab="")
#rgl.snapshot('~/Desktop/anoph.phase2/baye.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/baye.3D.nolabs.pdf', fmt = 'pdf')

#temp
all.sigs$vett<-c(rep("black", times=nrow(all.sigs)))
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.hum$chrom,r2vim.hum$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.prec$chrom,r2vim.prec$pos)]<-"red"
all.sigs$vett[paste(all.sigs$chrom,all.sigs$pos) %in% paste(r2vim.temp$chrom,r2vim.temp$pos)]<-"red"
table(all.sigs$vett)

plot3d(x=all.sigs$vim.hum,
       y=all.sigs$vim.prec,
       z=all.sigs$vim.temp,
       col = all.sigs$vett,
       xlab="", ylab="", zlab="")
#rgl.snapshot('~/Desktop/anoph.phase2/r2vim.3D.nolabs.png', fmt = 'png')
rgl.postscript('~/Desktop/anoph.phase2/r2vim.3D.nolabs.pdf', fmt = 'pdf')




