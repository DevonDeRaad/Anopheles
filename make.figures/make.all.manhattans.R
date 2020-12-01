#devtools::install_github('exaexa/scattermore')
library(scattermore)
library(gridExtra)
library(VennDiagram)
library(qvalue)
library(outliers)
library(e1071)
library(dplyr)
library(ggplot2)

#read in all sig values
all.sigs<-read.table(file="~/Desktop/anoph.phase2/all.sigs.adjusted.txt", header=T)
head(all.sigs)
#check correlations btwn methods
cor(all.sigs$temp.q,all.sigs$p_value_tmean)
cor(all.sigs[,3:11])

#make stacked manhattan plots
#make BPcum
nCHR <- length(unique(all.sigs$chrom))
all.sigs$BPcum <- NA
s <- 0
nbp <- c()
for (i in c("2R","2L","3R","3L","X")){
  nbp[i] <- max(all.sigs[all.sigs$chrom == i,]$pos)
  all.sigs[all.sigs$chrom == i,"BPcum"] <- all.sigs[all.sigs$chrom == i,"pos"] + s
  s <- s + nbp[i]
}

#check df
head(all.sigs)

#set plotting params
axis.set<-c()
f<-1
for (i in c("2R","2L","3R","3L","X")){
  axis.set[f]<-(max(all.sigs$BPcum[all.sigs$chrom == i]) + min(all.sigs$BPcum[all.sigs$chrom == i])) / 2
  f<-f+1
}
axis.set<-data.frame(chrom=c("2R","2L","3R","3L","X"),center=axis.set)

#
#hum
#

#make df of vetted outliers with genomic position
lfmm.hum <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.hum.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(lfmm.hum$chrom,lfmm.hum$pos),]
#plot lfmm hum
ylim <- abs(floor(log10(min(all.sigs$p_value_hannual)))) + 2 
#ggplot lfmm hum
lfmm.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(p_value_hannual), 
                         color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(p_value_hannual)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("LFMM humidity")

#
#make df of vetted outliers with genomic position
baye.hum <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.hum.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(baye.hum$chrom,baye.hum$pos),]
#plot baye hum
ylim <- abs(floor(log10(min(all.sigs$hum.q)))) + 2 
#ggplot baye hum
baye.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(hum.q), 
                     color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(hum.q)), 
             color='red') +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("BayeScEnv humidity")

#make df of vetted outliers with genomic position
r2vim.hum <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.hum.env.csv")
colnames(r2vim.hum)<-c("","chrom","pos")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(r2vim.hum$chrom,r2vim.hum$pos),]
#plot r2vim hum
ylim <- max(all.sigs$vim.hum) + 2 
#ggplot r2vim hum
r2vim.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.hum, 
                     color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=vim.hum), 
             color='red') +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "relative VIM") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("r2vim humidity")


#
#prec
#

#make df of vetted outliers with genomic position
lfmm.prec <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.prec.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(lfmm.prec$chrom,lfmm.prec$pos),]
#plot lfmm prec
ylim <- abs(floor(log10(min(all.sigs$p_value_hannual)))) + 2 
#ggplot lfmm prec
lfmm.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(p_value_pannual), 
                                    color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(p_value_pannual)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("LFMM precipitation")

#
#make df of vetted outliers with genomic position
baye.prec <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.prec.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(baye.prec$chrom,baye.prec$pos),]
#plot baye prec
ylim <- abs(floor(log10(min(all.sigs$prec.q)))) + 2 
#ggplot baye prec
baye.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(prec.q), 
                                    color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(prec.q)), 
             color='red') +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("BayeScEnv precipitation")

#make df of vetted outliers with genomic position
r2vim.prec <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.prec.env.csv")
colnames(r2vim.prec)<-c("","chrom","pos")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(r2vim.prec$chrom,r2vim.prec$pos),]
#plot r2vim prec
ylim <- max(all.sigs$vim.prec) + 2 
#ggplot r2vim prec
r2vim.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.prec, 
                                     color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=vim.prec), 
             color='red') +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "relative VIM") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("r2vim precipitation")


#
#temp
#

#make df of vetted outliers with genomic position
lfmm.temp <- read.csv("~/Downloads/ml.outputs/vetted.lfmm.temp.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(lfmm.temp$chrom,lfmm.temp$pos),]
#plot lfmm temp
ylim <- abs(floor(log10(min(all.sigs$p_value_hannual)))) + 2 
#ggplot lfmm temp
lfmm.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(p_value_tmean), 
                                     color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(p_value_tmean)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("LFMM temperature")

#
#make df of vetted outliers with genomic position
baye.temp <- read.csv("~/Downloads/ml.outputs/vetted.bayescenv.temp.env.csv")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(baye.temp$chrom,baye.temp$pos),]
#plot baye temp
ylim <- abs(floor(log10(min(all.sigs$temp.q)))) + 2 
#ggplot baye temp
baye.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(temp.q), 
                                     color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=-log10(temp.q)), 
             color='red') +
  geom_hline(yintercept = -log10(.01), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(q)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("BayeScEnv temperature")

#make df of vetted outliers with genomic position
r2vim.temp <- read.csv("~/Downloads/ml.outputs/vetted.r2vim.temp.env.csv")
colnames(r2vim.temp)<-c("","chrom","pos")
vet<-all.sigs[paste0(all.sigs$chrom,all.sigs$pos) %in% paste0(r2vim.temp$chrom,r2vim.temp$pos),]
#plot r2vim temp
ylim <- max(all.sigs$vim.temp) + 2 
#ggplot r2vim temp
r2vim.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.temp, 
                                      color = as.factor(chrom))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vet, 
             aes(x=BPcum, y=vim.temp), 
             color='red') +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "relative VIM") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
  #ggtitle("r2vim temperature")


#save all 9 plots stacked together
dev.off()
#grid.arrange(lfmm.hum.plot, baye.hum.plot, r2vim.hum.plot,
#             lfmm.prec.plot, baye.prec.plot, r2vim.prec.plot,
#             lfmm.temp.plot, baye.temp.plot, r2vim.temp.plot,
#             nrow = 9)

g <- arrangeGrob(lfmm.hum.plot, baye.hum.plot, r2vim.hum.plot,
             lfmm.prec.plot, baye.prec.plot, r2vim.prec.plot,
             lfmm.temp.plot, baye.temp.plot, r2vim.temp.plot,
             nrow = 9)

ggsave(g, filename="~/Desktop/anoph.phase2/all.manhattans.notitles.png", width = 8.5, height = 11, units ="in")
ggsave(g, filename="~/Desktop/anoph.phase2/all.manhattans.notitles.pdf", width = 8.5, height = 11, units ="in")


