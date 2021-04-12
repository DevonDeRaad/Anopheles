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
all.sigs<-read.csv("~/Downloads/all.sigs.csv")
head(all.sigs)
#read in vetted SNPs
vet <- read.csv("~/Desktop/anoph.phase2/vetted.outliers.csv")

#check correlations btwn methods
cor(all.sigs[,3:11])

#make stacked manhattan plots
#make BPcum
nCHR <- length(unique(all.sigs$Chromosome))
all.sigs$BPcum <- NA
s <- 0
nbp <- c()
for (i in c("2R","2L","3R","3L","X")){
  nbp[i] <- max(all.sigs[all.sigs$Chromosome == i,]$Position)
  all.sigs[all.sigs$Chromosome == i,"BPcum"] <- all.sigs[all.sigs$Chromosome == i,"Position"] + s
  s <- s + nbp[i]
}

#check df
head(all.sigs)

#set plotting params
axis.set<-c()
f<-1
for (i in c("2R","2L","3R","3L","X")){
  axis.set[f]<-(max(all.sigs$BPcum[all.sigs$Chromosome == i]) + min(all.sigs$BPcum[all.sigs$Chromosome == i])) / 2
  f<-f+1
}
axis.set<-data.frame(chrom=c("2R","2L","3R","3L","X"),center=axis.set)

#
#hum
#####
#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "lfmm.hum" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot lfmm hum
ylim <- abs(floor(log10(min(all.sigs$lfmm.hum)))) + 2 
#ggplot lfmm hum
lfmm.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(lfmm.hum), 
                         color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(lfmm.hum)), 
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


#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "baye.hum" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot baye hum
ylim <- abs(floor(log10(min(all.sigs$baye.hum)))) + 2 
#ggplot baye hum
baye.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(baye.hum), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(baye.hum)), 
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
#ggtitle("baye humidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "vim.hum" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot vim hum
ylim <- abs(floor(log10(min(all.sigs$vim.hum)))) + 2 
#ggplot vim hum
vim.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.hum, 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=vim.hum), 
             color='red', alpha =.75) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
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
#ggtitle("vim humidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "cors.hum" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot cors hum
ylim <- abs(floor(log10(min(all.sigs$cors.hum)))) + 2 
#ggplot cors hum
cors.hum.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(cors.hum), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(cors.hum)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.0001), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#ggtitle("cors humidity")

#####
#prec
#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "lfmm.prec" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot lfmm prec
ylim <- abs(floor(log10(min(all.sigs$lfmm.prec)))) + 2 
#ggplot lfmm prec
lfmm.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(lfmm.prec), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(lfmm.prec)), 
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
#ggtitle("LFMM precidity")


#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "baye.prec" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot baye prec
ylim <- abs(floor(log10(min(all.sigs$baye.prec)))) + 2 
#ggplot baye prec
baye.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(baye.prec), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(baye.prec)), 
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
#ggtitle("baye precidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "vim.prec" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot vim prec
ylim <- abs(floor(log10(min(all.sigs$vim.prec)))) + 2 
#ggplot vim prec
vim.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.prec, 
                                   color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=vim.prec), 
             color='red', alpha =.75) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
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
#ggtitle("vim precidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "cors.prec" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot cors prec
ylim <- abs(floor(log10(min(all.sigs$cors.prec)))) + 2 
#ggplot cors prec
cors.prec.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(cors.prec), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(cors.prec)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.0001), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#ggtitle("cors precidity")

#####
#temp
#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "lfmm.temp" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot lfmm temp
ylim <- abs(floor(log10(min(all.sigs$lfmm.temp)))) + 2 
#ggplot lfmm temp
lfmm.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(lfmm.temp), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(lfmm.temp)), 
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
#ggtitle("LFMM tempidity")


#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "baye.temp" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot baye temp
ylim <- abs(floor(log10(min(all.sigs$baye.temp)))) + 2 
#ggplot baye temp
baye.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(baye.temp), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(baye.temp)), 
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
#ggtitle("baye tempidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "vim.temp" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot vim temp
ylim <- abs(floor(log10(min(all.sigs$vim.temp)))) + 2 
#ggplot vim temp
vim.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.temp, 
                                   color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=vim.temp), 
             color='red', alpha =.75) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
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
#ggtitle("vim tempidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "cors.temp" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot cors temp
ylim <- abs(floor(log10(min(all.sigs$cors.temp)))) + 2 
#ggplot cors temp
cors.temp.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(cors.temp), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(cors.temp)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.0001), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#ggtitle("cors tempidity")

#####
#random
#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "lfmm.rand" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot lfmm rand
ylim <- abs(floor(log10(min(all.sigs$lfmm.rand)))) + 2 
#ggplot lfmm rand
lfmm.rand.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(lfmm.rand), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(lfmm.rand)), 
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
#ggtitle("LFMM randidity")


#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "baye.rand" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot baye rand
ylim <- abs(floor(log10(min(all.sigs$baye.rand)))) + 2 
#ggplot baye rand
baye.rand.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(baye.rand), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(baye.rand)), 
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
#ggtitle("baye randidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "vim.rand" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot vim rand
ylim <- abs(floor(log10(min(all.sigs$vim.rand)))) + 2 
#ggplot vim rand
vim.rand.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.rand, 
                                   color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=vim.rand), 
             color='red', alpha =.75) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
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
#ggtitle("vim randidity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "cors.rand" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot cors rand
ylim <- abs(floor(log10(min(all.sigs$cors.rand)))) + 2 
#ggplot cors rand
cors.rand.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(cors.rand), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(cors.rand)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.0001), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#ggtitle("cors randidity")

#####
#corr
#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "lfmm.corr" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot lfmm corr
ylim <- abs(floor(log10(min(all.sigs$lfmm.corr)))) + 2 
#ggplot lfmm corr
lfmm.corr.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(lfmm.corr), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(lfmm.corr)), 
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
#ggtitle("LFMM corridity")


#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "baye.corr" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot baye corr
ylim <- abs(floor(log10(min(all.sigs$baye.corr)))) + 2 
#ggplot baye corr
baye.corr.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(baye.corr), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(baye.corr)), 
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
#ggtitle("baye corridity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "vim.corr" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot vim corr
ylim <- abs(floor(log10(min(all.sigs$vim.corr)))) + 2 
#ggplot vim corr
vim.corr.plot<-ggplot(all.sigs, aes(x = BPcum, y = vim.corr, 
                                   color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=vim.corr), 
             color='red', alpha =.75) +
  geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
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
#ggtitle("vim corridity")

#make df of vetted outliers with genomic position
vetted<-vet[vet$test == "cors.corr" & vet$vet == "pass",]
vettee<-all.sigs[paste(all.sigs$Chromosome,all.sigs$Position) %in% paste(vetted$chrom,vetted$pos),]
#plot cors corr
ylim <- abs(floor(log10(min(all.sigs$cors.corr)))) + 2 
#ggplot cors corr
cors.corr.plot<-ggplot(all.sigs, aes(x = BPcum, y = -log10(cors.corr), 
                                    color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_point(data=vettee, 
             aes(x=BPcum, y=-log10(cors.corr)), 
             color='red', alpha =.75) +
  geom_hline(yintercept = -log10(.0001), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("black","darkgrey","black","darkgrey","darkgrey")) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#ggtitle("cors corridity")


#save all 10 plots stacked together twice
dev.off()
#grid.arrange(lfmm.hum.plot, baye.hum.plot, r2vim.hum.plot,
#             lfmm.prec.plot, baye.prec.plot, r2vim.prec.plot,
#             lfmm.temp.plot, baye.temp.plot, r2vim.temp.plot,
#             nrow = 9)

g <- arrangeGrob(lfmm.hum.plot, baye.hum.plot, vim.hum.plot, cors.hum.plot,
                 lfmm.prec.plot, baye.prec.plot, vim.prec.plot, cors.prec.plot,
                 lfmm.temp.plot, baye.temp.plot, nrow = 10)

ggsave(g, filename="~/Desktop/anoph.phase2/all.manhattans.1.png", width = 8.5, height = 11, units ="in")

h <- arrangeGrob(vim.temp.plot, cors.temp.plot, lfmm.corr.plot, baye.corr.plot,
                 vim.corr.plot, cors.corr.plot, lfmm.rand.plot, baye.rand.plot,
                 vim.rand.plot, cors.rand.plot, nrow = 10)

ggsave(h, filename="~/Desktop/anoph.phase2/all.manhattans.2.png", width = 8.5, height = 11, units ="in")






