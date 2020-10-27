
#script to combine significance measures from all 3 methods
#lfmm p values come from Marlon and get converted here
#bayescenv q values are combined in the script '/Users/devder/Desktop/anoph.phase2/bayescenv.results/combine.bayescenv.qvals.R'
#r2vim vims are combined in the script '~/Desktop/anoph.phase2/r2vim/combine.rel.vim.mins.R'
#final output file for ML classification output here:

library(qvalue)
library(dplyr)
library(ggplot2)
library(gridExtra)

#read in bayescan qvals
bayescan<-read.csv("~/Desktop/anoph.phase2/bayescenv.results/baye.qvalues.output.csv")
bayescan[1:10,]
table(bayescan$temp.q)[1:10]
#convert 0 qvals to 2e-05
bayescan[bayescan == 0]<-2e-05
table(bayescan$temp.q)[1:10]


#read in lfmm pvals
lfmm<-read.csv("~/Desktop/anoph.phase2/lfmm/whole_genome_LFMM_pvalues.csv")
lfmm[1:10,]

#convert lfmm p values to q values
#visual inspection to make sure this is working the way we hope
q<-qvalue(lfmm$p_value_pannual)
hist(lfmm$p_value_pannual)
hist(lfmm$p_value_pannual[lfmm$p_value_pannual < .01])
hist(q$qvalues)
hist(q$qvalues[q$qvalues < .01])

#looks good
pannual<-qvalue(lfmm$p_value_pannual)
hannual<-qvalue(lfmm$p_value_hannual)
tmean<-qvalue(lfmm$p_value_tmean)

#convert
lfmm$p_value_pannual<-pannual$qvalues
lfmm$p_value_hannual<-hannual$qvalues
lfmm$p_value_tmean<-tmean$qvalues

colnames(lfmm)[1:2]<-c("chrom","pos")
lfmm[1:10,]

#read in r2vim vim values
r2vim<-read.csv("~/Desktop/anoph.phase2/r2vim/all.vars.rel.vim.mins.csv")
head(r2vim)
colnames(r2vim)[2]<-"pos"
dim(r2vim)
dim(lfmm)
dim(bayescan)
#merge
df<-merge(bayescan, lfmm, by=c("chrom","pos"))
head(df)
table(df$pos %in% r2vim$pos)
dff<-merge(df, r2vim, by=c("chrom","pos"))
head(dff)
#write it out
write.table(dff, "~/Desktop/anoph.phase2/all.sigs.adjusted.txt", quote = F, row.names = F)

