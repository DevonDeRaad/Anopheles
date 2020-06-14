##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: LFMM precipitation, humidity, temperature 
# Author(s): Devon A DeRaad, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 04/03/2020 (dd/mm/yyyy)
##################################

#packages
#library(lfmm)
## to load lfmm functions because lfmm package installation was not possible
devtools::load_all("C:/Users/m538c005/Documents/R/Extra_lfmm/lfmm-master")
##

#bring in marlons dataset
phase1.env <- read.csv("FST_LFMM_BayeScan/individuals_variables_localities.csv")
phase2.env <- read.csv("FST_LFMM_BayeScan/localities_all_vairables_categorical.csv")

#subset from all phase2 variables to only phase1 individuals
phase1.all.env <- match_df(phase2.env, phase1.env, on = "ox_code")
#this df will give us the ordered precip and hum vectors needed to run lfmm for all phase1 individuals

#################
#Running lfmm2 with precipitation

#variable
X <- phase1.all.env$pannual

#2L
#read in reduced data matrix from phase1
#now we simply read in our matrix
Y.2L <- data.matrix(read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.2L.geno.matrix.012", 
                             header = FALSE, row.names = 1, sep = "\t"))

#test Y to see how many PCs to retain
y.pc <- prcomp(Y.2L) #will run for 15 minutes
spc.y2L <- summary(y.pc)
plot(spc.y2L$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
abline(v = 4, col = "red", lty = 2)

#The ridge_lfmm function is applied to Y (genotype) and X (phenotypes, SD = 1).
## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = Y.2L, X = X, K = 4)

#The ridge_lfmm function returns an object that contains the latent variable score matrix U,
#the latent variable loading matrix U, and B the effect sizes for all SNPs.
#The result can be used to perform an association study:
pv <- lfmm_test(Y = Y.2L, X = X, lfmm = mod.lfmm, calibrate = "gif")

#pull out ordered p values
pvalues <- pv$calibrated.pvalue

#bring in dataframe with matrix positions for 2L
matrix.pos.2L <- read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.2L.geno.matrix.012.pos", 
                        sep = "\t", header = F)
#binding pvalues to positions
matrix.pos.2L$pvalue <- pvalues[,1]

#2R
Y.2R <- data.matrix(read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.2R.geno.matrix.012", 
                             header = FALSE, row.names = 1, sep = "\t"))
y.pc <- prcomp(Y.2R); spc.y2R <- summary(y.pc)
plot(spc.y2R$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
abline(v = 5, col = "red", lty = 2)
mod.lfmm <- lfmm_ridge(Y = Y.2R, X = X, K = 5)
pv <- lfmm_test(Y = Y.2R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.2R <- read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.2R.geno.matrix.012.pos", 
                        sep = "\t", header = F)
matrix.pos.2R$pvalue <- pvalues[,1]

#3L
Y.3L <- data.matrix(read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.3L.geno.matrix.012", 
                             header = FALSE, row.names = 1, sep = "\t"))
y.pc <- prcomp(Y.3L); spc.y3L <- summary(y.pc)
plot(spc.y3L$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
abline(v = 5, col = "red", lty = 2)
mod.lfmm <- lfmm_ridge(Y = Y.3L, X = X, K = 5)
pv <- lfmm_test(Y = Y.3L, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3L <- read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.3L.geno.matrix.012.pos", 
                        sep = "\t", header = F)
matrix.pos.3L$pvalue <- pvalues[,1]

#3R
Y.3R <- data.matrix(read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.3R.geno.matrix.012", 
                             header = FALSE, row.names = 1, sep = "\t"))
y.pc <- prcomp(Y.3R); spc.y3R <- summary(y.pc)
plot(spc.y3R$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
abline(v = 5, col = "red", lty = 2)
mod.lfmm <- lfmm_ridge(Y = Y.3R, X = X, K = 5)
pv <- lfmm_test(Y = Y.3R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3R <- read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.3R.geno.matrix.012.pos", 
                        sep = "\t", header = F)
matrix.pos.3R$pvalue <- pvalues[,1]

#X
Y.X <- data.matrix(read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.X.geno.matrix.012", 
                            header = FALSE, row.names = 1, sep = "\t"))
y.pc <- prcomp(Y.X); spc.yx <- summary(y.pc)
plot(spc.yx$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
abline(v = 5, col = "red", lty = 2)
mod.lfmm <- lfmm_ridge(Y = Y.X, X = X, K = 5)
pv <- lfmm_test(Y = Y.X, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.X <- read.csv("FST_LFMM_BayeScan/ag1000g.phase1.ar3.pass.biallelic.maf05.X.geno.matrix.012.pos", 
                       sep = "\t", header = F)
matrix.pos.X$pvalue <- pvalues[,1]

#####
#Critical pvalues of LFMMs (adding adjusted p values as columns for each chromosome)
matrix.pos.2L$fdr <- p.adjust(matrix.pos.2L$pvalue, method = "fdr")
matrix.pos.2L$bonferroni <- p.adjust(matrix.pos.2L$pvalue, method = "bonferroni")

matrix.pos.2R$fdr <- p.adjust(matrix.pos.2R$pvalue, method = "fdr")
matrix.pos.2R$bonferroni <- p.adjust(matrix.pos.2R$pvalue, method = "bonferroni")

matrix.pos.3L$fdr <- p.adjust(matrix.pos.3L$pvalue, method = "fdr")
matrix.pos.3L$bonferroni <- p.adjust(matrix.pos.3L$pvalue, method = "bonferroni")

matrix.pos.3R$fdr <- p.adjust(matrix.pos.3R$pvalue, method = "fdr")
matrix.pos.3R$bonferroni <- p.adjust(matrix.pos.3R$pvalue, method = "bonferroni")

matrix.pos.X$fdr <- p.adjust(matrix.pos.X$pvalue, method = "fdr")
matrix.pos.X$bonferroni <- p.adjust(matrix.pos.X$pvalue, method = "bonferroni")

#####
#save initial data as rdata
save(phase1.env, phase2.env, phase1.all.env, Y.2L, Y.2R, Y.3L, Y.3R, Y.X, 
     file = "FST_LFMM_BayeScan/LFMM_initial_data.RData")

#####
#visualization 
fig.dir <- "FST_LFMM_BayeScan/Exploratory_LFMM_figures_pecipitation"
dir.create(fig.dir)

#2L
#The histogram of test significance values is expected to be flat, with a peak near zero
jpeg(paste0(fig.dir, "/2L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2L$pvalue, main = "2L LFMM p-values", xlab = "p-values")
hist(matrix.pos.2L$fdr, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2L$bonferroni, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

#A QQ-plot is displayed as follows.
jpeg(paste0(fig.dir, "/2L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2L$pvalue), rate = log(10)), -log10(matrix.pos.2L$pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

## Manhattan plot
jpeg(paste0(fig.dir, "/2L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2L$pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "2L Manhattan plot", col = "darkgrey")
crit.pval.fdr.2L <- -log10(max(matrix.pos.2L$pvalue[matrix.pos.2L$fdr <= 0.05]))
crit.pval.bon.2L <- -log10(max(matrix.pos.2L$pvalue[matrix.pos.2L$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2L, col = "blue", lty = 2)
abline(h = crit.pval.bon.2L, col = "red", lty = 2)
dev.off()

#2R
jpeg(paste0(fig.dir, "/2R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2R$pvalue, main = "2R LFMM p-values", xlab = "p-values")
hist(matrix.pos.2R$fdr, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2R$bonferroni, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/2R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2R$pvalue), rate = log(10)), -log10(matrix.pos.2R$pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/2R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2R$pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "2R Manhattan plot", col = "darkgrey")
crit.pval.fdr.2R <- -log10(max(matrix.pos.2R$pvalue[matrix.pos.2R$fdr <= 0.05]))
crit.pval.bon.2R <- -log10(max(matrix.pos.2R$pvalue[matrix.pos.2R$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2R, col = "blue", lty = 2)
abline(h = crit.pval.bon.2R, col = "red", lty = 2)
dev.off()

#3L
jpeg(paste0(fig.dir, "/3L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3L$pvalue, main = "3L LFMM p-values", xlab = "p-values")
hist(matrix.pos.3L$fdr, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3L$bonferroni, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3L$pvalue), rate = log(10)), -log10(matrix.pos.3L$pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3L$pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "3L Manhattan plot", col = "darkgrey")
crit.pval.fdr.3L <- -log10(max(matrix.pos.3L$pvalue[matrix.pos.3L$fdr <= 0.05]))
crit.pval.bon.3L <- -log10(max(matrix.pos.3L$pvalue[matrix.pos.3L$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3L, col = "blue", lty = 2)
abline(h = crit.pval.bon.3L, col = "red", lty = 2)
dev.off()

#3R
jpeg(paste0(fig.dir, "/3R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3R$pvalue, main = "3R LFMM p-values", xlab = "p-values")
hist(matrix.pos.3R$fdr, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3R$bonferroni, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3R$pvalue), rate = log(10)), -log10(matrix.pos.3R$pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3R$pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "3R Manhattan plot", col = "darkgrey")
crit.pval.fdr.3R <- -log10(max(matrix.pos.3R$pvalue[matrix.pos.3R$fdr <= 0.05]))
crit.pval.bon.3R <- -log10(max(matrix.pos.3R$pvalue[matrix.pos.3R$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3R, col = "blue", lty = 2)
abline(h = crit.pval.bon.3R, col = "red", lty = 2)
dev.off()

#X
jpeg(paste0(fig.dir, "/X_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.X$pvalue, main = "X LFMM p-values", xlab = "p-values")
hist(matrix.pos.X$fdr, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.X$bonferroni, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/X_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.X$pvalue), rate = log(10)), -log10(matrix.pos.X$pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "X LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/X_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.X$pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "X Manhattan plot", col = "darkgrey")
crit.pval.fdr.X <- -log10(max(matrix.pos.X$pvalue[matrix.pos.X$fdr <= 0.05]))
crit.pval.bon.X <- -log10(max(matrix.pos.X$pvalue[matrix.pos.X$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.X, col = "blue", lty = 2)
abline(h = crit.pval.bon.X, col = "red", lty = 2)
dev.off()

#####
#Selecting candidate loci
lfmm.candidate.loci.fdr.2L <- matrix.pos.2L[matrix.pos.2L$fdr <= 0.05,]
lfmm.candidate.loci.fdr.2R <- matrix.pos.2R[matrix.pos.2R$fdr <= 0.05,]
lfmm.candidate.loci.fdr.3L <- matrix.pos.3L[matrix.pos.3L$fdr <= 0.05,]
lfmm.candidate.loci.fdr.3R <- matrix.pos.3R[matrix.pos.3R$fdr <= 0.05,]
lfmm.candidate.loci.fdr.X <- matrix.pos.X[matrix.pos.X$fdr <= 0.05,]

lfmm.candidate.loci.bon.2L <- matrix.pos.2L[matrix.pos.2L$bonferroni <= 0.05,]
lfmm.candidate.loci.bon.2R <- matrix.pos.2R[matrix.pos.2R$bonferroni <= 0.05,]
lfmm.candidate.loci.bon.3L <- matrix.pos.3L[matrix.pos.3L$bonferroni <= 0.05,]
lfmm.candidate.loci.bon.3R <- matrix.pos.3R[matrix.pos.3R$bonferroni <= 0.05,]
lfmm.candidate.loci.bon.X <- matrix.pos.X[matrix.pos.X$bonferroni <= 0.05,]

#save candidate loci as matrices
save(lfmm.candidate.loci.fdr.2L, lfmm.candidate.loci.fdr.2R, lfmm.candidate.loci.fdr.3L,
     lfmm.candidate.loci.fdr.3R, lfmm.candidate.loci.fdr.X,
     lfmm.candidate.loci.bon.2L, lfmm.candidate.loci.bon.2R,
     lfmm.candidate.loci.bon.3L, lfmm.candidate.loci.bon.3R,
     lfmm.candidate.loci.bon.X, file = "FST_LFMM_BayeScan/LFMM_condidate_loci_precipitation.RData")


#################
#Running lfmm2 with humidity

#variable
X <- phase1.all.env$hannual

#2L
mod.lfmm <- lfmm_ridge(Y = Y.2L, X = X, K = 4)
pv <- lfmm_test(Y = Y.2L, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.2L$hum.pvalue <- pvalues[,1]

#2R
mod.lfmm <- lfmm_ridge(Y = Y.2R, X = X, K = 5)
pv <- lfmm_test(Y = Y.2R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.2R$hum.pvalue <- pvalues[,1]

#3L
mod.lfmm <- lfmm_ridge(Y = Y.3L, X = X, K = 5)
pv <- lfmm_test(Y = Y.3L, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3L$hum.pvalue <- pvalues[,1]

#3R
mod.lfmm <- lfmm_ridge(Y = Y.3R, X = X, K = 5)
pv <- lfmm_test(Y = Y.3R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3R$hum.pvalue <- pvalues[,1]

#X
mod.lfmm <- lfmm_ridge(Y = Y.X, X = X, K = 5)
pv <- lfmm_test(Y = Y.X, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.X$hum.pvalue <- pvalues[,1]

#####
#Critical pvalues of LFMMs
matrix.pos.2L$hum.fdr <- p.adjust(matrix.pos.2L$hum.pvalue, method = "fdr")
matrix.pos.2L$hum.bonferroni <- p.adjust(matrix.pos.2L$hum.pvalue, method = "bonferroni")

matrix.pos.2R$hum.fdr <- p.adjust(matrix.pos.2R$hum.pvalue, method = "fdr")
matrix.pos.2R$hum.bonferroni <- p.adjust(matrix.pos.2R$hum.pvalue, method = "bonferroni")

matrix.pos.3L$hum.fdr <- p.adjust(matrix.pos.3L$hum.pvalue, method = "fdr")
matrix.pos.3L$hum.bonferroni <- p.adjust(matrix.pos.3L$hum.pvalue, method = "bonferroni")

matrix.pos.3R$hum.fdr <- p.adjust(matrix.pos.3R$hum.pvalue, method = "fdr")
matrix.pos.3R$hum.bonferroni <- p.adjust(matrix.pos.3R$hum.pvalue, method = "bonferroni")

matrix.pos.X$hum.fdr <- p.adjust(matrix.pos.X$hum.pvalue, method = "fdr")
matrix.pos.X$hum.bonferroni <- p.adjust(matrix.pos.X$hum.pvalue, method = "bonferroni")

#####
#visualization 
fig.dir <- "FST_LFMM_BayeScan/Exploratory_LFMM_figures_humidity"
dir.create(fig.dir)

#2L
jpeg(paste0(fig.dir, "/2L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2L$hum.pvalue, main = "2L LFMM p-values", xlab = "p-values")
hist(matrix.pos.2L$hum.fdr, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2L$hum.bonferroni, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/2L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2L$hum.pvalue), rate = log(10)), -log10(matrix.pos.2L$hum.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

## Manhattan plot
jpeg(paste0(fig.dir, "/2L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2L$hum.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "2L Manhattan plot", col = "darkgrey")
crit.pval.fdr.2L <- -log10(max(matrix.pos.2L$hum.pvalue[matrix.pos.2L$hum.fdr <= 0.05]))
crit.pval.bon.2L <- -log10(max(matrix.pos.2L$hum.pvalue[matrix.pos.2L$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2L, col = "blue", lty = 2)
abline(h = crit.pval.bon.2L, col = "red", lty = 2)
dev.off()

#2R
jpeg(paste0(fig.dir, "/2R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2R$hum.pvalue, main = "2R LFMM p-values", xlab = "p-values")
hist(matrix.pos.2R$hum.fdr, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2R$hum.bonferroni, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/2R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2R$hum.pvalue), rate = log(10)), -log10(matrix.pos.2R$hum.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/2R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2R$hum.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "2R Manhattan plot", col = "darkgrey")
crit.pval.fdr.2R <- -log10(max(matrix.pos.2R$hum.pvalue[matrix.pos.2R$hum.fdr <= 0.05]))
crit.pval.bon.2R <- -log10(max(matrix.pos.2R$hum.pvalue[matrix.pos.2R$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2R, col = "blue", lty = 2)
abline(h = crit.pval.bon.2R, col = "red", lty = 2)
dev.off()

#3L
jpeg(paste0(fig.dir, "/3L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3L$hum.pvalue, main = "3L LFMM p-values", xlab = "p-values")
hist(matrix.pos.3L$hum.fdr, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3L$hum.bonferroni, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3L$hum.pvalue), rate = log(10)), -log10(matrix.pos.3L$hum.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3L$hum.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "3L Manhattan plot", col = "darkgrey")
crit.pval.fdr.3L <- -log10(max(matrix.pos.3L$hum.pvalue[matrix.pos.3L$hum.fdr <= 0.05]))
crit.pval.bon.3L <- -log10(max(matrix.pos.3L$hum.pvalue[matrix.pos.3L$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3L, col = "blue", lty = 2)
abline(h = crit.pval.bon.3L, col = "red", lty = 2)
dev.off()

#3R
jpeg(paste0(fig.dir, "/3R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3R$hum.pvalue, main = "3R LFMM p-values", xlab = "p-values")
hist(matrix.pos.3R$hum.fdr, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3R$hum.bonferroni, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3R$hum.pvalue), rate = log(10)), -log10(matrix.pos.3R$hum.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3R$hum.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "3R Manhattan plot", col = "darkgrey")
crit.pval.fdr.3R <- -log10(max(matrix.pos.3R$hum.pvalue[matrix.pos.3R$hum.fdr <= 0.05]))
crit.pval.bon.3R <- -log10(max(matrix.pos.3R$hum.pvalue[matrix.pos.3R$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3R, col = "blue", lty = 2)
abline(h = crit.pval.bon.3R, col = "red", lty = 2)
dev.off()

#X
jpeg(paste0(fig.dir, "/X_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.X$hum.pvalue, main = "X LFMM p-values", xlab = "p-values")
hist(matrix.pos.X$hum.fdr, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.X$hum.bonferroni, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/X_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.X$hum.pvalue), rate = log(10)), -log10(matrix.pos.X$hum.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "X LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/X_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.X$hum.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "X Manhattan plot", col = "darkgrey")
crit.pval.fdr.X <- -log10(max(matrix.pos.X$hum.pvalue[matrix.pos.X$hum.fdr <= 0.05]))
crit.pval.bon.X <- -log10(max(matrix.pos.X$hum.pvalue[matrix.pos.X$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.X, col = "blue", lty = 2)
abline(h = crit.pval.bon.X, col = "red", lty = 2)
dev.off()

#####
#Selecting candidate loci
lfmm.candidate.loci.fdr.hum.2L <- matrix.pos.2L[matrix.pos.2L$hum.fdr <= 0.05,]
lfmm.candidate.loci.fdr.hum.2R <- matrix.pos.2R[matrix.pos.2R$hum.fdr <= 0.05,]
lfmm.candidate.loci.fdr.hum.3L <- matrix.pos.3L[matrix.pos.3L$hum.fdr <= 0.05,]
lfmm.candidate.loci.fdr.hum.3R <- matrix.pos.3R[matrix.pos.3R$hum.fdr <= 0.05,]
lfmm.candidate.loci.fdr.hum.X <- matrix.pos.X[matrix.pos.X$hum.fdr <= 0.05,]

lfmm.candidate.loci.bon.hum.2L <- matrix.pos.2L[matrix.pos.2L$hum.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.hum.2R <- matrix.pos.2R[matrix.pos.2R$hum.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.hum.3L <- matrix.pos.3L[matrix.pos.3L$hum.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.hum.3R <- matrix.pos.3R[matrix.pos.3R$hum.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.hum.X <- matrix.pos.X[matrix.pos.X$hum.bonferroni <= 0.05,]

#save candidate loci as matrices
save(lfmm.candidate.loci.fdr.hum.2L, lfmm.candidate.loci.fdr.hum.2R, lfmm.candidate.loci.fdr.hum.3L,
     lfmm.candidate.loci.fdr.hum.3R, lfmm.candidate.loci.fdr.hum.X,
     lfmm.candidate.loci.bon.hum.2L, lfmm.candidate.loci.bon.hum.2R,
     lfmm.candidate.loci.bon.hum.3L, lfmm.candidate.loci.bon.hum.3R,
     lfmm.candidate.loci.bon.hum.X, file = "FST_LFMM_BayeScan/LFMM_condidate_loci_humidity.RData")


#################
#Running lfmm2 with temperature

#variable
X <- phase1.all.env$tmean

#2L
mod.lfmm <- lfmm_ridge(Y = Y.2L, X = X, K = 4)
pv <- lfmm_test(Y = Y.2L, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.2L$tem.pvalue <- pvalues[,1]

#2R
mod.lfmm <- lfmm_ridge(Y = Y.2R, X = X, K = 5)
pv <- lfmm_test(Y = Y.2R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.2R$tem.pvalue <- pvalues[,1]

#3L
mod.lfmm <- lfmm_ridge(Y = Y.3L, X = X, K = 5)
pv <- lfmm_test(Y = Y.3L, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3L$tem.pvalue <- pvalues[,1]

#3R
mod.lfmm <- lfmm_ridge(Y = Y.3R, X = X, K = 5)
pv <- lfmm_test(Y = Y.3R, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.3R$tem.pvalue <- pvalues[,1]

#X
mod.lfmm <- lfmm_ridge(Y = Y.X, X = X, K = 5)
pv <- lfmm_test(Y = Y.X, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
matrix.pos.X$tem.pvalue <- pvalues[,1]

#####
#Critical pvalues of LFMMs
matrix.pos.2L$tem.fdr <- p.adjust(matrix.pos.2L$tem.pvalue, method = "fdr")
matrix.pos.2L$tem.bonferroni <- p.adjust(matrix.pos.2L$tem.pvalue, method = "bonferroni")

matrix.pos.2R$tem.fdr <- p.adjust(matrix.pos.2R$tem.pvalue, method = "fdr")
matrix.pos.2R$tem.bonferroni <- p.adjust(matrix.pos.2R$tem.pvalue, method = "bonferroni")

matrix.pos.3L$tem.fdr <- p.adjust(matrix.pos.3L$tem.pvalue, method = "fdr")
matrix.pos.3L$tem.bonferroni <- p.adjust(matrix.pos.3L$tem.pvalue, method = "bonferroni")

matrix.pos.3R$tem.fdr <- p.adjust(matrix.pos.3R$tem.pvalue, method = "fdr")
matrix.pos.3R$tem.bonferroni <- p.adjust(matrix.pos.3R$tem.pvalue, method = "bonferroni")

matrix.pos.X$tem.fdr <- p.adjust(matrix.pos.X$tem.pvalue, method = "fdr")
matrix.pos.X$tem.bonferroni <- p.adjust(matrix.pos.X$tem.pvalue, method = "bonferroni")

#####
#visualization 
fig.dir <- "FST_LFMM_BayeScan/Exploratory_LFMM_figures_temperature"
dir.create(fig.dir)

#2L
jpeg(paste0(fig.dir, "/2L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2L$tem.pvalue, main = "2L LFMM p-values", xlab = "p-values")
hist(matrix.pos.2L$tem.fdr, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2L$tem.bonferroni, main = "2L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/2L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2L$tem.pvalue), rate = log(10)), -log10(matrix.pos.2L$tem.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/2L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2L$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "2L Manhattan plot", col = "darkgrey")
crit.pval.fdr.2L <- -log10(max(matrix.pos.2L$tem.pvalue[matrix.pos.2L$tem.fdr <= 0.05]))
crit.pval.bon.2L <- -log10(max(matrix.pos.2L$tem.pvalue[matrix.pos.2L$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2L, col = "blue", lty = 2)
abline(h = crit.pval.bon.2L, col = "red", lty = 2)
dev.off()

#2R
jpeg(paste0(fig.dir, "/2R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.2R$tem.pvalue, main = "2R LFMM p-values", xlab = "p-values")
hist(matrix.pos.2R$tem.fdr, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.2R$tem.bonferroni, main = "2R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/2R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.2R$tem.pvalue), rate = log(10)), -log10(matrix.pos.2R$tem.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "2R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/2R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.2R$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "2R Manhattan plot", col = "darkgrey")
crit.pval.fdr.2R <- -log10(max(matrix.pos.2R$tem.pvalue[matrix.pos.2R$tem.fdr <= 0.05]))
crit.pval.bon.2R <- -log10(max(matrix.pos.2R$tem.pvalue[matrix.pos.2R$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.2R, col = "blue", lty = 2)
abline(h = crit.pval.bon.2R, col = "red", lty = 2)
dev.off()

#3L
jpeg(paste0(fig.dir, "/3L_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3L$tem.pvalue, main = "3L LFMM p-values", xlab = "p-values")
hist(matrix.pos.3L$tem.fdr, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3L$tem.bonferroni, main = "3L LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3L_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3L$tem.pvalue), rate = log(10)), -log10(matrix.pos.3L$tem.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3L LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3L_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3L$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "3L Manhattan plot", col = "darkgrey")
crit.pval.fdr.3L <- -log10(max(matrix.pos.3L$tem.pvalue[matrix.pos.3L$tem.fdr <= 0.05]))
crit.pval.bon.3L <- -log10(max(matrix.pos.3L$tem.pvalue[matrix.pos.3L$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3L, col = "blue", lty = 2)
abline(h = crit.pval.bon.3L, col = "red", lty = 2)
dev.off()

#3R
jpeg(paste0(fig.dir, "/3R_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.3R$tem.pvalue, main = "3R LFMM p-values", xlab = "p-values")
hist(matrix.pos.3R$tem.fdr, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.3R$tem.bonferroni, main = "3R LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/3R_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.3R$tem.pvalue), rate = log(10)), -log10(matrix.pos.3R$tem.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "3R LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/3R_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.3R$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "3R Manhattan plot", col = "darkgrey")
crit.pval.fdr.3R <- -log10(max(matrix.pos.3R$tem.pvalue[matrix.pos.3R$tem.fdr <= 0.05]))
crit.pval.bon.3R <- -log10(max(matrix.pos.3R$tem.pvalue[matrix.pos.3R$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.3R, col = "blue", lty = 2)
abline(h = crit.pval.bon.3R, col = "red", lty = 2)
dev.off()

#X
jpeg(paste0(fig.dir, "/X_hist_pvalues.jpg"), width = 166, height = 200, units = "mm", res = 300)
par(mfrow = c(3, 1), cex = 0.7)
hist(matrix.pos.X$tem.pvalue, main = "X LFMM p-values", xlab = "p-values")
hist(matrix.pos.X$tem.fdr, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (FDR)")
hist(matrix.pos.X$tem.bonferroni, main = "X LFMM adjusted p-values", xlab = "Adjusted p-values (Bonferroni)")
dev.off()

jpeg(paste0(fig.dir, "/X_LFMM_QQ-plot.jpg"), width = 166, height = 150, units = "mm", res = 300)
qqplot(rexp(length(matrix.pos.X$tem.pvalue), rate = log(10)), -log10(matrix.pos.X$tem.pvalue), 
       xlab = "Expected quantile", ylab = "-Log10 p-values", main = "X LFMM QQ-plot",
       pch = 19, cex = .4); abline(0,1)
dev.off()

jpeg(paste0(fig.dir, "/X_Manhattan_plot.jpg"), width = 166, height = 90, units = "mm", res = 300)
par(cex = 0.7);plot(-log10(matrix.pos.X$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
                    main = "X Manhattan plot", col = "darkgrey")
crit.pval.fdr.X <- -log10(max(matrix.pos.X$tem.pvalue[matrix.pos.X$tem.fdr <= 0.05]))
crit.pval.bon.X <- -log10(max(matrix.pos.X$tem.pvalue[matrix.pos.X$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.X, col = "blue", lty = 2)
abline(h = crit.pval.bon.X, col = "red", lty = 2)
dev.off()

#####
#Selecting candidate loci
lfmm.candidate.loci.fdr.tem.2L <- matrix.pos.2L[matrix.pos.2L$tem.fdr <= 0.05,]
lfmm.candidate.loci.fdr.tem.2R <- matrix.pos.2R[matrix.pos.2R$tem.fdr <= 0.05,]
lfmm.candidate.loci.fdr.tem.3L <- matrix.pos.3L[matrix.pos.3L$tem.fdr <= 0.05,]
lfmm.candidate.loci.fdr.tem.3R <- matrix.pos.3R[matrix.pos.3R$tem.fdr <= 0.05,]
lfmm.candidate.loci.fdr.tem.X <- matrix.pos.X[matrix.pos.X$tem.fdr <= 0.05,]

lfmm.candidate.loci.bon.tem.2L <- matrix.pos.2L[matrix.pos.2L$tem.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.tem.2R <- matrix.pos.2R[matrix.pos.2R$tem.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.tem.3L <- matrix.pos.3L[matrix.pos.3L$tem.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.tem.3R <- matrix.pos.3R[matrix.pos.3R$tem.bonferroni <= 0.05,]
lfmm.candidate.loci.bon.tem.X <- matrix.pos.X[matrix.pos.X$tem.bonferroni <= 0.05,]

#save candidate loci as matrices
save(lfmm.candidate.loci.fdr.tem.2L, lfmm.candidate.loci.fdr.tem.2R, lfmm.candidate.loci.fdr.tem.3L,
     lfmm.candidate.loci.fdr.tem.3R, lfmm.candidate.loci.fdr.tem.X,
     lfmm.candidate.loci.bon.tem.2L, lfmm.candidate.loci.bon.tem.2R,
     lfmm.candidate.loci.bon.tem.3L, lfmm.candidate.loci.bon.tem.3R,
     lfmm.candidate.loci.bon.tem.X, file = "FST_LFMM_BayeScan/LFMM_condidate_loci_temperature.RData")


#################
#whole genome
matrix.pos.whole.gen <- rbind(matrix.pos.2R, matrix.pos.2L, matrix.pos.3R,
                              matrix.pos.3L, matrix.pos.X)

#p adjustment
matrix.pos.whole.gen$fdr <- p.adjust(matrix.pos.whole.gen$pvalue, method = "fdr")
matrix.pos.whole.gen$bonferroni <- p.adjust(matrix.pos.whole.gen$pvalue, method = "bonferroni")

#saving result matrices as RData
save(matrix.pos.2R, matrix.pos.2L, matrix.pos.3R, matrix.pos.3L, matrix.pos.X, 
     matrix.pos.whole.gen, file = "FST_LFMM_BayeScan/LFMM_position_matrices.RData")

#Selecting candidate loci
sel.lfmm.pre.fdr <- matrix.pos.whole.gen[matrix.pos.whole.gen$fdr <= 0.05, 1:4]
sel.lfmm.pre.bon <- matrix.pos.whole.gen[matrix.pos.whole.gen$bonferroni <= 0.05, c(1:3, 5)]
sel.lfmm.pre.T100 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$pvalue), ], 100)[, 1:3]
sel.lfmm.pre.T10 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$pvalue), ], 10)[, 1:3]

sel.lfmm.hum.fdr <- matrix.pos.whole.gen[matrix.pos.whole.gen$hum.fdr <= 0.05, 1:4]
sel.lfmm.hum.bon <- matrix.pos.whole.gen[matrix.pos.whole.gen$hum.bonferroni <= 0.05, c(1:3, 5)]
sel.lfmm.hum.T100 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$hum.pvalue), ], 100)[, 1:3]
sel.lfmm.hum.T10 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$hum.pvalue), ], 10)[, 1:3]

sel.lfmm.tem.fdr <- matrix.pos.whole.gen[matrix.pos.whole.gen$tem.fdr <= 0.05, 1:4]
sel.lfmm.tem.bon <- matrix.pos.whole.gen[matrix.pos.whole.gen$tem.bonferroni <= 0.05, c(1:3, 5)]
sel.lfmm.tem.T100 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$tem.pvalue), ], 100)[, 1:3]
sel.lfmm.tem.T10 <- head(matrix.pos.whole.gen[order(matrix.pos.whole.gen$tem.pvalue), ], 10)[, 1:3]

#save candidate loci as matrices for whole genome and distinct variables
save(sel.lfmm.pre.fdr, sel.lfmm.pre.bon, sel.lfmm.pre.T100, sel.lfmm.pre.T10, 
     file = "FST_LFMM_BayeScan/LFMM_condidate_loci_precipitation_wgen.RData")

save(sel.lfmm.hum.fdr, sel.lfmm.hum.bon, sel.lfmm.hum.T100, sel.lfmm.hum.T10, 
     file = "FST_LFMM_BayeScan/LFMM_condidate_loci_humidity_wgen.RData")

save(sel.lfmm.tem.fdr, sel.lfmm.tem.bon, sel.lfmm.tem.T100, sel.lfmm.tem.T10, 
     file = "FST_LFMM_BayeScan/LFMM_condidate_loci_temperature_wgen.RData")


#################
#LFMM summary figures
fdir <- "FST_LFMM_BayeScan/Summary_figures_LFMM"
dir.create(fdir)

#####
#relationship p-values precipitation and p-values temperature and p-values humity
jpeg(paste0(fdir, "/Whole_genome_pvalue_comparison_Prec_Hum.jpg"), width = 166, height = 150, units = "mm", res = 300)
par(cex = 0.7, mar = c(4.5, 4.5, 0.5, 0.5))
plot(-log10(matrix.pos.whole.gen$pvalue), -log10(matrix.pos.whole.gen$hum.pvalue), 
     xlab = "-Log10 p-values (precipitation)", ylab = "-Log10 p-values (humidity)", col = "darkgrey")
abline(v = 20, col = "red", lty = 2)
abline(h = 20, col = "red", lty = 2)
dev.off()

jpeg(paste0(fdir, "/Whole_genome_pvalue_comparison_Prec_Temp.jpg"), width = 166, height = 150, units = "mm", res = 300)
par(cex = 0.7, mar = c(4.5, 4.5, 0.5, 0.5))
plot(-log10(matrix.pos.whole.gen$pvalue), -log10(matrix.pos.whole.gen$tem.pvalue), 
     xlab = "-Log10 p-values (precipitation)", ylab = "-Log10 p-values (temperature)", col = "darkgrey")
abline(v = 20, col = "red", lty = 2)
abline(h = 20, col = "red", lty = 2)
dev.off()

jpeg(paste0(fdir, "/Whole_genome_pvalue_comparison_Hum_Temp.jpg"), width = 166, height = 150, units = "mm", res = 300)
par(cex = 0.7, mar = c(4.5, 4.5, 0.5, 0.5))
plot(-log10(matrix.pos.whole.gen$hum.pvalue), -log10(matrix.pos.whole.gen$tem.pvalue), 
     xlab = "-Log10 p-values (humidity)", ylab = "-Log10 p-values (temperature)", col = "darkgrey")
abline(v = 20, col = "red", lty = 2)
abline(h = 20, col = "red", lty = 2)
dev.off()

#####
#LFMM Manhattan plot all variables
##chromosom id place and colors
n2R <- nrow(matrix.pos.2R); n2L <- nrow(matrix.pos.2L); n3R <- nrow(matrix.pos.3R)
n3L <- nrow(matrix.pos.3L); nX <- nrow(matrix.pos.X)
chrom.place <- rbind(c(1, n2R / 2), c(1, n2R + (n2L / 2)), c(1, n2R + n2L + (n3R / 2)),
                     c(1, n2R + n2L + n3R + (n3L / 2)), c(1, n2R + n2L + n3R + n3L + (nX / 2)))

cols <- c(rep("gray45", nrow(matrix.pos.2R)), rep("gray70", nrow(matrix.pos.2L)),
          rep("gray45", nrow(matrix.pos.3R)), rep("gray70", nrow(matrix.pos.3L)),
          rep("gray45", nrow(matrix.pos.X)))

##plot
jpeg(paste0(fdir, "/Whole_genome_Manhattan_plot_all_hres.jpg"), width = 145, height = 210, units = "mm", res = 600)
par(cex = 0.5, mfrow = c(3, 1), mar = c(4, 4.5, 0.25, 0.5))
plot(-log10(matrix.pos.whole.gen$pvalue), pch = 19, cex = .4, xlab = "", ylab = "-Log10 p-values", 
     main = "", col = cols)
crit.pval.fdr.whole.gen <- -log10(max(matrix.pos.whole.gen$pvalue[matrix.pos.whole.gen$fdr <= 0.05]))
crit.pval.bon.whole.gen <- -log10(max(matrix.pos.whole.gen$pvalue[matrix.pos.whole.gen$bonferroni <= 0.05]))
abline(h = crit.pval.fdr.whole.gen, col = "blue", lty = 2); abline(h = crit.pval.bon.whole.gen, col = "red", lty = 2)
abline(h = 10, col = "purple", lty = 2); abline(h = 15, col = "darkgreen", lty = 2); abline(h = 20, col = "green", lty = 2)
text(chrom.place[, 2], chrom.place[, 1], labels = c("2R", "2L", "3R", "3L", "X"), col = "white")
legend("topleft", legend = "Precipitation", bty = "n")

plot(-log10(matrix.pos.whole.gen$hum.pvalue), pch = 19, cex = .4, xlab = "", ylab = "-Log10 p-values", 
     main = "", col = cols)
crit.pval.fdr.whole.gen <- -log10(max(matrix.pos.whole.gen$hum.pvalue[matrix.pos.whole.gen$hum.fdr <= 0.05]))
crit.pval.bon.whole.gen <- -log10(max(matrix.pos.whole.gen$hum.pvalue[matrix.pos.whole.gen$hum.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.whole.gen, col = "blue", lty = 2); abline(h = crit.pval.bon.whole.gen, col = "red", lty = 2)
abline(h = 10, col = "purple", lty = 2); abline(h = 15, col = "darkgreen", lty = 2); abline(h = 20, col = "green", lty = 2)
text(chrom.place[, 2], chrom.place[, 1], labels = c("2R", "2L", "3R", "3L", "X"), col = "white")
legend("topleft", legend = "Humidity", bty = "n")

plot(-log10(matrix.pos.whole.gen$tem.pvalue), pch = 19, cex = .4, xlab = "SNP", ylab = "-Log10 p-values", 
     main = "", col = cols)
crit.pval.fdr.whole.gen <- -log10(max(matrix.pos.whole.gen$tem.pvalue[matrix.pos.whole.gen$tem.fdr <= 0.05]))
crit.pval.bon.whole.gen <- -log10(max(matrix.pos.whole.gen$tem.pvalue[matrix.pos.whole.gen$tem.bonferroni <= 0.05]))
abline(h = crit.pval.fdr.whole.gen, col = "blue", lty = 2); abline(h = crit.pval.bon.whole.gen, col = "red", lty = 2)
abline(h = 10, col = "purple", lty = 2); abline(h = 15, col = "darkgreen", lty = 2); abline(h = 20, col = "green", lty = 2)
text(chrom.place[, 2], chrom.place[, 1], labels = c("2R", "2L", "3R", "3L", "X"), col = "white")
legend("topleft", legend = "Temperature", bty = "n")
dev.off()