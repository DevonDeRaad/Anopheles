##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: LFMM precipitation, humidity, temperature 
# Author(s): Devon A DeRaad, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 04/03/2020 (dd/mm/yyyy)
##################################

#packages
library(lfmm)
library(scales)

#Data
#bring in marlons dataset
phase1_all_env <- read.csv("LFMM/phase1_all_vairables_last.csv")

#variables
colnames(phase1_all_env)
Xvars <- c("p_annual", "h_annual", "t_mean", "autocorrelated", "random")

#chromosomes
chroms <- c("2R", "2L", "3R", "3L", "X") 

#genomic data
chro_matrices <- paste0("LFMM/", c("zip.maf05.ld100.matrix.2R.012",
                                   "zip.maf05.ld100.matrix.2L.012",
                                   "zip.maf05.ld100.matrix.3R.012",
                                   "zip.maf05.ld100.matrix.3L.012",
                                   "zip.maf05.ld100.matrix.X.012"))

chro_matrices <- lapply(chro_matrices, function(x) {
  data.matrix(read.csv(x, header = FALSE, row.names = 1, sep = "\t"))
})


#finding appropriate number of factors
dir.create("LFMM/latent_factor_exploration")
dir.create("LFMM/PCA_latent_factor_decision")

##populations
pops <- as.factor(paste(phase1_all_env$longitude, phase1_all_env$latitude))

col_pal <- rainbow(length(levels(pops)))
cols <- col_pal[pops]

#running lfmm in loop for all variables
k_lfmm <- lapply(1:length(Xvars), function(i) {
  X <- phase1_all_env[, Xvars[i]]
  
  #running for all chomosomes in loop
  all_chroms <- lapply(1:length(chro_matrices), function(j) {
    #pca to check variable contribution
    if (i == 1) {
      y.pc <- prcomp(chro_matrices[[j]]) #will run for 15 minutes
      spc.y <- summary(y.pc)
      png(paste0("LFMM/PCA_latent_factor_decision/Variance_PCA_", chroms[j], ".png"), 
          width = 4, height = 3, units = "in", res = 300)
      par(cex = 0.75, mar = c(4, 4, 1, 1))
      plot(spc.y$importance[2, 1:20], xlab = "PC", ylab = "Variance explained")
      dev.off()
    }
    
    #The ridge_lfmm function is applied to Y (genotype) and X (phenotypes, SD = 1). to find k factors
    mod_lfmm <- lfmm_ridge(Y = chro_matrices[[j]], X = X, K = 13)
    
    png(paste0("LFMM/latent_factor_exploration/", Xvars[i], "_", chroms[j], ".png"), 
        width = 7, height = 9.33, units = "in", res = 300)
    par(mfrow = c(4, 3), cex = 0.5, mar = c(4, 4, 1, 1))
    for (k in 1:12) {
      plot(mod_lfmm$U[, k:(k+1)], pch = 16, col = alpha(cols, 0.6), 
           xlab = paste("factor", k), ylab =  paste("factor", k+1))
    }
    dev.off()
    
  })
})

#check resultant plots and define number of factors per chromosome
#k factors for each chromosome 
kfactors <- c(5, 4, 5, 5, 5)

#running final lfmm in loop for all variables
all_lfmm <- lapply(1:length(Xvars), function(i) {
  X <- phase1_all_env[, Xvars[i]]
  
  #running for all chomosomes in loop
  all_chroms <- lapply(1:length(chro_matrices), function(j) {
    #The ridge_lfmm function is applied to Y (genotype) and X (phenotypes, SD = 1).
    ## Fit an LFMM, i.e, compute B, U, V estimates 
    mod_lfmm <- lfmm_ridge(Y = chro_matrices[[j]], X = X, K = kfactors[j])
    
    #The ridge_lfmm function returns an object that contains the latent variable score matrix U,
    pv <- lfmm_test(Y = chro_matrices[[j]], X = X, lfmm = mod_lfmm, calibrate = "gif")
    
    #pull out ordered p values
    pvalues <- pv$calibrated.pvalue[, 1]
  })
  
  #whole genome pvalues
  unlist(all_chroms)
})

pvals <- do.call(cbind, all_lfmm)
colnames(pvals) <- paste0("p_value_", Xvars)


#position matrices
pos_matrices <- paste0("LFMM/", c("maf05.ld100.matrix.2R.012.pos",
                                  "maf05.ld100.matrix.2L.012.pos",
                                  "maf05.ld100.matrix.3R.012.pos",
                                  "maf05.ld100.matrix.3L.012.pos",
                                  "maf05.ld100.matrix.X.012.pos"))

pos_matrices <- lapply(pos_matrices, read.csv,  sep = "\t", header = F)

wg_positions <- do.call(rbind, pos_matrices)

# combining results with position matrices
wg_lfmm <- cbind(wg_positions, pvals)
colnames(wg_lfmm)[1:2] <- c("Chromosome", "Position")
head(wg_lfmm)

#saving results
save(wg_lfmm, pos_matrices, all_lfmm, file = "LFMM/LFMM_complete_results.RData") 

write.csv(wg_lfmm, "LFMM/whole_genome_LFMM_pvalues.csv", row.names = FALSE) 
