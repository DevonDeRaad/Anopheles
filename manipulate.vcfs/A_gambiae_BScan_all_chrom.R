##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: BayeScan candidate loci selection 
# Author(s): Devon A DeRaad, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 05/03/2020 (dd/mm/yyyy)
##################################

#pacakges

# files and other variables
bsfiles <- paste0("FST_LFMM_BayeScan/", c("precip", "hum", "temp"), ".baye.txt")
chr <- c("2R", "2L", "3R", "3L", "X")

# position of SNPs
posfiles <- list.files("FST_LFMM_BayeScan/", pattern = "pos$", full.names = T)[c(2, 1, 4, 3, 5)]

pos <- lapply(posfiles, read.table, sep = "\t")
posi <- do.call(rbind, pos)
colnames(posi) <- c("Chrom", "pos")

# candidate loci selection
byescan_results_corrections <- lapply(1:length(bsfiles), function(x) {
  ## reading file
  bs <- read.table(bsfiles[x], sep = " ", header = T)
  colnames(bs)[ncol(bs)] <- "Chrom"
  
  ## reorganizing according to chromosoms
  bs <- bs[order(match(as.character(bs$X), chr)), ]
  bs <- cbind(bs, pos = posi$pos)
  
  ## correction
  bs$fdr <- p.adjust(bs$qval_g, method = "fdr")
  bs$bonferroni <- p.adjust(bs$qval_g, method = "bonferroni")
  
  bs
})


## selection precipitation
bs <- byescan_results_corrections[[1]]
sel.bsca.pre.fdr <- bs[bs$fdr <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.pre.bon <- bs[bs$bonferroni <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.pre.T10 <- tail(sel.bsca.pre.fdr[order(sel.bsca.pre.fdr$fst), ], 10)

## saving results
save(sel.bsca.pre.fdr, sel.bsca.pre.bon, sel.bsca.pre.T10, 
     file = "FST_LFMM_BayeScan/BayeScan_candidate_loci_precipitation.Rdata")


## selection humidity
bs <- byescan_results_corrections[[2]]
sel.bsca.hum.fdr <- bs[bs$fdr <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.hum.bon <- bs[bs$bonferroni <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.hum.T10 <- tail(sel.bsca.hum.fdr[order(sel.bsca.hum.fdr$fst), ], 10)

## saving results
save(sel.bsca.hum.fdr, sel.bsca.hum.bon, sel.bsca.hum.T10, 
     file = "FST_LFMM_BayeScan/BayeScan_candidate_loci_humidity.Rdata")


## selection temperature
bs <- byescan_results_corrections[[3]]
sel.bsca.tem.fdr <- bs[bs$fdr <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.tem.bon <- bs[bs$bonferroni <= 0.05, c("Chrom", "pos", "qval_g", "fst")]
sel.bsca.tem.T10 <- tail(sel.bsca.tem.fdr[order(sel.bsca.tem.fdr$fst), ], 10)

## saving results
save(sel.bsca.tem.fdr, sel.bsca.tem.bon, sel.bsca.tem.T10, 
     file = "FST_LFMM_BayeScan/BayeScan_candidate_loci_temperature.Rdata")

names(byescan_results_corrections) <- c("prec", "hum", "temp")
save(byescan_results_corrections, file = "FST_LFMM_BayeScan/BayeScan_results_corrections.RData")
