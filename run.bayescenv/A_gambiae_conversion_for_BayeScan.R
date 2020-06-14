##################################
# Project: Anopheles gambiae ecological niche and genomics 
# Process: Converting SNP matrix into geste format with genotype info 
# Author(s): Devon A DeRaad, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 01/04/2020 (dd/mm/yyyy)
##################################

#define a function which will accept the snp matrix and output a geste matrix for a given population
geste.mat <- function(x) {
  #save each of the four columns needed for geste file as vectors and then cbind them into a dataframe
  #calculate the number of haplotypes (2 * # of non-missing genotypes)
  haplotypes<-(colSums(x >=0))*2
  #calculate how many alleles at each site (we have already filtered to only bi-allelic sites, so just rep 2)
  alleles<-rep(max(x), times = ncol(x))
  #calculate non-reference haplotypes, 1 means het, 2 means hom alt, so add the # of 1 genotypes to 2* # of 2 genotypes
  count.alt.haplotypes<-(colSums(x == 1) + (colSums(x ==2)*2))
  #calculate ref haplotypes, 0 means hom ref, 1 means het, so 2* number of 0 sites plus 1* number of het sites
  count.ref.haplotypes<-((colSums(x == 0)*2) + colSums(x ==1))
  #build df
  df.x<-as.data.frame(cbind(haplotypes, alleles, count.ref.haplotypes, count.alt.haplotypes))
  #fix rownames
  rownames(df.x) <- (1:length(x[1,]))
  
  return(df.x)
}

#bring in dataframe with locality and environmental data for each individual
phase1.alldata <- read.csv(file = "FST_LFMM_BayeScan/phase1.allvariables.csv")

#genetic matrice names
chrom_names <- list.files("FST_LFMM_BayeScan/", pattern = "geno.matrix.012$", full.names = TRUE)
chroms <- c("2L", "2R", "3L", "3R", "X")

#new folder for geste files
dir.create("Geste_files")

#geste file names
gfiles <- paste0("FST_LFMM_BayeScan/Geste_files/", chroms,".geste.txt")

#latitudes indicating distinct populations
lats <- unique(phase1.alldata$latitude)

#running in loop
pm <- lapply(1:length(chroms), function(x) {
  #bring in chrom SNP matrix
  chrom <- data.matrix(read.csv(chrom_names[x], header = FALSE, row.names = 1, sep = "\t"))
  
  #create the file header
  tc <- paste0("[loci]=", ncol(chrom), "\n") # "[loci]=_\n","\n","[populations]=14\n"
  cat(tc, "\n", "[populations]=14\n", file = gfiles[x], sep = "")
  
  #add pops to file
  rd <- lapply(1:length(lats), function(y) {
    #subet matrix to only individuals (rows) from pop1
    snp.mat <- chrom[phase1.alldata$latitude == lats[y],]
    #turn snp matrix into geste matrix
    ges <- geste.mat(snp.mat)
    #cat the appropriate '[pop]=_' identifier as a header
    cat("\n", paste0("[pop]=", y, "\n"), file = gfiles[x], sep = "", append = T)
    #append the geste matrix to the population header
    write.table(ges, file = gfiles[x], col.names = F, quote = F, append = T)
  })
  cat(x, "of", length(chroms), "finished\n")
})