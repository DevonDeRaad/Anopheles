#!/bin/sh
#
#SBATCH --job-name=bi.allele              # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase2     # Set working d$
#SBATCH --mem-per-cpu=10gb            # memory requested
#SBATCH --time=8000

#this script will subset the phase2 vcfs to only new samples not included in phase1
#then subset the SNPs in them to only SNPs in phase1
#then convert the filtered vcfs into matrices for allele freq counting and plotting in R

#pull out only unique phase2 individuals from vcfs
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2L.vcf --keep phase2.uniques.txt --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2L
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2R.vcf --keep phase2.uniques.txt --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2R
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3L.vcf --keep phase2.uniques.txt --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3L
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3R.vcf --keep phase2.uniques.txt --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3R
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.X.vcf --keep phase2.uniques.txt --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.X


#pull out only SNPs present in our phase1 maf filtered dataset
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2L.recode.vcf --positions /home/d669d153/scratch/anoph.phase1/ag1000g.phase1.ar3.pass.biallelic.maf05.2L.geno.matrix.012.pos --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2L
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.2R.recode.vcf --positions /home/d669d153/scratch/anoph.phase1/ag1000g.phase1.ar3.pass.biallelic.maf05.2R.geno.matrix.012.pos --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2R
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3L.recode.vcf --positions /home/d669d153/scratch/anoph.phase1/ag1000g.phase1.ar3.pass.biallelic.maf05.3L.geno.matrix.012.pos --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3L
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.3R.recode.vcf --positions /home/d669d153/scratch/anoph.phase1/ag1000g.phase1.ar3.pass.biallelic.maf05.3R.geno.matrix.012.pos --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3R
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.X.recode.vcf --positions /home/d669d153/scratch/anoph.phase1/ag1000g.phase1.ar3.pass.biallelic.maf05.X.geno.matrix.012.pos --recode --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.X


#convert the filtered vcfs into matrices for allele frequency counting
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2L.recode.vcf --012 --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2L.geno.matrix
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2R.recode.vcf --012 --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.2R.geno.matrix
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3L.recode.vcf --012 --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3L.geno.matrix
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3R.recode.vcf --012 --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.3R.geno.matrix
vcftools --vcf /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.X.recode.vcf --012 --out /home/d669d153/scratch/anoph.phase2/ag1000g.phase2.ar1.pass.biallelic.phase1pos.X.geno.matrix
