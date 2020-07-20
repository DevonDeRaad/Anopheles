#!/bin/sh
#
#SBATCH --job-name=bayescenv.matrix             # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=25             # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1      # Set working d$
#SBATCH --mem-per-cpu=2gb            # memory requested
#SBATCH --time=50000


#run bayescan on the snp matrix
#/home/d669d153/work/BayeScan2.1/binaries/BayeScan2.1_linux64bit -snp all_chroms_geno_matrix.012 -o matrixbayescanOD10 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10
#module list
#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv X.geste -env phase1.zhum.txt -o bayescan/X.hum -threads 10

/home/d669d153/work/bayescenv-new/source/bayescenv X.geste -env phase1.zprecip.txt -o bayescan/X.precip -threads 10

/home/d669d153/work/bayescenv-new/source/bayescenv X.geste -env phase1.ztemp.txt -o bayescan/X.temp -threads 10



#!/bin/sh
#
#SBATCH --job-name=bayescenv.matrix             # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=25             # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1      # Set working d$
#SBATCH --mem-per-cpu=2gb            # memory requested
#SBATCH --time=50000


#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv 2L.geste.txt -env phase1.zhum.txt -o bayescan/2L.hum -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 2L.geste.txt -env phase1.zprecip.txt -o bayescan/2L.precip -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 2L.geste.txt -env phase1.ztemp.txt -o bayescan/2L.temp -threads 25




#!/bin/sh
#
#SBATCH --job-name=bayescenv.matrix             # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=25             # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1      # Set working d$
#SBATCH --mem-per-cpu=2gb            # memory requested
#SBATCH --time=50000


#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv 2R.geste.txt -env phase1.zhum.txt -o bayescan/2R.hum -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 2R.geste.txt -env phase1.zprecip.txt -o bayescan/2R.precip -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 2R.geste.txt -env phase1.ztemp.txt -o bayescan/2R.temp -threads 25



#!/bin/sh
#
#SBATCH --job-name=bayescenv.matrix             # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=25             # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1      # Set working d$
#SBATCH --mem-per-cpu=2gb            # memory requested
#SBATCH --time=50000


#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv 3L.geste.txt -env phase1.zhum.txt -o bayescan/3L.hum -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 3L.geste.txt -env phase1.zprecip.txt -o bayescan/3L.precip -threads 25

/home/d669d153/work/bayescenv-new/source/bayescenv 3L.geste.txt -env phase1.ztemp.txt -o bayescan/3Ltemp -threads 25






