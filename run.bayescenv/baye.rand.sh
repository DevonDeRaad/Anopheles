#!/bin/sh
#
#SBATCH --job-name=bayescenv.anopheles.rand             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --array=1-5
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1/bayescenv     # Set working d$
#SBATCH --mem-per-cpu=1g            # memory requested
#SBATCH --time=20000

#I renamed each geste file according to the following conventions to fit the 'array' variable 1-5
#2R.geste.txt > 1.geste.txt
#2L.geste.txt > 2.geste.txt
#3R.geste.txt > 3.geste.txt
#3L.geste.txt > 4.geste.txt
#X.geste.txt > 5.geste.txt

#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv ${SLURM_ARRAY_TASK_ID}.geste.txt -env phase1.zrand.txt -o ${SLURM_ARRAY_TASK_ID}.rand -threads 15


