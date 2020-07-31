#!/bin/sh
#
#SBATCH --job-name=bayescenv.anopheles.temp             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --array=1-5
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1     # Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=20000

#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv ${SLURM_ARRAY_TASK_ID}.geste.txt -env phase1.ztemp.txt -o ${SLURM_ARRAY_TASK_ID}.temp -threads 15


#!/bin/sh
#
#SBATCH --job-name=bayescenv.anopheles.hum             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --array=1-5
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1     # Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=20000

#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv ${SLURM_ARRAY_TASK_ID}.geste.txt -env phase1.zhum.txt -o ${SLURM_ARRAY_TASK_ID}.hum -threads 15


#!/bin/sh
#
#SBATCH --job-name=bayescenv.anopheles.precip             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --array=1-5
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/anoph.phase1     # Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=20000

#run bayescenv on snp matrix
/home/d669d153/work/bayescenv-new/source/bayescenv ${SLURM_ARRAY_TASK_ID}.geste.txt -env phase1.zprecip.txt -o ${SLURM_ARRAY_TASK_ID}.precip -threads 15