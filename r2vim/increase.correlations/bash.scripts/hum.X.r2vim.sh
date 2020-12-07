#!/bin/sh
#
#SBATCH --job-name=hum.X.r2vim             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/r2vim        # Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=10000

module load R

R -e "Sys.setenv(RSTUDIO_PANDOC='/panfs/pfs.local/work/bi/bin/pandoc/bin');  rmarkdown::render('hum.X.r2vim.Rmd',output_file='hum.X.r2vim.results.html')"

