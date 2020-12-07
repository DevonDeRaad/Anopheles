#!/bin/sh
#
#SBATCH --job-name=prec.3R.r2vim             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=10               # CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/r2vim        # Set working d$
#SBATCH --mem-per-cpu=3gb            # memory requested
#SBATCH --time=50000

module load R
R -e "Sys.setenv(RSTUDIO_PANDOC='/panfs/pfs.local/work/bi/bin/pandoc/bin');  rmarkdown::render('prec.3R.r2vim.Rmd',output_file='prec.3R.r2vim.results.html')"

