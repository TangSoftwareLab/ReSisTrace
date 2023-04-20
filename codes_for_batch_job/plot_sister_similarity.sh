#!/bin/bash
#SBATCH --job-name=r_plot
#SBATCH --account=project_2003482
#SBATCH --output=log/plot_sister_similarity/output_%j.txt
#SBATCH --error=log/plot_sister_similarity/errors_%j.txt
#SBATCH --partition=hugemem_longrun
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1200G

# Load r-env-singularity
module load r-env-singularity/4.1.1

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/<project>" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save plot_sister_similarity.R
seff $SLURM_JOBID
