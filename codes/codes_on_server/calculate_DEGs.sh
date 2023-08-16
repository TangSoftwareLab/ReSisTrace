#!/bin/bash
#SBATCH --job-name=deg_+samp
#SBATCH --account=project_2003482
#SBATCH --output=log/calculate_DEGs/+samp_output_%j.txt
#SBATCH --error=log/calculate_DEGs/+samp_errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=350G

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
srun singularity_wrapper exec Rscript --no-save calculate_DEGs.R +samp
seff $SLURM_JOBID

