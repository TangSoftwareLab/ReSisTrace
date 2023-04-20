#!/bin/bash
#SBATCH --job-name=conserve_gene
#SBATCH --account=project_2003482
#SBATCH --output=log/calculate_conserve_gene/output_%j.txt
#SBATCH --error=log/calculate_conserve_gene/errors_%j.txt
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
srun singularity_wrapper exec Rscript --no-save calculate_conserve_gene.R
seff $SLURM_JOBID
