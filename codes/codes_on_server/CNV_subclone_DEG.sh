#!/bin/bash -l
#SBATCH --job-name=DEG_+subclone
#SBATCH --account=project_2003482
#SBATCH --output=log/CNV_subclone_DEG/+subclone_output_%j.txt
#SBATCH --error=log/CNV_subclone_DEG/+subclone_errors_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
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
srun singularity_wrapper exec Rscript --no-save CNV_subclone_DEG.R +subclone
seff $SLURM_JOBID
