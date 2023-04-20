#!/bin/bash
#SBATCH --job-name=scvelo
#SBATCH --account=project_2003482
#SBATCH --output=log/scvelo/output_%j.txt
#SBATCH --error=log/scvelo/errors_%j.txt
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=300000

# Load  bioconda
export PROJAPPL=/projappl/project_2003482
module load bioconda
source activate scvelo

# Run the py script
python3 runscvelo.py
