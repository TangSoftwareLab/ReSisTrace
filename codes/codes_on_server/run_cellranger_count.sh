#!/bin/bash
#SBATCH --job-name=+samp
#SBATCH --export=ALL
#SBATCH --account=project_2003482
#SBATCH --output=log/+samp_output_%j.txt
#SBATCH --error=log/+samp_errors_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=hugemem
#SBATCH --nodes=1 --ntasks-per-node=40
#SBATCH --signal=2
#SBATCH --mem=1200G

export PATH=/scratch/project_2003482/cellranger/cellranger-5.0.1:$PATH

cellranger count --id=+samp \
    --fastqs=+dir+samp \
    --sample=+samp \
    --transcriptome=/scratch/project_2003482/genome_with_lineage_label/GRCH38_AnnV25_label/ \
    --jobmode=local --localcores=40

seff $SLURM_JOBID
