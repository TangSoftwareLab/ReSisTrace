#!/bin/bash -l
#SBATCH --job-name=+samp
#SBATCH --account=project_2003482
#SBATCH --output=log/umi_tools/+samp_output_%j.txt
#SBATCH --error=log/umi_tools/+samp_errors_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1

export PROJAPPL=/projappl/project_2003482
module load bioconda
conda activate umi_tools

python3 umitools_cross_sample.py /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/extract_lineage_label/sorted/+samp/ /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/umi_tools/+samp/ /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/code/log/umi_tools/python_log/+samp.txt "directional"

seff $SLURM_JOBID
