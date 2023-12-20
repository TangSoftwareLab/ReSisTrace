#!/bin/bash -l
#SBATCH --job-name=+samp
#SBATCH --account=project_2003482
#SBATCH --output=log/extract_lineage_label/+samp_output_%j.txt
#SBATCH --error=log/extract_lineage_label/+samp_errors_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1

module load biokit

samtools view /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/cellranger/+samp/outs/possorted_genome_bam.bam pBA439_UMI_20 | awk '{{for (i=1; i<=NF; ++i)      {{if($i ~ "^CB:Z:"){{print $10,",",substr($12,6),",",substr($i,6), ",", substr($(i+3),6)}}}}}}' > /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/extract_lineage_label/+samp.txt

seff $SLURM_JOBID
