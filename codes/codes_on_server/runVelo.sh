#!/bin/bash
#SBATCH --job-name=vt+samp
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --account=project_2003482
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH -o +logdir+samp_out_%j.txt
#SBATCH -e +logdir+samp_err_%j.txt
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail

export PROJAPPL=/projappl/project_2003482
module load bioconda
source activate velocyto
module load samtools/1.7

samtools sort -l 7 -m 2048M -t CB -O BAM -@ 3 -o /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/cellranger/+samp/outs/cellsorted_possorted_genome_bam.bam /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/cellranger/+samp/outs/possorted_genome_bam.bam

velocyto run10x -vvv -m /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/raw/mm10_rmsk.gtf +dir/+samp/ /scratch/project_2003482/genome_with_lineage_label/GRCH38_AnnV25_label/genes/genes.gtf
