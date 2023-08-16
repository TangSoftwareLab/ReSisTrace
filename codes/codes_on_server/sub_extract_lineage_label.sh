dir="/scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/raw/LT-AML-Kuramochi-031121/"
samples=`ls -1 ${dir}`
echo ${samples[@]/FASTQC}

for sample in ${samples[@]/FASTQC}
do
    sed s/+samp/${sample}/g run_extract_lineage_label.sh | sed s@+dir@${dir}@g | sbatch
done
