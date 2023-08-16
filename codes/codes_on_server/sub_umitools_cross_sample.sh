dir="/scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/extract_lineage_label/sorted/"

samples=`ls -1 ${dir}`
echo ${samples}

for sample in ${samples}
do
    sed s/+samp/${sample}/g run_umitools_cross_sample.sh | sbatch
done
