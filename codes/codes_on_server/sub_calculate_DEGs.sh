dir="/scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/umi_tools/"

samples=`ls -1 ${dir}`
echo ${samples}

for S in ${samples}
do
  sed -e "s/+samp/${S}/g" ./calculate_DEGs.sh | sbatch
done
