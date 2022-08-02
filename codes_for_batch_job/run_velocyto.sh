inputdir="/scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/cellranger/"
logdir="./log/velocyto_log/"
samples=`ls -1 ${inputdir} | grep -v "log" | grep -v "sh$"`

mkdir ${logdir} 

for sample in ${samples}
do
  echo ${sample}
  sed s/+samp/${sample}/g /scratch/project_2003482/2021.11.26_KURAMOCHI_AML_lineage_tracing/code/runVelo.sh | sed s@+dir@${inputdir}@g | sed s@+logdir@${logdir}@g | sbatch
done
