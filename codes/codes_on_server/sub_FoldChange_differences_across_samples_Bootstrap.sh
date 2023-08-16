bootstrapTimes=(10 100 1000 10000) 

for t in ${bootstrapTimes[@]}
do
  sed -e "s/+iter/${t}/g" ./FoldChange_differences_across_samples_Bootstrap.sh | sbatch
done
