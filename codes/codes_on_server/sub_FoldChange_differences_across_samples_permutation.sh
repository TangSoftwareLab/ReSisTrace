permutationTimes=(10 100 1000 10000) 

for t in ${permutationTimes[@]}
do
  sed -e "s/+iter/${t}/g" ./FoldChange_differences_across_samples_permutation.sh | sbatch
done
