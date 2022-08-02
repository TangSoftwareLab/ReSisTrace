subclones=(A B C D E) 

for t in ${subclones[@]}
do
  sed -e "s/+subclone/${t}/g" ./CNV_subclone_DEG.sh | sbatch
done
