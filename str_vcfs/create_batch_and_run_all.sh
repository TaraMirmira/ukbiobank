for i in {1..21}; do
 echo $i
 bash generate_sbatch.sh $i
 sbatch batch_job.sh
 rm batch_job.sh
done

