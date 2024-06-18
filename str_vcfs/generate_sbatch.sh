chr=$1

echo "#!/usr/bin/env bash" >> batch_job.sh


echo '#SBATCH --partition=ind-shared' >> batch_job.sh
echo '#SBATCH --account=ddp268' >> batch_job.sh
echo "#SBATCH --job-name=pgen$chr" >> batch_job.sh
echo '#SBATCH --nodes=1' >> batch_job.sh
echo '#SBATCH --mem=128G' >> batch_job.sh
echo '#SBATCH --ntasks-per-node=1' >> batch_job.sh
echo '#SBATCH --time=05:00:00' >> batch_job.sh
echo "#SBATCH --output=\"/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs/pgen$chr.%j.%N.out\"" >> batch_job.sh
echo '#SBATCH --mail-type=FAIL' >> batch_job.sh
echo '#SBATCH --mail-user=tmirmira' >> batch_job.sh
echo '#SBATCH --export=ALL' >> batch_job.sh

echo 'source ~/.bashrc' >> batch_job.sh

echo 'source /etc/profile.d/modules.sh' >> batch_job.sh
echo 'conda activate env_popcorn2' >> batch_job.sh

echo 'cd /expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs' >> batch_job.sh
#echo "bash fix_vcf.sh $chr" >> batch_job.sh
#echo "bash create_pgen.sh $chr" >> batch_job.sh
#echo "bash sorted_pgens.sh $chr" >> batch_job.sh
echo "bash merge_pgen.sh $chr" >> batch_job.sh
