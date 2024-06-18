#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp268
#SBATCH --job-name=ids
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --output="/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs/ids.%j.%N.out"
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tmirmira
#SBATCH --export=ALL
source ~/.bashrc
source /etc/profile.d/modules.sh
conda activate env_popcorn2
cd /expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs
bash fix_dup_ids_all.sh
