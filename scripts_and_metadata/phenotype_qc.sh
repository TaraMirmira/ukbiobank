infile=$1
phenotype=$2
samples=$3
outfile=$4

scripts=/expanse/projects/gymreklab/tmirmira/ukbiobank/scripts_and_metadata
data=/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes

echo $data/$infile

python $scripts/check_covar.py $samples $data/$infile temp
python $scripts/normalize_data.py temp $phenotype $data/$outfile
rm temp


