#bash phenotype_qc.sh creatinine_phen_and_covars.txt creatinine ../samples/AFR_unrelated_samples.txt AFR_norm_creatinine_phen_covars.txt
#bash phenotype_qc.sh creatinine_phen_and_covars.txt creatinine ../samples/EUR_unrelated_samples.txt EUR_norm_creatinine_phen_covars.txt



infile=creatinine_phen_and_covars.txt
phenotype=egfr

scripts=/expanse/projects/gymreklab/tmirmira/ukbiobank/scripts_and_metadata
data=/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes

echo $data/$infile

#AFR
python $scripts/check_covar.py ../samples/AFR_unrelated_samples.txt $data/$infile temp
python $scripts/transform_creatinine.py temp AFR $data/AFR_transf_creatinine.txt
python $scripts/normalize_data.py $data/AFR_transf_creatinine.txt $phenotype $data/AFR_norm_egfr_phen_covars.txt
rm temp


#EUR
python $scripts/check_covar.py ../samples/EUR_unrelated_samples.txt $data/$infile temp
python $scripts/transform_creatinine.py temp EUR $data/EUR_transf_creatinine.txt
python $scripts/normalize_data.py $data/EUR_transf_creatinine.txt $phenotype $data/EUR_norm_egfr_phen_covars.txt
rm temp
