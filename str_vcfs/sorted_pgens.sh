chr=$1
  
#for strs
plink2 --pfile chr${chr}strs --keep intersect_all_samples.txt --indiv-sort f intersect_all_samples.txt --out tomerge/chr${chr}strs_sorted --make-pgen

#for snps
plink2 --pfile /expanse/projects/gymreklab/CAST/ukbb/array_imputed/pfile_converted/chr${chr} --keep intersect_all_samples_fid_iid.txt --out tomerge/chr${chr}snps_sorted --make-pgen

