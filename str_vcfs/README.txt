instructions:

snp pgens located at /expanse/projects/gymreklab/CAST/ukbb/array_imputed/pfile_converted/

starting input is imputed STRs from Johnathan located at /expanse/projects/gymreklab/jmargoli/ukbiobank/str_imputed/runs/first_pass/vcfs/annotated_strs/chr$chr.vcf.gz

run: fix_vcf.sh i for i in 1 to 22


next create STR pgens

run: create_pgen.sh i for i in 1 to 22

next get the intersection of all samples across all chromosomes

run: python intersect_all_chrs_samples.py

The output of the above line will match the order of the snp pgens. This means we just need to reorder
the str pgens (done in this way because str pgens are smaller, don't want to have to reorder snp pgens)




