instructions:

snp pgens located at /expanse/projects/gymreklab/CAST/ukbb/array_imputed/pfile_converted/

starting input is imputed STRs from Johnathan located at /expanse/projects/gymreklab/jmargoli/ukbiobank/str_imputed/runs/first_pass/vcfs/annotated_strs/chr$chr.vcf.gz

run: fix_vcf.sh i for i in 1 to 22

this will create chr$i.vcf in tmirmira/ukbiobank/str_vcfs

next create STR pgens

run: create_pgen.sh i for i in 1 to 22

this will create chr$istrs.pgen/psam/pvar in tmirmira/ukbiobank/str_vcfs

next get the intersection of all samples across all chromosomes

run: python intersect_all_chrs_samples.py

The output of the above line will match the order of the snp pgens. This means we just need to reorder
the str pgens (done in this way because str pgens are smaller, don't want to have to reorder snp pgens)

next: reorder str pgens to match sample order of snp pgens and restrict snps to samples in intersection (see above)

run: sorted_pgens.sh i for i in 1 to 22
^may need to run the lines for snps vs strs separately due to memory requirements being much larger for snps (can require 128GB for snps esp chr 1, chr 2, etc)

At this point, will have snp pgens and str pgens ready for merging in directory called tomerge

next merge snp and str pgens

run: merge_pgen.sh i for i in 1 to 22

merged pgen written to directory called merged

last step: fix dup ids by constructing rsids that are rsid+ref+alt

run: bash fix_dup_ids_all.sh (this will modify the pvar files 


The strs have been condensed into dosages. The original min and max len of each str can be found in chr${i}minmax in tmirmira/ukbiobank/str_vcfs



