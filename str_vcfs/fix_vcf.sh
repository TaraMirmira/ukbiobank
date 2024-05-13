chr=$1


orig_vcf=/expanse/projects/gymreklab/jmargoli/ukbiobank/str_imputed/runs/first_pass/vcfs/annotated_strs/chr$chr.vcf.gz
new_vcf_loc=/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs

bcftools query -l $orig_vcf > $new_vcf_loc/chr${chr}_sample_ids

awk -F_ '{print $1}' chr${chr}_sample_ids > chr${chr}_correct_sample_ids
bcftools reheader -s chr${chr}_correct_sample_ids $orig_vcf -o $new_vcf_loc/temp_chr$chr.vcf.gz



bcftools view --header-only $new_vcf_loc/temp_chr$chr.vcf.gz > chr${chr}_orig_header
awk 'NR==2 {print "##command=HipSTR-v0.6.2"} {print}' chr${chr}_orig_header > chr${chr}_new_header

bcftools reheader -h chr${chr}_new_header $new_vcf_loc/temp_chr$chr.vcf.gz -o $new_vcf_loc/chr$chr.vcf.gz

tabix -p vcf chr$chr.vcf.gz
rm $new_vcf_loc/temp_chr$chr.vcf.gz 
rm chr${chr}_correct_sample_ids
rm chr${chr}_orig_header
rm chr${chr}_new_header
