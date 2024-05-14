chr=$1

tomerge=/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs/tomerge
merged=/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs/merged
snps=${tomerge}/chr${chr}snps_sorted
strs=${tomerge}/chr${chr}strs_sorted
outpre=$merged/chr${chr}_snpstrs

a=$(wc -l < ${snps}.pvar)
variantct_snps=$((a - 1))
variantct_strs=$(grep -Ev '^#' ${strs}.pvar | wc -l)
#variantct=$((variantct_strs + variantct_sps - 1))

a=$(wc -l < ${snps}.psam)
num_samples=$((a - 1))

#(echo -e "#IID\tSEX" && ( zcat $vcf | grep -m1 -E '^#CHROM' | cut -f10- | sed 's/\t/\tNA\n/g;s/$/\tNA/' )) > chr${chr}strs.psam

cp $snps.psam $outpre.psam

#cp $snps.pvar $outpre.pvar
#awk '!/^#/ { print >> "'"$outpre"'.pvar" }' "$strs.pvar"
awk '!/^##/' $strs.pvar > temp$chr.pvar


python batch_merge_pgen.py $snps.pgen $strs.pgen $snps.pvar temp$chr.pvar $variantct_snps $variantct_strs $outpre.pgen $outpre.pvar $num_samples
rm temp$chr.pvar
