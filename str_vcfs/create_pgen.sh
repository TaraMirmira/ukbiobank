chr=$1

vcf=chr${chr}.vcf.gz

#zcat $vcf | awk -F $'\t' -v 'OFS=\t' '!/^#/ {$4="A"; $5="T";} /;PERIOD=/ {print $1, $2, $3, $4, $5;}' > chr${chr}strs.pvar

bcftools view -h $vcf > chr${chr}strs.pvar

zcat $vcf | awk -F $'\t' -v OFS='\t' '!/^#/ && /;PERIOD=/ { $4="A"; $5="T"; print $1, $2, $3, $4, $5; }' >> chr${chr}strs.pvar


variantct=$(grep -Ev '^#' chr${chr}strs.pvar | wc -l)

(echo -e "#IID\tSEX" && ( zcat $vcf | grep -m1 -E '^#CHROM' | cut -f10- | sed 's/\t/\tNA\n/g;s/$/\tNA/' )) > chr${chr}strs.psam

python batch_add_str_dosages.py $vcf chr${chr}strs.pgen chr${chr}minmax $variantct
