for i in {6..20}; do
 echo $i
 #careful running following line, commented out to prevent accidental runs but needed during original creation of merged snp/str datva
 #mv merged/chr${i}_snpstrs.pvar merged/chr${i}_snpstrs_orig.pvar
 python fix_dup_ids.py merged/chr${i}_snpstrs_orig.pvar merged/chr${i}_snpstrs.pvar
done
