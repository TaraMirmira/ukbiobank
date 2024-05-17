for i in {6..20}; do
 echo $i
 python fix_dup_ids.py merged/chr${i}_snpstrs_orig.pvar merged/chr${i}_snpstrs.pvar
done
