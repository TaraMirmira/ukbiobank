import pandas as pd
import sys

outfile = sys.argv[1]

snp_pgen_dir="/expanse/projects/gymreklab/CAST/ukbb/array_imputed/pfile_converted/"
str_pgen_dir="/expanse/projects/gymreklab/tmirmira/ukbiobank/str_vcfs/"

#confirmed the psam files are identical across snp pgens and identical across str pgens so looking
#at chr 1 for each is sufficient

snp_samples = pd.read_csv(snp_pgen_dir+"chr"+str(1)+".psam", sep="\t")['IID'].values

str_samples = set(pd.read_csv(str_pgen_dir+"chr"+str(1)+"strs.psam", sep="\t")['#IID'].values)


#done specifically to maintain sampels in the order they are in the snp files
#the str files are smaller so they can be reordered using plink
with open(outfile, 'w') as f:
    f.write("IID\n")
    for item in snp_samples:
        if item in str_samples:
            f.write("%s\n" % item)



