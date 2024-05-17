import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]


inpvar = pd.read_csv(infile, sep="\t")

def add_alt_ref(rsid, ref, alt):
    ref_alt = ref + "_" + alt
    if rsid.endswith(ref_alt):
        return rsid
    else:
        return rsid+"_"+ref+"_"+alt

inpvar['ID'] = inpvar.apply(lambda x: add_alt_ref(x.ID, x.REF, x.ALT), axis=1)

inpvar.to_csv(outfile, sep="\t", index=False)


