import os
import pandas as pd
import numpy as np

data = "/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes/"
f_ldl_data = data+"ldl_cholesterol_direct_30780.txt"
f_statin = data+"statin_use_covariate.txt"
f_aliquot = data+"ldl_cholesterol_direct_aliquot_30782.txt"
f_age = data+"ALL_UKB_age.phen"
f_sex = data+"ALL_UKB_sex.phen"
f_pcs = data+"ALL_UKB_pcs.phen"


ldl = pd.read_csv(f_ldl_data, sep="\t", header=None, skiprows=[0])
ldl.drop(columns=[0, 4], inplace=True)
ldl.rename(columns={1:"eid",2:"ldl1",3:"ldl2"}, inplace=True)

#first drop samples that are na in both ldl1 and ldl2
ldl = ldl[~(ldl["ldl1"].isnull() & ldl["ldl2"].isnull())]

statin = pd.read_csv(f_statin, sep=" ", header=None)
statin.rename(columns={0:"eid", 1:"statin1", 2:"statin2"}, inplace=True)

aliquot = pd.read_csv(f_aliquot, sep="\t", skiprows=[0], header=None)
aliquot.drop(columns=[0, 4], inplace=True)
aliquot.rename(columns={1:"eid",2:"a1",3:"a2"}, inplace=True)


ldl_statin = ldl.merge(statin, on="eid", how="inner")
ldl_statin_aliquot = ldl_statin.merge(aliquot, on="eid", how="inner")


#From J's paper: "preferentially choosing the measurement at the initial 
#assessment when measurements were taken at both visits"

#if ldl1 != NA, keep ldl1 and statin = 0 if statin1 = 0
#elif ldl2 != NA, keep ldl2 and statin = 0 if statin 2 = 0
#else should not happen

def pick_ldl_and_statin_val(e, l1, l2, s1, s2, a1, a2):
    if not np.isnan(l1):
        return [e, l1, 0, s1, a1] #0 for using first LDL measurement
    elif not np.isnan(l2):
        return [e, l2, 1, s2, a2] #1 for using second LDL measurement
    else:
        print("should not happen")
        exit(0)

lsa_condensed = ldl_statin_aliquot.apply(lambda x: pick_ldl_and_statin_val(x.eid, x.ldl1, x.ldl2, x.statin1, x.statin2, x.a1, x.a2), axis=1, result_type='expand')

lsa_condensed.rename(columns={0:"eid", 1:"ldl",2:"repeat",3:"statin", 4:"aliquot"}, inplace=True)

lsa_condensed = lsa_condensed[lsa_condensed['aliquot'].notnull()]

#one hot encode aliquot
lsa_condensed = pd.get_dummies(lsa_condensed, columns=['aliquot'], prefix='', prefix_sep='')

#add age
age = pd.read_csv(f_age, header = None, sep=" ")
age.drop(columns=0, inplace=True)
age.rename(columns={1:"eid", 2:"age"}, inplace=True)
age= age[age['age'].notnull()]

#add sex
sex = pd.read_csv(f_sex, header = None, sep=" ")
sex.drop(columns=0, inplace=True)
sex.rename(columns={1:"eid", 2:"sex"}, inplace=True)
sex = sex[sex['sex'].notnull()]

#add pcs
pcs = pd.read_csv(f_pcs, header = None, sep=" ")
pcs.drop(columns=0, inplace=True)
pcs.rename(columns={1:"eid", 2:"pc1", 3:"pc2",4:"pc3",5:"pc4",6:"pc5",7:"pc6",8:"pc7",9:"pc8",10:"pc9",11:"pc10"}, inplace=True)
pcs.dropna(inplace=True)

ldl_phen_covars = lsa_condensed.merge(age, on="eid", how="inner")
ldl_phen_covars = ldl_phen_covars.merge(sex, on="eid", how="inner")
ldl_phen_covars = ldl_phen_covars.merge(pcs, on="eid", how="inner")
#for associatr pipeline
#np.save(data+'ldl_phen_and_covars.npy', ldl_phen_covars.values, allow_pickle=False)

#for other stuff
colnames = ["ID", "LDL", "repeat", "statin", "a1", "a2","a3","a4", "age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10"]
ldl_phen_covars.columns = colnames
max_ones = 0
max_ones_col = ""
for a in ["a1", "a2", "a3", "a4"]:
    num_ones = np.sum(ldl_phen_covars[a])
    if num_ones > max_ones:
        max_ones = num_ones
        max_ones_col = a
ldl_phen_covars.drop(columns=max_ones_col, inplace=True)
#for one hot encoding for lin reg, can't include all, have to drop one col otherwise will get 
#problems with multicollinearity e.g. corr mat could not be inverted (VIF_INFINITE) with plink
#glm. One col will be perfectly predicted by rest, so need to remove it

ldl_phen_covars.to_csv(data+"ldl_phen_and_covars.txt", sep=" ", index=False)

