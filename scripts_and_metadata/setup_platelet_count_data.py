import os
import pandas as pd
import numpy as np

data = "/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes/"
f_platelet_data = data+"platelet_count_30080.txt"
f_device = data+"platelet_count_device_id_30083.txt"
f_age = data+"ALL_UKB_age.phen"
f_sex = data+"ALL_UKB_sex.phen"
f_pcs = data+"ALL_UKB_pcs.phen"


platelet = pd.read_csv(f_platelet_data, sep="\t", header=None, skiprows=[0])
platelet.drop(columns=[0, 5], inplace=True)
platelet.rename(columns={1:"eid",2:"platelet1",3:"platelet2",4:"platelet3"}, inplace=True)

#first drop samples that are na in all cols
platelet = platelet[~(platelet["platelet1"].isnull() & platelet["platelet2"].isnull() & platelet["platelet3"].isnull())]

device = pd.read_csv(f_device, sep="\t", header=None)
device.drop(columns=[0, 5], inplace=True)
device.rename(columns={1:"eid", 2:"device1", 3:"device2",4:"device3"}, inplace=True)

platelet = platelet.merge(device, on="eid", how="inner")


#From J's paper: "preferentially choosing the measurement at the initial 
#assessment when measurements were taken at both visits"

def pick_platelet_and_device_val(e, p1, p2, p3, d1, d2, d3):
    if not np.isnan(p1):
        return [e, p1, 0, d1] #0 for using first platelet measurement
    elif not np.isnan(l2):
        return [e, p2, 1, d2] #1 for using second platelet measurement
    elif not np.isnan(p3):
        return [e, p3, 2, d3] #2 for using third platelet measurement
    else:
        print("should not happen")
        exit(0)

platelet_condensed = platelet.apply(lambda x: pick_platelet_and_device_val(x.eid, x.platelet1, x.platelet2, x.platelet3, x.device1, x.device2, x.device3), axis=1, result_type='expand')

print(platelet_condensed)

exit(0)


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

ldl_phen_covars.to_csv(data+"platelet_count_phen_and_covars.txt", sep=" ", index=False)

