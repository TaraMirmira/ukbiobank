import os
import pandas as pd
import numpy as np

data = "/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes/"
f_creatinine_data = data+"creatinine_30700.txt"
f_aliquot = data+"creatinine_aliquot_30702.txt"
f_age = data+"ALL_UKB_age.phen"
f_sex = data+"ALL_UKB_sex.phen"
f_pcs = data+"ALL_UKB_pcs.phen"


with open(f_creatinine_data, 'r') as file:
    header = file.readline().strip().split('\t')
    lines = file.readlines()
    creatinine_data = [line.lstrip('\t').strip().split('\t') for line in lines]

creatinine = pd.DataFrame(creatinine_data, columns=header)
creatinine.rename(columns={"30700-0.0":"creatinine1","30700-1.0":"creatinine2"}, inplace=True)
creatinine['eid'] = creatinine['eid'].astype(int)

#first drop samples that are na in all cols
creatinine = creatinine[~(creatinine["creatinine1"].isnull() & creatinine["creatinine2"].isnull())]

with open(f_aliquot, 'r') as file:
    header = file.readline().strip().split('\t')
    lines = file.readlines()
    aliquot_data = [line.lstrip('\t').strip().split('\t') for line in lines]

# Create DataFrame
aliquot= pd.DataFrame(aliquot_data, columns=header)
aliquot.rename(columns={"30702-0.0":"aliquot1", "30702-1.0":"aliquot2"}, inplace=True)
aliquot['eid'] = aliquot['eid'].astype(int)

creatinine= creatinine.merge(aliquot, on="eid", how="inner")


#From J's paper: "preferentially choosing the measurement at the initial 
#assessment when measurements were taken at both visits"

def pick_creatinine_and_aliquot_val(e, c1, c2, a1, a2):
    if not pd.isna(c1):
        return [e, c1, 0, a1] #0 for using first platelet measurement
    elif not pd.isna(c2):
        return [e, c2, 1, a2] #1 for using second platelet measurement
    else:
        print("should not happen")
        exit(0)


ca_condensed = creatinine.apply(lambda x: pick_creatinine_and_aliquot_val(x.eid, x.creatinine1, x.creatinine2, x.aliquot1, x.aliquot2), axis=1, result_type='expand')


ca_condensed.rename(columns={0:"eid", 1:"creatinine",2:"repeat",3:"aliquot"}, inplace=True)

ca_condensed = ca_condensed[ca_condensed['aliquot'].notnull()]

#one hot encode device and repeat
ca_condensed = pd.get_dummies(ca_condensed, columns=['aliquot'], prefix='a', prefix_sep='', dtype=int, drop_first=True)
ca_condensed = pd.get_dummies(ca_condensed, columns=['repeat'], prefix='', prefix_sep='', dtype=int, drop_first=True)

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

ca_phen_covars = ca_condensed.merge(age, on="eid", how="inner")
ca_phen_covars = ca_phen_covars.merge(sex, on="eid", how="inner")
ca_phen_covars = ca_phen_covars.merge(pcs, on="eid", how="inner")


colnames = ["ID", "creatinine", "a1", "a2", "a3", "a4", "a5", "age", "sex", "pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10"]
ca_phen_covars.columns = colnames
#max_ones = 0
#max_ones_col = ""
#for a in ["d1", "d2", "d3", "d4"]:
#    num_ones = np.sum(pc_phen_covars[a])
#    if num_ones > max_ones:
#        max_ones = num_ones
#        max_ones_col = a
#pc_phen_covars.drop(columns=max_ones_col, inplace=True)
#for one hot encoding for lin reg, can't include all, have to drop one col otherwise will get 
#problems with multicollinearity e.g. corr mat could not be inverted (VIF_INFINITE) with plink
#glm. One col will be perfectly predicted by rest, so need to remove it

ca_phen_covars.to_csv(data+"creatinine_phen_and_covars.txt", sep=" ", index=False)

