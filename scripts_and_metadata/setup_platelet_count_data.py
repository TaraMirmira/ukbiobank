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


with open(f_device, 'r') as file:
    # Read the header
    header = file.readline().strip().split('\t')
    # Read the remaining lines
    lines = file.readlines()
    # Strip leading tabs and split by tab
    device_data = [line.lstrip('\t').strip().split('\t') for line in lines]

# Create DataFrame
device = pd.DataFrame(device_data, columns=header)
device.rename(columns={"30083-0.0":"device1", "30083-1.0":"device2","30083-2.0":"device3"}, inplace=True)
device['eid'] = device['eid'].astype(int)

platelet = platelet.merge(device, on="eid", how="inner")


#From J's paper: "preferentially choosing the measurement at the initial 
#assessment when measurements were taken at both visits"

def pick_platelet_and_device_val(e, p1, p2, p3, d1, d2, d3):
    if not pd.isna(p1):
        return [e, p1, 0, d1] #0 for using first platelet measurement
    elif not pd.isna(p2):
        return [e, p2, 1, d2] #1 for using second platelet measurement
    elif not pd.isna(p3):
        return [e, p3, 2, d3] #2 for using third platelet measurement
    else:
        print("should not happen")
        exit(0)

pcd_condensed = platelet.apply(lambda x: pick_platelet_and_device_val(x.eid, x.platelet1, x.platelet2, x.platelet3, x.device1, x.device2, x.device3), axis=1, result_type='expand')


pcd_condensed.rename(columns={0:"eid", 1:"platelet_count",2:"repeat",3:"device"}, inplace=True)

pcd_condensed = pcd_condensed[pcd_condensed['device'].notnull()]

#one hot encode device and repeat
pcd_condensed = pd.get_dummies(pcd_condensed, columns=['device'], prefix='', prefix_sep='')
pcd_condensed = pd.get_dummies(pcd_condensed, columns=['repeat'], prefix='', prefix_sep='')

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

pc_phen_covars = pcd_condensed.merge(age, on="eid", how="inner")
pc_phen_covars = pc_phen_covars.merge(sex, on="eid", how="inner")
pc_phen_covars = pc_phen_covars.merge(pcs, on="eid", how="inner")


colnames = ["ID", "platelet_count", "d1", "d2", "d3", "d4", "r1", "r2", "r3", "age", "sex", "pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10"]
pc_phen_covars.columns = colnames
max_ones = 0
max_ones_col = ""
for a in ["d1", "d2", "d3", "d4"]:
    num_ones = np.sum(pc_phen_covars[a])
    if num_ones > max_ones:
        max_ones = num_ones
        max_ones_col = a
pc_phen_covars.drop(columns=max_ones_col, inplace=True)
#for one hot encoding for lin reg, can't include all, have to drop one col otherwise will get 
#problems with multicollinearity e.g. corr mat could not be inverted (VIF_INFINITE) with plink
#glm. One col will be perfectly predicted by rest, so need to remove it

#repeat for repeat
max_ones = 0
max_ones_col = ""
for a in ["r1", "r2", "r3"]:
    num_ones = np.sum(pc_phen_covars[a])
    if num_ones > max_ones:
        max_ones = num_ones
        max_ones_col = a
pc_phen_covars.drop(columns=max_ones_col, inplace=True)

pc_phen_covars.to_csv(data+"platelet_count_phen_and_covars.txt", sep=" ", index=False)

