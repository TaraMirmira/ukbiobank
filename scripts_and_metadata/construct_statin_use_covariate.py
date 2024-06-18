import pandas as pd
import numpy as np

df = pd.read_csv("/expanse/projects/gymreklab/tmirmira/ukbiobank/phenotypes/ukb46782_medications_field20003.csv")
x = np.loadtxt("/expanse/projects/gymreklab/tmirmira/ukbiobank/scripts_and_metadata/statin_encoding.txt", dtype=str)
x = list(x)
x = [float(a) for a in x]

instance0 = df.filter(regex='^20003-0.*|'+'eid')
instance1 = df.filter(regex='^20003-1.*|'+'eid')


df_statin_instance0 = instance0.isin(x)
df_statin_instance1 = instance1.isin(x)


#if any statins were measured at the first instance
statin_use0 = df_statin_instance0.apply(lambda row: row.any(), axis=1)
#if any statins were measured at the second instance
statin_use1 = df_statin_instance1.apply(lambda row: row.any(), axis=1)

covar = pd.DataFrame({'eid':df['eid'],'statin_use0':statin_use0, 'statin_use1':statin_use1})
covar = covar.astype('int')

#for use with plink
covar.to_csv("statin_use_covariate.txt", sep=" ", header=False, index=False)

#for use with Jonathan's pipeline
np.save('statin_use_covariate.npy', covar.values, allow_pickle=False)
