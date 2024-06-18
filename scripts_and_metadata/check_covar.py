import numpy as np
import pandas as pd
import sys

f_samples = sys.argv[1]
f_phen_covars = sys.argv[2]
outfile = sys.argv[3]

samples = pd.read_csv(f_samples)
f_phen_covars = pd.read_csv(f_phen_covars, sep=" ")

df = f_phen_covars.merge(samples, on="ID", how="inner")

for c in df.columns:
    if c != "ID":
        if df[c].nunique() == 2: #cat covars wil be one hot encoded
            #is categorical covar
            n1 = df[df[c]==0].shape[0]
            n2 = df[df[c]==1].shape[0]
            if float(n1)/float(df.shape[0]) < 0.001:
                df = df[df[c] == 1]
            elif float(n2)/float(df.shape[0]) < 0.001:
                df = df[df[c] == 0]

df.to_csv(outfile, sep=" ", index=False)



