import pandas as pd
import scipy.stats as stats
import sys



def Inverse_Quantile_Normalization(M):

    #After quantile normalization of samples, standardize expression of each gene
    M = M.transpose()
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    Q = Q.transpose()
    return Q


f_phen = sys.argv[1]
colname = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(f_phen, sep=" ")
df[colname] = Inverse_Quantile_Normalization(df[[colname]])
df['ID'] = df['ID'].astype('int')
df.insert(0, 'FID', df.iloc[:, 0])
df.rename(columns={'ID':'IID'}, inplace=True)
df.to_csv(outfile, sep=" ", index=False)
