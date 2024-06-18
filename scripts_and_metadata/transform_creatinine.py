import pandas as pd
import sys

infile = sys.argv[1]
pop = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(infile, sep=" ")
isblack = (pop == "AFR")

def egfr_calc(scr, ismale, isblack, age):
    k = 0.9 if ismale else 0.7
    a = -0.411 if ismale else -0.429
    egfr = 141 * (min(scr / k, 1)**a) * (max(scr / k, 1)**-1.209) * (0.993**age)
    if not ismale:
        egfr *= 1.1018
    if isblack:
        egfr *= 1.159
    return egfr

egfr = df.apply(lambda x: egfr_calc(x.creatinine, x.sex, isblack, x.age), axis = 1)

df['creatinine'] = egfr
df.rename(columns={'creatinine':'egfr'}, inplace=True)
df.to_csv(outfile, sep=" ", index=False)


