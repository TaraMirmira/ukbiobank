import sys
from pgenlib import PgenWriter
from pgenlib import PgenReader
import numpy as np


snps_infile = sys.argv[1]
strs_infile = sys.argv[2]
outfile = sys.argv[3]
snp_ct = int(sys.argv[4])
str_ct = int(sys.argv[5])
variant_ct = snp_ct + str_ct
num_samples = int(sys.argv[6])


pgen_writer = PgenWriter(bytes(str(outfile), "utf8"), 
        num_samples, 
        variant_ct=variant_ct, 
        dosage_present=True)


dosages = np.empty(num_samples, dtype=np.float32)

reader1 = PgenReader(bytes(str(snps_infile), "utf8"))
for i in range(0, snp_ct):
    reader1.read_dosages(i, dosages)
    #pgen_writer.append_dosages_batch(dosages.reshape((1,dosages.shape[0])))
    pgen_writer.append_dosages(dosages)
reader1.close()


reader2 = PgenReader(bytes(str(strs_infile), "utf8"))
for i in range(0, str_ct):
    reader2.read_dosages(i, dosages)
    #pgen_writer.append_dosages_batch(dosages.reshape((1,dosages.shape[0])))
    pgen_writer.append_dosages(dosages)
reader2.close()

pgen_writer.close()
