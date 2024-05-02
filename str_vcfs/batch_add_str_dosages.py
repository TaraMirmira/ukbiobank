import trtools.utils.tr_harmonizer as trh
import numpy as np
import sys
import cyvcf2
from cyvcf2 import VCF
from pgenlib import PgenWriter
import math


infile = sys.argv[1]
outfile = sys.argv[2]
min_max_len_file = sys.argv[3]
variant_ct = int(sys.argv[4])


vcf_reader = VCF(infile)
harmonizer = trh.TRRecordHarmonizer(vcf_reader)

pgen_writer = PgenWriter(bytes(str(outfile), "utf8"), 
        len(vcf_reader.samples), 
        variant_ct=variant_ct, 
        dosage_present=True)


min_max_len = []

batch_size = 1000
num_batches = math.ceil(variant_ct/batch_size)


#I think it safe to assume num variants will always be larger than batch size so we have at least one batch
dosages = np.empty((batch_size, len(vcf_reader.samples)), dtype=np.float32)
num_variants_processed_batch = 0
num_variants_processed = 0
for idx, trrecord in enumerate(harmonizer):
    ap1 = trrecord.vcfrecord.format('AP1')
    ref1 = 1 - np.sum(ap1, axis=1)
    ref1 = np.clip(ref1, 0, 1) #if negative due to rounding errors, cutoff at 0
    ap2 = trrecord.vcfrecord.format('AP2')
    ref2 = 1 - np.sum(ap2, axis=1)
    ref2 = np.clip(ref2, 0, 1)

    lens = np.array(trrecord.alt_allele_lengths)
    ref_len = trrecord.ref_allele_length
    min_len = min(min(lens), ref_len)
    max_len = max(max(lens), ref_len)
    lens = lens.reshape((lens.shape[0],1))

    h1_probs = np.dot(ap1, lens)
    h1_probs = np.clip(h1_probs, 0, max(lens))
    h2_probs = np.dot(ap2, lens)
    h2_probs = np.clip(h2_probs, 0, max(lens))

    ref1_prob = ref1 * ref_len
    ref2_prob = ref2 * ref_len
    ref1_prob = ref1_prob.reshape(ref1_prob.shape[0], 1)
    ref2_prob = ref2_prob.reshape(ref2_prob.shape[0], 1)

    d = h1_probs + h2_probs + ref1_prob + ref2_prob
    #d = np.round(d, 2)
    d_transf = (d - 2*min_len) / (2*max_len - 2*min_len) 
    if (np.any(d_transf >= 1.1) or np.any(d_transf <= -0.1)):
        assert(False)
    d_transf = np.clip(d_transf, 0, 1)
    #assert(np.all(d_transf<= 1) and np.all(d_transf >=0))

    dosages[num_variants_processed_batch] = d_transf.flatten()

    str_id = trrecord.record_id
    min_max_len.append([str_id, min_len, max_len])

    num_variants_processed += 1
    num_variants_processed_batch += 1
    if ((num_variants_processed % batch_size == 0) or (num_variants_processed == variant_ct)):
        print("done with a batch")
        print("num variants processed = ", num_variants_processed)
        next_batch = min(batch_size, variant_ct - num_variants_processed)
        pgen_writer.append_dosages_batch(dosages[:num_variants_processed_batch])
        num_variants_processed_batch = 0 #reset


np.savetxt(min_max_len_file, np.array(min_max_len, dtype=str), fmt='%s')
assert(num_variants_processed == variant_ct)
pgen_writer.close()
vcf_reader.close()
