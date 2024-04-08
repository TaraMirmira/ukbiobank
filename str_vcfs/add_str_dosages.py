import trtools.utils.tr_harmonizer as trh
import numpy as np
import sys
from cyvcf2 import VCF


infile = sys.argv[1]
outfile = sys.argv[2]
min_max_len_file = sys.argv[3]


vcf_reader = VCF(infile)
harmonizer = trh.TRRecordHarmonizer(vcf_reader)

# Update FORMAT field to include GT:DS
vcf_reader.add_format_to_header({'ID': 'DS', 'Description': 'Dosage', 'Type': 'Float', 'Number': '1'})
vcf_reader.add_format_to_header({'ID': 'GT:DS', 'Description': 'Genotype:Dosage', 'Type': 'String', 'Number': '2'})


vcf_writer = cyvcf2.Writer(outfile, invcf)

def get_gt(d):
    if d < 0.5:
        return "0/0"
    else:
        return "0/1"

v_get_gt = np.vectoriz(get_gt)

min_max_len = []

for trrecord in harmonizer:
    ap1 = trrecord.vcfrecord.format('AP1')
    ref1 = 1 - np.sum(ap1, axis=1)
    ref1 = np.clip(ref1, 0, np.inf) #if negative due to rounding errors, cutoff at 0
    ap2 = trrecord.vcfrecord.format('AP2')
    ref2 = 1 - np.sum(ap2, axis=1)
    ref2 = np.clip(ref2, 0, np.inf)

    lens = np.array(trrecord.alt_allele_lengths)
    ref_len = trrecord.ref_allele_length
    min_len = min(min(lens), ref_len)
    max_len = max(max(lens), ref_len)
    lens = lens.reshape((lens.shape[0],1))

    h1_probs = np.dot(ap1, lens)
    h2_probs = np.dot(ap2, lens)
    ref1_prob = ref1 * ref_len
    ref2_prob = ref2 * ref_len
    ref1_prob = ref1_prob.reshape(ref1_prob.shape[0], 1)
    ref2_prob = ref2_prob.reshape(ref2_prob.shape[0], 1)
    d = h1_probs + h2_probs + ref1_prob + ref2_prob
    d = np.round(d, 2)
    d_transf = (d - 2*min_len) / (2*max_len - 2*min_len) 
    assert(np.all(d_transf<= 1) and np.all(d_transf >=0))
    print(d_transf.shape)

    genotypes = v_get_gt(d_transf) 

    trrecord.vcfrecord.FORMAT = "GT:DS"

    dosage_str = ":".join([f"{gt}:{ds:.2f}" for gt, ds in zip(genotypes, d_transf)])
    trrecord.vcfrecord.genotypes = [dosage_str]

    # Write variant to output VCF file
    vcf_writer.write_record(trrecord.vcfrecord)

    str_id = trrecord.record_id
    min_max_len.append([str_id, min_len, max_len])

    #break

np.savetxt(min_max_len_file, np.array(min_max_len))
