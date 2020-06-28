import gzip, time, sys
import numpy as np
import pandas as pd
from tqdm import tqdm

# define equivalent alleles
equiv = dict()
equiv["AC"] = set(["TG", "AC", "TC", "AG"])
equiv["AG"] = set(["TC", "AG", "TG", "AC"])
equiv["CA"] = set(["GT", "CA", "GA", "CT"])
equiv["CT"] = set(["GA", "CT", "GT", "CA"])
equiv["TC"] = set(["AG", "TC", "AC", "TG"])
equiv["TG"] = set(["AC", "TG", "AG", "TC"])
equiv["GA"] = set(["CT", "GA", "CA", "GT"])
equiv["GT"] = set(["CA", "GT", "CT", "GA"])

# define reversed alleles
reverse = dict()
reverse["AC"] = set(["GT", "CA", "CT", "GA"])
reverse["AG"] = set(["CT", "GA", "GT", "CA"])
reverse["CA"] = set(["TG", "AC", "AG", "TC"])
reverse["CT"] = set(["AG", "TC", "TG", "AC"])
reverse["TC"] = set(["GA", "CT", "CA", "GT"])
reverse["TG"] = set(["CA", "GT", "GA", "CT"])
reverse["GA"] = set(["TC", "AG", "AC", "TG"])
reverse["GT"] = set(["AC", "TG", "TC", "AG"])

# define strand ambiguous snp
ambiguous = set(["AT", "CG", "TA", "GC"])

# load legend
legend_fnm = sys.argv[1]
legend = pd.read_csv(legend_fnm, header=None, delim_whitespace=True)
legend.columns = ['CHR', 'SNP', 'CM', 'BP', 'A2_ref', 'A1_ref']

# filter legend
legend.drop_duplicates(subset='SNP', keep=False, inplace=True)
legend.drop_duplicates(subset='BP', keep=False, inplace=True)
legend.reset_index(drop=True, inplace=True)

# load sumstats
sumstats_fnm = sys.argv[2]
sumstats = pd.read_csv(sumstats_fnm, delim_whitespace=True)
sumstats.rename(columns={'#CHR': 'CHR', 'POS': 'BP'}, inplace=True)

# get z-score
sumstats['Z'] = sumstats['all_inv_var_meta_beta'] / sumstats['all_inv_var_meta_sebeta']

# merge sumstats and legend
sumstats = sumstats.merge(legend, on=['CHR', 'BP'])
sumstats.sort_values(by=['CHR', 'BP'], inplace=True)

# filter and flip alleles
nsnp = sumstats.shape[0]

ss_BP = sumstats['BP'].values
ss_Z = sumstats['Z'].values; ss_beta = sumstats['all_inv_var_meta_beta'].values
ss_A1 = sumstats['ALT'].values; ss_A2 = sumstats['REF'].values

ref_BP = sumstats['BP'].values
ref_A1 = sumstats['A1_ref'].values; ref_A2 = sumstats['A2_ref'].values

drop_idx = []
for idx in range(nsnp):

    # not bi-allelic
    if len(ss_A1[idx]) > 1 or len(ss_A2[idx]) > 1:
        drop_idx.append(idx)
        continue

    # different BP
    if ss_BP[idx] != ref_BP[idx]:
        drop_idx.append(idx)
        continue

    # get alleles
    ss_alleles  = ss_A1[idx]  + ss_A2[idx]
    ref_alleles = ref_A1[idx] + ref_A2[idx]

    # strand ambiguous
    if (ss_alleles in ambiguous) or (ref_alleles in ambiguous):
        drop_idx.append(idx)
        continue

    # equivalent
    if ss_alleles in equiv[ref_alleles]:
        ss_A1[idx] = ref_A1[idx]
        continue

    # reversed
    if ss_alleles in reverse[ref_alleles]:
        ss_A1[idx]    = ref_A1[idx]
        ss_A2[idx]    = ref_A2[idx]
        ss_Z[idx]     = -1.0 * ss_Z[idx]
        ss_beta[idx]  = -1.0 * ss_beta[idx]
        continue

sumstats.drop(drop_idx, inplace=True)

# get sample size
tot_ncase = 3199.0
tot_nctrl = 897488.0
tot_n = 4.0 / (1.0/tot_ncase + 1.0/tot_nctrl)

sumstats['N'] = tot_n

# save to file
sumstats = sumstats[['SNP_y', 'CHR', 'BP', 'ALT', 'REF', 'Z', 'N',
                     'all_inv_var_meta_beta', 'all_inv_var_meta_sebeta']]
sumstats.columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N', 'BETA', 'SE']
sumstats.to_csv(sys.argv[3], sep='\t', index=False)
