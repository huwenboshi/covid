import sys

import numpy as np
import pandas as pd

thres = 1e-5
window_size = 500000

sumstats1_fnm = sys.argv[1]
out_fnm = sys.argv[2]

sumstats1 = pd.read_csv(sumstats1_fnm, delim_whitespace=True)

signif1_bp = sumstats1[sumstats1['P'] < thres][['SNP', 'CHR', 'BP']]

shared_snp = signif1_bp['SNP']
shared_bp = signif1_bp[signif1_bp['SNP'].isin(shared_snp)][['CHR', 'BP']]

shared_bp['START'] = np.fmax(shared_bp['BP']-window_size, 0)
shared_bp['STOP'] = np.fmax(shared_bp['BP']+window_size, 0)

shared_bp['CHR'] = 'chr' + shared_bp['CHR'].astype(str)
shared_bp['START'] = shared_bp['START'].astype(int)
shared_bp['STOP'] = shared_bp['STOP'].astype(int)

shared_bp = shared_bp[['CHR', 'START', 'STOP']]
shared_bp.to_csv(out_fnm, sep='\t', index=False, header=False)
