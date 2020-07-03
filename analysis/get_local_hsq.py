import argparse, math, sys, logging, time, gzip

import numpy as np
import pandas as pd

from pandas_plink import read_plink 

# main function
def main():
   
    # maf threshold
    maf_thres = 0.05

    # get command line
    args = get_command_line()

    # load regions
    regions = pd.read_table(args.regions, header=None)
    regions.columns = ['CHR', 'START', 'STOP']
    regions['CHR'] = regions['CHR'].str.replace('chr', '').astype(int)
    all_chrom = pd.unique(regions['CHR'])

    # load legend maf file
    maf = load_frqfile(args.bfile, 1, 22)
    legend = load_legend(args.bfile, 1, 22)
    legmaf = []
    for i in range(22):
        tmp_df = legend[i].merge(maf[i], on=['CHR', 'SNP'])
        legmaf.append(tmp_df)

    # load sumstat
    sumstat = pd.read_table(args.sumstat)

    # iterate through regions
    out = []
    nregions = regions.shape[0]
    for i in range(nregions):
       
        # extract region sumstat
        chrom = regions.iloc[i]['CHR']
        start = regions.iloc[i]['START']
        stop = regions.iloc[i]['STOP']
        sumstat_region = sumstat[(sumstat['CHR']==chrom) &
                                 (sumstat['BP']>=start) &
                                 (sumstat['BP']<=stop)]
       
        # extract region maf
        legmaf_chrom = legmaf[chrom-1]
        legmaf_region = legmaf_chrom[(legmaf_chrom['BP']>=start) &
                             (legmaf_chrom['BP']<=stop) &
                             (legmaf_chrom['MAF']>=maf_thres) &
                             (legmaf_chrom['SNP'].isin(sumstat_region['SNP']))]
      
        # extract region zscore
        zscore_region = sumstat_region[sumstat_region['SNP']\
            .isin(legmaf_region['SNP'])]['Z'].values
        nsample_region = sumstat_region[sumstat_region['SNP']\
            .isin(legmaf_region['SNP'])]['N'].values
        beta_region = zscore_region / np.sqrt(nsample_region)
        nsample = np.mean(nsample_region)

        # load genotype data
        bim, fam, bed = read_plink(args.bfile+str(chrom), verbose=False)
        geno = bed[legmaf_region.index.values,:].compute()
        geno = standardize(geno)
        
        # estimate local hsq
        ld = np.corrcoef(geno)
        hsq, hsq_se = estimate_localhsq(ld, beta_region, nsample)
        
        # append to output
        out.append([chrom, start, stop, hsq, hsq_se])

    # write to output
    out = pd.DataFrame(out)
    out.columns = ['CHR', 'START', 'STOP', 'HSQ', 'HSQ_SE']
    out.to_csv(args.out, index=False, sep='\t')

# load minor allele frequencies
def load_frqfile(frqfile_fnm, start_chrom, stop_chrom):
    
    all_frq = []
    for i in range(start_chrom, stop_chrom+1):
        all_frq.append(pd.read_table('{}{}.frq'.format(frqfile_fnm,i),
            delim_whitespace=True))
    
    return all_frq

# load minor allele frequencies
def load_legend(frqfile_fnm, start_chrom, stop_chrom):
    
    all_leg = []
    for i in range(start_chrom, stop_chrom+1):
        tmp_df = pd.read_table('{}{}.bim'.format(frqfile_fnm,i),
            delim_whitespace=True, header=None)
        tmp_df.columns = ['CHR', 'SNP', 'CM', 'BP', 'A2', 'A1']
        all_leg.append(tmp_df)
    
    return all_leg

# estimate local hsq
def estimate_localhsq(ld, beta, nsample, k=50):
  
    # get initial estimate
    ld_w, ld_v = eig_decomp(ld)
    hsq = 0.0
    for i in range(k):
        hsq += 1/ld_w[i] * np.square(np.dot(ld_v[:,i], beta))

    # correct for bias
    hsq = (nsample * hsq - k) / (nsample - k)
    hsq_var = 4/np.square(nsample-k)*hsq + 2*k/np.square(nsample-k)
    hsq_se = np.sqrt(hsq_var)
    
    return hsq, hsq_se

# get eigenvalue decomposition
def eig_decomp(ld):
    ld_w, ld_v = np.linalg.eigh(ld)
    idx = ld_w.argsort()[::-1]
    ld_w = ld_w[idx]
    ld_v = ld_v[:,idx]
    return (ld_w.T, ld_v)

# standardize columns
def standardize(mat):
    
    mat -= mat.mean(axis=0)
    mat /= (mat.std(axis=0) + 1e-8)
    mat = np.nan_to_num(mat)

    return mat

# get command line input
def get_command_line():

    parser = argparse.ArgumentParser(description='Estimate local hsq')

    parser.add_argument('--bfile', dest='bfile', type=str,
                        help='Genotype file', required=True)

    parser.add_argument('--regions', dest='regions', type=str,
                        help='Region file', required=True)

    parser.add_argument('--sumstat', dest='sumstat', type=str,
                        help='Sumstat file', required=True)

    parser.add_argument('--out', dest='out', type=str,
                        help='Output file name', required=True)

    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
