#!/usr/bin/python
# (c) 2018-2023 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from scipy import stats

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

# valid letters
valid_letters = set(['A', 'C', 'G', 'T'])

# main function
def main():

    all_pop = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    all_pop = ['AFR', 'AMR']

    # get command line argument and initialize log
    args = get_command_line()

    # parse command line 
    argmap = check_command_line(args)

    # load gwas sumstat
    sumstat = pd.read_table(argmap['sumstat'],
        delim_whitespace=True, nrows=100000)
    sumstat.drop(['SNP'], axis=1, inplace=True)
    sumstat.rename(columns={'#CHR':'CHR', 'POS':'BP', 'REF':'A2', 'ALT':'A1',
        'rsid':'SNP'}, inplace=True)

    # find significant hit
    sumstat = sumstat[sumstat['all_inv_var_meta_p']<argmap['thres']]

    # get all relevant chromosomes
    all_chrom = pd.unique(sumstat['CHR'])

    # iterate through pop and chrom
    pop_sumstat = dict()
    
    # iterate through pop
    for pop in all_pop:
        
        # iterate through chrom
        sumstat_chrom = []
        for chrom in all_chrom:
            frq_f = '{}/{}/1000G.{}.QC.{}.frqx'.format(argmap['frqdir'],
                pop, pop, chrom)
            frq = pd.read_table(frq_f, sep='\t')
            frq.rename(columns={'A2':'A1', 'A1':'A2'}, inplace=True)
            sumstat_chrom.append(sumstat.merge(frq, on=['SNP']))
        
        # concatenate across chromosomes
        sumstat_tmp = pd.concat(sumstat_chrom, ignore_index=True)
        homa1 = sumstat_tmp['C(HOM A2)'] # plink encoding is reversed
        het = sumstat_tmp['C(HET)']
        homa2 = sumstat_tmp['C(HOM A1)'] # plink encoding is reversed
        tot = homa1 + het + homa2

        # get frequency of the effect allele
        nsnp = sumstat_tmp.shape[0]
        freq = np.zeros(nsnp)
        ss_a1a2 = (sumstat_tmp['A1_x'] + sumstat_tmp['A2_x']).values
        ref_a1a2 = (sumstat_tmp['A1_y'] + sumstat_tmp['A2_y']).values
        for i in range(nsnp):
            if ss_a1a2[i] in equiv[ref_a1a2[i]]:
                freq[i] = (homa1[i] + het[i]) / tot[i]
            elif ss_a1a2[i] in reverse[ref_a1a2[i]]:
                freq[i] = (homa2[i] + het[i]) / tot[i]
            else:
                freq[i] = np.nan

        # update freq
        sumstat_tmp['{}_A1_FREQ'.format(pop)] = freq
       
        # select columns to output
        sumstat_tmp = sumstat_tmp[['SNP', 'CHR_x', 'BP', 'A1_x', 'A2_x',
            'all_inv_var_meta_beta', 'all_inv_var_meta_sebeta',
            'all_inv_var_meta_p', '{}_A1_FREQ'.format(pop)]]
        sumstat_tmp.rename(columns={'CHR_x':'CHR', 'A1_x':'A1', 'A2_x':'A2'},
            inplace=True)

        # update pop sumstat
        pop_sumstat[pop] = sumstat_tmp

    # create output
    sumstat_out = pop_sumstat[all_pop[0]]
    for i in range(1, len(all_pop)):
        sumstat_out = sumstat_out.merge(pop_sumstat[all_pop[i]],
            how='outer', on=['SNP'])

    # output columns
    out_cols = ['SNP', 'CHR_x', 'BP_x', 'A1_x', 'A2_x',
        'all_inv_var_meta_beta_x', 'all_inv_var_meta_sebeta_x',
        'all_inv_var_meta_p_x']
    for pop in all_pop:
        out_cols.append('{}_A1_FREQ'.format(pop))

    # get final output
    sumstat_out = sumstat_out[out_cols]
    sumstat_out.rename(columns={'CHR_x':'CHR', 'BP_x':'BP', 'A1_x':'A1',
        'A2_x':'A2', 'all_inv_var_meta_beta_x':'BETA',
        'all_inv_var_meta_sebeta_x':'SE', 'all_inv_var_meta_p_x':'P'},
        inplace=True)
    print(sumstat_out.to_markdown())


# check command line
def check_command_line(args):

    # parse command line
    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    return argmap

# get command line input
def get_command_line():
    
    # create the help document
    parser = argparse.ArgumentParser(description='Get allele frequencies'
        'of COVID-19 GWAS associations in different populations')

    parser.add_argument('--sumstat', dest='sumstat', type=str,
        required=False, help='Summary statistics file')

    parser.add_argument('--frqdir', dest='frqdir', type=str,
        required=False, help='Prefix of the frequency files')

    parser.add_argument('--thres', dest='thres', type=float,
        required=False, help='Threshold for significance')

    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
