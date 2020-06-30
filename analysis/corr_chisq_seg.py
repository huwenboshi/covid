#!/usr/bin/python
# (c) 2018-2023 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from scipy import stats

# main function
def main():

    # get command line argument and initialize log
    args = get_command_line()

    # parse command line 
    argmap = check_command_line(args)

    # load frequency file
    start_chrom, stop_chrom = argmap['use-chrom']
    freq = load_frqfile(args.frqfile, start_chrom, stop_chrom)
    freq = freq[freq['MAF'] > 0.05].reset_index(drop=True)

    # load seg annotations
    annot = load_annot(argmap['seg-annot'], start_chrom, stop_chrom)
    annot = annot.merge(freq, on=['SNP'])
    annot = annot.reset_index(drop=True)

    # load gwas sumstat
    sumstat = pd.read_table(argmap['sumstat'], delim_whitespace=True)
    sumstat['CHISQ'] = np.square(sumstat['Z'])

    # merge with annotation
    sumstat_annot = sumstat.merge(annot, on=['SNP'])
    sumstat_annot = sumstat_annot.reset_index(drop=True)
    
    # get rankr and se
    chisq = sumstat_annot['CHISQ'].values
    seg = sumstat_annot['{}.bed'.format(argmap['tissue'])].values
    rankr, rankr_se = cor_chisq_seg(chisq, seg)
   
    # print out result
    print(argmap['tissue'], rankr, rankr_se)

# get rankr and se
def cor_chisq_seg(chisq, seg):

    # get overall estimate
    rankr = stats.spearmanr(chisq, seg)[0]

    # get jackknife pseudo values
    nsnp = chisq.shape[0]
    nblock = 200
    ps_rankr = []
    blocks = create_block(0, nsnp-1, nblock)
    for blk in blocks:
        chisq_blk = np.delete(chisq, blk)
        seg_blk = np.delete(seg, blk)
        ps_rankr.append(stats.spearmanr(chisq_blk, seg_blk)[0])
    ps_rankr = np.array(ps_rankr)

    # get standard error
    mean_ps_rankr = np.mean(ps_rankr, axis=0)
    diffsq = np.square(ps_rankr - mean_ps_rankr)
    rankr_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0)) 

    return rankr, rankr_se

# load annotation
def load_annot(prefix, start_chrom, end_chrom):

    annot = pd.DataFrame()
    for i in range(start_chrom, end_chrom+1):
        filename = '{}{}.annot.gz'.format(prefix, i)
        tmp = pd.read_table(filename, delim_whitespace=True)
        annot = pd.concat([annot, tmp], ignore_index=True)

    return annot

# intersect snp
def intersect_snp_estimation(annot_snps, freq, min_maf):
    annot_idx = np.where(freq['MAF']>min_maf)[0]
    return annot_idx

# load minor allele frequency file
def load_frqfile(frqfile_fnm, start_chrom, stop_chrom):
    
    all_frq = []
    for i in range(start_chrom, stop_chrom+1):
        all_frq.append(pd.read_table('{}{}.frq'.format(frqfile_fnm,i),
            delim_whitespace=True))
    
    all_frq = pd.concat(all_frq, axis=0, ignore_index=True)

    return all_frq

# create blocks for block jackknife
def create_block(start_idx, stop_idx, nblock):

    block_size = int(np.ceil(float(stop_idx-start_idx+1)/float(nblock)))
    cuts = range(start_idx, stop_idx, block_size)
    blocks = []
    
    for i in range(len(cuts)-1):
        start = cuts[i]
        stop = cuts[i+1]
        blocks.append(np.arange(start, stop))
        
    blocks.append(np.arange(cuts[len(cuts)-1], stop_idx+1))

    return blocks

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
    parser = argparse.ArgumentParser(description='Get correlation between'
        'chisq statistics and SEG annotation')

    parser.add_argument('--tissue', dest='tissue', type=str,
        required=False, help='Tissue name')

    parser.add_argument('--sumstat', dest='sumstat', type=str,
        required=False, help='Baseline annotation file')

    parser.add_argument('--seg-annot', dest='seg-annot', type=str,
        required=False, help='SEG annotation file')

    parser.add_argument('--use-chrom', dest='use-chrom', type=int,
        required=False, default=[1, 22], nargs=2,
        help='Specify which chromosomes to use')

    parser.add_argument('--frqfile', dest='frqfile', type=str,
        required=False, help='Prefix of the frequency files')

    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
